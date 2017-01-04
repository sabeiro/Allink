/*
 * imf2.c - calculate eigenvalues of modified Seifert-Langer theory
 * (C) Copryight 2011 Martin Hoemberg <mhoembe@gwdg.de>
 */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <assert.h>

#define SQR(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))
#define MAX(X,Y)        (((X) > (Y)) ? (X) : (Y))
#define MIN(X,Y)        (((X) < (Y)) ? (X) : (Y))


/* parameters of the modified Seifert-Langer theory */
struct param
{
	double b;		/* intermonolayer friction */
	double km;		/* monolayer compressibility */
	double km1;		/* dilational viscosity */
	double kappa;		/* bending rigidity */
	double phi;		/* bilayer mass density */
	double d;		/* dist between neutral surface and midplane */
	double gamma;		/* tension */
	double eta;		/* solvent viscosity */
	double mu;		/* monolayer surface viscosity */
};

static int cmp(const void *a, const void *b)
{
	double x = *(double*)a;
	double y = *(double*)b;

	if (x < y) return -1;
	if (x > y) return 1;
	return 0;
}

static void findroots(complex double a, complex double b, complex double c,
							complex double e[])
{
	int i;
	complex double p = (b - a * a / 3.) / 3.;
	complex double q = (2. * a * a * a / 27. - a * b / 3. + c) / 2.;
	complex double D = creal(p * p * p + q * q); /* MODEL SPECIFIC!!! */
	double eps;

	assert(abs(cimag(p)) < 1e-10);
	assert(abs(creal(q)) < 1e-10);
	assert(abs(cimag(p * p * p + q * q)) < 1e-10);

	complex double upr = -q + csqrt(D);
	double mod_upr = pow(cabs(upr), 1. / 3.);
	double arg_upr = carg(upr);
	complex double umr = -q - csqrt(D);
	double mod_umr = pow(cabs(umr), 1. / 3.);
	double arg_umr = carg(umr);

	complex double rp = .5 * (-1. + I * sqrt(3.));
	complex double rm = .5 * (-1. - I * sqrt(3.));
	complex double up = mod_upr * cexp(I * arg_upr / 3.);

	for (eps = 1e-30; eps < 1e-6; eps *= 10.) {
		complex double um[3];
		double sort[3];

		for (i = 0; i < 3; i++) {
			um[i] = mod_umr * cexp(I*(2. * M_PI * i + arg_umr) / 3.);
			sort[i] = cabs((um[i] * up + p) / p);
		}
		qsort(sort, 3, sizeof(*sort), cmp);

		for (i = 0; i < 3; i++) {
			double test = cabs((um[i] * up + p) / p);
			if (test == sort[0] && test < eps) {
				e[0] = up + um[i] - a / 3.;
				e[1] = rp * up + rm * um[i] - a / 3.;
				e[2] = rm * up + rp * um[i] - a / 3.;
				return;
			}
		}
	}
	fprintf(stderr, "This should never happen!\n");
	fprintf(stderr, "up=%lg%+lg*I p=%lg%+lg*I q=%lg%+lg*I D=%lg\n",
		creal(up), cimag(up), creal(p), cimag(p),
		creal(q), cimag(q), creal(D));
}

static void identifyroots3(double d[], const complex double e[])
{
	int i, recnt = 0;
	double re, im;
	complex double e2[3];
	
	for (i = 0; i < 3; i++) {
		re = fabs(creal(I * e[i]));
		im = fabs(cimag(I * e[i]));
		if (im < 1e-12) im = 0.;
		e2[i] = re + I * im;

		if (im == 0.) recnt++;
	}

	if (recnt == 1) {
		for (i = 0; i < 3; i++) {
			if (cimag(e2[i]) == 0.) {
				d[0] = creal(e2[i]);
			} else {
				d[1] = creal(e2[i]);
				d[2] = cimag(e2[i]);
			}
		}
	} else if (recnt == 2) {
		fprintf(stderr, "Two real dispersion branches dont exist.\n");
	} else {
		d[0] = MIN(MIN(creal(e2[0]), creal(e2[1])), creal(e2[2]));
		d[1] = MAX(MAX(creal(e2[0]), creal(e2[1])), creal(e2[2]));
		d[2] = 0.;
		
	}
	d[3] = recnt;
}

static void calc_dispersion3(struct param *p, double q, complex double e[])
{
	double q2 = q * q;
	double q4 = q2 * q2;
	double tkappa = p->kappa + 2. * SQR(p->d) * p->km;
	double denom = p->phi * (p->km1 * q2 + p->mu * q2 + 2. * p->eta * q + 2. * p->b);
	
	/* prefactor omega^0 */
	complex double c = I * p->km * q4 * (p->kappa * q2 + p->gamma) / denom;

	/* prefactor omega^1 */
	complex double b = -q2 * (q4 * p->kappa * p->km1 + q4 * tkappa * p->mu
		+ 2. * q2 * q * tkappa * p->eta
		+ 2. * q2 * tkappa * p->b + q2 * p->gamma * (p->km1 + p->mu)
		+ 2. * q * p->gamma * p->eta + 4. * p->km * q * p->eta
		+ 2. * p->gamma * p->b) / denom;

	/* prefactor omega^2 */
	complex double a = -I * q * (q4 * q * 2. * SQR(p->d) * p->km1 * p->mu
			+ 4. * SQR(p->d) * q4 * p->km1 * p->eta
			+ 4. * q2 * q * SQR(p->d) * p->km1 * p->b
			+ 4. * p->eta * q2 * p->km1 + 4. * p->eta * q2 * p->mu
			+ 8. * SQR(p->eta) * q + p->phi * p->km * q
			+ 8. * p->eta * p->b) / denom;

	/* calculate roots */
	findroots(a, b, c, e);
}

static void identifyroots2(double d[], const complex double e[])
{
	int i, recnt = 0;
	complex double e2[2];
	
	for (i = 0; i < 2; i++) {
		double re = fabs(creal(I * e[i]));
		double im = fabs(cimag(I * e[i]));
		if (im < 1e-12) im = 0.;
		e2[i] = re + I * im;
		if (im == 0.) recnt++;
	}

	switch (recnt) {
	case 0:
		d[0] = 0.;
		d[1] = creal(e2[0]);
		d[2] = cimag(e2[0]);
		break;
	case 2:
		d[0] = MIN(creal(e2[0]), creal(e2[1]));
		d[1] = MAX(creal(e2[0]), creal(e2[1]));
		d[2] = 0.;
		break;
	default:
		fprintf(stderr, "Invalid number of real dispersions\n");
	}
}

/* simpler dispersion relation for the case when Phi = 0 */
static void calc_dispersion2(struct param *p, double q, complex double e[])
{
	double denom = 2. * p->eta * q * q * (p->km1 + p->mu) +
		4. * p->eta * p->b +
		q*q*q*q*q * SQR(p->d) * p->km1 * p->mu +
		2. * q*q*q*q * SQR(p->d) * p->km1 * p->eta +
		2. * q*q*q * SQR(p->d) * p->km1 * p->b +
		4. * SQR(p->eta) * q;
	double tkappa = p->kappa + 2. * SQR(p->d) * p->km;
	double tmp1 = p->kappa * SQR(q) + p->gamma;
	double tmp2 = tkappa * SQR(q) + p->gamma;
//	printf("tkappa=%lg tmp1=%lg tmp2=%lg denom=%lg\n", tkappa, tmp1, tmp2, denom);

	/* prefactor omega^0 */
	double b = -.5 * q*q*q * p->km * tmp1 / denom;

	/* prefactor omega^1 */
	complex double a = -.5 * I * q * (
			4. * p->eta * q * p->km +
			2. * p->b * tmp2 +
			2. * q * p->eta * tmp2 +
			q*q * p->mu * tmp2 +
			q*q * p->km1 * tmp1) / denom;
//	printf("a=%lg%+lg b=%lg%+lg\n", creal(a), cimag(a), creal(b), cimag(b));

	e[0] = -.5 * a + csqrt(.25 * a * a - b);
	e[1] = -.5 * a - csqrt(.25 * a * a - b);
}

int main(int argc, char **argv)
{
	struct param p;
	double q;
	complex double e[3];
	double d[4];
	int i;
	
	if (argc < 10) {
		fprintf(stderr, "Syntax: %s [b] [mu] [km] [km1] [kappa] [phi] [d] "
			"[gamma] [eta]\n\n",
			argv[0]);
		return EXIT_FAILURE;
	}	

	/* setup */	
	p.b = atof(argv[1]); // 1.7
	p.mu = atof(argv[2]); // 66.2
	p.km = atof(argv[3]); // 1.04
	p.km1 = atof(argv[4]); // 11
	p.kappa = atof(argv[5]); // 18
	p.phi = atof(argv[6]); // 28.56
	p.d = atof(argv[7]); // 2.0
	p.gamma = atof(argv[8]); // 0.1
	p.eta = atof(argv[9]); // 0.0

	printf("#");
	for (i = 0; i < argc; i++)
		printf(" %s", argv[i]);
	printf("\n");
	printf("# b=%lg mu=%lg km=%lg km1=%lg kappa=%lg phi=%lg d=%lg gamma=%lg eta=%lg\n",
		p.b, p.mu, p.km, p.km1, p.kappa, p.phi, p.d, p.gamma, p.eta);

	double x = (p.kappa + 2. * SQR(p.d) * p.km) / p.kappa;
	double y = p.km1 / p.mu;
	printf("# q1*=%lg q2*=%lg\n",
		sqrt(2. * (4. * SQR(p.eta) - p.gamma * p.phi) / (p.phi	* (p.kappa + 2. * SQR(p.d) * p.km))),
		pow(SQR(p.km) * p.phi / (p.kappa * SQR(p.km1)) * (y + 1.) / SQR(x - 1.) * (x + y + sqrt(SQR(x-y)+4*y)), 1./4.));

	/* calculation */
	for (q = 1e-6; q < 1e3; q *= 1.01) {
		if (p.phi > 0.) {
			calc_dispersion3(&p, q, e);
			identifyroots3(d, e);
		} else {
			calc_dispersion2(&p, q, e);
			identifyroots2(d, e);
		}
	
		printf("%lg %lg %lg %lg %lg %lg %lg", q,
				creal(I * e[0]), cimag(I * e[0]),
				creal(I * e[1]), cimag(I * e[1]),
				creal(I * e[2]), cimag(I * e[2]));
		printf(" %lg %lg %lg %lg\n", d[0], d[1], d[2], d[3]);
	}

	return EXIT_SUCCESS;
}

