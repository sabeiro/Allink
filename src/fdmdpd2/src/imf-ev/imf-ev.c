/*
 * imf-ev.c - calculate eigenvalues of modified Seifert-Langer theory
 * (C) Copryight 2011 Martin Hoemberg <mhoembe@gwdg.de>
 *
 */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <gsl/gsl_eigen.h>
#include <string.h>

#define SQR(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))

/* parameters of the modified Seifert-Langer theory */
struct param
{
	double b;		/* intermonolayer friction */
	double eta;		/* surface viscosity */
	double km;		/* monolayer compressibility */
	double kappa;		/* bending rigidity */
	double mPhi0;		/* lipid mass density */
	double d;		/* dist between neutral surface and midplane */
	double gamma;		/* tension */
	double eta_w;		/* solvent viscosity */
	double Gamma;		/* normal friction coefficient */
};

static int calc_discriminant(gsl_matrix *A)
{
	int i, j;
	double x[9];

	for (i = 0; i < 3; i++) for (j = 0; j < 3; j++) {
		x[3 * i + j] = gsl_matrix_get(A, i, j);
	}

	double a = -1. * (x[8] + x[4] + x[0]);
	double b = -1. * (x[6] * x[2] + x[7] * x[5] - x[8] * x[4] - x[8] * x[0] + x[3] * x[1] - x[4] * x[0]);
	double c = x[7] * x[5] * x[0] + x[6] * x[2] * x[4] - x[7] * x[3] * x[2] - x[6] * x[1] * x[5]
		+ x[8] * x[3] * x[1] - x[8] * x[4] * x[0];
	
	double P = (b - a * a / 3.) / 3.;
	double Q = (2. * CUBE(a) / 27. - a * b / 3. + c) / 2.;

	double D = CUBE(P) + SQR(Q);
	if (D > 0.)
		return 1;
	else if (D < 0.)
		return -1;
	else if (D == 0. && Q == 0.)
		return 3;
	else
		return 2;
}

static void findroots(complex double a, complex double b, complex double c,
							complex double e[])
{
	int i;
	complex double p = (b - a * a / 3.) / 3.;
	complex double q = (2. * a * a * a / 27. - a * b / 3. + c) / 2.;
	complex double D = p * p * p + q * q;

	complex double upr = -q + csqrt(D);
	double mod_upr = pow(cabs(upr), 1. / 3.);
	double arg_upr = carg(upr);
	complex double umr = -q - csqrt(D);
	double mod_umr = pow(cabs(umr), 1. / 3.);
	double arg_umr = carg(umr);

	complex double rp = .5 * (-1. + I * sqrt(3.));
	complex double rm = .5 * (-1. - I * sqrt(3.));

	complex double up = mod_upr * cexp(I * arg_upr / 3.);
	for (i = 0; i < 3; i++) {
		complex double um = mod_umr * cexp(I*(2. * M_PI * i + arg_umr) / 3.);
		if (cabs((um * up + p) / p) < 1e-7) {
			e[0] = up + um - a / 3.;
			e[1] = rp * up + rm * um - a / 3.;
			e[2] = rm * up + rp * um - a / 3.;
			return;
		}
	}
	fprintf(stderr, "This should never happen!\n");
	fprintf(stderr, "up=%lg%+lg*I p=%lg%+lg*I q=%lg%+lg*I D=%lg%+lg+I\n",
		creal(up), cimag(up), creal(p), cimag(p),
		creal(q), cimag(q), creal(D), cimag(D));
}

static complex double mget(gsl_matrix_complex *A, int i, int j)
{
	gsl_complex g = gsl_matrix_complex_get(A, i, j);
	return GSL_REAL(g) + I * GSL_IMAG(g);
}

static void calc_ev(gsl_matrix_complex *A, complex double e[])
{
	complex double m[9];
       	m[0] = mget(A, 0, 0); /* a */
	m[1] = mget(A, 0, 1); /* b */
	m[2] = mget(A, 0, 2); /* c */
	m[3] = mget(A, 1, 0); /* d */
	m[4] = mget(A, 1, 1); /* e */
	m[5] = mget(A, 1, 2); /* f */
	m[6] = mget(A, 2, 0); /* g */
	m[7] = mget(A, 2, 1); /* h */
	m[8] = mget(A, 2, 2); /* i */

	complex double a = -1. * (m[0] + m[4] + m[8]);
	complex double b = m[0] * m[8] - m[2] * m[6] - m[1] * m[3] + m[4] * m[8] + m[0] * m[4] - m[5] * m[7];
	complex double c = m[6] * m[2] * m[4] - m[0] * m[4] * m[8] - m[2] * m[3] * m[7] +
	       	m[1] * m[3] * m[8] + m[0] * m[5] * m[7] - m[1] * m[6] * m[5];

	findroots(a, b, c, e);
}

static void mset(gsl_matrix_complex *A, int i, int j, complex double z)
{
	gsl_complex y;
	GSL_REAL(y) = creal(z);
	GSL_IMAG(y) = cimag(z);
	gsl_matrix_complex_set(A, i, j, y);
}

static int calc_dispersion(struct param *p, gsl_eigen_nonsymm_workspace *w,
					complex double dispersion[], int *D, double q)
{
	gsl_matrix_complex *A = gsl_matrix_complex_calloc(3, 3);
//	double omega = sqrt(p->kappa * SQR(SQR(q)) / p->mPhi0);
//	double tau = (2. * p->b + p->eta * SQR(q)) / (p->km * SQR(q));
//	double P = 2. * p->km / (SQR(q) * p->kappa);
//	double eps = 2. * SQR(p->d) * p->km / p->kappa;

	double tkappa = p->kappa + 2 * SQR(p->d) * p->km;
	double q2 = SQR(q);
	double q4 = SQR(q2);
	mset(A, 0, 1, -1.);
	complex double factor = (p->Gamma == 0.) ? 1. : 1. + (tkappa * q4 + p->gamma * q2) / p->Gamma;
	mset(A, 1, 0, (tkappa * q4 + p->gamma * q2) / factor / p->mPhi0);
	mset(A, 1, 1, 4. * p->eta_w * q / p->mPhi0);
	mset(A, 1, 2, -2. * p->km * p->d * q2 / p->mPhi0);
	mset(A, 2, 0, -1. * p->km * p->d * q4 / (2. * p->b + 2. * p->eta_w * q + p->eta * q2));
	mset(A, 2, 2, p->km * q2 / (2. * p->b + 2. * p->eta_w * q + p->eta * q2));

	/* calculate discriminant to check number of complex eigenvalues */
//	*D = calc_discriminant(A);

//	gsl_matrix_set(A, 0, 2, -1.);
//	gsl_matrix_set(A, 1, 0, -eps / (p->d * P * tau));
//	gsl_matrix_set(A, 1, 1, 1. / tau);
//	gsl_matrix_set(A, 2, 0, (1. + eps) * SQR(omega));
//	gsl_matrix_set(A, 2, 1, -SQR(omega) * P * p->d);

	calc_ev(A, dispersion);

//	ret = gsl_eigen_nonsymm_complex(A, ev, w);
//	if (ret) {
//		fprintf(stderr, "Calculation of eigenvalues failed. Only %d "
//			"eigenvalues converged.\n\n", (int)w->n_evals);
//		return 1;
//	}
//	dispersion[0] = dispersion[1] = dispersion[2] = -1.;
//	for (i = 0; i < 3; i++) {
//		gsl_complex e = gsl_vector_complex_get(ev, i);
//		dispersion[i] = GSL_REAL(e) + I * GSL_IMAG(e);
//		if (GSL_IMAG(e) == 0.) {
//			dispersion[0] = GSL_REAL(e);
//		} else if (GSL_IMAG(e) > 0.) {
//			dispersion[1] = GSL_REAL(e);
//			dispersion[2] = GSL_IMAG(e);
//		}
//	}
//
//	for (i = 0; i < 3; i++) {
//		if (dispersion[i] < 0.) {
//			printf("dispersion[%d] < 0! calculation failed.\n", i);
//			return 1;
//		}
//	}

	gsl_matrix_complex_free(A);
//	gsl_vector_complex_free(ev);
	return 0;
}	

int main(int argc, char **argv)
{
	struct param p;
	gsl_eigen_nonsymm_workspace *w;
	complex double dispersion[3]; /* gamma1, gamma2, omegaB */
	double q;

	if (argc < 10) {
		fprintf(stderr, "Syntax: %s [b] [eta] [km] [kappa] [phi] [d] "
			"[gamma] [eta_w] [Gamma]\n\n",
			argv[0]);
		return EXIT_FAILURE;
	}	

	/* setup */	
	w = gsl_eigen_nonsymm_alloc(3);
	gsl_eigen_nonsymm_params(0, 1, w);
	p.b = atof(argv[1]); // 1.7
	p.eta = atof(argv[2]); // 66.2
	p.km = atof(argv[3]); // 1.04
	p.kappa = atof(argv[4]); // 18
	p.mPhi0 = atof(argv[5]); // 28.56
	p.d = atof(argv[6]); // 2.0
	p.gamma = atof(argv[7]); // 0.1
	p.eta_w = atof(argv[8]); // 0.0
	p.Gamma = atof(argv[9]); // 0.0

	printf("# b=%lg eta=%lg km=%lg kappa=%lg mPhi0=%lg d=%lg gamma=%lg\n",
		p.b, p.eta, p.km, p.kappa, p.mPhi0, p.d, p.gamma);

	/* calculation */
	for (q = 1e-8; q < 1e3; q *= 1.01) {
//	for (q = 1.; q < 1e11; q *= 1.01) {
		int D;
		calc_dispersion(&p, w, dispersion, &D, q);
		printf("%lg %d %lg %lg %lg %lg %lg %lg\n", q, D,
				creal(dispersion[0]), cimag(dispersion[0]),
				creal(dispersion[1]), cimag(dispersion[1]),
				creal(dispersion[2]), cimag(dispersion[2]));
	}

	/* cleanup */
	gsl_eigen_nonsymm_free(w);
	return EXIT_SUCCESS;
}

