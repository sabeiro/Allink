/*
 * solver.c - fit modified Seifert-Langer theory to measured data points
 * (C) Copryight 2011 Martin Hoemberg <mhoembe@gwdg.de>
 *
 */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mkl_rci.h>
#include <complex.h>
#include <limits.h>
#include <assert.h>

#define SQR(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))
#define MAX(a,b) ((a)>(b) ? (a) : (b))
#define MIN(a,b) ((a)>(b) ? (b) : (a))
#define ARRAY_SIZE(arr) (sizeof(arr) / sizeof((arr)[0]))

static int n_gamma1 = 0;	/* number of data points for gamma1 */
static int n_gamma2 = 0;	/* number of data points for gamma2 */
static int n_omega2 = 0;	/* number of data points for omega2 */
static double *x_gamma1 = NULL;	/* data points for gamma1 */
static double *x_gamma2 = NULL;	/* data points for gamma2 */
static double *x_omega2 = NULL;	/* data points for omega2 */
static double *q_gamma1 = NULL;	/* q's for gamma1 */
static double *q_gamma2 = NULL;	/* q's for gamma2 */
static double *q_omega2 = NULL;	/* q's for omega2 */

#define FIT_KM		0
#define FIT_MU		1
#define FIT_D		2
#define FIT_KM1		3
#define FIT_KAPPA	4
#define FIT_PHI		5
#define FIT_B		6
#define FIT_GAMMA	7
#define FIT_ETA		8

#define FIT_MAX		9
#define FIT_PARAMS	4

#define ERRCHECK(x)							\
	do {								\
		int __res = (x);					\
	       	if (__res == TR_INVALID_OPTION) {			\
			fprintf(stderr, "dtrnlspbc_init() TR_INVALID_"	\
			"OPTION - there was an error in the input para" \
			"meters.\n\n");					\
			exit(1);					\
		} else if (__res == TR_OUT_OF_MEMORY) {			\
			fprintf(stderr, "dtrnlspbc_init() TR_OUT_OF_M"	\
			"EMORY - there was a memory error.\n\n");	\
			exit(1);					\
		}							\
	} while (0);							



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

	for (eps = 1e-15; eps < 1e-4; eps *= 10.) {
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
		d[2] = 1000.;
		
	}
}

static void calc_dispersion3(double x[], double dispersion[], double q)
{
	complex double e[3];
	double q2 = q * q;
	double q4 = q2 * q2;
	double tkappa = x[FIT_KAPPA] + 2. * SQR(x[FIT_D]) * x[FIT_KM];
	double denom = x[FIT_PHI] * (x[FIT_KM1] * q2 + x[FIT_MU] * q2 + 2. * x[FIT_ETA] * q + 2. * x[FIT_B]);
	
	/* prefactor omega^0 */
	complex double c = I * x[FIT_KM] * q4 * (x[FIT_KAPPA] * q2 + x[FIT_GAMMA]) / denom;

	/* prefactor omega^1 */
	complex double b = -q2 * (q4 * x[FIT_KAPPA] * x[FIT_KM1] + q4 * tkappa * x[FIT_MU]
		+ 2. * q2 * q * tkappa * x[FIT_ETA]
		+ 2. * q2 * tkappa * x[FIT_B] + q2 * x[FIT_GAMMA] * (x[FIT_KM1] + x[FIT_MU])
		+ 2. * q * x[FIT_GAMMA] * x[FIT_ETA] + 4. * x[FIT_KM] * q * x[FIT_ETA]
		+ 2. * x[FIT_GAMMA] * x[FIT_B]) / denom;

	/* prefactor omega^2 */
	complex double a = -I * q * (q4 * q * 2. * SQR(x[FIT_D]) * x[FIT_KM1] * x[FIT_MU]
			+ 4. * SQR(x[FIT_D]) * q4 * x[FIT_KM1] * x[FIT_ETA]
			+ 4. * q2 * q * SQR(x[FIT_D]) * x[FIT_KM1] * x[FIT_B]
			+ 4. * x[FIT_ETA] * q2 * x[FIT_KM1] + 4. * x[FIT_ETA] * q2 * x[FIT_MU]
			+ 8. * SQR(x[FIT_ETA]) * q + x[FIT_PHI] * x[FIT_KM] * q
			+ 8. * x[FIT_ETA] * x[FIT_B]) / denom;

	/* calculate roots */
	findroots(a, b, c, e);
	identifyroots3(dispersion, e);
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
static void calc_dispersion2(double x[], double dispersion[], double q)
{
	complex double e[2];
	double denom = 2. * x[FIT_ETA] * q * q * (x[FIT_KM1] + x[FIT_MU]) +
		4. * x[FIT_ETA] * x[FIT_B] +
		q*q*q*q*q * SQR(x[FIT_D]) * x[FIT_KM1] * x[FIT_MU] +
		2. * q*q*q*q * SQR(x[FIT_D]) * x[FIT_KM1] * x[FIT_ETA] +
		2. * q*q*q * SQR(x[FIT_D]) * x[FIT_KM1] * x[FIT_B] +
		4. * SQR(x[FIT_ETA]) * q;
	double tkappa = x[FIT_KAPPA] + 2. * SQR(x[FIT_D]) * x[FIT_KM];
	double tmp1 = x[FIT_KAPPA] * SQR(q) + x[FIT_GAMMA];
	double tmp2 = tkappa * SQR(q) + x[FIT_GAMMA];

	/* prefactor omega^0 */
	double b = -.5 * q*q*q * x[FIT_KM] * tmp1 / denom;

	/* prefactor omega^1 */
	complex double a = -.5 * I * q * (
			4. * x[FIT_ETA] * q * x[FIT_KM] +
			2. * x[FIT_B] * tmp2 +
			2. * q * x[FIT_ETA] * tmp2 +
			q*q * x[FIT_MU] * tmp2 +
			q*q * x[FIT_KM1] * tmp1) / denom;

	e[0] = -.5 * a + csqrt(.25 * a * a - b);
	e[1] = -.5 * a - csqrt(.25 * a * a - b);

	identifyroots2(dispersion, e);
}

static void calc_dispersion(double x[], double dispersion[], double q)
{
	if (x[FIT_PHI] == 0.) {
		calc_dispersion2(x, dispersion, q);
	} else {
		calc_dispersion3(x, dispersion, q);
	}
}


/* calculate functional */
static void fcn(int *m, int *n, double x[], double f[])
{
	int i;
	double dispersion[3];
	(void) m; /* silence compiler warnings about unused parameters */
	(void) n;

	for (i = 0; i < n_gamma1; i++) {
		calc_dispersion(x, dispersion, q_gamma1[i]);
		f[i] = log(dispersion[0] / x_gamma1[i]);
//		printf("1. f[%d]=%lg\n", i, f[i]);
	}
	for (i = 0; i < n_gamma2; i++) {
		calc_dispersion(x, dispersion, q_gamma2[i]);
		f[n_gamma1 + i] = log(dispersion[1] / x_gamma2[i]);
//		printf("2. f[%d]=%lg\n", n_gamma1 + i, f[n_gamma1 + i]);
	}
	for (i = 0; i < n_omega2; i++) {
		calc_dispersion(x, dispersion, q_omega2[i]);
		f[n_gamma1 + n_gamma2 + i] = log(dispersion[2] / x_omega2[i]);
//		printf("3. f[%d]=%lg\n", n_gamma1 + n_gamma2 + i, f[n_gamma1 + n_gamma2 + i]);
	}
}

static int read_file(char *fn, int *n, double **_x, double **_q)
{
	int cnt = 0;

	FILE *FH = fopen(fn, "r");
	if (FH == NULL) {
		fprintf(stderr, "Couldn't open '%s' for reading.\n", fn);
		exit(1);
	}

	/* count number of entries */
	while (!feof(FH)) {
		char line[LINE_MAX];
		char *s = fgets(line, LINE_MAX, FH);

		if (s) {
			s = strchr(line, '#');
			if (s == line) continue;
			cnt++;
		}
	}
	printf("%d data points found in file '%s'.\n", cnt, fn);

	(*n) = cnt;
	(*_x) = calloc(cnt, sizeof(**_x));
	(*_q) = calloc(cnt, sizeof(**_q));
	rewind(FH);

	cnt = 0;
	while (!feof(FH)) {
		char line[LINE_MAX];
		char *s = fgets(line, LINE_MAX, FH);
		double q, x;

		if (s) {
			s = strchr(line, '#');
			if (s == line) continue;

			if (sscanf(line, "%lg %lg\n", &q, &x) != 2) {
				fprintf(stderr, "Error while reading line "
						"'%s'\n", line);
				exit(1);
			}
			(*_x)[cnt] = x;
//			(*_q)[cnt] = sqrt(q);
			(*_q)[cnt] = q;
			cnt++;
		}
	}

	fclose(FH);
	return 0;
}

static void solve(int n, double x[], double lower[], double upper[])
{
	_TRNSPBC_HANDLE_t handle = 0L;
	int m = n_gamma1 + n_gamma2 + n_omega2;
	double eps[6] = {1e-12, 1e-12, 1e-12, 1e-12, 1e-12, 1e-12};
	double fjac_eps = 1e-8;
	int iter1 = 5000;
	int iter2 = 100;
	double rs = 0.;
	int req;
	double fvec[m];
	double fjac[m * n];
	int step = 0;
	
	ERRCHECK(dtrnlspbc_init(&handle, &n, &m, x, lower, upper, eps, &iter1, &iter2, &rs));

	/* initialize Jacobian and functional vector */
	ERRCHECK(djacobi(fcn, &n, &m, fjac, x, &fjac_eps));
	fcn(&m, &n, x, fvec);
	do {
		printf("Step %d: ", ++step);
		ERRCHECK(dtrnlspbc_solve(&handle, fvec, fjac, &req));

		switch (req) {
		case 2:
			printf("Calc Jacobian\n");
			ERRCHECK(djacobi(fcn, &n, &m, fjac, x, &fjac_eps));
			break;
		case 1:
			printf("Calc Functional\n");
			fcn(&m, &n, x, fvec);
			break;
		case 0:
			printf("Successful iteration\n");
			break;
		case -1:
			printf("Number of iterations exceeded\n");
			break;
		case -2:
			printf("Delta < eps(1)\n");
			break;
		case -3:
			printf("||F(x)||_2 < eps(2)\n");
			break;
		case -4:
			printf("The Jacobian is singular (decrease eps(3)?)\n");
			break;
		case -5:
			printf("||s||_2 < eps(4)\n");
			break;
		case -6:
			printf("||F(x)||_2 - ||F(x) - J(x)s||_2 < eps(5)\n");
			break;
		default:
			printf("This should never happen.\n");
		}
	} while (req >= 0);

	ERRCHECK(dtrnlspbc_delete(&handle));
}

int main(int argc, char **argv)
{
	double x[FIT_MAX], lower[FIT_MAX], upper[FIT_MAX];
	char str[FIT_MAX][20];
	int i, cnt = 1;

	assert(FIT_PARAMS <= FIT_MAX);

	if (argc < 12) {
		fprintf(stderr, "Syntax: %s [b] [mu] [d] [kappa] [phi] [km]"
		" [km1] [Sigma] [gamma1.txt] [gamma2.txt] [omega2.txt]\n\n",
		argv[0]);
		return EXIT_FAILURE;
	}

	/* intermonolayer friction */
	x[FIT_B] = atof(argv[cnt++]);
	lower[FIT_B] = 0.1;
	upper[FIT_B] = 500.;
	strcpy(str[FIT_B], "b");

	/* surface shear viscosity */
	x[FIT_MU] = atof(argv[cnt++]);
	lower[FIT_MU] = 0.;
	upper[FIT_MU] = 500.;
	strcpy(str[FIT_MU], "mu");

	/* distance between neutral surface and bilayer midplane */
	x[FIT_D] = atof(argv[cnt++]);
	lower[FIT_D] = 0.1;
	upper[FIT_D] = 10.;
	strcpy(str[FIT_D], "d");

	/* unnormalized bending rigidity */
	x[FIT_KAPPA] = atof(argv[cnt++]);
	lower[FIT_KAPPA] = 1.;
	upper[FIT_KAPPA] = 50.;
	strcpy(str[FIT_KAPPA], "kappa");

	/* bilayer mass density */	
	x[FIT_PHI] = atof(argv[cnt++]);
	lower[FIT_PHI] = 0.;
	upper[FIT_PHI] = 400.;
	strcpy(str[FIT_PHI], "phi");

	/* monolayer compression modulus */
	x[FIT_KM] = atof(argv[cnt++]);
	lower[FIT_KM] = 0.1;
	upper[FIT_KM] = 100.;
	strcpy(str[FIT_KM], "km");

	/* monolayer dilational viscosity */	
	x[FIT_KM1] = atof(argv[cnt++]);
	lower[FIT_KM1] = 0.1;
	upper[FIT_KM1] = 400.;
	strcpy(str[FIT_KM1], "km1");
	
	/* surface tension */
	x[FIT_GAMMA] = atof(argv[cnt++]);
	lower[FIT_GAMMA] = -20.;
	upper[FIT_GAMMA] = 20.;
	strcpy(str[FIT_GAMMA], "gamma");

	/* solvent shear viscosity */
	x[FIT_ETA] = 0.;
	lower[FIT_ETA] = 0.;
	upper[FIT_ETA] = 0.;
	strcpy(str[FIT_ETA], "eta");
	
	/* check given initial values */
	for (i = 0; i < FIT_MAX; i++) {
		if (x[i] < lower[i]) {
			printf("%s=%lg is too small!\n", str[i], x[i]);
			return EXIT_FAILURE;
		}
		if (x[i] > upper[i]) {
			printf("%s=%lg is too large!\n", str[i], x[i]);
			return EXIT_FAILURE;
		}
	}

	/* read data files */
	read_file(argv[cnt++], &n_gamma1, &x_gamma1, &q_gamma1);
	read_file(argv[cnt++], &n_gamma2, &x_gamma2, &q_gamma2);
	read_file(argv[cnt++], &n_omega2, &x_omega2, &q_omega2);

	printf("Fitting parameters:");
	for (i = 0; i < FIT_PARAMS; i++)
		printf(" %s", str[i]);
	printf("\n");

	/* solve for the solutions */
	solve(FIT_PARAMS, x, lower, upper);

	/* print results */
	printf("Results:\n=====================\n");
	for (i = 0; i < FIT_MAX; i++)
		printf("%s = %lg\n", str[i], x[i]);

	/* [b] [mu] [km] [km1] [kappa] [phi] [d] [gamma] [eta] */
	printf("../imf2 %lg %lg %lg %lg %lg %lg %lg %lg %lg\n",
		x[FIT_B], x[FIT_MU], x[FIT_KM], x[FIT_KM1], x[FIT_KAPPA],
		x[FIT_PHI], x[FIT_D], x[FIT_GAMMA], x[FIT_ETA]);
	return EXIT_SUCCESS;
}


