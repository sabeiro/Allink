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
#include <gsl/gsl_eigen.h>
#include <limits.h>
#include <assert.h>

//#define FIXED_VISCOSITY

#define SQR(x) ((x)*(x))
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

/* global variables for computing the dispersion relations */
static gsl_matrix *A = NULL;
static gsl_eigen_nonsymm_workspace *w = NULL;
static gsl_vector_complex *ev = NULL;

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


/*
 * The four parameters are contained in x[]:
 * x[0] = 2 d^2 km / kappa \
 * x[1] = b / km           | fit parameters
 * x[2] = eta / km         /
 * x[3] = kappa / mPhi0 (fixed value)
 */
static void calc_dispersion(double x[], double dispersion[], double q)
{
	double re[3], im[3], q2 = SQR(q);
	int i;
	
	gsl_matrix_set_zero(A);
	gsl_matrix_set(A, 0, 2, -1.);
	gsl_matrix_set(A, 1, 0, -1. * SQR(q2) / (2. * x[1] + q2 * x[2]));
	gsl_matrix_set(A, 1, 1, q2 / (2. * x[1] + q2 * x[2]));
	gsl_matrix_set(A, 2, 0, (1. + x[0]) * x[3] * SQR(q2));
	gsl_matrix_set(A, 2, 1, -1. * x[0] * q2 * x[3]);

	if (gsl_eigen_nonsymm(A, ev, w)) {
		fprintf(stderr, "Calculation of eigenvalues failed. Only %d "
			"eigenvalues converged.\n\n", (int)w->n_evals);
		exit(1);
	}
	for (i = 0; i < 3; i++) {
		dispersion[i] = -1.;
		re[i] = GSL_REAL(gsl_vector_complex_get(ev, i)); 
		im[i] = GSL_IMAG(gsl_vector_complex_get(ev, i)); 
	}

	if (im[0] == 0. && im[1] == 0. && im[2] == 0.) {
		/* three real eigenvalues */
		fprintf(stderr, "Warning: 3 real eigenvalues.\n");
		dispersion[0] = MAX(MAX(re[0], re[1]), re[2]);
		dispersion[1] = MIN(MIN(re[0], re[1]), re[2]);
		dispersion[2] = 1000.;
	} else {
		/* one real, two complex conjugated eigenvalues */
		for (i = 0; i < 3; i++) {
			if (im[i] == 0.) {
				dispersion[0] = re[i];
			} else if (fabs(im[i]) > 0.) {
				dispersion[1] = re[i];
				dispersion[2] = fabs(im[i]);
			}
		}
	}

	for (i = 0; i < 3; i++) {
		if (dispersion[i] < 0.) {
			fprintf(stderr, "dispersion[%d]<0! calculation failed."
			"\ng1=%lg g2=%lg w2=%lg\new=%lg+I*%lg %lg+I*%lg %lg+I*"
			"%lg\n",i, dispersion[0], dispersion[1], dispersion[2],
		       	re[0], im[0], re[1], im[1], re[2], im[2]);
			exit(1);
		}
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
	}
	for (i = 0; i < n_gamma2; i++) {
		calc_dispersion(x, dispersion, q_gamma2[i]);
		f[n_gamma1 + i] = log(dispersion[1] / x_gamma2[i]);
	}
	for (i = 0; i < n_omega2; i++) {
		calc_dispersion(x, dispersion, q_omega2[i]);
		f[n_gamma1 + n_gamma2 + i] = log(dispersion[2] / x_omega2[i]);
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
	
	A = gsl_matrix_calloc(3, 3);
	ev = gsl_vector_complex_alloc(3);
	w = gsl_eigen_nonsymm_alloc(3);
	gsl_eigen_nonsymm_params(0, 1, w);
	
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
	gsl_matrix_free(A);
	gsl_vector_complex_free(ev);
	gsl_eigen_nonsymm_free(w);
}

int main(int argc, char **argv)
{
	if (argc < 10) {
		fprintf(stderr, "Syntax: %s [b] [eta] [km] [d] [kappa] [mPhi0]"
		" [gamma1.txt] [gamma2.txt] [omega2.txt]\n\n", argv[0]);
		return EXIT_FAILURE;
	}

	double b = atof(argv[1]);     /* b */
	double eta = atof(argv[2]);   /* eta */
	double km = atof(argv[3]);    /* km */
	double d = atof(argv[4]);     /* d */
	double kappa = atof(argv[5]); /* kappa */
	double mPhi0 = atof(argv[6]); /* mPhi0 */

	double x[4] = { 2. * SQR(d) * km / kappa, b / km,
		eta / km, kappa / mPhi0 };
	double lower[4] = { 0., 0., 0., 0. };
	double upper[4] = {100., 20., 20., 10.};

	read_file(argv[7], &n_gamma1, &x_gamma1, &q_gamma1);
	read_file(argv[8], &n_gamma2, &x_gamma2, &q_gamma2);
	read_file(argv[9], &n_omega2, &x_omega2, &q_omega2);

	solve(3, x, lower, upper);

	d = sqrt(x[0] * kappa / 2. / km);
	b = x[1] * km;
	eta = x[2] * km;
	
	printf("Results:\n=====================\n");
	printf("A=%lg B=%lg C=%lg D=%lg\n", x[0], x[1], x[2], x[3]);
	printf("b = %lg\neta = %lg\nkm = %lg\nd = %lg\nkappa = %lg\n",
					b, eta, km, d, kappa);
	printf("%lg %lg %lg %lg %lg %lg\n", b, eta, km, kappa, mPhi0, d);

	return EXIT_SUCCESS;
}

