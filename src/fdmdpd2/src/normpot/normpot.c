/*
 * normpot.c - estimation ofthe normal potential between two coms
 * (C) Copyright 2010 by Martin Hoemberg <mhoembe@gwdg.de>
 *
 * Syntax: ./normpot [L] [profile.txt]
 */

#define _GNU_SOURCE
#include "common.h"
#include <complex.h>
#include <fenv.h>
#include <fftw3.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_blas.h>

static complex double *fft_out;
static double *fft_in;
fftw_plan p_fwd;
fftw_plan p_bwd;


struct data
{
	double *phi;		/* normal profile */
	complex double *tphi;	/* FFT of phi */
	gsl_vector *phi_vec;	/* \hat\phi / N */
	double dx;		/* binwidth */
	double phi_av;		/* number of particles */
	int n;			/* # bins in the profile */
	int p;			/* # degrees of freedom */
	int max;		/* max number of bins */
};

/* ------------------------------------------------------------------------- */

static int load_profile(struct data *d, const char *fn)
{
	FILE *FH = fopen(fn, "r");
	int i;

	if (FH == NULL)
		return error(ENOENT, "Couldn't open '%s' for reading", fn);

	d->phi = NULL;
	d->n = 0;
	d->max = 0;

	while (!feof(FH)) {
		double x;
		char buf[LINE_MAX];

		if (fgets(buf, sizeof(buf), FH) == NULL) break;
		if (strchr(buf, '#') == buf) continue;

		if (d->n >= d->max) {
			d->max += 0x1000;
			d->phi = realloc(d->phi, d->max * sizeof(*d->phi));
			if (d->phi == NULL)
				return novm("d->phi");
		}
		sscanf(buf, "%lg %lg", &x, &d->phi[d->n++]);
		if (d->n == 1)
			d->dx -= x;
		if (d->n == 2)
			d->dx += x;
	}

	/* rescaling so that the phis are dimensionless */
	for (i = 0; i < d->n; i++) {
		d->phi[i] *= d->dx;
		d->phi_av += d->phi[i];
	}
	d->phi_av /= d->n;
	debug("%d bins, dx=%lg, phi_av=%lg", d->n, d->dx, d->phi_av);

	fclose(FH);
	return 0;
}

static int fft_profile(struct data *d)
{
	int k; 

	d->tphi = calloc(d->n / 2 + 1, sizeof(*d->tphi));
	if (d->tphi == NULL) return novm("d->tphi");

	/* setup FFT */	
	fft_out = fftw_malloc((d->n / 2 + 1) * sizeof(*fft_out));
	if (fft_out == NULL) return novm("fft_out");
	fft_in = fftw_malloc(d->n * sizeof(*fft_in));
	if (fft_in == NULL) return novm("fft_in");
	p_fwd = fftw_plan_dft_r2c_1d(d->n, fft_in, fft_out, FFTW_MEASURE);
	p_bwd = fftw_plan_dft_c2r_1d(d->n, fft_out, fft_in, FFTW_MEASURE);

	/* calculate FFT of density profile */
	memset(fft_out, 0, (d->n / 2 + 1) * sizeof(*fft_out));
	memcpy(fft_in, d->phi, d->n * sizeof(*d->phi));
	fftw_execute(p_fwd);
	for (k = 0; k < d->n / 2 + 1; k++)
		d->tphi[k] = fft_out[k];

	assert(d->n % 2 == 0);
	d->p = MIN(d->p, d->n / 2);
	debug("Using d->p=%d dofs", d->p);

	d->phi_vec = gsl_vector_alloc(d->n);
	for (k = 0; k < d->n; k++)
		gsl_vector_set(d->phi_vec, k, d->phi[k] / d->phi_av / d->n);

	return 0;
}

static int setup(struct data *d, int argc, char **argv)
{
	int ret;

	if (argc < 3) {
		fprintf(stderr, "Syntax: %s [p] [profile.txt]\n\n", argv[0]);
		return EINVAL;
	}

	d->p = atoi(*(++argv));
	if (d->p <= 0) return error(EDOM, "0 < p <= n/2");
		
	ret = load_profile(d, *(++argv));
	if (ret) return ret;

	ret = fft_profile(d);
	if (ret) return ret;

	return 0;
}

static void cleanup(struct data *d)
{
	free(d->phi);
	free(d->tphi);
	fftw_free(fft_in);
	fftw_free(fft_out);
	fftw_destroy_plan(p_fwd);
	fftw_destroy_plan(p_bwd);
	gsl_vector_free(d->phi_vec);
}

/* ------------------------------------------------------------------------- */

static void calc_AB(struct data *d, const gsl_vector *x, gsl_vector *A, double *B)
{
	int i;

	fft_out[0] = 0.;
	for (i = 0; i < d->p; i++) {
		int k = i + 1;
		fft_out[k] = gsl_vector_get(x, i) * d->tphi[k];
	}
	for (i = d->p; i < d->n / 2; i++) {
		int k = i + 1;
		fft_out[k] = 0.;
	}
	
	fftw_execute(p_bwd);

	for (i = 0, *B = 0.; i < d->n; i++) {
		assert(fft_in[i] < 1e5);
		double r = exp(-1. * fft_in[i] / d->n);
		gsl_vector_set(A, i, r);
		(*B) += r;
	}
	assert(isnormal(*B));
}

static int eq_f(const gsl_vector *x, void *params, gsl_vector *f)
{
	struct data *d = params;
	double B;

	calc_AB(d, x, f, &B);

	gsl_vector_scale(f, 1. / B);
	gsl_vector_sub(f, d->phi_vec);

	return GSL_SUCCESS;
}

static int eq_df(const gsl_vector *x, void *params, gsl_matrix *J)
{
	struct data *d = params;
	double B;
	int i, j;

	gsl_vector *A = gsl_vector_alloc(d->n);
	gsl_matrix *T = gsl_matrix_alloc(d->n, d->n);
	gsl_matrix *H = gsl_matrix_alloc(d->n, d->p);

	calc_AB(d, x, A, &B);

	/* calculate df / dJ, n x n matrix */
	for (i = 0; i < d->n; i++) {
		for (j = 0; j < d->n; j++) {
			double kd = (i == j) ? 1. : 0.;
			double r = gsl_vector_get(A, i) / B *
					(gsl_vector_get(A, j) / B - kd);
			gsl_matrix_set(T, i, j, r);
		}
	}

	/* calculate dJ / du, n x p matrix */
	gsl_matrix_set_zero(H);
	for (i = 0; i < d->n; i++) {
		complex double r;
		double re;
		for (j = 0; j < d->p; j++) {
			int k = j + 1;
			r = cexp((2. * M_PI * I * i * k) / d->n) * d->tphi[k];
			re = creal(r) / d->n;
			gsl_matrix_set(H, i, j, re);
		}
	}

	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., T, H, 0., J);

	gsl_vector_free(A);
	gsl_matrix_free(T);
	gsl_matrix_free(H);
	return GSL_SUCCESS;
}

static int eq_fdf(const gsl_vector *x, void *params, gsl_vector *f, gsl_matrix *J)
{
	eq_f(x, params, f);
	eq_df(x, params, J);
	return GSL_SUCCESS;
}

static int print_results(struct data *d, gsl_vector *x)
{
	int i;

	fft_out[0] = 0.;
	for (i = 0; i < d->p; i++) {
		int k = i + 1;
		fft_out[k] = gsl_vector_get(x, i);
	}
	for (i = d->p; i < d->n / 2; i++) {
		int k = i + 1;
		fft_out[k] = 0.;
	}
	fftw_execute(p_bwd);

	double u_int[d->n];
	memcpy(u_int, fft_in, d->n * sizeof(*fft_in));

	gsl_vector *A = gsl_vector_alloc(d->n);
	double B;

	calc_AB(d, x, A, &B);

	for (i = 0; i < d->n; i++) {
		double y = i * d->dx;
		double v = u_int[i] / d->n;
		double u = fft_in[i] / d->n;
		double f1 = d->phi[i];
		double f2 = d->phi_av * d->n * gsl_vector_get(A, i) / B;
		printf("%lg %lg %lg %lg %lg\n", y, u, v, f1, f2);
	}

	gsl_vector_free(A);
	return 0;
}

static int print_state(int iter, gsl_multifit_fdfsolver *s)
{
	fprintf(stderr, "iter=%d, |f|=%lg\n", iter, gsl_blas_dnrm2(s->f));
	return 0;	
}

static int find_potential(struct data *d)
{
	const gsl_multifit_fdfsolver_type *T = gsl_multifit_fdfsolver_lmsder;
	gsl_multifit_fdfsolver *s;
	gsl_multifit_function_fdf f;
	gsl_vector *initial = gsl_vector_alloc(d->p);
	int iter = 0, status;

	f.f = eq_f;
	f.df = eq_df;
	f.fdf = eq_fdf;
	f.n = d->n;
	f.p = d->p;
	f.params = d;

	gsl_vector_set_all(initial, 1.);
	s = gsl_multifit_fdfsolver_alloc(T, d->n, d->p);

	gsl_multifit_fdfsolver_set(s, &f, initial);

	do {
		iter++;
		status = gsl_multifit_fdfsolver_iterate(s);
		fprintf(stderr, "status = %s\n", gsl_strerror(status));
		print_state(iter, s);
		if (status) break;
		status = gsl_multifit_test_delta(s->dx, s->x, 1e-2, 1e-2);
	} while (status == GSL_CONTINUE && iter < 1000);
	fprintf(stderr, "status = %s\n", gsl_strerror(status));

	print_results(d, s->x);
	gsl_multifit_fdfsolver_free(s);
	gsl_vector_free(initial);
	return 0;
}

/* ------------------------------------------------------------------------- */

int main(int argc, char **argv)
{
	int ret;
	struct data d;

	memset(&d, 0, sizeof(d));

	feenableexcept(FE_UNDERFLOW | FE_OVERFLOW | FE_INVALID | FE_DIVBYZERO);

	ret = setup(&d, argc, argv);
	if (ret) goto error;

	ret = find_potential(&d);
	if (ret) goto error;

	cleanup(&d);
	fftw_cleanup();

	return EXIT_SUCCESS;
error:
	return EXIT_FAILURE;
}

