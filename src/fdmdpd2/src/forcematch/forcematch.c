/*
 * forcematch.c - force matching for the coms
 * (C) Copyright 2010 by Martin Hoemberg <mhoembe@gwdg.de>
 *
 * This program uses the multiscale coarse-graining (MSCG) method developed by
 * Noid et al, J. Chem. Phys. 128, 244114 (2008),
 * Noid et al, J. Chem. Phys. 128, 244115 (2008)
 *
 * Basis Splines (B-splines) are used as basis functions, since they form
 * a complete set of basis vectors with compact support and they have a
 * continuous derivative. They are already used by
 *
 * Lu and Voth, J. Phys. Chem. B 113, 1501-1510 (2009)
 *
 *
 * Syntax: ./forcematch [config.txt] [profile.txt | - ]
 */

#define DIM 2

#define _GNU_SOURCE
#include "common.h"
#include <omp.h>
#include <mkl.h>
#include <gsl/gsl_bspline.h>

struct data
{
	double *L;		/* boxlengths */
	double *x_cg;		/* positions of cg-sites */
	double *G;		/* huge matrix, in column-major order! */
	double *f;		/* forces from "atomistic" model */
	double *x;		/* estimated model parameters */
	short *leaflet;		/* +1: upper leaflet; -1: lower leaflet */ 
	double chisq;		/* final chi-squared */
	int N;			/* # cg sites */
	int nt;			/* # loaded timesteps */
	int k;			/* degree of the b-splines */
	int ntmax;		/* max. # of timesteps to load at once */
	int nttotal;		/* total amount of processed timesteps */
	int ncoeffs;		/* number of fitting coefficients for each l */
	double rmin;		/* lower cutoff of the potentials */
	double rmax;		/* upper cutoff */
};

/* ------------------------------------------------------------------------- */

static int load_coms(struct data *d, FILE *FH)
{
	int i, j;
	int upper = 0, lower = 0;
	double L[DIM];

	d->nt = 0;
	do {
		char buf[LINE_MAX], *s;
		double t, tmp;
		int N;

		do {
			s = fgets(buf, LINE_MAX, FH);
			if (s == NULL) break;
		} while (strncmp(s, "# L=", 4) != 0);
		if (s == NULL) break;

#if (DIM == 3)
		if (sscanf(buf, "# L=%lg %lg %lg N=%d t=%lg",
				       	&L[0], &L[1], &L[2], &N, &t) != 5)
			return error(EINVAL, "Line '%s' not understood", buf);
#else 
		if (sscanf(buf, "# L=%lg %lg %lg N=%d t=%lg", &tmp,
				       	&L[0], &L[1], &N, &t) != 5)
			return error(EINVAL, "Line '%s' not understood", buf);
#endif
		if (d->N < 0) { /* this is the first timestep at all! */
			d->N = N;

			d->G = calloc((size_t)DIM * d->ntmax * d->N *
				2 * d->ncoeffs, sizeof(*d->G));
			if (d->G == NULL) return novm("data->G");

			d->f = calloc((size_t)DIM * d->ntmax * d->N, sizeof(*d->f));
			if (d->f == NULL) return novm("data->f");

			d->x_cg = calloc((size_t)DIM * N * d->ntmax, sizeof(*d->x_cg));
			if (d->x_cg == NULL) return novm("d->x_cg");
			d->leaflet = calloc(N * d->ntmax, sizeof(*d->leaflet));
			if (d->leaflet == NULL) return novm("d->leaflet");

			d->L = calloc((size_t)DIM * d->ntmax, sizeof(*d->L));
			if (d->L == NULL) return novm("d->L");
		}
		if (d->N != N) 
			return error(EINVAL, "d->N=%d N=%d mismatch", d->N, N);
		for (i = 0; i < DIM; i++)
			d->L[DIM * d->nt + i] = L[i];
#if (DIM == 3)
		debug("load_coms: nt=%d t=%lg L=%lg %lg %lg N=%d",
				d->nt, t, L[0], L[1], L[2], N);
#else
		debug("load_coms: nt=%d t=%lg L=%lg %lg N=%d",
				d->nt, t, L[0], L[1], N);
#endif

		for (i = 0; i < N; i++) {
			double x[DIM], f[DIM], v[DIM], o;
			if (fgets(buf, LINE_MAX, FH) == NULL)
				return error(EINVAL, "premature end of file");
		
#if (DIM == 3)	
			/* x, v, r, f */	
			if (sscanf(buf, "%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg", 
				&x[0], &x[1], &x[2], &v[0], &v[1], &v[2],
				&o, &tmp, &tmp, &f[0], &f[1], &f[2]) != 12)
				return error(EINVAL, "Line '%s' not "
							"understood", buf);
#else
			/* x, v, r, f */	
			if (sscanf(buf, "%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg", 
				&tmp, &x[0], &x[1], &tmp, &v[0], &v[1],
				&o, &tmp, &tmp, &tmp, &f[0], &f[1]) != 12)
				return error(EINVAL, "Line '%s' not "
							"understood", buf);
#endif
			for (j = 0; j < DIM; j++) {
				while (x[j] >= L[j]) x[j] -= L[j];
				while (x[j] < 0.) x[j] += L[j];
				d->x_cg[DIM * (d->nt * N + i) + j] = x[j];
				d->f[DIM * (d->nt * N + i) + j] = f[j];
			}

			/* depending on the monolayer, the normal coordinate
			 * is quenched to +/- 1. */
			d->leaflet[d->nt * N + i] = (o > 0. ) ? 1 : -1;
			if (o > 0.) {
				upper++;
			} else {
				lower++;
			}
		}
		
		d->nt++;
		d->nttotal++;
	} while (d->nt < d->ntmax);

	debug("load_coms: nt=%d timesteps, N=%d cg-sites, upper=%d, lower=%d",
		       d->nt, d->N, upper, lower);
	return 0;
}

/* ------------------------------------------------------------------------- */

static void insert_force(struct data *d, int i, int j, int k,
		gsl_matrix *dB, gsl_bspline_workspace *bw,
	       	gsl_bspline_deriv_workspace *bdw, double dbar[])
{
	size_t rows = DIM * d->nt * d->N;
	double r, delta[DIM];
	int l, m, sameleaflet;

	/* calculate the distance between both lipids */
	{
		double *x_cg = d->x_cg + i * d->N * DIM;

		for (l = 0; l < DIM; l++) {
			double L = d->L[DIM * i + l];
			delta[l] = x_cg[DIM * k + l] - x_cg[DIM * j + l];
			while (delta[l] >= .5 * L) delta[l] -= L;
			while (delta[l] < -.5 * L) delta[l] += L;
		}
		sameleaflet = (d->leaflet[i * d->N + j] == d->leaflet[i * d->N + k]);

#if (DIM == 3)
		if (!sameleaflet) {
			for (l = 0; l < DIM; l++) {
				delta[l] -= dbar[l] * d->leaflet[i * d->N + j];
			}
		}

		r = sqrt(SQR(delta[0]) + SQR(delta[1]) + SQR(delta[2]));
#else
		r = hypot(delta[0], delta[1]);
#endif
	}
	if (r < d->rmin || r > d->rmax) return;
	
	/* iterate over radial coefficients */
	size_t istart, iend;
	gsl_bspline_deriv_eval_nonzero(r, 2, dB, &istart, &iend, bw, bdw);

	for (m = 0; m < d->k; m++) {

		double Fr = -1. * gsl_matrix_get(dB, m, 1);
		size_t col = m + istart + ((sameleaflet) ? 0 : d->ncoeffs);
		size_t row1 = DIM * (i * d->N + j);
		size_t row2 = DIM * (i * d->N + k);

		if (r > 0.) {
			for (l = 0; l < DIM; l++) {
				/* write force matrix in COLUMN-MAJOR order */
				d->G[row1 + col * rows + l] -= Fr * delta[l]/r;
				d->G[row2 + col * rows + l] += Fr * delta[l]/r;
			}
		}
	}
}

#if (DIM == 3)
static void add_dbar(struct data *d, int i, int j, int k, double dbar[],
								int *count)
{
	int l;
	
	if (d->leaflet[i * d->N + j] != d->leaflet[i * d->N + k]) {
		double *x_cg = d->x_cg + i * d->N * DIM;

		for (l = 0; l < DIM; l++) {
			double delta = 0., L = d->L[DIM * i + l];

			/* since leaflet is either +/- 1 it can be used to
			 * get the correct sign of dbar here */
			delta += d->leaflet[i * d->N + j] * x_cg[DIM * j + l];
			delta += d->leaflet[i * d->N + k] * x_cg[DIM * k + l];

			while (delta >= .5 * L) delta -= L;
			while (delta < -.5 * L) delta += L;

			dbar[l] += delta;
		}
		(*count)++;
	}
}
#endif

static int setup_forces(struct data *d)
{
	int i, count = 0;
	double dbar[3] = {0., 0., 0.};
		
	memset(d->G, 0, (size_t)DIM * d->ntmax * d->N * 2 * d->ncoeffs *
								sizeof(*d->G));

	#pragma omp parallel shared(d,dbar,count)
	{
		int nbreak = d->ncoeffs + 2 - d->k;
		gsl_matrix *dB = gsl_matrix_alloc(d->k, 3);
		gsl_bspline_workspace *bw = gsl_bspline_alloc(d->k, nbreak);
		gsl_bspline_deriv_workspace *bdw = gsl_bspline_deriv_alloc(d->k);
		gsl_bspline_knots_uniform(d->rmin, d->rmax, bw);
#if (DIM == 3)
		int count2 = 0;
		double dbar2[3] = {0., 0., 0.};
		
		#pragma omp for schedule(dynamic) nowait
		for (i = 0; i < d->nt; i++) {
			int j, k;

			#pragma omp critical
			{
				debug("calculating dbar for timestep %d on "
					"thread %d.", i, omp_get_thread_num());
			}

			for (j = 0; j < d->N; j++) {
				for (k = j + 1; k < d->N; k++) {
					add_dbar(d, i, j, k, dbar2, &count2);
				}
			}
		}

		#pragma omp atomic
		count += count2;

		for (i = 0; i < 3; i++) {
			#pragma omp atomic
			dbar[i] += dbar2[i];
		}
		
		#pragma omp barrier

		#pragma omp single
		{
			for (i = 0; i < 3; i++) {
				dbar[i] /= count;
			}
		}
#endif
		#pragma omp for schedule(dynamic) nowait
		for (i = 0; i < d->nt; i++) {
			int j, k;

			#pragma omp critical
			{
				debug("finding pairs for timestep %d on thread"
					" %d.", i, omp_get_thread_num());
			}

			for (j = 0; j < d->N; j++) {
				for (k = j + 1; k < d->N; k++) {
					insert_force(d, i, j, k, dB, bw, bdw, dbar);
				}
			}
		}

		gsl_bspline_free(bw);
		gsl_bspline_deriv_free(bdw);
		gsl_matrix_free(dB);
	}
#if (DIM == 3)
	printf("# dbar=%lg %lg %lg\n", dbar[0], dbar[1], dbar[2]);
#endif
	return 0;
}

/* ------------------------------------------------------------------------- */

static int leastsquares(struct data *d)
{
	MKL_INT m = DIM * d->nt * d->N;
	MKL_INT n = 2 * d->ncoeffs;
	MKL_INT p = 2;
	MKL_INT lwork = -1LL, info = 0LL;
	int l;
	double opt = 0, *work = NULL, D[p];
	double x[n];
	double chisq;

	/* setup matrices in column-major (FORTRAN) order */
	double *B = calloc(p * n, sizeof(*B));
	if (B == NULL) return novm("B");

	for (l = 0; l < p; l++) {
		B[l + ((l + 1) * d->ncoeffs - 1) * p] = 1.;
		D[l] = 0.;
	}

	/* dry run to get the optimal work size */
	dgglse(&m, &n, &p, d->G, &m, B, &p, d->f, D, x, &opt, &lwork, &info);
	lwork = opt;
	debug("Optimal size reported by dgglse() is %lld", lwork);
	if (info != 0LL)
		return error(EINVAL, "dgglse() failed with %lld", info);

	work = calloc(lwork, sizeof(*work));
	if (work == NULL) return novm("work");

	/* serious run */
	info = 0LL;
	dgglse(&m, &n, &p, d->G, &m, B, &p, d->f, D, x, work, &lwork, &info);
	if (info != 0LL)
		return error(EINVAL, "dgglse() failed with %lld", info);

	/* calculate chi-squared */
	for (l = n - p, chisq = 0.; l < m; l++)
		chisq += SQR(d->f[l]);
	chisq /= m;
	debug("leastsquares: chisq=%lg nt=%d", chisq, d->nt);
	d->chisq += chisq * d->nt;

	/* test for normality */
	for (l = 0; l < n; l++) {
		if (!isfinite(x[l])) {
			double o = x[l];
			x[l] = 0.;
			fprintf(stderr, "x[%d]: %lg -> %lg.\n", l, o, x[l]);
		}
		d->x[l] += d->nt * x[l];
	}

	free(work);
	free(B);
	return 0;
}

/* ------------------------------------------------------------------------- */

static int print_results(struct data *d)
{
	int i, l;
	double r;

	d->chisq /= d->nttotal;
	printf("# chisq = %e\n# r/|z|\tf\tdf/dr\tg\tdg/dz\n", d->chisq);
	
	int nbreak = d->ncoeffs + 2 - d->k;
	gsl_matrix *dB = gsl_matrix_alloc(d->ncoeffs, 3);
	gsl_bspline_workspace *bw = gsl_bspline_alloc(d->k, nbreak);
	gsl_bspline_deriv_workspace *bdw = gsl_bspline_deriv_alloc(d->k);
	gsl_bspline_knots_uniform(d->rmin, d->rmax, bw);

	for (r = d->rmin; r < d->rmax; r += .01) {
		printf("%lg", r);
		gsl_bspline_deriv_eval(r, 2, dB, bw, bdw);
		gsl_vector_const_view vu = gsl_matrix_const_column(dB, 0);
		gsl_vector_const_view vf = gsl_matrix_const_column(dB, 1);
		for (l = 0; l < 2; l++) {
			double u = 0., du = 0.;
			for (i = 0; i < d->ncoeffs; i++) {
				u += d->x[l * d->ncoeffs + i] / d->nttotal * gsl_vector_get(&vu.vector, i);
				du += d->x[l * d->ncoeffs + i] / d->nttotal * gsl_vector_get(&vf.vector, i);
			}
			printf("\t%lg\t%lg", u, du);
		}
		printf("\n");
	}

	printf("## # k=%d ncoeffs=%d rc=%lg dim=%d\n", d->k, d->ncoeffs, d->rmax, DIM);
	for (i = 0; i < d->ncoeffs; i++) {
		printf("## %d\t%lg\t%lg\n", i, d->x[i] / d->nttotal, d->x[i + d->ncoeffs] / d->nttotal);
	}

	gsl_matrix_free(dB);
	gsl_bspline_free(bw);
	gsl_bspline_deriv_free(bdw);

	return 0;
}

/* ------------------------------------------------------------------------- */

static int load_config(struct data *d, const char *fn)
{
	FILE *FH = fopen(fn, "r");
	
	memset(d, 0, sizeof(*d));

	if (FH == NULL)
		return error(ENOENT, "Couldn't open '%s' for reading", fn);

	/* 1st line: order of b-spline functions, typically 4 or 6 */
	if (fscanf(FH, "%d ", &d->k) != 1) return error(EINVAL, "k");
	/* 2rd line: number of fitting coefficients for each l */
	if (fscanf(FH, "%d ", &d->ncoeffs) != 1) return error(EINVAL, "ncoeffs");
	/* 3th line: rmin, rmax / cutoffs of the interactions */
	if (fscanf(FH, "%lg %lg ", &d->rmin, &d->rmax) != 2)
				       return error(EINVAL, "rmin/rmax");
	/* 4th line: maximum number of timesteps for block averaging */
	if (fscanf(FH, "%d ", &d->ntmax) != 1) return error(EINVAL, "ntmax");
		
	debug("k=%d, ncoeffs=%d, rmin=%lg, rmax=%lg, ntmax=%d",
		d->k, d->ncoeffs, d->rmin, d->rmax, d->ntmax);

	if (d->k < 3 || d->k > 10)
		return error(EDOM, "Invalid k=%d specified", d->k);
	if (d->ncoeffs < 5 || d->ncoeffs > 100000)
		return error(EDOM, "Invalid ncoeffs=%d specified", d->ncoeffs);
	if (d->rmin < 0. || d->rmin >= d->rmax)
		return error(EDOM, "Invalid rmin=%lg specified", d->rmin);
	if (d->rmax < d->rmin || d->rmax >= 100.)
		return error(EDOM, "Invalid rmax=%lg specified", d->rmax);
	
	d->x = calloc(2 * d->ncoeffs, sizeof(*d->x));
	if (d->x == NULL) return novm("d->x");
	d->N = -1;

	fclose(FH);
	return 0;
}

/* ------------------------------------------------------------------------- */

int main(int argc, char **argv)
{
	int ret;
	struct data d;
	FILE *FH;
	const char *fn;
	
	if (argc < 3) {
		fprintf(stderr, "Syntax: %s [ config.txt ] [ coms.dat | - ]\n",
								argv[0]);
		goto error;
	}

	fprintf(stderr, "Using OpenMP with %d threads on %d CPUs!\n",
				omp_get_max_threads(), omp_get_num_procs());

	ret = load_config(&d, *(++argv));
	if (ret) return ret;

	fn = *(++argv);
	FH = (strcmp(fn, "-") == 0) ? stdin : fopen(fn, "r");
	if (FH == NULL)
		return error(ENOENT, "Couldn't open '%s' for reading", fn);

	do {
		ret = load_coms(&d, FH);
		if (ret) return ret;
	
		if (d.nt > 0) {
			ret = setup_forces(&d);
			if (ret) return ret;

			ret = leastsquares(&d);
			if (ret) goto error;
		}
	} while (d.nt == d.ntmax);

	ret = print_results(&d);
	if (ret) return ret;

	/* cleanup */
	if (FH != stdin) fclose(FH);
	free(d.G);
	free(d.f);
	free(d.x_cg);
	free(d.x);

	fprintf(stderr, "Success.\n\n");
	return EXIT_SUCCESS;
error:
	fprintf(stderr, "Failed.\n\n");
	return EXIT_FAILURE;
}

