/*
 * dynamics.c - In-situ calculation of various dynamical quantities
 * Copyright (C) 2008-2010 Martin Hoemberg <mhoembe@gwdg.de>
 *
 * This file is part of FDMDPD2.
 *
 * FDMDPD2 is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * FDMDPD2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with FDMDPD2.  If not, see <http://www.gnu.org/licenses/>.
 */

/* $Id: dynamics.c 317 2011-08-01 07:57:35Z hoemberg $ */

#include "fdmdpd2.h"
#include <fftw3.h>

/*
 * NOTE: There are two different mapping schemes of points to the grid.
 * Currently it is only possible to select one of them at compile-time.
 */
#define BSPLINES

enum mode
{
	LOAD,		/* load data */
	SAVE,		/* save data to disk */
};

/* --------------------------- common variables ---------------------------- */
static int inuse = 0;			/* is the whole thing enabled? */
static long comdyn = 0;			/* use com coordinates? */
static int blocks;			/* max. number of blocks */
static int elements;			/* max. number of elements per block */
static int oldxv_enabled = 0;		/* do we need the old coordinates? */
static VEC2 *oldxv = NULL;		/* stored bead positions */
static int *blocklen = NULL;		/* lengths of each block */
static int *c = NULL;			/* counter for each block / element */
static int mcount = 0;			/* # how often measure() was called */
static int cur_blocks = 1;		/* current number of blocks */
static double delta_a = 0.;		/* time between two calls */
static char *iobuf = NULL;		/* I/O buffer */
static size_t iobufsize;		/* size for I/O buffer */
static struct beads *bcom = NULL;	/* beads structure for coms */
static int nbl_enabled = 0;		/* do we need a neighbor list? */
static VEC2 *nbl_pos = NULL;		/* tmp. buffer for old positions */
static long dump_current = 0;		/* dump current data everytime? */
static char *cur_iobuf = NULL;		/* I/O buffer for current dumps */
static FILE *dumpFH = NULL;		/* file handle for current dumps */

/* -------------------------------- FFT variables -------------------------- */
static int fft_enabled = 0;		/* enable fft routines? */
static int fft_K[2];			/* number of bins in x direction */
static long fft_n = 4;			/* degree of b-spline: 2, 4, 6 */
static int fft_rsize;			/* size in real space: K[0] * K[1] */
static int fft_qsize;			/* in q-space: K[0] * (K[1] / 2 + 1) */
static fftw_plan fft_plan;		/* FFTW plan */
static complex double *fft_tmp;		/* temporary FFT buffer */
static double *fft_Q;			/* real space (input buffer) */
static complex double *fft_b;		/* prefactors for b-spline expansion */
static complex double *fft_S;		/* k-space (output buffer) */
static double iboxlen[2];		/* inv. tang. box lengths */

/* --------------------- intermediate structure function ------------------- */
static int isf_enabled = 0;		/* enable measurement of isf? */
static complex double *isf_fftrho;	/* stored FFTs of the density */
static complex double *isf_Fs;		/* dynamic self structure factor */
static complex double *isf_Fc;		/* dynamic coll. structure factor */
static double *isf_Fc0;			/* static coll. structure factor */
static double *isf_r2 = NULL;		/* beads' mean-square displacement */
static double *isf_r4 = NULL;		/* for non-gaussianity parameter */

/* ------------------ fourpoint (fp) correlation functions ----------------- */
static int fp_enabled = 0;		/* enable measurement of fp? */
static double fp_a = 1.0;		/* cutoff parameter "a" */
static struct nblist fp_nbl;		/* neighbor list */
static double *fp_Q = NULL;		/* <Q(t)> and higher moements */
static double *fp_S4 = NULL;		/* S_4^ol: SS/SD/DD after another */

/* --------------- velocity autocorrelation functions (vacf)---------------- */
static int vacf_enabled = 0;		/* enable measurement of vacf? */
static VEC *vacf_c = NULL;		/* correlations in three dimensions */
#ifdef MCOM
/* ----------------- force autocorrelation functions (facf)----------------- */
static int facf_enabled = 0;		/* enable measurement of facf? */
static double *facf_g = NULL;		/* g(t,r), see below for definition */
static double *facf_h = NULL;		/* h(t,r), see below for definition */
static int *facf_c = NULL;		/* counters for g(t,r) and h(t,r) */
static double facf_dr = 0.1;		/* radial binwidth */
static int facf_rbins;			/* radial bins */
static VEC *facf_oldf = NULL;		/* memory of all forces */
static struct nblist facf_nbl;		/* neighbor list */
#endif
/* ---------- autocorrelation function of undulation modes (uacf) ---------- */
static int uacf_enabled = 0;		/* enable measurement of uacf? */
static complex double *uacf_fft0;	/* current FFT buffer */
static complex double *uacf_fft;	/* stored FFTs of the height */
static complex double *uacf_S;		/* dynamic structure factor */
static double *uacf_S0;			/* static structure factor */
static double *uacf_Qnorm;		/* normalization buffer */

/* ---------- autocorrelation function of particle current (jacf) ---------- */
/* the correlation function is a 3x3 rank-two tensor 
 * the _s_ quantities refer to the sum over the two leaflets, whereas the _d_
 * quantities are the difference between the two leaflets (upper - lower)/2 */
static int jacf_enabled = 0;		/* enable measurement of jacf ? */
static complex double *jacf_s_fft[2];	/* stored FFTs of the velocity field */
static complex double *jacf_d_fft[2];
static complex double *jacf_s_fft0[2];	/* current FFT of the velocity field */
static complex double *jacf_d_fft0[2];
static complex double *jacf_s_C[4];	/* dynamic correlation function */
static complex double *jacf_d_C[4];
static double *jacf_s_C0[4];		/* static correlation function */
static double *jacf_d_C0[4];

/* ---------------------- global (extern) declarations --------------------- */
cfg_opt_t dyn_opts[] = {
	CFG_SIMPLE_INT("n", &fft_n),
	CFG_SIMPLE_INT("com", &comdyn),
	CFG_SIMPLE_INT("dump_current", &dump_current),
	CFG_END()
};

/* ------------------------------------------------------------------------- */

#ifdef BSPLINES
/*
 * METHOD #1: Cardinal B-Spline interpolation
 *
 * This method is taken from the SPME algorithm used in electrostatics. [1]
 * High values of K, i.e. 64, 128, 256 are needed. The order of the b-splines
 * can be selected in the config-file (recommended: n=4). However, numerical
 * artifacts are expected at the length scale of the b-splines, as well as at
 * the end of the first Brillouin zone. Another issue is the occurence
 * of anisotropy effects at the length scale of the b-splines.
 *
 * [1] U Essmann et al., JCP 103, 8577 (1995)
 */ 

static double (*fft_Mn)(const double u);/* cardinal b-spline function */

static __attribute__((const)) double bspline_M2(const double u)
{
	double x = fabs(u - 1.);
	return (x < 1.) ? 1. - x : 0.;
}

static __attribute__((const)) double bspline_M4(const double u)
{
	double x = fabs(u - 2.);
	if (x < 1.) {
		return (.5 * x - 1.) * x * x + 2./3.;
	} else if (x < 2.) {
		return ((-1./6. * x + 1.) * x - 2.) * x + 4./3.;
	} else {
		return 0.;
	}
}

static __attribute__((const)) double bspline_M6(const double u)
{
	double x = fabs(u - 3.);
	if (x < 1.) {
		return ((-1./12. * x + .25) * x * x -.5) * x * x + 11./20.;
	} else if (x < 2.) {
		return ((((1./24. * x - 3./8.) * x + 5./4.) * x - 7./4.) * x +
							5./8.) * x + 17./40.;
	} else if (x < 3.) {
		return ((((-1./120. * x + 1./8.) * x - 3./4.) * x + 9./4.) * x
						- 27./8.) * x + 81./40.;
	} else {
		return 0.;
	}
}

/*
 * calculates the complex weighting factors from the Cardinal b-Splines for
 * one component of one k-vector. This is eq. (4.4) in the article.
 */
static complex double fft_field_bm(const int m, const int K)
{
	complex double num, denom;
	int k;

	if (((fft_n & 1) == 1) && (abs(2 * m) == K))
		return 0.;

	num = cexp(-2. * M_PI * I * (fft_n - 1.) * m / K);
		
	for (k = 0, denom = 0.; k <= fft_n - 2; k++)
		denom += fft_Mn(k + 1) *
			cexp(-2. * M_PI * I * m * k / K);
	return num / denom;
}

/* calculation of the constant coefficient matrix */
static void fft_field_calc_b(void)
{
	complex double b1, b2;
	int i, j; /* array indeces */
	int m1, m2; /* Miller's indices */
	
	for (i = 0; i < fft_K[0]; i++) {
		m1 = (i <= fft_K[0] / 2) ? i : (i - fft_K[0]);
		b1 = fft_field_bm(m1, fft_K[0]);
		for (j = 0; j < fft_K[1] / 2 + 1; j++) {
			m2 = j;
			b2 = fft_field_bm(m2, fft_K[1]);
			fft_b[i * (fft_K[1] / 2 + 1) + j] = b1 * b2;
		}
	}
}

/* add a vector to the input field Q of the FFT */
static void add2Q(const double x[static 2], const double q, double Q[])
{
	double u, M[2][fft_n];
	int i, j, fu[2];

	for (i = 0; i < 2; i++) {
		u = x[i] * iboxlen[i];
		u -= (u >= 1.) ? 1. : 0.;
		u += (u < 0.) ? 1. : 0.;
		assert(u >= 0. && u < 1.);
		u *= fft_K[i];
		fu[i] = (int)floor(u);

		for (j = 0; j < fft_n; j++) {
			M[i][j] = fft_Mn(u - (fu[i] - j));
		}
	}

	for (i = 0; i < fft_n; i++) {
		int k0 = fu[0] - i;
		k0 += (k0 < 0) ? fft_K[0] : 0;
		for (j = 0; j < fft_n; j++) {
			int k1 = fu[1] - j;
			k1 += (k1 < 0) ? fft_K[1] : 0;
			Q[k0 * fft_K[1] + k1] += q * M[0][i] * M[1][j];
		}
	}
}

/* ------------------------------------------------------------------------- */

#else
/* 
 * METHOD #2: Standard, next gridpoint interpolation
 *
 * This method is easier and perhaps faster. Typical discretizations are
 * smaller than before, i.e. K=16,32 have proven to work correctly. The
 * structure factor (sinc-expressions) of this mapping function is discussed
 * in the Appendix of
 *
 * I Cooke and M Deserno, JCP 123, 224710 (2005)
 */

static void add2Q(const double x[static 2], const double q, double Q[])
{
	int i, fu[2];

	for (i = 0; i < 2; i++) {
		fu[i] = (int)floor(x[i] * iboxlen[i] * fft_K[i] + .5);
		fu[i] %= fft_K[i];
		assert(fu[i] >= 0 && fu[i] < fft_K[i]);
	}
	Q[fu[0] * fft_K[1] + fu[1]] += q;
}

static __attribute__((const)) double sinc(const double x)
{
	return (x == 0.) ? 1. : sin(M_PI * x) / (M_PI * x);
}

static void fft_field_calc_b(void)
{
	int i, j;
	double m0, m1, b0, b1;

	for (i = 0; i < fft_K[0]; i++) {
		m0 = (i <= fft_K[0] / 2) ? i : (i - fft_K[0]);
		b0 = sinc(m0 / fft_K[0]);
		assert(isfinite(b0));
		for (j = 0; j < fft_K[1] / 2 + 1; j++) {
			m1 = j;
			b1 = sinc(m1 / fft_K[1]);
			assert(isfinite(b1));
			fft_b[i * (fft_K[1] / 2 + 1) + j] = 1. / (b0 * b1);
		}
	}
}
#endif

/* ------------------------------------------------------------------------- */

static void fft_add1(const double x[static 2])
{
	add2Q(x, 1., fft_Q);
}

static void fft_addQ(const double x[static 2], const double q)
{
	add2Q(x, q, fft_Q);
}

static void fft_alloc(void)
{
	/* field to calculate FFT */
	fft_b = calloc(fft_qsize, sizeof(*fft_b));
	if (fft_b == NULL) novm("fft_b");
	fft_Q = fftw_malloc(fft_rsize * sizeof(*fft_Q));
	if (fft_Q == NULL) novm("fft_Q");
	fft_S = fftw_malloc(fft_qsize * sizeof(*fft_S));
	if (fft_S == NULL) novm("fft_S");

	fft_plan = fftw_plan_dft_r2c_2d(fft_K[0], fft_K[1], fft_Q, fft_S,
								FFTW_PATIENT);
	fft_tmp = calloc(fft_qsize, sizeof(*fft_tmp));
	if (fft_tmp == NULL) novm("fft_tmp");
}

static void fft_free(void)
{
	if (!fft_enabled) return;

	free(fft_tmp);
	fftw_destroy_plan(fft_plan);
	fftw_free(fft_S);
	fftw_free(fft_Q);
	free(fft_b);
}

static void fft_loadsave_header(FILE *FH, struct beads *restrict b,
							const enum mode mode)
{
	if (!fft_enabled) return;

	if (mode == LOAD) {
		if (fscanf(FH, "# K=%d %d\n", &fft_K[0], &fft_K[1]) != 2)
			fatal(EINVAL, "fft_loadsave_header");
		
		debug("isf: K=%d %d n=%ld", fft_K[0], fft_K[1], fft_n);
		
		fft_rsize = fft_K[0] * fft_K[1];
		fft_qsize = fft_K[0] * (fft_K[1] / 2 + 1);

#ifdef BSPLINES
		if (ismaster)
			printf("Using B-Splines interpolation.\n");

		if (fft_n == 2) {
			fft_Mn = bspline_M2;
		} else if (fft_n == 4) {
			fft_Mn = bspline_M4;
		} else if (fft_n == 6) {
			fft_Mn = bspline_M6;
		} else {
			fatal(EINVAL, "Mn=%ld is impossible", fft_n);
		}
#else
		if (ismaster)
			printf("Using Standard interpolation.\n");
#endif	
		fft_alloc();
		fft_field_calc_b();
	} else {
		if (ismaster)
			fprintf(FH, "# K=%d %d\n", fft_K[0], fft_K[1]);
	}
}

/* ------------------------------------------------------------------------- */

static void oldxv_loadsave(FILE *FH, struct beads *restrict b,
					const int offset, const enum mode mode)
{
	if (!oldxv_enabled) return;

	PASSPORT *pp = b->passport + col_index * b->groupsize;
	VEC2 *buf, *xv = oldxv + offset * b->groupsize;
	int i, d;
	
	buf = calloc(b->nN, sizeof(*buf));
	if (buf == NULL) novm("buf");
		
	switch (mode) {
	case LOAD:
		for (i = 0; i < b->nN; i++) {
			if (fscanf(FH, "%lf %lf %lf %lf %lf %lf\n",
				&buf[i][0], &buf[i][1], &buf[i][2],
				&buf[i][3], &buf[i][4], &buf[i][5]) != 6)
				fatal(EINVAL, "Error reading pos from "
								"line %d", i);
		}

		for (i = 0; i < b->groupsize; i++) {
			int id = GET_ID(pp[i]);
			if (!GET_EXISTS(pp[i])) continue;
			
			for (d = 0; d < (int)ARRAY_SIZE(*buf); d++)
				xv[i][d] = buf[id][d];
		}
		
		break;
	case SAVE:
		for (i = 0; i < b->groupsize; i++) {
			int id = GET_ID(pp[i]);
			if (!GET_EXISTS(pp[i])) continue;

			for (d = 0; d < (int)ARRAY_SIZE(*buf); d++)
				buf[id][d] = xv[i][d];
		}
		if (FH) {
			MPI_Reduce(MPI_IN_PLACE, buf, ARRAY_SIZE(*buf) * b->nN,
					MPI_DOUBLE, MPI_SUM, 0, comm_grid);
			for (i = 0; i < b->nN; i++) {
				fprintf(FH, "%.12lg %.12lg %.12lg %.12lg "
					"%.12lg %.12lg\n",
				       	buf[i][0], buf[i][1], buf[i][2],
					buf[i][3], buf[i][4], buf[i][5]);
			}
		} else {
			MPI_Reduce(buf, NULL, ARRAY_SIZE(*buf) * b->nN,
					MPI_DOUBLE, MPI_SUM, 0, comm_grid);
		}

		break;
	}

	free(buf);
}

static void oldxv_insert(struct beads *restrict b)
{
	if (!oldxv_enabled) return;

	int i, d;

	for (i = 0; i < b->groupsize; i++) {
		int j = b->groupsize * col_index + i;
		for (d = 0; d < 3; d++) {
			oldxv[i][d] = b->xv[j][d] + b->nx[j][d] * b->l[d];
			oldxv[i][d + 3] = b->xv[j][d + 3];
		}
	}
}

/* FIXME: in the last block the program will shift out of memory! */
static void oldxv_shift(struct beads *restrict b, const int block)
{
	if (!oldxv_enabled) return;

	VEC2 *bl = oldxv + block * elements * b->groupsize;
	memmove(bl + b->groupsize, bl, elements * b->groupsize * sizeof(*bl));
}


/* ------------------------------------------------------------------------- */

static void isf_alloc(void)
{
	size_t be = blocks * elements;

	/* block buffers */
	isf_fftrho = calloc(be * fft_qsize, sizeof(*isf_fftrho));
	if (isf_fftrho == NULL) novm("isf_fftrho");

	/* mean-square displacement, non-Gaussianity */
	isf_r2 = calloc(be, sizeof(*isf_r2));
	if (isf_r2 == NULL) novm("isf_r2");
	isf_r4 = calloc(be, sizeof(*isf_r4));
	if (isf_r4 == NULL) novm("isf_r4");
	
	/* dynamic structure factor result storage */
	isf_Fs = calloc(be * fft_qsize, sizeof(*isf_Fs));
	if (isf_Fs == NULL) novm("isf_Fs");
	isf_Fc = calloc(be * fft_qsize, sizeof(*isf_Fc));
	if (isf_Fc == NULL) novm("isf_Fc");

	/* static structure factor result storage */
	isf_Fc0 = calloc(fft_qsize, sizeof(*isf_Fc0));
	if (isf_Fc0 == NULL) novm("isf_Fc0");
}

static void isf_free(void)
{
	if (!isf_enabled) return;

	free(isf_fftrho);
	free(isf_Fs);
	free(isf_Fc);
	free(isf_Fc0);
	free(isf_r2);
	free(isf_r4);
}

static void isf_loadsave_header(const enum mode mode)
{
	if (!isf_enabled) return;
	if (mode == LOAD) isf_alloc();
}

static void isf_loadsave_r2r4(FILE *FH, const int offset, const enum mode mode)
{
	if (!isf_enabled) return;

	double r2, r4;

	switch (mode) {
	case LOAD:
		if (fscanf(FH, "%lf %lf\n", &r2, &r4) != 2)
			fatal(EIO, "Reading r2r4 offset %d", offset);
		if (ismaster) {
			isf_r2[offset] = r2;
			isf_r4[offset] = r4;
		}
		break;
	case SAVE:
		MPI_Reduce(&isf_r2[offset], &r2, 1, MPI_DOUBLE, MPI_SUM, 0,
								comm_grid);
		MPI_Reduce(&isf_r4[offset], &r4, 1, MPI_DOUBLE, MPI_SUM, 0,
								comm_grid);
		if (ismaster) fprintf(FH, "%.12lg %.12lg\n", r2, r4);
		break;
	}
}

static void isf_loadsave_fc0(FILE *FH, const enum mode mode)
{
	if (!isf_enabled) return;
	
	int i;
	double *fc0;

	fc0 = calloc(fft_qsize, sizeof(*fc0));
	if (fc0 == NULL) novm("fc0");

	switch (mode) {
	case LOAD:
		for (i = 0; i < fft_qsize; i++) {
			if (fscanf(FH, "%lf\n", &fc0[i]) != 1)
				fatal(EIO, "reading fc0 line %d", i);
		}
		if (ismaster)
			memcpy(isf_Fc0, fc0, fft_qsize * sizeof(*fc0));
		break;
	case SAVE:
//		MPI_Reduce(isf_Fc0, fc0, fft_qsize, MPI_DOUBLE, MPI_SUM,
//		       					0, comm_grid);
		if (ismaster) {
			for (i = 0; i < fft_qsize; i++)
				fprintf(FH, "%.12lg\n", isf_Fc0[i]);
		}
		break;
	}

	free(fc0);
}

static void isf_loadsave_fft(FILE *FH, const int offset, const enum mode mode)
{
	if (!isf_enabled) return;

	complex double *fftrho, *Fs, *Fc;
	int i;

	switch (mode) {
	case LOAD:
		fftrho = isf_fftrho + offset * fft_qsize;
		Fs = isf_Fs + offset * fft_qsize;
		Fc = isf_Fc + offset * fft_qsize;
		
		for (i = 0; i < fft_qsize; i++) {
			double re1, re2, re3, im1, im2, im3;
			if (fscanf(FH, "%lf %lf %lf %lf %lf %lf\n",
				&re1, &im1, &re2, &im2, &re3, &im3) != 6)
				fatal(EINVAL, "Error while reading line"
								" %d", i);
			fftrho[i] = re1 + im1 * I;
			if (ismaster) {
				Fs[i] = re2 + im2 * I;
				Fc[i] = re3 + im3 * I;
			}
		}
		break;
	case SAVE:
		fftrho = isf_fftrho + offset * fft_qsize;
		Fs = calloc(fft_qsize, sizeof(*Fs));
		if (Fs == NULL) novm("Fs");
		Fc = isf_Fc + offset * fft_qsize;
		
		MPI_Reduce(isf_Fs + offset * fft_qsize, Fs, 2 * fft_qsize,
					MPI_DOUBLE, MPI_SUM, 0, comm_grid);
		if (ismaster) {
			for (i = 0; i < fft_qsize; i++)
				fprintf(FH, "%.12lg %.12lg %.12lg %.12lg "
							"%.12lg %.12lg\n",
					creal(fftrho[i]), cimag(fftrho[i]),
					creal(Fs[i]), cimag(Fs[i]),
					creal(Fc[i]), cimag(Fc[i]));
		}

		free(Fs);
		break;
	}
}

/* ------------------------------------------------------------------------- */

/* measure the collective intermediate scattering function */
static void isfc_fft_density(struct beads *restrict b)
{
	if (!isf_enabled) return;
	
	VEC2 *xv = b->xv + b->groupsize * col_index;
	PASSPORT *pp = b->passport + col_index * b->groupsize;
	int i;
	double x[2];
	
	memset(fft_Q, 0, fft_rsize * sizeof(*fft_Q));

	for (i = 0; i < b->groupsize; i++) {
		if (!GET_EXISTS(pp[i])) continue;
		x[0] = xv[i][TANG1];
		x[1] = xv[i][TANG2];
		fft_add1(x);
	}

	fftw_execute(fft_plan);

	for (i = 0; i < fft_qsize; i++)
		fft_tmp[i] = fft_b[i] * fft_S[i];

	MPI_Allreduce(MPI_IN_PLACE, fft_tmp, 2 * fft_qsize,
					MPI_DOUBLE, MPI_SUM, comm_grid);
}

static void isfc_dump_current(FILE *FH)
{
	if (!isf_enabled) return;
	int i;
			
	for (i = 0; i < fft_qsize; i++)
		fprintf(FH, "%.12lg %.12lg\n", creal(fft_tmp[i]),
							cimag(fft_tmp[i]));
}

static void isfc_measure(const int offset)
{
	if (!isf_enabled) return;

	int i;
	complex double *fc = isf_Fc + offset * fft_qsize;
	complex double *old_rhofft = isf_fftrho + offset * fft_qsize;
	
	for (i = 0; i < fft_qsize; i++)
		fc[i] += fft_tmp[i] * old_rhofft[i];
}

static void isfc_insert(void)
{
	if (!isf_enabled) return;

	int i;

	for (i = 0; i < fft_qsize; i++)
		isf_fftrho[i] = conj(fft_tmp[i]);
}

/* FIXME: in the last block the program will shift out of memory! */
static void isfc_shift(const int block)
{
	if (!isf_enabled) return;

	complex double *bl = isf_fftrho + block * elements * fft_qsize;
	memmove(bl + fft_qsize, bl, elements * fft_qsize * sizeof(*bl));
}

static void isfc_static(void)
{
	if (!isf_enabled) return;

	int i;

	for (i = 0; i < fft_qsize; i++) {
		complex double S = fft_tmp[i];
		isf_Fc0[i] += creal(S * conj(S));
	}
}

/* self intermediate scattering function and mean-square displacement */
static void isfs_measure(struct beads *restrict b, const int offset)
{
	if (!isf_enabled) return;

	complex double *fs = isf_Fs + offset * fft_qsize;
	VEC2 *pos = b->xv + b->groupsize * col_index;
	VEC2 *old_pos = oldxv + offset * b->groupsize;
	PASSPORT *pp = b->passport + col_index * b->groupsize;
	VEC *nx = b->nx + b->groupsize * col_index;
	int i, j, d;
	
	memset(fft_Q, 0, fft_rsize * sizeof(*fft_Q));

	for (i = 0; i < b->groupsize; i++) {
		double tmp[2];
		double r2 = 0.;
		
		if (!GET_EXISTS(pp[i])) continue;
		//if (GET_TYPE(pp[i]) == 0) continue; /* only A */

		for (j = 0; j < 2; j++) {
			d = (j == 0) ? TANG1 : TANG2;
			tmp[j] = pos[i][d] + nx[i][d] * b->l[d] - old_pos[i][d];
			r2 += SQR(tmp[j]);
			while (tmp[j] >= b->l[d]) tmp[j] -= b->l[d];
			while (tmp[j] < 0.) tmp[j] += b->l[d];
		}
		fft_add1(tmp);

		isf_r2[offset] += r2;
		isf_r4[offset] += SQR(r2);
	}

	fftw_execute(fft_plan);

	for (i = 0; i < fft_qsize; i++)
		fs[i] += fft_b[i] * fft_S[i];
}

/* ------------------------------------------------------------------------- */

/*
 * autocorrelation function of the undulation power spectrum (UACF).
 * 
 * The static part is power spectrum of the undulations which readily yields
 * the bending rigidity of the bilayer. From the decay of the UACF one
 * can extract the inter-monolayer friction coefficient "b" by means of the
 * Seifert-Langer theory.
 *
 * References:
 *   SA Shkulipa, WK den Otter, WJ Briels, JCP 125, 234905 (2006)
 *   U Seifert and SA Langer, EPL 23, 71 (1993)
 *   U Seifert and SA Langer, Biophys.Chem. 49, 13 (1994)
 *
 */

static void uacf_alloc(void)
{
	size_t be = blocks * elements;

	/* block buffers */
	uacf_fft = calloc(be * fft_qsize, sizeof(*uacf_fft));
	if (uacf_fft == NULL) novm("uacf_fft");
	
	/* current FFT buffer */
	uacf_fft0 = calloc(fft_qsize, sizeof(*uacf_fft0));
	if (uacf_fft0 == NULL) novm("uacf_fft0");

	/* dynamic structure factor result storage */
	uacf_S = calloc(be * fft_qsize, sizeof(*uacf_S));
	if (uacf_S == NULL) novm("uacf_S");

	/* static structure factor result storage */
	uacf_S0 = calloc(fft_qsize, sizeof(*uacf_S0));
	if (uacf_S0 == NULL) novm("uacf_S0");

	/* temporary buffer for normalization purpose */
	uacf_Qnorm = calloc(fft_rsize, sizeof(*uacf_Qnorm));
	if (uacf_Qnorm == NULL) novm("uacf_Qnorm");
}

static void uacf_free(void)
{
	if (!uacf_enabled) return;

	free(uacf_fft);
	free(uacf_fft0);
	free(uacf_S);
	free(uacf_S0);
	free(uacf_Qnorm);
}

static void uacf_loadsave_header(const enum mode mode)
{
	if (!uacf_enabled) return;
	if (mode == LOAD) uacf_alloc();
}

static void uacf_loadsave_S0(FILE *FH, const enum mode mode)
{
	if (!uacf_enabled) return;
	
	int i;
	double *fc0;

	fc0 = calloc(fft_qsize, sizeof(*fc0));
	if (fc0 == NULL) novm("fc0");

	switch (mode) {
	case LOAD:
		for (i = 0; i < fft_qsize; i++) {
			if (fscanf(FH, "%lf\n", &fc0[i]) != 1)
				fatal(EIO, "reading uacf-S0 line %d",i);
		}
		if (ismaster)
			memcpy(uacf_S0, fc0, fft_qsize * sizeof(*fc0));
		break;
	case SAVE:
		if (ismaster) {
			for (i = 0; i < fft_qsize; i++)
				fprintf(FH, "%.12lg\n", uacf_S0[i]);
		}
		break;
	}

	free(fc0);
}

static void uacf_loadsave_fft(FILE *FH, const int offset, const enum mode mode)
{
	if (!uacf_enabled) return;

	complex double *fft = uacf_fft + offset * fft_qsize;
	complex double *S = uacf_S + offset * fft_qsize;
	int i;

	switch (mode) {
	case LOAD:
		for (i = 0; i < fft_qsize; i++) {
			double re1, re2, im1, im2;
			if (fscanf(FH, "%lf %lf %lf %lf\n",
						&re1, &im1, &re2, &im2) != 4)
				fatal(EINVAL, "Error while reading line"
								" %d", i);
			fft[i] = re1 + im1 * I;
			if (ismaster) {
				S[i] = re2 + im2 * I;
			}
		}
		break;
	case SAVE:
		if (ismaster) {
			for (i = 0; i < fft_qsize; i++)
				fprintf(FH, "%.12lg %.12lg %.12lg %.12lg\n",
					creal(fft[i]), cimag(fft[i]),
					creal(S[i]), cimag(S[i]));
		}

		break;
	}
}

static void uacf_dump_current(FILE *FH)
{
	if (!uacf_enabled) return;
	int i;
			
	for (i = 0; i < fft_qsize; i++)
		fprintf(FH, "%lg %lg\n", creal(uacf_fft0[i]),
							cimag(uacf_fft0[i]));
}

static void uacf_calc_fft(struct beads *restrict b)
{
	if (!uacf_enabled) return;
	
	VEC2 *xv = b->xv + b->groupsize * col_index;
	PASSPORT *pp = b->passport + col_index * b->groupsize;
	int i;
	double mean, mcnt;
	
	memset(fft_Q, 0, fft_rsize * sizeof(*fft_Q));
	memset(uacf_Qnorm, 0, fft_rsize * sizeof(*uacf_Qnorm));

	for (i = 0; i < b->groupsize; i++) {
		double x[2] = {xv[i][TANG1], xv[i][TANG2]};

		if (!GET_EXISTS(pp[i])) continue;
		if (GET_TYPE(pp[i]) == 1) continue; /* only A */

		fft_addQ(x, xv[i][NORMAL]);
		add2Q(x, 1., uacf_Qnorm);
	}
	
	MPI_Allreduce(MPI_IN_PLACE, fft_Q, fft_rsize,
					MPI_DOUBLE, MPI_SUM, comm_grid);
	MPI_Allreduce(MPI_IN_PLACE, uacf_Qnorm, fft_rsize,
					MPI_DOUBLE, MPI_SUM, comm_grid);
	for (i = 0, mean = 0.; i < fft_rsize; i++) {
		double a = uacf_Qnorm[i];
		fft_Q[i] = (a > 0.) ? fft_Q[i] / a : 0.;
		mcnt += (a > 0.) ? 1. : 0.;
		mean += fft_Q[i];
	}
	mean /= mcnt;
	for (i = 0; i < fft_rsize; i++)
		fft_Q[i] -= mean;

	fftw_execute(fft_plan);

	for (i = 0; i < fft_qsize; i++)
		uacf_fft0[i] = fft_b[i] * fft_S[i];
}

static void uacf_measure(const int offset)
{
	if (!uacf_enabled) return;

	int i;
	complex double *S = uacf_S + offset * fft_qsize;
	complex double *old_fft = uacf_fft + offset * fft_qsize;
	
	for (i = 0; i < fft_qsize; i++)
		S[i] += uacf_fft0[i] * old_fft[i];
}

static void uacf_insert(void)
{
	if (!uacf_enabled) return;

	int i;

	for (i = 0; i < fft_qsize; i++)
		uacf_fft[i] = conj(uacf_fft0[i]);
}

/* FIXME: in the last block the program will shift out of memory! */
static void uacf_shift(const int block)
{
	if (!uacf_enabled) return;

	complex double *bl = uacf_fft + block * elements * fft_qsize;
	memmove(bl + fft_qsize, bl, elements * fft_qsize * sizeof(*bl));
}

static void uacf_static(void)
{
	if (!uacf_enabled) return;

	int i;

	for (i = 0; i < fft_qsize; i++) {
		complex double S = uacf_fft0[i];
		uacf_S0[i] += creal(S * conj(S));
	}
}

/* ------------------------------------------------------------------------- */

/* 4-point functions */

static void fp_loadsave_Q(FILE *FH, const int offset, const enum mode mode)
{
	if (!fp_enabled) return;

	if (mode == LOAD) {
		double Q[5];
		if (fscanf(FH, "%lf %lf %lf %lf %lf\n", &Q[0], &Q[1], &Q[2],
							&Q[3], &Q[4]) != 5)
			fatal(EIO, "Reading Q at offset %d", offset);
		if (ismaster) {
			fp_Q[5 * offset + 0] = Q[0];
			fp_Q[5 * offset + 1] = Q[1];
			fp_Q[5 * offset + 2] = Q[2];
			fp_Q[5 * offset + 3] = Q[3];
			fp_Q[5 * offset + 4] = Q[4];
		}
	} else { /* SAVE */
		if (ismaster)
			fprintf(FH, "%.12lg %.12lg %.12lg %.12lg %.12lg\n",
			       	fp_Q[5 * offset + 0],
			       	fp_Q[5 * offset + 1],
			       	fp_Q[5 * offset + 2],
			       	fp_Q[5 * offset + 3],
			       	fp_Q[5 * offset + 4]);
	}
}

static void fp_loadsave_S4(FILE *FH, const int offset, const enum mode mode)
{
	if (!fp_enabled) return;
	
	int i;
	
	if (mode == LOAD) {
		double *s4 = fp_S4 + offset * 3 * fft_qsize;

		for (i = 0; i < fft_qsize; i++) {
			double x, y, z;

			if (fscanf(FH, "%lf %lf %lf\n", &x, &y, &z) != 3)
				fatal(EINVAL, "Error while reading line"
								" %d", i);
			if (ismaster) {
				s4[3 * i + 0] = x;
				s4[3 * i + 1] = y;
				s4[3 * i + 2] = z;
			}
		}
	} else { /* SAVE */
		double *s4 = calloc(3 * fft_qsize, sizeof(*s4));

		MPI_Reduce(fp_S4 + offset * 3 * fft_qsize, s4, 3 * fft_qsize,
					MPI_DOUBLE, MPI_SUM, 0, comm_grid);
		if (ismaster) {
			for (i = 0; i < fft_qsize; i++) {
				fprintf(FH, "%.12lg %.12lg %.12lg\n",
			       		s4[3 * i + 0],
					s4[3 * i + 1],
					s4[3 * i + 2]);
			}
		}
		free(s4);
	}
}

static void fp_init(struct beads *restrict b)
{
	size_t be = blocks * elements;
	
	nblist_alloc(b, &fp_nbl, groups_per_row * b->groupsize,
			b->passport, b->xv, groups_per_col * b->groupsize,
			b->passport + groups_per_row * b->groupsize,
			nbl_pos, fp_a, NBL_NOFLAGS);

	fp_Q = calloc(5 * be, sizeof(*fp_Q));
	if (fp_Q == NULL) novm("fp_Q");

	fp_S4 = calloc(3 * be * fft_qsize, sizeof(*fp_S4));
	if (fp_S4 == NULL) novm("fp_S4");
}

static void fp_loadsave_header(FILE *FH, struct beads *restrict b,
						const enum mode mode)
{
	if (!fp_enabled) return;
	
	if (mode == LOAD) {
		if (fscanf(FH, "# a=%lg\n", &fp_a) != 1)
			fatal(EIO, "fp_loadsave_header");
		fp_init(b);
	} else {
		if (ismaster)
			fprintf(FH, "# a=%lg\n", fp_a);
	}
}

static void fp_free(void)
{
	if (!fp_enabled) return;
	
	nblist_free(&fp_nbl);
	free(fp_Q);
	free(fp_S4);
}

static void nbl_prep_old_coords(struct beads *restrict b, const int offset)
{
	int i, d;

	for (i = 0; i < b->groupsize; i++) {
		for (d = 0; d < 3; d++) {
			double tmp = oldxv[offset * b->groupsize + i][d];
			while (tmp >= b->l[d]) tmp -= b->l[d];
			while (tmp < 0.) tmp += b->l[d];
			nbl_pos[row_index * b->groupsize + i][d] = tmp;
		}
	}

	MPI_Allgather(MPI_IN_PLACE, 0, MPI_DOUBLE, nbl_pos,
		ARRAY_SIZE(*nbl_pos) * b->groupsize, MPI_DOUBLE, comm_col);
}

static void fp_S4ol(struct beads *restrict b, double w[], const int offset)
{
	VEC2 *xv = oldxv + offset * b->groupsize;
	PASSPORT *pp = b->passport + col_index * b->groupsize;
	int i, j, d;

	MPI_Allreduce(MPI_IN_PLACE, w, b->groupsize * groups_per_col,
						MPI_DOUBLE, MPI_SUM, comm_col);
	memset(fft_Q, 0, fft_rsize * sizeof(*fft_Q));
	for (i = 0; i < b->groupsize; i++) {
		double tmp[2];
		if (!GET_EXISTS(pp[i])) continue;

		for (j = 0; j < 2; j++) {
			d = (j == 0) ? TANG1 : TANG2;
			tmp[j] = xv[i][d];
			while (tmp[j] >= b->l[d]) tmp[j] -= b->l[d];
			while (tmp[j] < 0.) tmp[j] += b->l[d];
		}
		fft_addQ(tmp, w[row_index * b->groupsize + i]);
	}

	fftw_execute(fft_plan);
	for (i = 0; i < fft_qsize; i++)
		fft_S[i] *= fft_b[i];
}

static void fp_measure(struct beads *restrict b, const int offset)
{
	if (!fp_enabled) return;

	int i, d;
	double Q[2];
	double a2 = SQR(fp_a);
	double *wS, *wD;

	nblist_rebuild(&fp_nbl); 

	memset(Q, 0, ARRAY_SIZE(Q) * sizeof(*Q));

	wS = calloc(2 * b->groupsize * groups_per_col, sizeof(*wS));
	wD = wS + b->groupsize * groups_per_col;
	if (wS == NULL) novm("wS");

	for (i = 0; i < fp_nbl.count; i++) {
		int p_row = GET_POS(fp_nbl.passport[2 * i + 0]);
		int p_col = GET_POS(fp_nbl.passport[2 * i + 1])
				- groups_per_row * b->groupsize;
		int i_row = GET_ID(fp_nbl.passport[2 * i + 0]);
		int i_col = GET_ID(fp_nbl.passport[2 * i + 1]);
		double r2 = 0.;

		for (d = 0; d < 3; d++) {
			double tmp = b->xv[p_row][d] - nbl_pos[p_col][d];
			while (tmp >= .5 * b->l[d]) tmp -= b->l[d];
			while (tmp < -.5 * b->l[d]) tmp += b->l[d];
			r2 += SQR(tmp);
		}

		if (r2 < a2) {
			if (i_row == i_col) { /* self part */
				Q[0] += 1.;
				wS[p_col] += 1.;
			} else { /* distinct part */
				Q[1] += 1.;
				wD[p_col] += 1.;
			}
		}
	}

	fp_S4ol(b, wS, offset);
	for (i = 0; i < fft_qsize; i++) {
		fp_S4[3 * (offset * fft_qsize + i) + 0] += creal(fft_S[i] * conj(fft_S[i]));
		fft_tmp[i] = fft_S[i]; 
	}
		
	fp_S4ol(b, wD, offset);
	for (i = 0; i < fft_qsize; i++) {
		fp_S4[3 * (offset * fft_qsize + i) + 1] += creal(fft_tmp[i] * conj(fft_S[i]));
		fp_S4[3 * (offset * fft_qsize + i) + 2] += creal(fft_S[i] * conj(fft_S[i]));
	}

	MPI_Reduce(ismaster ? MPI_IN_PLACE : Q, ismaster ? Q : NULL,
			ARRAY_SIZE(Q), MPI_DOUBLE, MPI_SUM, 0, comm_grid);

	if (ismaster) {
		fp_Q[5 * offset + 0] += Q[0]; /* <Q_S> */
		fp_Q[5 * offset + 1] += Q[1]; /* <Q_D> */
		fp_Q[5 * offset + 2] += Q[0] * Q[0]; /* <Q_S * Q_S> */
		fp_Q[5 * offset + 3] += Q[0] * Q[1]; /* <Q_S * Q_D> */
		fp_Q[5 * offset + 4] += Q[1] * Q[1]; /* <Q_D * Q_D> */
	}

	free(wS);
}

/* ------------------------------------------------------------------------- */

static void vacf_alloc(void)
{
	/* autocorrelation functions in three dimensions */
	vacf_c = calloc(blocks * elements, sizeof(*vacf_c));
	if (vacf_c == NULL) novm("vacf_c");
}

static void vacf_free(void)
{
	if (!vacf_enabled) return;
	free(vacf_c);
}

static void vacf_loadsave_header(const enum mode mode)
{
	if (!vacf_enabled) return;
	if (mode == LOAD) vacf_alloc();
}

static void vacf_loadsave(FILE *FH, const int offset, const enum mode mode)
{
	if (!vacf_enabled) return;
	VEC cx;

	switch (mode) {
	case LOAD:
		if (fscanf(FH, "%lf %lf %lf\n", &cx[0], &cx[1], &cx[2]) != 3)
			fatal(EIO, "Reading vacf offset %d", offset);
		if (ismaster) {
			vacf_c[offset][0] = cx[0];
			vacf_c[offset][1] = cx[1];
			vacf_c[offset][2] = cx[2];
		}
		break;
	case SAVE:
		MPI_Reduce(vacf_c[offset], cx, ARRAY_SIZE(cx), MPI_DOUBLE,
							MPI_SUM, 0, comm_grid);
		if (ismaster) fprintf(FH, "%.12lg %.12lg %.12lg\n",
							cx[0], cx[1], cx[2]);
		break;
	}
}

static void vacf_measure(struct beads *restrict b, const int offset)
{
	if (!vacf_enabled) return;

	VEC2 *pos = b->xv + b->groupsize * col_index;
	VEC2 *old_pos = oldxv + offset * b->groupsize;
	PASSPORT *pp = b->passport + col_index * b->groupsize;
	int i, d;
	
	for (i = 0; i < b->groupsize; i++) {
		if (!GET_EXISTS(pp[i])) continue;

		for (d = 0; d < 3; d++) {
			vacf_c[offset][d] += pos[i][3 + d] * old_pos[i][3 + d];
		}
	}
}

/* ------------------------------------------------------------------------- */

#ifdef MCOM
/* force difference ("delta-f") autocorrelation function */

static void facf_init(struct beads *restrict b)
{
	double rc = mcom_get_rc(b->mcom);
	facf_rbins = (int)ceil(rc / facf_dr) + 1;
	
	facf_g = calloc(blocks * elements * facf_rbins, sizeof(*facf_g));
	if (facf_g == NULL) novm("facf_g");
	facf_h = calloc(blocks * elements * facf_rbins, sizeof(*facf_h));
	if (facf_h == NULL) novm("facf_h");
	facf_c = calloc(blocks * elements * facf_rbins, sizeof(*facf_c));
	if (facf_c == NULL) novm("facf_c");

	/* we don't want any artificial friction in the com-system */
	memset(b->gamma, 0, sizeof(b->gamma));
	memset(b->gamma_p, 0, sizeof(b->gamma_p));
			
	facf_oldf = calloc(blocks * elements * b->groupsize, sizeof(*facf_oldf));
	if (facf_oldf == NULL) novm("facf_oldf");

	nblist_alloc(b, &facf_nbl, groups_per_row * b->groupsize,
			b->passport, b->xv, groups_per_col * b->groupsize,
			b->passport + groups_per_row * b->groupsize,
			nbl_pos, rc, NBL_NOFLAGS);
}

static void facf_free(void)
{
	if (!facf_enabled) return;
	free(facf_g);
	free(facf_h);
	free(facf_c);
	free(facf_oldf);
	nblist_free(&facf_nbl);
}

static void facf_loadsave_header(FILE *FH, struct beads *restrict b,
							const enum mode mode)
{
	if (!facf_enabled) return;
	
	if (!comdyn)
		fatal(EINVAL, "comdyn must be enabled for facf()");
	assert(b == bcom);

	if (mode == LOAD) {
		int tmp;

		mcom_init(FH, b);
		mcom_init_nblist(b, b->mcom);
	
		/* the 2nd argument is ignored at this place and is only needed
		 * for the extract.pl program. */	
		if (fscanf(FH, "# facf_dr=%lg rbins=%d\n", &facf_dr, &tmp) != 2)
			fatal(EIO, "facf_dr");
			
		facf_init(b);
	} else {
		if (ismaster) {
			mcom_print_header(FH, b->mcom);
			fprintf(FH, "# facf_dr=%lg rbins=%d\n", facf_dr,
								facf_rbins);
		}
	}
}

static void facf_loadsave(FILE *FH, struct beads *restrict b,
					const int offset, const enum mode mode)
{
	if (!facf_enabled) return;

	PASSPORT *pp = b->passport + col_index * b->groupsize;
	VEC *buf, *oldf = facf_oldf + offset * b->groupsize;
	int i, d, c2[facf_rbins];
	double g[facf_rbins], h[facf_rbins];
	
	buf = calloc(b->nN, sizeof(*buf));
	if (buf == NULL) novm("buf");
		
	switch (mode) {
	case LOAD:
		for (i = 0; i < facf_rbins; i++) {
			if (fscanf(FH, "%lf %lf %d\n", &g[i], &h[i], &c2[i]) != 3)
				fatal(EINVAL, "Error reading ghc from "
								"line %d", i);
		}
		if (ismaster) {
			for (i = 0; i < facf_rbins; i++) {
				facf_g[offset * facf_rbins + i] = g[i];
				facf_h[offset * facf_rbins + i] = h[i];
				facf_c[offset * facf_rbins + i] = c2[i];
			}
		}
		for (i = 0; i < b->nN; i++) {
			if (fscanf(FH, "%lf %lf %lf\n",
				&buf[i][0], &buf[i][1], &buf[i][2]) != 3)
				fatal(EINVAL, "Error reading pos from "
								"line %d", i);
		}
		for (i = 0; i < b->groupsize; i++) {
			int id = GET_ID(pp[i]);
			if (!GET_EXISTS(pp[i])) continue;
			
			for (d = 0; d < (int)ARRAY_SIZE(*buf); d++)
				oldf[i][d] = buf[id][d];
		}
		break;
	case SAVE:
		MPI_Reduce(facf_g + offset * facf_rbins, g, facf_rbins,
					MPI_DOUBLE, MPI_SUM, 0, comm_grid);
		MPI_Reduce(facf_h + offset * facf_rbins, h, facf_rbins,
					MPI_DOUBLE, MPI_SUM, 0, comm_grid);
		MPI_Reduce(facf_c + offset * facf_rbins, c2, facf_rbins,
					MPI_INT, MPI_SUM, 0, comm_grid);
		if (ismaster) {
			for (i = 0; i < facf_rbins; i++)
				fprintf(FH, "%.12lg %.12lg %d\n",
							g[i], h[i], c2[i]);
		}
		for (i = 0; i < b->groupsize; i++) {
			int id = GET_ID(pp[i]);
			if (!GET_EXISTS(pp[i])) continue;

			for (d = 0; d < (int)ARRAY_SIZE(*buf); d++)
				buf[id][d] = oldf[i][d];
		}
		MPI_Reduce(ismaster ? MPI_IN_PLACE : buf, ismaster ? buf : NULL,
					ARRAY_SIZE(*buf) * b->nN, MPI_DOUBLE,
					MPI_SUM, 0, comm_grid);
		if (ismaster) {
			for (i = 0; i < b->nN; i++) {
				fprintf(FH, "%.12lg %.12lg %.12lg\n",
				       	buf[i][0], buf[i][1], buf[i][2]);
			}
		}

		break;
	}

	free(buf);
}

static void facf_measure(struct beads *restrict b, const int offset)
{
	if (!facf_enabled) return;

	VEC *old_f = facf_oldf + offset * b->groupsize, *fbuf, dr;
	double r2, rc2 = SQR(mcom_get_rc(b->mcom));
	int i, d;

	nblist_rebuild(&facf_nbl); 
	
	fbuf = calloc(groups_per_col * b->groupsize, sizeof(*fbuf));
	if (fbuf == NULL) novm("fbuf");

	MPI_Allgather(old_f, ARRAY_SIZE(*old_f) * b->groupsize, MPI_DOUBLE,
		fbuf, ARRAY_SIZE(*fbuf) * b->groupsize, MPI_DOUBLE, comm_col);

	for (i = 0; i < facf_nbl.count; i++) {
		int i_row = GET_ID(facf_nbl.passport[2 * i + 0]);
		int i_col = GET_ID(facf_nbl.passport[2 * i + 1]);
		if (i_row == i_col) continue;
		
		int p_row = GET_POS(facf_nbl.passport[2 * i + 0]);
		int p_col = GET_POS(facf_nbl.passport[2 * i + 1])
					- groups_per_row * b->groupsize;

		for (d = 0, r2 = 0.; d < 3; d++) {
			dr[d] = b->xv[p_row][d] - nbl_pos[p_col][d];
			while (dr[d] >= .5 * b->l[d]) dr[d] -= b->l[d];
			while (dr[d] < -.5 * b->l[d]) dr[d] += b->l[d];
			r2 += SQR(dr[d]);
		}

		if (r2 < rc2) {
			double sp = 0., sp1 = 0., sp2 = 0.;
			for (d = 0; d < 3; d++) {
				sp += b->f[p_row][d] * fbuf[p_col][d];
				sp1 += b->f[p_row][d] * dr[d];
				sp2 += fbuf[p_col][d] * dr[d];
			}
			int rbin = (int)floor(sqrt(r2) / facf_dr);
			assert(rbin >= 0 && rbin < facf_rbins);

			facf_g[offset * facf_rbins + rbin] += sp;
			facf_h[offset * facf_rbins + rbin] += sp1 * sp2 / r2;
			facf_c[offset * facf_rbins + rbin]++;
		}
	}

	free(fbuf);
}

static void facf_insert(struct beads *restrict b)
{
	if (!facf_enabled) return;

	int i, d;

	for (i = 0; i < b->groupsize; i++) {
		int j = b->groupsize * col_index + i;
		for (d = 0; d < 3; d++)
			facf_oldf[i][d] = b->f[j][d];
	}
}

/* FIXME: in the last block the program will shift out of memory! */
static void facf_shift(struct beads *restrict b, const int block)
{
	if (!facf_enabled) return;

	VEC *bl = facf_oldf + block * elements * b->groupsize;
	memmove(bl + b->groupsize, bl, elements * b->groupsize * sizeof(*bl));
}

#endif

/* ------------------------------------------------------------------------- */

/* particle current autocorrelation function (jacf) */
/* See Hansen/McDonald eqs. (7.4.6), (7.4.7), and (7.4.24) for definitions */  

static void jacf_alloc(void)
{
	size_t be = blocks * elements;
	int i;

	for (i = 0; i < ARRAY_SIZE(jacf_s_fft); i++) {
		/* block buffers */
		jacf_s_fft[i] = calloc(be * fft_qsize, sizeof(**jacf_s_fft));
		if (jacf_s_fft[i] == NULL) novm("jacf_s_fft[]");
		jacf_d_fft[i] = calloc(be * fft_qsize, sizeof(**jacf_d_fft));
		if (jacf_d_fft[i] == NULL) novm("jacf_d_fft[]");
		
		jacf_s_fft0[i] = calloc(fft_qsize, sizeof(**jacf_s_fft0));
		if (jacf_s_fft0[i] == NULL) novm("jacf_s_fft0[]");
		jacf_d_fft0[i] = calloc(fft_qsize, sizeof(**jacf_d_fft0));
		if (jacf_d_fft0[i] == NULL) novm("jacf_d_fft0[]");
	}

	for (i = 0; i < ARRAY_SIZE(jacf_s_C); i++) {
		/* dynamic correlation function storage */
		jacf_s_C[i] = calloc(be * fft_qsize, sizeof(**jacf_s_C));
		if (jacf_s_C[i] == NULL) novm("jacf_s_C[]");
		jacf_d_C[i] = calloc(be * fft_qsize, sizeof(**jacf_d_C));
		if (jacf_d_C[i] == NULL) novm("jacf_d_C[]");
	
		/* static correlation function */
		jacf_s_C0[i] = calloc(fft_qsize, sizeof(**jacf_s_C0));
		if (jacf_s_C0[i] == NULL) novm("jacf_s_C0[]");
		jacf_d_C0[i] = calloc(fft_qsize, sizeof(**jacf_d_C0));
		if (jacf_d_C0[i] == NULL) novm("jacf_d_C0[]");
	}
}

static void jacf_free(void)
{
	if (!jacf_enabled) return;
	
	int i;
	for (i = 0; i < ARRAY_SIZE(jacf_s_fft); i++) {
		free(jacf_s_fft[i]);
		free(jacf_d_fft[i]);
		free(jacf_s_fft0[i]);
		free(jacf_d_fft0[i]);
	}
	for (i = 0; i < ARRAY_SIZE(jacf_s_C); i++) {
		free(jacf_s_C[i]);
		free(jacf_d_C[i]);
		free(jacf_s_C0[i]);
		free(jacf_d_C0[i]);
	}
}

static void jacf_loadsave_header(const enum mode mode)
{
	if (!jacf_enabled) return;
	if (mode == LOAD) jacf_alloc();
}

static void jacf_loadsave_C0(FILE *FH, const enum mode mode)
{
	if (!jacf_enabled) return;
	
	double t1, t2;
	int i, j;

	switch (mode) {
	case LOAD:
		for (i = 0; i < fft_qsize; i++) {
			for (j = 0; j < ARRAY_SIZE(jacf_s_C0); j++) {
				if (fscanf(FH, "%lf %lf ", &t1, &t2) != 2)
					fatal(EIO, "reading C0 line %d "
							"element %d", i, j);
				jacf_s_C0[j][i] = (ismaster) ? t1 : 0.;
				jacf_d_C0[j][i] = (ismaster) ? t2 : 0.;
			}
		}
		break;
	case SAVE:
		if (ismaster) {
			for (i = 0; i < fft_qsize; i++) {
				for (j = 0; j < ARRAY_SIZE(jacf_s_C0); j++) {
					fprintf(FH, "%.12lg %.12lg ",
				       	jacf_s_C0[j][i], jacf_d_C0[j][i]);
				}
				fprintf(FH, "\n");
			}
		}
		break;
	}
}

static void jacf_loadsave_fft(FILE *FH, const int offset, const enum mode mode)
{
	if (!jacf_enabled) return;

	int i, j, o = offset * fft_qsize;
	double re1, im1, re2, im2;

	switch (mode) {
	case LOAD:
		for (i = 0; i < fft_qsize; i++) {
			for (j = 0; j < ARRAY_SIZE(jacf_s_fft); j++) {
				if (fscanf(FH, "%lf %lf %lf %lf ", &re1,
						       	&im1, &re2, &im2) != 4)
					fatal(EINVAL, "Error while read"
					"ing jacf_fft line %d element %d",i,j);
				jacf_s_fft[j][i + o] = re1 + im1 * I;
				jacf_d_fft[j][i + o] = re2 + im2 * I;
			}
			for (j = 0; j < ARRAY_SIZE(jacf_s_C); j++) {
				if (fscanf(FH, "%lf %lf %lf %lf ",
						&re1, &im1, &re2, &im2) != 4)
					fatal(EINVAL, "Error while read"
					"ing jacf_C line %d element %d", i, j);
				jacf_s_C[j][i+o] = ismaster ? re1 + im1*I : 0.;
				jacf_d_C[j][i+o] = ismaster ? re2 + im2*I : 0.;
			}
		}
		break;
	case SAVE:
		if (ismaster) {
			for (i = 0; i < fft_qsize; i++) {
				for (j = 0; j < ARRAY_SIZE(jacf_s_fft); j++) {
					fprintf(FH, "%.12lg %.12lg %.12lg %.12lg ",
						creal(jacf_s_fft[j][i + o]),
						cimag(jacf_s_fft[j][i + o]),
						creal(jacf_d_fft[j][i + o]),
						cimag(jacf_d_fft[j][i + o]));
				}
				for (j = 0; j < ARRAY_SIZE(jacf_s_C); j++) {
					fprintf(FH, "%.12lg %.12lg %.12lg %.12lg ",
						creal(jacf_s_C[j][i + o]),
						cimag(jacf_s_C[j][i + o]),
						creal(jacf_d_C[j][i + o]),
						cimag(jacf_d_C[j][i + o]));

				}
				fprintf(FH, "\n");
			}
		}
		break;
	}
}

static void jacf_dump_current(FILE *FH)
{
	if (!jacf_enabled) return;
	int i, j;
			
	for (i = 0; i < fft_qsize; i++) {
		for (j = 0; j < ARRAY_SIZE(jacf_s_fft0); j++) {
			fprintf(FH, "%.12lg %.12lg %.12lg %.12lg ",
			creal(jacf_s_fft0[j][i]), cimag(jacf_s_fft0[j][i]),
			creal(jacf_d_fft0[j][i]), cimag(jacf_d_fft0[j][i]));
		}
		fprintf(FH, "\n");
	}
}

static void jacf_do_fft(struct beads *restrict b)
{
	if (!jacf_enabled) return;

	int dir[ARRAY_SIZE(jacf_s_fft)] = { TANG1, TANG2 };
	
	VEC2 *xv = b->xv + b->groupsize * col_index;
	PASSPORT *pp = b->passport + col_index * b->groupsize;
	int i, j, n, idx;
	double x[2];
	
	for (n = 0; n < ARRAY_SIZE(jacf_s_fft); n++) {
		/* sum of the two monolayers */
		memset(fft_Q, 0, fft_rsize * sizeof(*fft_Q));
		for (i = 0; i < b->groupsize; i++) {
			if (!GET_EXISTS(pp[i])) continue;
			x[0] = xv[i][TANG1];
			x[1] = xv[i][TANG2];
			fft_addQ(x, xv[i][3 + dir[n]]);
		}
		fftw_execute(fft_plan);

		for (i = 0; i < fft_qsize; i++)
			jacf_s_fft0[n][i] = fft_b[i] * fft_S[i];

		MPI_Allreduce(MPI_IN_PLACE, jacf_s_fft0[n], 2 * fft_qsize,
					MPI_DOUBLE, MPI_SUM, comm_grid);

		/* difference between the two monolayers */
		memset(fft_Q, 0, fft_rsize * sizeof(*fft_Q));
		idx = 0;
		for (i = 0; i < b->local_n; i++) {
			double sgn = (b->local_n_b[i] == 0) ? 1. : -1.; 
			for (j = 0; j < b->N[b->local_n_b[i]]; j++, idx++) {
				if (!GET_EXISTS(pp[idx])) continue;
				x[0] = xv[idx][TANG1];
				x[1] = xv[idx][TANG2];
				fft_addQ(x, sgn * xv[idx][3 + dir[n]]);
			}
		}
		fftw_execute(fft_plan);

		for (i = 0; i < fft_qsize; i++)
			jacf_d_fft0[n][i] = fft_b[i] * fft_S[i];

		MPI_Allreduce(MPI_IN_PLACE, jacf_d_fft0[n], 2 * fft_qsize,
					MPI_DOUBLE, MPI_SUM, comm_grid);
	}
}

static void jacf_measure(const int offset)
{
	if (!jacf_enabled) return;

	static const int d = ARRAY_SIZE(jacf_s_fft0);
	int i, j, k;
	for (i = 0; i < d; i++) {
		for (j = 0; j < d; j++) {
			complex double *sC = jacf_s_C[d * i + j] + offset * fft_qsize;
			complex double *dC = jacf_d_C[d * i + j] + offset * fft_qsize;
			complex double *old_s_fft = jacf_s_fft[j] + offset * fft_qsize;
			complex double *old_d_fft = jacf_d_fft[j] + offset * fft_qsize;

			for (k = 0; k < fft_qsize; k++) {
				sC[k] += jacf_s_fft0[i][k] * old_s_fft[k];
				dC[k] += jacf_d_fft0[i][k] * old_d_fft[k];
			}
		}
	}
}

static void jacf_insert(void)
{
	if (!jacf_enabled) return;

	int i, j;

	for (i = 0; i < ARRAY_SIZE(jacf_s_fft); i++) {
		for (j = 0; j < fft_qsize; j++) {
			jacf_s_fft[i][j] = conj(jacf_s_fft0[i][j]);
			jacf_d_fft[i][j] = conj(jacf_d_fft0[i][j]);
		}
	}
}

/* FIXME: in the last block the program will shift out of memory! */
static void jacf_shift(const int block)
{
	if (!jacf_enabled) return;

	int i;

	for (i = 0; i < ARRAY_SIZE(jacf_s_fft); i++) {
		complex double *bl = jacf_s_fft[i] + block * elements * fft_qsize;
		memmove(bl + fft_qsize, bl, elements * fft_qsize * sizeof(*bl));
		bl = jacf_d_fft[i] + block * elements * fft_qsize;
		memmove(bl + fft_qsize, bl, elements * fft_qsize * sizeof(*bl));
	}
}

static void jacf_static(void)
{
	if (!jacf_enabled) return;

	static const int d = ARRAY_SIZE(jacf_s_fft0);
	int i, j, k;

	for (i = 0; i < d; i++) {
		for (j = 0; j < d; j++) {
			for (k = 0; k < fft_qsize; k++) {
				complex double Si = jacf_s_fft0[i][k];
				complex double Sj = jacf_s_fft0[j][k];
				jacf_s_C0[d * i + j][k] += creal(Si * conj(Sj));
				Si = jacf_d_fft0[i][k];
				Sj = jacf_d_fft0[j][k];
				jacf_d_C0[d * i + j][k] += creal(Si * conj(Sj));
			}
		}
	}
}

/* ------------------------------------------------------------------------- */

static void init_blocks(void)
{
	int i, j, k;
	for (k = 1; k < mcount; k++) {
		cur_blocks = 1;
		for (i = k / elements; i > 0; i /= elements)
			cur_blocks++;
		for (i = 0; i < cur_blocks; i++) {
			if (k % (int)(pow(elements + 1, i)) == 0) {
				blocklen[i]++;
				int len = MIN(blocklen[i], elements);
				for (j = 0; j < len; j++)
					c[i * elements + j]++;
			}

		}
	}
}

static void dyn_loadsave_header(FILE *FH, struct beads *restrict b,
							const enum mode mode)
{
	char obs[512], *s;
	int tmp_nN;
	double dt;
	
	if (mode == LOAD) {
		if (fscanf(FH, "# blocks=%d elements=%d count=%d nN=%d dt=%lg "
					"obs=%s\n", &blocks, &elements,
					&mcount, &tmp_nN, &dt, obs) != 6)
			fatal(EINVAL, "header");

		if (b->nN != tmp_nN)
			fatal(EINVAL, "Particle number mismatch %d!=%d",
								b->nN, tmp_nN);
		for (s = obs; *s; s++) {
			switch (*s) {
			case 'i':
				isf_enabled = 1;
				fft_enabled = 1;
				oldxv_enabled = 1;
				break;
			case 'f':
				fp_enabled = 1;
				oldxv_enabled = 1;
				nbl_enabled = 1;
				break;
			case 'v':
				vacf_enabled = 1;
				oldxv_enabled = 1;
				break;
#ifdef MCOM
			case 'F':
				facf_enabled = 1;
				oldxv_enabled = 1;
				nbl_enabled = 1;
				break;
#endif
			case 'u':
				uacf_enabled = 1;
				fft_enabled = 1;
				break;
			case 'j':
				jacf_enabled = 1;
				fft_enabled = 1;
				break;
			default:
				fatal(EINVAL, "Observable '%c' is unknown",*s);
			}
		}
		
		size_t be = blocks * elements;
		c = calloc(be, sizeof(*c));
		if (c == NULL) novm("c");
		blocklen = calloc(blocks, sizeof(*blocklen));
		if (blocklen == NULL) novm("blocklen");
		if (oldxv_enabled) {
			oldxv = calloc(be * b->groupsize, sizeof(*oldxv));
			if (oldxv == NULL) novm("oldxv");
		}

		if (nbl_enabled) {
			nbl_pos = calloc((size_t)groups_per_col * b->groupsize,
							sizeof(*nbl_pos));
			if (nbl_pos == NULL) novm("nbl_pos");
		}

		init_blocks();
	} else {
		memset(obs, 0, sizeof(obs));
		if (isf_enabled) strcat(obs, "i");
		if (fp_enabled) strcat(obs, "f");
		if (vacf_enabled) strcat(obs, "v");
#ifdef MCOM
		if (facf_enabled) strcat(obs, "F");
#endif
		if (uacf_enabled) strcat(obs, "u");
		if (jacf_enabled) strcat(obs, "j");
		
		if (ismaster)
			fprintf(FH, "# blocks=%d elements=%d count=%d nN=%d "
					"dt=%lg obs=%s\n", blocks, elements,
					mcount, b->nN, delta_a, obs);
	}
}

static void dyn_loadsave(FILE *FH, struct beads *restrict b,
						const enum mode mode)
{
	int total, block, elem, i, offset;
	char buf[512];

	/* 1. read/write configuration */
	dyn_loadsave_header(FH, b, mode);
	fft_loadsave_header(FH, b, mode);
	isf_loadsave_header(mode);
	fp_loadsave_header(FH, b, mode);
	vacf_loadsave_header(mode);
#ifdef MCOM
	facf_loadsave_header(FH, b, mode);
#endif
	uacf_loadsave_header(mode);
	jacf_loadsave_header(mode);

	if (mcount == 0) return;

	isf_loadsave_fc0(FH, mode);
	uacf_loadsave_S0(FH, mode);
	jacf_loadsave_C0(FH, mode);

	/* how many elements in total do we have to read? */
	for (i = 0, total = 0; i < blocks; i++)
		total += MIN(blocklen[i], elements);

	/*
	 * in the very last element, no measurement has been performed, yet.
	 * However, there are already coordinates stored. That's why we have to
	 * read one more element.
	 */
	total++;

	for (i = 0; i < total; i++) {
		block = i / elements;
		elem = i % elements;
		
		if (mode == LOAD) {
			fgets(buf, sizeof(buf), FH);
		} else {
			if (ismaster)
				fprintf(FH, "# block %d element %d\n",
								block, elem);
		}
		offset = block * elements + elem;
		
		oldxv_loadsave(FH, b, offset, mode);

		isf_loadsave_r2r4(FH, offset, mode);
		isf_loadsave_fft(FH, offset, mode);

		fp_loadsave_Q(FH, offset, mode);
		fp_loadsave_S4(FH, offset, mode);

		vacf_loadsave(FH, offset, mode);
#ifdef MCOM
		facf_loadsave(FH, b, offset, mode);
#endif
		uacf_loadsave_fft(FH, offset, mode);
		jacf_loadsave_fft(FH, offset, mode);
	}
}

/* ------------------------------------------------------------------------- */

static int get_first_index(void)
{
	int rank, size;
	int i, ret, *lipids;

	MPI_Comm_rank(comm_grid, &rank);
	MPI_Comm_size(comm_grid, &size);

	lipids = alloca(size * sizeof(*lipids));
	lipids[rank] = bcom->local_n;
	
	MPI_Allgather(MPI_IN_PLACE, 0, MPI_INT, lipids, 1, MPI_INT, comm_grid);
	for (i = 0, ret = 0; i < rank; i++)
		ret += lipids[i];
	return ret;
}

/* TODO: this method is inconsistent with the one from mcom.c */
static void calc_cm(const struct beads *restrict b)
{
	int i, j, d, idx, i0;
	VEC *nx = b->nx + col_index * b->groupsize;
	VEC2 *xv = b->xv + col_index * b->groupsize;
	
	/* find index for first molecule of this process */
	i0 = get_first_index();

	for (i = 0, idx = 0; i < b->local_n; i++) {
		VEC2 cm = {0., 0., 0., 0., 0., 0.};
		for (j = 0; j < b->N[b->local_n_b[i]]; j++, idx++) {
			for (d = 0; d < 3; d++)
				cm[d] += xv[idx][d] + nx[idx][d] * b->l[d];
			for (d = 3; d < 6; d++)
				cm[d] += xv[idx][d];
		}
		for (d = 0; d < (int)ARRAY_SIZE(cm); d++)
			cm[d] /= b->N[b->local_n_b[i]];
		set_bead(bcom, cm, i0 + i, col_index * bcom->groupsize + i, 0);
	}
}

static void single_bead(const struct beads *restrict b, const int index2)
{
	int i, d, idx, i0, N;
	VEC *nx = b->nx + col_index * b->groupsize;
	VEC2 *xv = b->xv + col_index * b->groupsize;
	
	/* find index for first molecule of this process */
	i0 = get_first_index();

	for (i = 0, idx = index2; i < b->local_n; i++) {
		VEC2 tmp;
		N = b->N[b->local_n_b[i]];
		if (likely(N > index2)) {
			for (d = 0; d < 3; d++)
				tmp[d] = xv[idx][d] + nx[idx][d] * b->l[d];
			for (d = 3; d < 6; d++)
				tmp[d] = xv[idx][d];
			/* FIXME: type is incorrect, but that doesn't matter */
			set_bead(bcom, tmp, i0 + i,
					col_index * bcom->groupsize + i, 0);
		}
		idx += N;
	}
}

/* ------------------------------------------------------------------------- */

static void comdyn_update(const struct beads *restrict b)
{
	bcom->local_n = b->local_n;
	bcom->time = b->time;
	bcom->step = b->step;
	memcpy(bcom->l, b->l, sizeof(b->l));

	/* clear passports */
	memset(bcom->passport, 0, (groups_per_row + groups_per_col) *
			       bcom->groupsize * sizeof(*bcom->passport));
	if (comdyn >= 1000) {
		single_bead(b, comdyn - 1000);
	} else {
		calc_cm(b);
	}
	bcast_beads(bcom);
}

static void comdyn_load(const struct beads *restrict b)
{
	int sum_n = 0, i, size, *n;

	bcom = calloc(1, sizeof(*bcom));
	if (bcom == NULL) novm("bcom");
	
	bcom->blocks = b->blocks;
	alloc_blocks(bcom);

	MPI_Comm_size(comm_grid, &size);
	n = alloca(size * sizeof(*n));
	
	for (i = 0; i < bcom->blocks; i++) {
		bcom->n[i] = b->n[i];
		bcom->N[i] = 1;
		sum_n += b->n[i];
	}
	MPI_Allgather(&sum_n, 1, MPI_INT, n, 1, MPI_INT, comm_grid);
	for (i = 0; i < size; i++) {
		if (n[i] > bcom->groupsize)
			bcom->groupsize = n[i];
	}
	
	alloc_beads(bcom);
	
	memcpy(bcom->local_n_b, b->local_n_b, b->blocks*sizeof(*b->local_n_b));
	memcpy(bcom->l, b->l, sizeof(b->l));
	bcom->nN = sum_n;
	bcom->dt = b->dt;

	comdyn_update(b);
}

static void comdyn_free()
{
	if (!comdyn) return;

	free_beads(bcom);
	free(bcom);
}

/* ------------------------------------------------------------------------- */

void dyn_init(double _delta_a, size_t bufsize)
{
	inuse = 1;
	delta_a = _delta_a;

	iobufsize = bufsize;
	iobuf = malloc(iobufsize);
	if (iobuf == NULL) novm("iobuf");

	if (comdyn && ismaster)
		printf("COM dynamics enabled, comdyn=%ld.\n", comdyn);

	if (dump_current) {
		cur_iobuf = malloc(iobufsize);
		if (cur_iobuf == NULL) novm("cur_iobuf");
		dumpFH = fopen("dyn_dump.dat", "w");
		if (dumpFH == NULL) fatal(EIO, "Couldn't open "
						"'dyn_dump.dat' for reading");
		setbuffer(dumpFH, cur_iobuf, iobufsize);
	}
}

void dyn_load(const char *fn, struct beads *restrict b)
{
	FILE *FH;
	
	if (!inuse) return;

	FH = fopen(fn, "r");
	if (FH == NULL)
		fatal(ENOENT, "Can't open '%s' for reading", fn);
	setbuffer(FH, iobuf, iobufsize);

	if (comdyn) {
		comdyn_load(b);
		dyn_loadsave(FH, bcom, LOAD);
	} else {
		dyn_loadsave(FH, b, LOAD);
	}

	fclose(FH);
}

void dyn_save(const char *fn, struct beads *restrict b)
{
	FILE *FH = NULL;

	if (!inuse) return;
	
	if (ismaster) {
		FH = fopen(fn, "w");
		if (FH == NULL)
			fatal(ENOENT, "Can't open '%s' for writing", fn);
		setbuffer(FH, iobuf, iobufsize);
	}

	if (comdyn) {
		comdyn_update(b);
		dyn_loadsave(FH, bcom, SAVE);
	} else {
		dyn_loadsave(FH, b, SAVE);
	}
	
	if (ismaster) {
		fclose(FH);
		printf("File '%s' successfully written to disk.\n", fn);
	}
}

static void __dyn_measure(struct beads *restrict b)
{
	int i, j;

	if (!inuse) return;
	
	iboxlen[0] = 1. / b->l[TANG1];
	iboxlen[1] = 1. / b->l[TANG2];

	isfc_fft_density(b);
	uacf_calc_fft(b);
	jacf_do_fft(b);
	if (ismaster) {
		isfc_static();
		uacf_static();
		jacf_static();

		if (dump_current) {
			fprintf(dumpFH, "# time=%lg\n", b->time);
			//isfc_dump_current(dumpFH);
			uacf_dump_current(dumpFH);
			//jacf_dump_current(dumpFH);
		}
	}

	while (mcount / (int)pow(elements, cur_blocks) > 0)
		cur_blocks++;

	for (i = cur_blocks - 1; i >= 0 && mcount > 0; i--) {
		if (mcount % (int)(pow(elements + 1, i)) == 0) {
			blocklen[i]++;
			int len = MIN(blocklen[i], elements);
			for (j = 0; j < len; j++) {
				int offset = i * elements + j;
	
				if (nbl_enabled)
					nbl_prep_old_coords(b, offset);

				isfc_measure(offset);
				isfs_measure(b, offset);
				fp_measure(b, offset);
				vacf_measure(b, offset);
#ifdef MCOM
				facf_measure(b, offset);
#endif
				uacf_measure(offset);
				jacf_measure(offset);
				c[offset]++;
			}

			/* 2) update blocks */
			oldxv_shift(b, i);
			isfc_shift(i);
#ifdef MCOM
			facf_shift(b, i);
#endif
			uacf_shift(i);
			jacf_shift(i);
		}
	}
	/* put current configuration to the front */
	oldxv_insert(b);
	isfc_insert();
#ifdef MCOM
	facf_insert(b);
#endif
	uacf_insert();
	jacf_insert();
	mcount++;
}

#ifdef MCOM
/*
 * Communicate forces and sum up the intra- and intermolecular contributions.
 * After this function this process has full knowledge of the forces acting
 * on its local monomers (but not about non-local monomers).
 */
static void add_forces(const struct beads *restrict b)
{
	int i, d, d0;
	VEC *f_row = b->f + b->groupsize * col_index;
	VEC *f_col = b->f + b->groupsize * (groups_per_row + row_index);

	MPI_Allreduce(MPI_IN_PLACE, b->f,
				3 * groups_per_row * b->groupsize,
				MPI_DOUBLE, MPI_SUM, comm_row);
	MPI_Allreduce(MPI_IN_PLACE, b->f + b->groupsize * groups_per_row,
				3 * groups_per_col * b->groupsize,
				MPI_DOUBLE, MPI_SUM, comm_col);

	/* check for dimensionality, FIXME: Q&D */
	d0 = (mcom_get_diml(b->mcom) == 3) ? 0 : 1;

	/* We subtract f_intra to obtain the force difference between the
	 * MCOM forces and the atomistic forces. */
	for (i = 0; i < b->groupsize; i++) {
		assert(d0 == 0 || (f_row[i][0] == 0. && f_col[i][0] == 0.));
		for (d = d0; d < 3; d++)
			f_row[i][d] += f_col[i][d] - b->f_intra[i][d];
	}

	MPI_Allgather(MPI_IN_PLACE, 0, MPI_DOUBLE, b->f, 3 * b->groupsize,
							MPI_DOUBLE, comm_row);
}
#endif

void dyn_measure(struct beads *restrict b)
{
	if (comdyn) {
		comdyn_update(b);
#ifdef MCOM
		if (facf_enabled) {
			memset(bcom->f, 0, (groups_per_row + groups_per_col) *
					bcom->groupsize * sizeof(*bcom->f));
			mcom_calc_forces(bcom, bcom->mcom, NBL_INVALID);

			assert(bcom->groupsize >= b->local_n);
			mcom_calc_resulting_forces(b, bcom->f_intra);
			add_forces(bcom);
		}
#endif
		__dyn_measure(bcom);
	} else {
		__dyn_measure(b);
	}
}

void dyn_free(void)
{
	if (!inuse) return;

	free(c);
	free(blocklen);
	free(iobuf);
	if (oldxv_enabled)
		free(oldxv);
	if (nbl_enabled)
		free(nbl_pos);
	if (dump_current) {
		fclose(dumpFH);
		free(cur_iobuf);
	}
	
	isf_free();
	fp_free();
	vacf_free();
	fft_free();
#ifdef MCOM
	facf_free();
#endif
	uacf_free();
	jacf_free();
	comdyn_free();
}

