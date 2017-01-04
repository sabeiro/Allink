/*
 * friction.c - calculate friction coefficients for coarse-grained interactions
 * (C) Copyright 2010 by Martin Hoemberg <mhoembe@gwdg.de>
 *
 * Syntax: ./friction [ fm.txt ] [ comforces.dat | - ]
 *
 * "fm.txt" is the output of forcematch and "comforces.dat" is the output
 * of comforces. "-" specifies that these data are read from stdin.
 */

/* $Id$ */

#include "common.h"
#include <gsl/gsl_bspline.h>

/* ------------------------------------------------------------------------- */

typedef double VEC[3];
typedef double VEC2[6];

static int ncoeffs = 0;			/* # b-spline coefficients */
static double *coeffs = NULL;		/* b-spline coefficients */
static gsl_matrix *dB = NULL;		/* b-spline evaluation matrix */
static double rc = 0.;			/* interaction cutoff */
static double rc2 = 0.;			/* rc^2 */
static gsl_bspline_workspace *bw;
static gsl_bspline_deriv_workspace *bdw;
static double L[3];			/* box lengths */
static int nmols = -1;			/* # molecules */
static VEC2 *com;			/* com coordinates and velocities */
static VEC2 *f_at_cg;			/* total atomistic & cg forces */
static VEC *re;				/* com coordinates and velocities */
static double dr2_max = 0.;		/* max. distance encountered */
static int dim = 0;			/* 2 or 3 */
static double systime = 0.;

static double *gpara_num = NULL;
static double *gpara_num2 = NULL;
static double *gpara_denom = NULL;
static double *gpara_denom2 = NULL;
static double *gperp_num = NULL;
static double *gperp_num2 = NULL;
static double *gperp_denom = NULL;
static double *gperp_denom2 = NULL;
static int *gcount = NULL;

static double rbw = 0.05;
static int rbins = 0;

/* ------------------------------------------------------------------------- */

/*
 * open file containing the coefficients for the b-splines
 */
static int load_fm(const char *fn)
{
	FILE *FH = fopen(fn, "r");
	char buf[LINE_MAX];
	int i = 0, k, d;

	debug("Trying to open '%s'", fn);	
	if (FH == NULL)
		return error(ENOENT, "Couldn't open '%s' for reading", fn);

	if (get_data_line(FH, "# k=%d ncoeffs=%d rc=%lg dim=%d",
				&k, &ncoeffs, &rc, &d) != 4)
		return error(EINVAL, "Couldn't understand mcom header");
	
	rc2 = SQR(rc);
	debug("k=%d ncoeffs=%d rc=%lg dim=%d", k, ncoeffs, rc, d);

	/* allocate histograms */
	rbins = (int)ceil(rc / rbw) + 1;

	gpara_num = calloc(rbins, sizeof(*gpara_num));
	if (gpara_num == NULL) return novm("gpara_num");
	gpara_num2 = calloc(rbins, sizeof(*gpara_num2));
	if (gpara_num2 == NULL) return novm("gpara_num2");
	gpara_denom = calloc(rbins, sizeof(*gpara_denom));
	if (gpara_denom == NULL) return novm("gpara_denom");
	gpara_denom2 = calloc(rbins, sizeof(*gpara_denom2));
	if (gpara_denom2 == NULL) return novm("gpara_denom2");
	
	gperp_num = calloc(rbins, sizeof(*gperp_num));
	if (gperp_num == NULL) return novm("gperp_num");
	gperp_num2 = calloc(rbins, sizeof(*gperp_num2));
	if (gperp_num2 == NULL) return novm("gperp_num2");
	gperp_denom = calloc(rbins, sizeof(*gperp_denom));
	if (gperp_denom == NULL) return novm("gperp_denom");
	gperp_denom2 = calloc(rbins, sizeof(*gperp_denom2));
	if (gperp_denom2 == NULL) return novm("gperp_denom2");
	
	gcount = calloc(rbins, sizeof(*gcount));
	if (gcount == NULL) return novm("gcount");
	
	bw = gsl_bspline_alloc(k, ncoeffs - k + 2);
	gsl_bspline_knots_uniform(0., rc, bw);
	bdw = gsl_bspline_deriv_alloc(k);

	dB = gsl_matrix_alloc(k, 2);

	coeffs = calloc(2 * ncoeffs, sizeof(*coeffs));
	if (coeffs == NULL) return novm("coeffs");

	switch (d) {
	case 2:
		debug("mcom: setting 2d interactions");
		dim = 2;
		break;
	case 3:
		debug("mcom: setting 3d interactions");
		dim = 3;
		break;
	default:
		return error(EINVAL, "mcom dimension d=%d", d);
	}

	while (fgets(buf, sizeof(buf), FH) != NULL && i < ncoeffs) {
		int tmp;

		if (strchr(buf, '#') == buf) continue;
		if (sscanf(buf, "%d %lg %lg", &tmp, coeffs + i,
					coeffs + ncoeffs + i) != 3)
			return error(EINVAL, "Error reading coeff %d", i);
		if (tmp != i) return error(EINVAL, "Error reading '%s'", buf);
		i++;
	}
	fclose(FH);

	if (i < ncoeffs)
		return error(EINVAL, "Only %d coefficients found, but expected"
				" %d in file '%s'", i, ncoeffs, fn);
	return 0;
}

static int load_comforces(FILE *FH)
{
	int n, i;
	char buf[LINE_MAX];

	if (fgets(buf, sizeof(buf), FH) == NULL) return 1;
	
	if (sscanf(buf, "# L=%lg %lg %lg N=%d t=%lg\n", &L[0], &L[1], &L[2],
							&n, &systime) != 5)
		return error(EINVAL, "load_comforces header");
	fprintf(stderr, "calculating t=%lg n=%d\n", systime, n);

	if (nmols < 0) {
		nmols = n;
		com = calloc(n, sizeof(*com));
		if (com == NULL) return novm("load_comforces com");
		re = calloc(n, sizeof(*re));
		if (re == NULL) return novm("load_comforces re");
		f_at_cg = calloc(n, sizeof(*f_at_cg));
		if (f_at_cg == NULL) return novm("load_comforces f_at_cg");
	} else if (nmols != n) {
		return error(EINVAL, "nmols != n (%d, %d)", nmols, n);
	}

	memset(f_at_cg, 0, nmols * sizeof(*f_at_cg));

	for (i = 0; i < n; i++) {

		if (fgets(buf, sizeof(buf), FH) == NULL)
			return error(EINVAL, "premature eof");
		if (sscanf(buf, "%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg",
				&com[i][0], &com[i][1], &com[i][2],
				&com[i][3], &com[i][4], &com[i][5],
				&re[i][0], &re[i][1], &re[i][2],
				&f_at_cg[i][0], &f_at_cg[i][1], &f_at_cg[i][2]) != 12)
		       return error(EINVAL, "invalid line '%s'", buf);
	}
	
	return 0;
}

static void handle_ia(VEC2 drdv, VEC f_at, VEC f_cg, double r)
{
	double x, a = drdv[0] / r, b = drdv[1] / r, c = drdv[2] / r;
	int bin = (int)(r / rbw + .5);

	assert(bin >= 0 && bin < rbins);

	/* parallel component */
	x = (f_cg[0] - f_at[0]) * a + (f_cg[1] - f_at[1]) * b +
						(f_cg[2] - f_at[2]) * c;
//	printf("%lg ", x);
	gpara_num[bin] += x;
	gpara_num2[bin] += x * x;

	x = drdv[3] * a + drdv[4] * b + drdv[5] * c;
//	printf("%lg\n", x);
	gpara_denom[bin] += x;
	gpara_denom2[bin] += x * x;

	/* perpendicular components */
	VEC a1, a2;
	double q = 1. / hypot(b, c);
	a1[0] = 0.;
	a1[1] = -c * q;
	a1[2] = b * q;
	a2[0] = - 1. * (b * b + c * c) * q;
	a2[1] = a * b * q;
	a2[2] = a * c * q;
	x = (a1[0] + a2[0]) * f_at[0] + (a1[1] + a2[1]) * f_at[1] +
						(a1[2] + a2[2]) * f_at[2];
	gperp_num[bin] += x;
	gperp_num2[bin] += x * x;
	x = (a1[0] + a2[0]) * drdv[3] + (a1[1] + a2[1]) * drdv[4] +
						(a1[2] + a2[2]) * drdv[5];
	gperp_denom[bin] += x;
	gperp_denom2[bin] += x * x;
	
	/* counter */
	gcount[bin]++;
}

static int read_ia(FILE *FH)
{
	char buf[LINE_MAX];
	int id1, id2, d, i;
	VEC f_cg, f_at;
	VEC2 drdv;
	double dr2;
	
	do {
		if (fgets(buf, sizeof(buf), FH) == NULL)
			return error(EINVAL, "read_ia: Premature eof");
		if (strncmp(buf, "# --END--", 9) == 0)
			break;
		if (sscanf(buf, "%d %d %lg %lg %lg", &id1, &id2,
					&f_at[0], &f_at[1], &f_at[2]) != 5)
			return error(EINVAL, "read_ia: line '%s'", buf);

		/* calculate drdv */
		if (dim == 2) {
			f_at[0] = 0.;
			drdv[0] = 0.;
			drdv[3] = 0.;
			for (d = 1, dr2 = 0.; d < 3; d++) {
				drdv[d] = com[id2][d] - com[id1][d];
				while (drdv[d] >= .5 * L[d]) drdv[d] -= L[d];
				while (drdv[d] < -.5 * L[d]) drdv[d] += L[d];
				drdv[d + 3] = com[id2][d + 3] - com[id1][d + 3];
				dr2 += SQR(drdv[d]);
			}
		} else if (dim == 3) {
			for (d = 0, dr2 = 0.; d < 3; d++) {
				drdv[d] = com[id2][d] - com[id1][d];
				while (drdv[d] >= .5 * L[d]) drdv[d] -= L[d];
				while (drdv[d] < -.5 * L[d]) drdv[d] += L[d];
				drdv[d + 3] = com[id2][d + 3] - com[id1][d + 3];
				dr2 += SQR(drdv[d]);
			}
		}
		if (dr2 > dr2_max) dr2_max = dr2;

		if (dr2 < rc2) {
			double f = 0., r = sqrt(dr2);
			size_t is, ie;
			gsl_bspline_deriv_eval_nonzero(r, 1, dB, &is, &ie, bw, bdw);
			int o = is + ((re[id1][0] * re[id2][0] < 0.) ? ncoeffs : 0);
			for (i = 0; i < bw->k; i++)
				f -= coeffs[o + i] * gsl_matrix_get(dB, i, 1);

			f /= r;
			f_cg[0] = f * drdv[0];
			f_cg[1] = f * drdv[1];
			f_cg[2] = f * drdv[2];

			f_at_cg[id1][3] -= f_cg[0];
			f_at_cg[id1][4] -= f_cg[1];
			f_at_cg[id1][5] -= f_cg[2];
			f_at_cg[id2][3] += f_cg[0];
			f_at_cg[id2][4] += f_cg[1];
			f_at_cg[id2][5] += f_cg[2];
		
			handle_ia(drdv, f_at, f_cg, r);
		}
	} while (1);
	return 0;
}

static void output_df(void)
{
	int i;

	printf("# t=%lg n=%d L=%lg %lg %lg\n", systime, nmols, L[0], L[1], L[2]);
	for (i = 0; i < nmols; i++) {
		VEC f;
		f[0] = f_at_cg[i][0] - f_at_cg[i][3];
		f[1] = f_at_cg[i][1] - f_at_cg[i][4];
		f[2] = f_at_cg[i][2] - f_at_cg[i][5];
		printf("%lg %lg %lg %lg %lg %lg\n",
			com[i][0], com[i][1], com[i][2], f[0], f[1], f[2]);
	}
}

static void cleanup(void)
{
	gsl_bspline_free(bw);
	gsl_bspline_deriv_free(bdw);
	gsl_matrix_free(dB);
	free(coeffs);
	free(com);
	free(re);
	free(f_at_cg);
	free(gpara_num);
	free(gpara_num2);
	free(gpara_denom);
	free(gpara_denom2);
	free(gperp_num);
	free(gperp_num2);
	free(gperp_denom);
	free(gperp_denom2);
	free(gcount);
}

int main(int argc, char **argv)
{
	FILE *FH;
	int ret;

	if (argc < 3) {
		fprintf(stderr, "Syntax: %s [fm.txt] [comforces.dat|-]\n",
								argv[0]);
		return EXIT_FAILURE;
	}

	ret = load_fm(argv[1]);
	if (ret) return EXIT_FAILURE;

	if (strcmp(argv[2], "-") == 0) {
		FH = stdin;
	} else {
		FH = fopen(argv[2], "r");
		if (FH == NULL) {
			error(ENOENT, "Couldn't open %s for reading", argv[2]);
			return EXIT_FAILURE;
		}
	}

	while (load_comforces(FH) == 0) {
		ret = read_ia(FH);
		if (ret) return EXIT_FAILURE;
//		output_df();
	}

	fprintf(stderr, "dr_max = %lg\n", sqrt(dr2_max));

	int i;
	for (i = 0; i < rbins; i++) {
		if (gcount[i] == 0) continue;
		gpara_num[i] /= gcount[i];
		gpara_num2[i] /= gcount[i];
		gpara_denom[i] /= gcount[i];
		gpara_denom2[i] /= gcount[i];
		gperp_num[i] /= gcount[i];
		gperp_num2[i] /= gcount[i];
		gperp_denom[i] /= gcount[i];
		gperp_denom2[i] /= gcount[i];
		
		/* print parallel friction */
		double dn, dd, g;
		dn = sqrt((gpara_num2[i] - SQR(gpara_num[i])) / (gcount[i] - 1.));
		dd = sqrt((gpara_denom2[i] - SQR(gpara_denom[i])) / (gcount[i] - 1.));
		g = gpara_num[i] / gpara_denom[i];
		fprintf(stderr, "%d %lg %lg %lg", i, (i + .5) * rbw,
			g, fabs(g) * hypot(dn / gpara_num[i], dd / gpara_denom[i]));

		/* print perpendicular friction */
		dn = sqrt((gperp_num2[i] - SQR(gperp_num[i])) / (gcount[i] - 1.));
		dd = sqrt((gperp_denom2[i] - SQR(gperp_denom[i])) / (gcount[i] - 1.));
		g = gperp_num[i] / gperp_denom[i];
		fprintf(stderr, " %lg %lg\n", -g, fabs(g) * hypot(dn / gperp_num[i], dd / gperp_denom[i]));
	}

	if (FH != stdin) fclose(FH);
	cleanup();
	return EXIT_SUCCESS;
}

