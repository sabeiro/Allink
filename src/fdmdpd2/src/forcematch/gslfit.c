/*
 * gslfit.c - Produce B-Spline coefficients from existing interaction potential
 * (C) Copyright 2010 by Martin Hoemberg <mhoembe@gwdg.de>
 *
 * Syntax: ./gslfit [k] [ncoeffs] [rc] < input > output
 *
 */

#include "common.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include <assert.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_multifit.h>

#define POINTS_MAX 		10000

static int k;			/* degree of B-Splines */
static int ncoeffs;		/* number of coefficients */
static double rc;		/* cutoff of the potential */
static double r[POINTS_MAX];	/* r values */
static double U[POINTS_MAX];	/* U values */
static int count = 0;		/* # points */

static int read_data(FILE *FH)
{
	char buf[LINE_MAX];

	while (!feof(FH)) {
		if (fgets(buf, sizeof(buf), FH) == NULL) continue;
		if (buf[0] == '#') continue;
		assert(count < POINTS_MAX);
		if (sscanf(buf, "%lg %lg", &r[count], &U[count]) != 2)
			return error(EINVAL, "line '%s' not understood", buf);
		count++;
	}
	fprintf(stderr, "Read %d lines successfully.\n", count);
	return 0;
}

static void gslfit_data(void)
{
	int i, j;
	gsl_matrix *X = gsl_matrix_alloc(count, ncoeffs);
	gsl_vector *c = gsl_vector_alloc(ncoeffs);
	gsl_matrix *cov = gsl_matrix_alloc(ncoeffs, ncoeffs);
	gsl_multifit_linear_workspace *mw = gsl_multifit_linear_alloc(count,
								ncoeffs);
	gsl_bspline_workspace *bw = gsl_bspline_alloc(k, ncoeffs - 2);
	gsl_vector *B = gsl_vector_alloc(ncoeffs);
	double chisq = 0.;
	       
	gsl_bspline_knots_uniform(0., rc, bw);

	for (i = 0; i < count; i++) {
		gsl_bspline_eval(r[i], B, bw);
		for (j = 0; j < ncoeffs; j++) {
			gsl_matrix_set(X, i, j, gsl_vector_get(B, j));
		}
	}

	gsl_vector_view y = gsl_vector_view_array(U, count);
	gsl_multifit_linear(X, &y.vector, c, cov, &chisq, mw);
	
	for (i = 0; i < ncoeffs; i++)
		printf("%d %lg\n", i, gsl_vector_get(c, i));

	gsl_vector_free(c);       
	gsl_matrix_free(X);
	gsl_matrix_free(cov);
	gsl_multifit_linear_free(mw);
	gsl_bspline_free(bw);
	gsl_vector_free(B);
}

int main(int argc, char **argv)
{
	if (argc < 4) {
		fprintf(stderr, "Syntax: %s [k] [ncoeffs] [rc]\n\n",
								argv[0]);
		return EXIT_FAILURE;
	}

	k = atoi(argv[1]);
	ncoeffs = atoi(argv[2]);
	rc = atof(argv[3]);

	fprintf(stderr, "k=%d ncoeffs=%d rc=%lg\n", k, ncoeffs, rc);

	int ret = read_data(stdin);
	if (ret) return EXIT_FAILURE;

	gslfit_data();
	return EXIT_SUCCESS;
}

