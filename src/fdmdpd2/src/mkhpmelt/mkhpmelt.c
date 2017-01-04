/*
 * mkhpmelt.c - create homopolyer melt
 * (C) Copyright 2009 by M. Hoemberg <mhoembe@gwdg.de>
 *
 * This program creates an isotropic melt of homopolymers.
 */

#define _GNU_SOURCE

#include "rand.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* Box-Muller algorithm to create gaussian random numbers */
static double gaussian(struct rng *r, double m, double s)
{
	double u1 = rng_uniform(r);
	double u2 = rng_uniform(r);
	double z;

	do {
       		z = sqrt(-2. * log(1. - u1)) * cos(2. * M_PI * u2);
	} while (!isnormal(z));

	return m + s * z; 
}

static int print_hp(struct rng *r, double L[], int N, double Re)
{
	double x[3];
	double v[3]; 
	int i, d;
	double s = sqrt(Re * Re / (N - 1.) / 3.);

	for (d = 0; d < 3; d++)
		x[d] = L[d] * Re * rng_uniform(r);

	for (i = 0; i < N; i++) {

		for (d = 0; d < 3; d++)
			v[d] = gaussian(r, 0., 1.);

		printf("%lg %lg %lg %lg %lg %lg %d\n",
				x[0], x[1], x[2], v[0], v[1], v[2], 0);

		for (d = 0; d < 3; d++)
			x[d] += gaussian(r, 0., s);
	}

	return 0;
}

static int print_system(int argc, char **argv)
{
	double L[3];
	int n;
	int N;
	double kN;
	double rho_coex;
	double Re;
	int i;
	struct rng *r = rng_alloc(0);

	L[0] = atof(argv[1]);
	L[1] = atof(argv[2]);
	L[2] = atof(argv[3]);
	n = atoi(argv[4]);
	N = atoi(argv[5]);
	kN = atof(argv[6]);
	rho_coex = atof(argv[7]);
	Re = atof(argv[8]);

	double v = -2. * (kN + 3.) / rho_coex;
	double w = 1.5 * (kN + 2.) / (rho_coex * rho_coex);

	printf("# L=%lg %lg %lg n=%d N=%d t=0\n", L[0] * Re, L[1] * Re,
						       L[2] * Re, n, N);
	printf("# v=%lg 0 0 w=%lg 0 0 0\n", v, w);

	for (i = 0; i < n; i++)
		print_hp(r, L, N, Re);

	rng_free(r);

	return EXIT_SUCCESS;
}

int main(int argc, char **argv)
{
	if(argc < 9) {
		fprintf(stderr, "Syntax: mkhpmelt [Lx] [Ly] [Lz] [n] [N] "
				"[kappaN] [rho_coex] [R_e/sigma]\n\n");
		return EXIT_FAILURE;
	} else {
		print_system(argc, argv);
		return EXIT_SUCCESS;
	}
}

