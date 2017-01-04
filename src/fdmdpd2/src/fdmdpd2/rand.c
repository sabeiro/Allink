/**
 * rand.c - portable random number generators
 * Copyright (C) 2008 Martin Hoemberg <mhoembe@gwdg.de>
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

/* $Id: rand.c 50 2010-03-02 16:26:29Z hoemberg $ */

#define _GNU_SOURCE
#include <stdlib.h>
#include <math.h>
#include "rand.h"

/* ------------------------------------------------------------------------- */

#ifdef __IBMC__
#include <essl.h>

struct rng
{
	double vseed[10000];
	int iseed;
};

struct rng *rng_alloc(unsigned long seed)
{
	struct rng *rng = malloc(sizeof(*rng));
	rng->iseed = (int)seed;
	return rng;
}

void rng_free(struct rng *restrict rng)
{
	free(rng);
}

double rng_uniform(struct rng *restrict rng)
{
	double r;
	durxor(&rng->iseed, 1, &r, rng->vseed);
	return r; 
}

int rng_uniform_vector(struct rng *restrict rng, int n,
							double *restrict r)
{
	durxor(&rng->iseed, n, r, rng->vseed);
	return 0;
}

/* ------------------------------------------------------------------------- */

#elif defined __MKL__
#include <mkl_vsl.h>

struct rng
{
	VSLStreamStatePtr stream;
};

struct rng *rng_alloc(unsigned long seed)
{
	struct rng *rng = malloc(sizeof(*rng));
	vslNewStream(&rng->stream, VSL_BRNG_MT19937, seed);
	return rng;
}

void rng_free(struct rng *restrict rng)
{
	vslDeleteStream(&rng->stream);
	free(rng);
}

double rng_uniform(struct rng *restrict rng)
{
	double r;
	vdRngUniform(VSL_METHOD_DUNIFORM_STD, rng->stream, 1, &r, 0., 1.);
	return r; 
}

int rng_uniform_vector(struct rng *restrict rng, int n,
							double *restrict r)
{
	vdRngUniform(VSL_METHOD_DUNIFORM_STD, rng->stream, n, r, 0., 1.);
	return 0;
}

/* ------------------------------------------------------------------------- */

#elif defined __GSL__
#include <gsl/gsl_rng.h>

struct rng
{
	gsl_rng *r;
};

struct rng *rng_alloc(unsigned long seed)
{
	struct rng *rng = malloc(sizeof(*rng));
	rng->r = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(rng->r, seed);
	return rng; 
}

void rng_free(struct rng *restrict rng)
{
	gsl_rng_free(rng->r);
	free(rng); 
}

double rng_uniform(struct rng *restrict rng)
{
	return gsl_rng_uniform(rng->r);
}

int rng_uniform_vector(struct rng *restrict rng, int n,
							double *restrict r)
{
	int i;
	for (i = 0; i < n; i++)
		r[i] = gsl_rng_uniform(rng->r);
	return 0;
}

/* ------------------------------------------------------------------------- */

#elif defined _GNU_SOURCE

struct rng
{
	struct drand48_data buf;
};

struct rng *rng_alloc(unsigned long seed)
{
	struct rng *rng = malloc(sizeof(*rng));
	srand48_r(seed, &rng->buf);
	return rng;
}

void rng_free(struct rng *restrict rng)
{
	free(rng);
}

double rng_uniform(struct rng *restrict rng)
{
	double r;
	drand48_r(&rng->buf, &r);
	return r;
}

int rng_uniform_vector(struct rng *restrict rng, int n,
							double *restrict r)
{
	int i;
	for (i = 0; i < n; i++)
		drand48_r(&rng->buf, r + i);
	return 0;
}

/* ------------------------------------------------------------------------- */

#else

#error Cannot find random number generator!!

#endif

/* Box-Muller algorithm to create gaussian random numbers */
double rng_gaussian(struct rng *r, double m, double s)
{
	double u1 = rng_uniform(r);
	double u2 = rng_uniform(r);
	double z;
	
	do {
		z = sqrt(-2. * log(1. - u1)) * cos(2. * M_PI * u2);
	} while (!isnormal(z));

	return m + s * z;
}

