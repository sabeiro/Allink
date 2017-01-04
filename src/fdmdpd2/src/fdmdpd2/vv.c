/*
 * vv.c - Velocity Verlet integrators in various ensembles
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

/* $Id: vv.c 366 2013-08-08 14:19:08Z fuhrmans $ */

/*
 * This module implements several symplectic integration algorithms, which are
 * all based on the same general method described in the article:
 * A. Kolb and B. Duenweg, JCP 111, 4453 (1999)
 */ 

#include "fdmdpd2.h"
#include "rand.h"

/* ------------------------------------------------------------------------- */

/* Static Variables for the Langevin Piston */
static struct rng *vvrng = NULL;/* random numbers for Langevin piston */
static double P0 = 0.;		/* target pressure */
static double Q = 1e-4;		/* mass of the integrator DoF */
static double iQ = 0.;		/* inverse mass of the integrator DoF */
static double pgamma = 0.1;	/* friction constant */
static double psigma = 0.;	/* noise coefficient */

/* Static Variables for the NPtT ensemble */
static double A;		/* current area A(t) of the box */
static double pA = 0.;		/* area momentum */
static double rt = 1.;		/* constant ratio, L_y / L_z */

/* Static Variables for the NPtT ensemble with 2 DoFs */
static double L[2];		/* L_y(t) and L_z(t) */
static double pL[2] = {0., 0.};	/* conjugated momenta */

/* Static Variables for the NPT ensemble */
static double V;		/* current volume V(t) */
static double pV = 0.;		/* volume momentum */

/* Function pointers to the first and second part of the integration */
static void (*do_step1)(struct beads *restrict b);
static void (*do_step2)(struct beads *restrict b);

static int first = 0;		/* first bead to be integrated */
static int last = 0;		/* last bead to be integrated */

cfg_opt_t integrator_opts[] = {
	CFG_STR("ensemble", "NVT", CFGF_NONE),
	CFG_SIMPLE_FLOAT("Q", &Q),
	CFG_SIMPLE_FLOAT("gamma", &pgamma),
	CFG_SIMPLE_FLOAT("P0", &P0),
	CFG_END()
};

/* ------------------------------------------------------------------------- */

/* First part of NVT integration */
static void nvt_step1(struct beads *restrict b)
{
	int i, d;

	for (i = first; i < last; i++) {
		VEC x, dn;

		for (d = 0; d < 3; d++) {
			x[d] = b->xv[i][d] + b->xv[i][d + 3] * b->dt;
			dn[d] = floor(x[d] / b->l[d]);
			b->xv[i][d] = x[d] - dn[d] * b->l[d];
			b->nx[i][d] += dn[d];
		}
	}
}

/* Second part of NVT integration */
static void nvt_step2(struct beads *restrict b)
{
}

/* ------------------------------------------------------------------------- */

/* First part of NPT integration */ 
static void npt_step1(struct beads *restrict b)
{
	int i, d;
	double P, xi, l1, l2;
	TENSOR p;

	/* step 2 */
	get_pressure(b, &p);
	P = (p[0] + p[1] + p[2]) / 3.;
	xi = rng_uniform(vvrng) - .5; 
	pV += .5 * (P - P0) * b->dt;
	pV -= .5 * pgamma * pV * iQ * b->dt;
	pV += psigma * xi;

	/* step 3 */
	l1 = V;
	l2 = 1. / V;
	V += .5 * pV * iQ * b->dt;
	l1 /= V;
	l1 = pow(l1, 2. / 3.);
	
	/* step 5 */
	V += .5 * pV * iQ * b->dt;
	l2 *= V;
	l2 = pow(l2, 1. / 3.);
	b->l[0] = pow(V, 1. / 3.);
	b->l[1] = b->l[0];
	b->l[2] = b->l[0];

	/* step 4 and rescaling */
	for (i = first; i < last; i++) {
		VEC x, dn;

		for (d = 0; d < 3; d++) {
			x[d] = b->xv[i][d] + l1 * b->xv[i][d + 3] * b->dt;
			x[d] *= l2;
			b->xv[i][d + 3] /= l2;
			dn[d] = floor(x[d] / b->l[d]);
			b->xv[i][d] = x[d] - dn[d] * b->l[d];
			b->nx[i][d] += dn[d];
		}
	}
}

/* Second part of NPT integration */
static void npt_step2(struct beads *restrict b)
{
	double xi, P;
	TENSOR p;
	
	/* step 6 */
	get_pressure(b, &p);
	P = (p[0] + p[1] + p[2]) / 3.;
	xi = rng_uniform(vvrng) - .5;
	pV += .5 * V * (P - P0) * b->dt;
	pV -= .5 * pgamma * pV * iQ * b->dt;
	pV += psigma * xi;
}

/* ------------------------------------------------------------------------- */

/* First part of NPtT integration */ 
static void nptt_step1(struct beads *restrict b)
{
	int i, d;
	double Pt, xi, A1, A2;
	TENSOR p;

	/* step 2 */
	get_pressure(b, &p);
	Pt = .5 * (p[TANG1] + p[TANG2]);
	xi = rng_uniform(vvrng) - .5; 
	pA += .5 * b->l[NORMAL] * (Pt - P0) * b->dt;
	pA -= .5 * pgamma * pA * iQ * b->dt;
	pA += psigma * xi;

	/* step 3 */
	A1 = A;
	A2 = 1. / A;
	A += .5 * pA * iQ * b->dt;
	A1 /= A;
	
	/* step 5 */
	A += .5 * pA * iQ * b->dt;
	A2 *= A;
	b->l[TANG1] = sqrt(A * rt);
	b->l[TANG2] = sqrt(A / rt);

	/* step 4 and rescaling */
	for (i = first; i < last; i++) {
		VEC x, dn;
	
		x[NORMAL] = b->xv[i][NORMAL] + b->xv[i][NORMAL + 3] * b->dt;
		x[TANG1] = b->xv[i][TANG1] + A1 * b->xv[i][TANG1 + 3] * b->dt;
		x[TANG2] = b->xv[i][TANG2] + A1 * b->xv[i][TANG2 + 3] * b->dt;
		x[TANG1] *= sqrt(A2);
		x[TANG2] *= sqrt(A2);
		b->xv[i][TANG1 + 3] /= sqrt(A2);
		b->xv[i][TANG2 + 3] /= sqrt(A2);

		for (d = 0; d < 3; d++) {
			dn[d] = floor(x[d] / b->l[d]);
			b->xv[i][d] = x[d] - dn[d] * b->l[d];
			b->nx[i][d] += dn[d];
		}
	}
}

/* Second part of NPtT integration */
static void nptt_step2(struct beads *restrict b)
{
	double xi, Pt;
	TENSOR p;
	
	/* step 6 */
	get_pressure(b, &p);
	Pt = .5 * (p[TANG1] + p[TANG2]);
	xi = rng_uniform(vvrng) - .5;
	pA += .5 * b->l[NORMAL] * (Pt - P0) * b->dt;
	pA -= .5 * pgamma * pA * iQ * b->dt;
	pA += psigma * xi;
}
/* ------------------------------------------------------------------------- */

/* First part of NPxT integration */ 
static void npxt_step1(struct beads *restrict b)
{
	int i, d;
	double Pt, xi, A1, A2;
	TENSOR p;

	/* step 2 */
	get_pressure(b, &p);
	Pt = p[TANG1];
	xi = rng_uniform(vvrng) - .5; 
	pA += .5 * b->l[NORMAL] * (Pt - P0) * b->dt;
	pA -= .5 * pgamma * pA * iQ * b->dt;
	pA += psigma * xi;

	/* step 3 */
	A1 = A;
	A2 = 1. / A;
	A += .5 * pA * iQ * b->dt;
	A1 /= A;
	
	/* step 5 */
	A += .5 * pA * iQ * b->dt;
	A2 *= A;
	b->l[TANG1] = sqrt(A * rt);

	/* step 4 and rescaling */
	for (i = first; i < last; i++) {
		VEC x, dn;
	
		x[NORMAL] = b->xv[i][NORMAL] + b->xv[i][NORMAL + 3] * b->dt;
		x[TANG1] = b->xv[i][TANG1] + A1 * b->xv[i][TANG1 + 3] * b->dt;
		x[TANG2] = b->xv[i][TANG2] + b->xv[i][TANG2 + 3] * b->dt;
		x[TANG1] *= sqrt(A2);
		b->xv[i][TANG1 + 3] /= sqrt(A2);

		for (d = 0; d < 3; d++) {
			dn[d] = floor(x[d] / b->l[d]);
			b->xv[i][d] = x[d] - dn[d] * b->l[d];
			b->nx[i][d] += dn[d];
		}
	}
}

/* Second part of NPxT integration */
static void npxt_step2(struct beads *restrict b)
{
	double xi, Pt;
	TENSOR p;
	
	/* step 6 */
	get_pressure(b, &p);
	Pt = p[TANG1];
	xi = rng_uniform(vvrng) - .5;
	pA += .5 * b->l[NORMAL] * (Pt - P0) * b->dt;
	pA -= .5 * pgamma * pA * iQ * b->dt;
	pA += psigma * xi;
}

/* ------------------------------------------------------------------------- */

/* First part of anisotropic NPtT integration */
static void nptt2_step1(struct beads *restrict b)
{
	int i, d;
	double L1[2], L2[2], tmp;
	TENSOR p;

	/* step 2 */
	get_pressure(b, &p);
	tmp = pL[0];
	pL[0] += (b->l[NORMAL] * L[1] * (p[TANG1] - P0) +
				SQR(pL[1]) * iQ / CUBE(L[0])) * .5 * b->dt;
	pL[1] += (b->l[NORMAL] * L[0] * (p[TANG2] - P0) +
				SQR(tmp) * iQ / CUBE(L[1])) * .5 * b->dt;
	pL[0] -= .5 * pgamma * pL[0] * iQ * b->dt;
	pL[0] += psigma * (rng_uniform(vvrng) - .5);
	pL[1] -= .5 * pgamma * pL[1] * iQ * b->dt;
	pL[1] += psigma * (rng_uniform(vvrng) - .5);

	/* step 3 */
	L1[0] = L[0];
	L1[1] = L[1];
	L2[0] = 1. / L[0];
	L2[1] = 1. / L[1];
	tmp = L[0];
	L[0] += .5 * pL[0] * iQ * b->dt / SQR(L[1]);
	L[1] += .5 * pL[1] * iQ * b->dt / SQR(tmp);
	L1[0] /= L[0];
	L1[1] /= L[1];

	/* step 5 */
	tmp = L[0];
	L[0] += .5 * pL[0] * iQ * b->dt / SQR(L[1]);
	L[1] += .5 * pL[1] * iQ * b->dt / SQR(tmp);
	L2[0] *= L[0];
	L2[1] *= L[1];
	b->l[TANG1] = L[0];
	b->l[TANG2] = L[1];

	/* step 4 and rescaling */
	for (i = first; i < last; i++) {
		VEC y, dn;
	
		y[NORMAL] = b->xv[i][NORMAL] + b->xv[i][NORMAL + 3] * b->dt;
		y[TANG1] = b->xv[i][TANG1] + L1[0] * b->xv[i][TANG1 + 3] * b->dt;
		y[TANG2] = b->xv[i][TANG2] + L1[1] * b->xv[i][TANG2 + 3] * b->dt;
		y[TANG1] *= sqrt(L2[0]);
		y[TANG2] *= sqrt(L2[1]);
		b->xv[i][TANG1 + 3] /= sqrt(L2[0]);
		b->xv[i][TANG2 + 3] /= sqrt(L2[1]);

		for (d = 0; d < 3; d++) {
			dn[d] = floor(y[d] / b->l[d]);
			b->xv[i][d] = y[d] - dn[d] * b->l[d];
			b->nx[i][d] += dn[d];
		}
	}
}

/* Second part of anisotropic NPtT integration */
static void nptt2_step2(struct beads *restrict b)
{
	TENSOR p;
	double tmp;

	/* step 6 */
	get_pressure(b, &p);
	tmp = pL[0];
	pL[0] += (b->l[NORMAL] * L[1] * (p[TANG1] - P0) +
				SQR(pL[1]) * iQ / CUBE(L[0])) * .5 * b->dt;
	pL[1] += (b->l[NORMAL] * L[0] * (p[TANG2] - P0) +
				SQR(tmp) * iQ / CUBE(L[1])) * .5 * b->dt;
	pL[0] -= .5 * pgamma * pL[0] * iQ * b->dt;
	pL[0] += psigma * (rng_uniform(vvrng) - .5);
	pL[1] -= .5 * pgamma * pL[1] * iQ * b->dt;
	pL[1] += psigma * (rng_uniform(vvrng) - .5);
}

/* ------------------------------------------------------------------------- */

void vv_step1(struct beads *restrict b)
{
	int i, d;

	/* first step in all ensembles */
	for (i = first; i < last; i++) {
		for (d = 0; d < 3; d++) {
			assert(fabs(b->f[i][d]) < 50000.);
			assert(b->xv[i][d] >= 0. && b->xv[i][d] < b->l[d]);
			b->xv[i][d + 3] += .5 * b->dt * b->f[i][d];
		}
	}

	do_step1(b);
	debug("b->l: %lg %lg %lg", b->l[0], b->l[1], b->l[2]);
	
	MPI_Allgather(MPI_IN_PLACE, 0, MPI_DOUBLE, b->xv,
		       	6 * b->groupsize, MPI_DOUBLE, comm_row);
	MPI_Allgather(b->xv + first, 6 * b->groupsize, MPI_DOUBLE,
			b->xv + groups_per_row * b->groupsize,
			6 * b->groupsize, MPI_DOUBLE, comm_col);
}

void vv_step2(struct beads *restrict b)
{
	int i, d;

	do_step2(b);
	
	/* last step in all ensembles */
	b->v2max = 0.;
	for (i = first; i < last; i++) {
		b->v2[i] = 0.;
		for (d = 0; d < 3; d++) {
			b->xv[i][d + 3] += .5 * b->dt * b->f[i][d];
			b->v2[i] += SQR(b->xv[i][d + 3]);
		}	
		b->v2max = MAX(b->v2[i], b->v2max);
	}
	
	MPI_Allreduce(MPI_IN_PLACE, &b->v2max, 1, MPI_DOUBLE, MPI_MAX,
								comm_grid);
}

static void vv_free(void)
{
	rng_free(vvrng);
}

void vv_alloc(struct beads *restrict b, cfg_t *cfg)
{
	const char *en = cfg_getstr(cfg_getsec(cfg, "integrator"), "ensemble");

	if (strcasecmp(en, "NVT") == 0) {
		do_step1 = nvt_step1;
		do_step2 = nvt_step2;
		debug("NVT ensemble");
	} else if (strcasecmp(en, "NPT") == 0) {
		V = b->l[0] * b->l[1] * b->l[2];
		iQ = 1. / Q;
		do_step1 = npt_step1;
		do_step2 = npt_step2;
		debug("NPT ensemble: Q=%lg gamma=%lg P0=%lg", Q, pgamma, P0);
	} else if (strcasecmp(en, "NPtT") == 0) {
		A = b->l[TANG1] * b->l[TANG2];
		rt = b->l[TANG1] / b->l[TANG2];
		iQ = 1. / (Q * SQR(b->l[NORMAL]));
		do_step1 = nptt_step1;
		do_step2 = nptt_step2;
		debug("NPtT ensemble: Q=%lg gamma=%lg P0=%lg rt=%0.2f",
						       	Q, pgamma, P0, rt);
	} else if (strcasecmp(en, "NPtT2") == 0) {
		L[0] = b->l[TANG1];
		L[1] = b->l[TANG2];
		iQ = 1. / (Q * SQR(b->l[NORMAL]));
		do_step1 = nptt2_step1;
		do_step2 = nptt2_step2;
		debug("NPtT2 ensemble: Q=%lg gamma=%lg P0=%lg", Q, pgamma, P0);
	} else if (strcasecmp(en, "NPxT") == 0) {
		A = b->l[TANG1] * b->l[TANG2];
		rt = b->l[TANG1] / b->l[TANG2];
		L[0] = b->l[TANG1];
		L[1] = b->l[TANG2];
		iQ = 1. / (Q * SQR(b->l[NORMAL]));
		do_step1 = npxt_step1;
		do_step2 = npxt_step2;
		debug("NPxT ensemble: Q=%lg gamma=%lg P0=%lg", Q, pgamma, P0);
	} else {
		fatal(EINVAL, "The ensemble '%s' is unknown", en);
	}
	psigma = sqrt(6. * pgamma * b->dt);
	first = col_index * b->groupsize; 
	last = 	(col_index + 1) * b->groupsize;

	vvrng = rng_alloc(0);
	if (vvrng == NULL)
		novm("Couldn't init random number generator");

	atexit(vv_free);
}

