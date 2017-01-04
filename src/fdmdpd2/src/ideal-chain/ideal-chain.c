/*
 * ideal-chain.c - Consistent intramolecular interaction coefficients
 * (C) Copyright 2008-2010 Martin Hoemberg <mhoembe@gwdg.de>
 *
 * Syntax: ./ic [N] [kb] [l] [Re] [ks]
 *
 * Output: A value of "ks" which is consistent with N, kb, l and Re
 *
 * Method: bead displacements and reptation "slithering snake" moves
 */


#define _GNU_SOURCE

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define SQR(x) ((x)*(x))
typedef double VEC[3];

/* ------------------------------------------------------------------------- */

static int N;				/* beads per chain */
static double kb;			/* bond-angle stiffness constant */
static double l0;			/* equilibrium bond length */
static double Re;			/* rms end-to-end vector */
static double ks;			/* spring constant */
static FILE *VMDFH = NULL;		/* VMD(tm) configuration dump file */
static VEC *x = NULL;			/* array for bead coordinates */
static VEC *xb = NULL;			/* backup coordinates */
static double eold = 1e10;		/* current energy */

/* ------------------------------------------------------------------------- */

static int vmd_fopen(const char *fn)
{
	VMDFH = fopen(fn, "w");
	if (VMDFH == NULL) {
		fprintf(stderr, "Couldn't open '%s' for writing", fn);
		return 1;
	}
	
	fprintf(VMDFH, "atom 0:%d radius 1.0 name A\n", N - 1);
	fprintf(VMDFH, "bond 0::%d\n", N - 1);

	return 0;
}

static void vmd_fclose(void)
{
	fclose(VMDFH);
}

static int vmd_append(void)
{
	int i;
	VEC com = {0., 0., 0.};

	for (i = 0; i < N; i++) {
		com[0] += x[i][0];
		com[1] += x[i][1];
		com[2] += x[i][2];
	}
	com[0] /= N;
	com[1] /= N;
	com[2] /= N;
	
	fprintf(VMDFH, "\ntimestep\npbc 100.0 100.0 100.0\n");
	fprintf(VMDFH, "# kb=%lg ks=%lg l=%lg\n", kb, ks, l0);
	for (i = 0; i < N; i++)
		fprintf(VMDFH, "%lg %lg %lg\n",
				x[i][0] - com[0],
				x[i][1] - com[1],
				x[i][2] - com[2]);
	fflush(VMDFH);
	return 0;
}

/* ------------------------------------------------------------------------- */

static double ia_harmonic_spring(VEC xy[const restrict], 
				const int i, const int j)
{
	VEC dr;
	double dr2;
	double energy;
	double r;

	dr[0] = xy[j][0] - xy[i][0];
	dr[1] = xy[j][1] - xy[i][1];
	dr[2] = xy[j][2] - xy[i][2]; 
	dr2 = SQR(dr[0]) + SQR(dr[1]) + SQR(dr[2]);
	r = sqrt(dr2);

	energy = .5 * ks * SQR(r - l0);
	return energy;
}

/* bond angle potentials */
static double ia_bond_angle(VEC xy[const restrict],
			const int i, const int j, const int k)
{
	VEC d0, d1;
	double n0, n1;
	double ct;
	double energy;
			
	d0[0] = xy[j][0] - xy[i][0];
	d0[1] = xy[j][1] - xy[i][1];
	d0[2] = xy[j][2] - xy[i][2];
	n0 = sqrt(SQR(d0[0]) + SQR(d0[1]) + SQR(d0[2]));
	d0[0] /= n0;
	d0[1] /= n0;
	d0[2] /= n0;

	d1[0] = xy[k][0] - xy[j][0]; 
	d1[1] = xy[k][1] - xy[j][1]; 
	d1[2] = xy[k][2] - xy[j][2];
	n1 = sqrt(SQR(d1[0]) + SQR(d1[1]) + SQR(d1[2]));
	d1[0] /= n1;
	d1[1] /= n1;
	d1[2] /= n1;

	ct = d0[0] * d1[0] + d0[1] * d1[1] + d0[2] * d1[2];
	energy = kb * (1. - ct);

	return energy;
} 

static double energy(VEC coords[const restrict])
{
	double e = 0.;
	int i;

	for (i = 0; i < N - 1; i++)
		e += ia_harmonic_spring(coords, i, i + 1);
	for (i = 1; i < N - 1; i++)
		e += ia_bond_angle(coords, i - 1, i, i + 1);
	return e;
}

/* ------------------------------------------------------------------------- */

static int do_displacement_move(void)
{
	int i = N * drand48();
	double r, phi, costheta, sintheta;

	r = drand48() * 0.20;
	phi = 2. * M_PI * drand48();
	costheta = 2. * drand48() - 1.;
	sintheta = sqrt(1. - SQR(costheta));

	memcpy(xb, &x[i], sizeof(*x));

	x[i][0] += r * cos(phi) * sintheta;
	x[i][1] += r * sin(phi) * sintheta;
	x[i][2] += r * costheta;
	
	double enew = energy(x);

	if (enew < eold) goto accept;
	if (drand48() < exp(eold - enew)) goto accept;

	memcpy(&x[i], xb, sizeof(*x));

	return 0;
accept:
	eold = enew;
	return 1;
}

static double propose_r(void)
{
	double r, y, p;
	static double rmax = 3.;
	do {
		r = rmax * drand48();
		y = SQR(rmax) * drand48();
		p = SQR(r) * exp(-.5 * ks * SQR(r - l0));
	} while (y > p);
//	printf("proposed r=%lg\n", r);
	return r;
}

static int do_reptation_move(void)
{
	double enew, elastbond, eoldbond, rold;
	double r, phi, costheta, sintheta;
		
	r = propose_r();
	phi = 2. * M_PI * drand48();
	costheta = 2. * drand48() - 1.;
	sintheta = sqrt(1. - SQR(costheta));
	elastbond = .5 * ks * SQR(r - l0);
	
	if (drand48() < 0.5) { /* delete first monomer, put a new tail */
		memcpy(xb, x + 1, (N - 1) * sizeof(*x));
		rold = sqrt(SQR(x[0][0] - x[1][0]) + SQR(x[0][1] - x[1][1]) +
						SQR(x[0][2] - x[1][2]));
		xb[N - 1][0] = xb[N - 2][0] + r * cos(phi) * sintheta;
		xb[N - 1][1] = xb[N - 2][1] + r * sin(phi) * sintheta;
		xb[N - 1][2] = xb[N - 2][2] + r * costheta;
		
	} else { /* delete tail, put new head */
		memcpy(xb + 1, x, (N - 1) * sizeof(*x));
		rold = sqrt(SQR(x[N - 1][0] - x[N - 2][0]) +
					SQR(x[N - 1][1] - x[N - 2][1]) +
					SQR(x[N - 1][2] - x[N - 2][2]));
		xb[0][0] = xb[1][0] + r * cos(phi) * sintheta;
		xb[0][1] = xb[1][1] + r * sin(phi) * sintheta;
		xb[0][2] = xb[1][2] + r * costheta;
	}

	eoldbond = .5 * ks * SQR(rold - l0);
	enew = energy(xb);

	double q = elastbond - eoldbond + eold - enew;
	if (q > 0.) goto accept;
	if (drand48() < exp(q)) goto accept;

	return 0;
accept:
	memcpy(x, xb, N * sizeof(*x));
	eold = enew;
	return 1;
}

static double calculate_re(void)
{
	double re, re_sum = 0., re2_sum = 0.;
	int acc, reptation_acc = 0, displ_acc = 0;
	int reptation_rej = 0, displ_rej = 0;
	int re_count = 0;
	double blocksum = 0., blocksum2 = 0.;
	int blocksum_count = 0;
	double mw, sd;
	
	eold = energy(x);
	re = SQR(x[0][0] - x[N - 1][0]);
	re += SQR(x[0][1] - x[N - 1][1]);
	re += SQR(x[0][2] - x[N - 1][2]);
	re = sqrt(re);
			
	do {
		if (drand48() < .6) {
			acc = do_reptation_move();
			if (acc == 1)
				reptation_acc++;
			else
				reptation_rej++;
			
		} else {
			acc = do_displacement_move();
			if (acc == 1)
				displ_acc++;
			else
				displ_rej++;
		}

		if (acc == 1) {
			re = SQR(x[0][0] - x[N - 1][0]);
			re += SQR(x[0][1] - x[N - 1][1]);
			re += SQR(x[0][2] - x[N - 1][2]);
			re = sqrt(re);
		}

		re_sum += re;
		re2_sum += SQR(re);
		re_count++;

		if (re_count == 50000) {
			mw = re_sum / re_count;
			//sd = sqrt(re2_sum / re_count - SQR(mw));
			//printf("mw=%lg sd=%lg\n", mw, sd);

			blocksum += mw;
			blocksum2 += SQR(mw);
			blocksum_count++;

			re_count = 0;
			re_sum = 0.;
			re2_sum = 0.;

			if (blocksum_count > 20) {
				mw = blocksum / blocksum_count;
				sd = sqrt((blocksum2 / blocksum_count - SQR(mw)) / (blocksum_count - 1.));
				if (sd < 0.01) break;
			}
			vmd_append();
		}
	} while (1);

	mw = blocksum / blocksum_count;
	sd = sqrt((blocksum2 / blocksum_count - SQR(mw)) / (blocksum_count - 1.));

	printf("%lg %lg %lg %lg %lg\n", ks, mw, sd, 
		reptation_acc / (1. * reptation_acc + reptation_rej) * 100.,
		displ_acc / (1. * displ_acc + displ_rej) * 100.);
	
	return mw;
}

/* put in a (not very meaningful) initial conformation */
static void initial_coords(void)
{
	int i;

	for (i = 0; i < N; i++) {
		x[i][0] = i * 1.;
		x[i][1] = 0.;
		x[i][2] = 0.;
	}
}

int main(int argc, char **argv)
{
	double re;

	if (argc < 6) {
		fprintf(stderr, "Syntax: %s [N] [kb] [l] [Re] [ks]\n\n",
				       				argv[0]);
		return EXIT_FAILURE;
	}

	N = atoi(argv[1]);
	kb = atof(argv[2]);
	l0 = atof(argv[3]);
	Re = atof(argv[4]);
	ks = atof(argv[5]);
	printf("# N=%d kb=%lg ks=%lg l=%lg Re=%lg\n", N, kb, ks, l0, Re);

	srand48(0L);
	if (vmd_fopen("ic.vtf")) return EXIT_FAILURE;

	x = calloc(2 * N, sizeof(*x));
	xb = x + N;
	if (x == NULL) {
		fprintf(stderr, "novm: x / xb\n");
		return EXIT_FAILURE;
	}

	initial_coords();

	printf("# ks Re d(Re) acc(Rep) acc(Displ)\n");
	do {
		fprintf(stderr, "Using ks=%lg\n", ks);
		re = calculate_re();
		
		if ((re - Re) / Re > 0.01) {
			ks *= 1.10;
		} else if ((re - Re) / Re < -0.01) {
			ks *= 0.95;
		}

	} while (fabs((re - Re) / Re) > 0.01);

	fprintf(stderr, "Final result ks=%lg\n\n", ks);
	vmd_fclose();
	free(x);
	return EXIT_SUCCESS;
}

