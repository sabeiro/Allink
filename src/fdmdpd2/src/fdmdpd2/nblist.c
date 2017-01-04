/*
 * nblist.c - A Fast Pairlist-Construction Algorithm
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

/* $Id: nblist.c 296 2011-06-15 12:17:05Z hoemberg $ */

/* 
 * A Fast Pairlist-Construction Algorithm for Molecular Simulations
 * under Periodic Boundary Conditions 
 * T. N. Heinz and P. H. Huenenberger, J Comput Chem 25: 1474-1486, 2004.
 *
 * Most of the functions and variables have the same name as in the article,
 * so see detailed explanations there. The algorithm has been slightly
 * modified, so that it distinguishes between row- and column-particles.
 */

#include "fdmdpd2.h"

/* ------------------------------------------------------------------------- */

/* global */
static long K = 20;		/* max. number of particles per cell */
static int M;			/* number of cells */
static VEC l;			/* cell lengths */
static VEC l2;			/* square cell lengths */
static VEC il;			/* inverse cell lengths */
static int N[3];		/* number of cells in each direction */
static double dx = 0.5;		/* estimate cell lengths, used for init */
static int *Pc0 = NULL;		/* initial setup for cpa->Pc */

/* supported options in the configuration file */
cfg_opt_t nblist_opts[] = {
	CFG_SIMPLE_FLOAT("dx", &dx),
	CFG_SIMPLE_INT("K", &K),
	CFG_END()
};

/* ------------------------------------------------------------------------- */

/*
 * allocate memory for cell structure
 */
static void realloc_nblist(struct nblist *nbl)
{
	nbl->max_M = M;
	
	/* very crude approximation */
	nbl->Pm = realloc(nbl->Pm, .7 * M * sizeof(*nbl->Pm));
	if (unlikely(nbl->Pm == NULL)) novm("nbl->Pm");
       
	nbl->row.Pc = realloc(nbl->row.Pc, 2 * (M + 1) * sizeof(*nbl->row.Pc));
	if (unlikely(nbl->row.Pc == NULL)) novm("nbl->row.Pc");
	nbl->col.Pc = nbl->row.Pc + M + 1;
	
	nbl->row.C = realloc(nbl->row.C, 2 * (K * M) * sizeof(*nbl->row.C));
	if (unlikely(nbl->row.C == NULL)) novm("nbl->row.C (M=%d)", M);
	nbl->col.C = nbl->row.C + K * M;
}

/*
 * minimum image-function
 */
static int mi(int n, int NN)
{
	return n - NN * (int)((abs(n) + .5 * NN) / NN);
}

static int get_mask(int dm[const static 3], double rc2)
{
	double sum = 0.;
	int dn;

	dn = abs(mi(dm[0], N[0]));
	sum = (dn > 0) ? SQR(dn - 1) * l2[0] : 0.;
	if (sum > rc2) return 0;

	if (dm[0] == 0 || ((dm[1] == N[1] - 1) && (dm[2] == N[2] - 1)))
		dn = abs(mi(dm[1], N[1]));
	else
		dn = MIN(abs(mi(dm[1], N[1])), abs(mi(dm[1] + 1, N[1])));
	sum += (dn > 0) ? SQR(dn - 1) * l2[1] : 0.;
	if (sum > rc2) return 0;
		
	if ((dm[2] == N[2] - 1) || ((dm[0] == 0) && (dm[1] == 0)))
		dn = abs(mi(dm[2], N[2]));
	else
		dn = MIN(abs(mi(dm[2], N[2])), abs(mi(dm[2] + 1, N[2])));
	sum += (dn > 0) ? SQR(dn - 1) * l2[2] : 0.;
	if (sum > rc2) return 0;

	return 1;
}

/*
 * Calculate the mask pointers. The mask array is not stored to save memory.
 */ 
static void set_mask_pointer(struct nblist *nbl)
{
	int dm[3];
	int old_mask = 0;

	nbl->s = 0;
	for (dm[2] = 0; dm[2] < N[2]; dm[2]++) {
		for (dm[1] = 0; dm[1] < N[1]; dm[1]++) {
			for (dm[0] = 0; dm[0] < N[0]; dm[0]++) {

				int mask = get_mask(dm, nbl->rc2);
				int dmi = dm[0] + N[0] *
						(dm[1] + N[1] * dm[2]);

				if (old_mask == mask || dmi == 0) continue;

				if (mask == 1)
					nbl->Pm[2 * nbl->s + 0] = dmi;
				else
					nbl->Pm[2 * (nbl->s++) + 1] = dmi;
			       	
				old_mask = mask;
			}
		}
	}

	if (old_mask == 1)
		nbl->Pm[2 * (nbl->s++) + 1] = M;
}

/*
 * perform final cleanup
 */
static void cleanup(void)
{
	free(Pc0);
}

/*
 * recalculate grid size after the volume has changed
 */
static void recalc_cells(const struct beads *b)
{
	int i;
	int old = M;

	for (i = 0; i < 3; i++) {
		N[i] = (int)(b->l[i] / dx);
		l[i] = b->l[i] / N[i];
		l2[i] = SQR(l[i]);
		il[i] = 1. / l[i];
	}
	M = N[0] * N[1] * N[2];

	if (M > old) {
		debug("NBL: cells increased: %d --> %d", old, M);

		if (unlikely(Pc0 == NULL))
			atexit(cleanup);

		Pc0 = realloc(Pc0, (M + 1) * sizeof(*Pc0));
		if (Pc0 == NULL) novm("Pc0");
	
		for (i = 0; i < M + 1; i++)
			Pc0[i] = i * K;
	}
}

/*
 * recalculate cell structure, this function should be called initially, or
 * after the volume of the simulation box has changed.
 */
static void create_mask_pointer(const struct beads *b, struct nblist *nbl)
{
	recalc_cells(b);

	if (nbl->max_M < M)
		realloc_nblist(nbl);
	
	set_mask_pointer(nbl);
}

/* ------------------------------------------------------------------------- */

static void compact_cells(struct cpa *cpa)
{
	int i, first, last, dif;
#ifdef DEBUG
	int max = 0;
#endif

	for (i = 0, first = 0; i < M + 1; i++) {
		assert(cpa->Pc[i] < (i + 1) * K);
		last = cpa->Pc[i];
		dif = last - i * K;
		cpa->Pc[i] = first;
		if (dif == 0) continue;

		memmove(cpa->C + first, cpa->C + i * K, dif * sizeof(*cpa->C));
		first += dif;
#ifdef DEBUG
		max = MAX(max, dif);
#endif
	}
	debug("compact_cells: max used K=%d", max);
}

static void set_cells(struct cpa *cpa)
{
	int i, m_x, m_y, m_z, m, Pcm;

	memcpy(cpa->Pc, Pc0, (M + 1) * sizeof(*cpa->Pc));

	for (i = 0; i < cpa->n; i++) {
		if (!GET_EXISTS(cpa->pp[i])) continue;
		assert(!(cpa->xv[i][0] == 0. && cpa->xv[i][1] == 0. &&
							cpa->xv[i][2] == 0.));
		m_x = (int)floor(cpa->xv[i][0] * il[0]);
		assert(m_x >= 0 && m_x < N[0]);
		m_y = (int)floor(cpa->xv[i][1] * il[1]);
		assert(m_y >= 0 && m_y < N[1]);
		m_z = (int)floor(cpa->xv[i][2] * il[2]);
		assert(m_z >= 0 && m_z < N[2]);
		m = N[0] * (N[1] * m_z + m_y) + m_x;
		Pcm = cpa->Pc[m]++;
		cpa->C[Pcm] = cpa->pp[i];
	}
}

static void create_cell_pointer(struct nblist *nbl)
{
	set_cells(&nbl->row);
	compact_cells(&nbl->row);

	set_cells(&nbl->col);
	compact_cells(&nbl->col);
}

/* ------------------------------------------------------------------------- */

#ifdef NOTYET
static int calculate(int row, int col)
{
	int q = (rows >= cols) ? rows : cols;
	int j = (rows >= cols) ? col - row_index : row - col_index;
	j = (j + q) % q;

	if (j < empty_lines && j > 0) {
		return -1; /* don't calculate */
	} else if (q - j < empty_lines) {
		return 1; /* calculate */
	}
	return 0; /* let Plimpton decide */
}
#endif

/*
 * do this Plimpton checkerboarding to see whether or not this interaction
 * should be computed by this process.
 */
static int test_checkerboard(PASSPORT prow, PASSPORT pcol)
{
	int row_id = GET_ID(prow);
	int col_id = GET_ID(pcol);

	if (((row_id > col_id) && ((row_id + col_id) & 1) == 0) ||
		((row_id < col_id) && ((row_id + col_id) & 1) == 1) ||
		(row_id == col_id))
		return 0;
	return 1;
}

/*
 * append an interacting pair to the neighbor list. if there
 * is insufficient memory, the whole list is realloc'ed.
 */
static void append_nblist(struct nblist *nbl, PASSPORT prow, PASSPORT pcol)
{
	int c = nbl->count++;

	if (unlikely(nbl->count >= nbl->max)) {
		nbl->max += 0x200000;
		nbl->passport = realloc(nbl->passport,
					2 * nbl->max * sizeof(*nbl->passport));
		if (unlikely(nbl->passport == NULL))
			novm("nbl->passport");
		nbl->drdv = realloc(nbl->drdv, nbl->max * sizeof(*nbl->drdv));
		if (unlikely(nbl->drdv == NULL)) novm("nbl->drdv");
	}

	nbl->passport[2 * c + 0] = prow;
	nbl->passport[2 * c + 1] = pcol; 
}

/*
 * this version of loop_col honors Newton's third law, i.e. the interaction
 * of particle "i" with particle "j" is considered to be the same, and only
 * one of them is included in the pairlist.
 */
static void loop_col_w_n3(struct nblist *nbl, PASSPORT prow, int m1, int m2)
{
	int n2;

	for (n2 = nbl->col.Pc[m1]; n2 < nbl->col.Pc[m2]; n2++) {
		PASSPORT pcol = nbl->col.C[n2];
		if (test_checkerboard(prow, pcol)) {
			append_nblist(nbl, prow, pcol);
		}
	}
}

/*
 * this version of loop_col honors Newton's third law, like the one above.
 * However if two particles of species "C" meet, this interaction is ignored.
 * This functionality is useful for the phantom solvent.
 */
static void loop_col_w_n3_cc(struct nblist *nbl, PASSPORT prow, int m1, int m2)
{
	int n2;
	int is_c = (GET_TYPE(prow) == 2) ? 1 : 0;

	for (n2 = nbl->col.Pc[m1]; n2 < nbl->col.Pc[m2]; n2++) {
		PASSPORT pcol = nbl->col.C[n2];
		if (is_c && GET_TYPE(pcol) == 2) continue;
		if (test_checkerboard(prow, pcol)) {
			append_nblist(nbl, prow, pcol);
		}
	}
}

/*
 * this version of loop_col treats the interactions of the particles
 * "i" and "j" as separate things, and both are included in the pairlist.
 */
static void loop_col_wo_n3(struct nblist *nbl, PASSPORT prow, int m1, int m2)
{
	int n2;

	for (n2 = nbl->col.Pc[m1]; n2 < nbl->col.Pc[m2]; n2++) {
		PASSPORT pcol = nbl->col.C[n2];
		append_nblist(nbl, prow, pcol);
	}
}

/* ------------------------------------------------------------------------- */

/*
 * setup a concurrent neighbor list
 */
void nblist_alloc(const struct beads *b, struct nblist *nbl,
		int row_n, PASSPORT row_pp[const], VEC2 row_xv[const],
		int col_n, PASSPORT col_pp[const], VEC2 col_xv[const],
		double rc, enum nblist_flags flags)
{
	memset(nbl, 0, sizeof(*nbl));
	
	nbl->rc2 = SQR(rc);
	nbl->row.n = row_n;
	nbl->row.pp = row_pp;
	nbl->row.xv = row_xv;
	nbl->col.n = col_n;
	nbl->col.pp = col_pp;
	nbl->col.xv = col_xv;

	if (flags & NBL_NEWTON3) {
		nbl->loop_col = (flags & NBL_IGNORE_CC) ? loop_col_w_n3_cc :
								loop_col_w_n3;
	} else {
		nbl->loop_col = loop_col_wo_n3;
	}

	debug("NBL: dx=%lg K=%ld rc=%lg row_n=%d col_n=%d",
		dx, K, sqrt(nbl->rc2), nbl->row.n, nbl->col.n);

	create_mask_pointer(b, nbl);
}

/*
 * if the current state of the neighbor list is specified by the flags,
 * then the corresponding action is performed here.
 */
void nblist_update(struct beads *b, struct nblist *nbl, enum nblist_flags f)
{
	switch (f) {
	case NBL_VOLUME_CHANGED:
		create_mask_pointer(b, nbl);
	case NBL_INVALID: /* fall-through */
		nblist_rebuild(nbl);
	case NBL_VALID: /* fall-through */
		break;
	default:
		fatal(EINVAL, "nbl: this should never happen");
	}
}

/*
 * rebuild the whole neighbor list 
 */
void nblist_rebuild(struct nblist *nbl)
{
	int m, sl, n1;

	nbl->count = 0;
	create_cell_pointer(nbl);

	for (m = 0; m < M; m++) {
		for (n1 = nbl->row.Pc[m]; n1 < nbl->row.Pc[m + 1]; n1++) {
			PASSPORT prow = nbl->row.C[n1];

			nbl->loop_col(nbl, prow, m, m + 1);

			for (sl = 0; sl < nbl->s; sl++) {
				int m1 = m + nbl->Pm[2 * sl + 0];
				int m2 = m + nbl->Pm[2 * sl + 1];
				nbl->loop_col(nbl, prow,
						MIN(m1, M),
						MIN(m2, M));
				nbl->loop_col(nbl, prow,
						MAX(0, m1 - M),
						MAX(0, m2 - M));
			}
		}
	}
	debug("NBL: rebuilt neighbor list (%d entries)", nbl->count);
}

/*
 * free memory of a neighbor list
 */
void nblist_free(struct nblist *nbl)
{
	free(nbl->Pm);
	free(nbl->row.Pc);
	free(nbl->row.C);
	if (nbl->max > 0) {
		free(nbl->passport);
		free(nbl->drdv);
	}
}

