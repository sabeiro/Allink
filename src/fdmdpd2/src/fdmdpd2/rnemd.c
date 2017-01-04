/*
 * rnemd.c - Reverse NEMD method to apply shear to the system
 * Copyright (C) 2010 Martin Hoemberg <mhoembe@gwdg.de>
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

/* $Id: rnemd.c 296 2011-06-15 12:17:05Z hoemberg $ */

/* 
 * References:
 *   F. Müller-Plathe, PRE 59, 4894 (1999)
 *   T. J. Müller and F. Müller-Plathe, ChemPhysChem 10, 2305 (2009)
 *   C. M. Tenney and E. J. Maginn, JCP 132, 014103 (2010)
 */

#include "fdmdpd2.h"

struct rnemd_list
{
	double v;		/* particle velocity */
	int pos;		/* particle passport */
	int rank;		/* rank of this process */
};

/* --------------------------- common variables ---------------------------- */
static int inuse = 0;		/* is this module used? */
static double P_sum = 0.;	/* total sum of exchanged momentum */
static double t0 = 0.;		/* time interval */
static long N = 10;		/* # slabs */
static long dg = 1;		/* dimension of velocity gradient */
static long dd = 2;		/* dimension of velocity direction */
static long sw = 1;		/* velocity swaps for sweep */
static MPI_Datatype rnemd_type;
static int rank;
static int size;
static struct rnemd_list *list = NULL;
static int max = 0;

/* ------------------------------------------------------------------------- */

cfg_opt_t rnemd_opts[] = {
	CFG_SIMPLE_INT("slabs", &N),
	CFG_SIMPLE_INT("gradient", &dg),
	CFG_SIMPLE_INT("velocity", &dd),
	CFG_SIMPLE_INT("swaps", &sw),
	CFG_END()
};

/* ------------------------------------------------------------------------- */

void rnemd_init(struct beads *b)
{
	inuse = 1;

	if (N % 2 == 1) 
		fatal(EINVAL, "rnemd: N must be even");
	if (dd < 0 || dd >= 3 || dg < 0 || dg >= 3)
		fatal(EINVAL, "rnemd: gradient / velocity invalid");

	printf("rnemd: slabs=%ld swaps=%ld gradient=%ld velocity=%ld\n",
								N, sw, dg, dd);
	assert(N && sw);
	int blocklens[2] = {1, 2};
	struct rnemd_list q;
	ptrdiff_t indices[2] = {(ptrdiff_t)&q.v - (ptrdiff_t)&q,
		(ptrdiff_t)&q.pos - (ptrdiff_t)&q};
	MPI_Datatype old_types[2] = {MPI_DOUBLE, MPI_INT};
	MPI_Type_struct(ARRAY_SIZE(blocklens), blocklens, indices, old_types, &rnemd_type);
	MPI_Type_commit(&rnemd_type);
	
	MPI_Comm_rank(comm_grid, &rank);
	MPI_Comm_size(comm_grid, &size);
				
	max = MAX(0x10000, 2 * sw * size);
	list = calloc(max, sizeof(*list));
	if (list == NULL) novm("rnemd: list");

	t0 = b->time;
}

void rnemd_free(void)
{
	if (inuse) {
		MPI_Type_free(&rnemd_type);
		if (list) free(list);
	}
}

static int rnemd_cmp(const void *__a, const void *__b)
{
	struct rnemd_list *a = (struct rnemd_list *)__a;
	struct rnemd_list *b = (struct rnemd_list *)__b;
	return (int)floor((a->v - b->v) * 10000.);
}

/* convention: 
 * the most positive velocities are withdrawn from slab #0 and added to #N/2.
 */
void rnemd_sweep(struct beads *b)
{
	if (!inuse) return;

	double idx = N / b->l[dg];
	int i, c = 0;
	int o = col_index * b->groupsize;
	int o2 = (row_index + groups_per_row) * b->groupsize;
	struct rnemd_list list2[2 * sw];

	/* create hitlist of this process */
	for (i = 0; i < b->groupsize; i++) {
		int k = (int)floor(b->xv[o + i][dg] * idx);
		double v = b->xv[o + i][dd + 3];

		if ((k == 0 && v > 0.) || (k == N / 2 && v < 0.)) {
			if (c >= max) {
				max += 0x10000;
				list = realloc(list, max * sizeof(*list));
				if (list == NULL) novm("rnemd: list");
			}
			list[c].v = v;
			list[c].pos = i;
			list[c++].rank = rank;
		}
	}
	debug("rnemd: %d beads in both slabs", c);
	qsort(list, c, sizeof(*list), rnemd_cmp);

	/* join hitlists of all processes */
	memset(list2, 0, 2 * sw * sizeof(*list2));
	for (i = 0; i < MIN(sw, c / 2); i++) {
		list2[i] = list[i];
		list2[i + sw] = list[c - i - 1];
	}
	MPI_Allgather(list2, 2 * sw, rnemd_type, list, 2 * sw, rnemd_type, comm_grid);
	qsort(list, 2 * sw * size, sizeof(*list), rnemd_cmp);

	/* now replace the velocities */
	for (i = 0; i < sw; i++) {
		int ie = 2 * size * sw - 1 - i;

		if (rank == list[i].rank) {
			assert(b->xv[o + list[i].pos][dd + 3] * list[ie].v < 0.);
			b->xv[o + list[i].pos][dd + 3] = list[ie].v;
			b->xv[o2 + list[i].pos][dd + 3] = list[ie].v;
			P_sum += list[ie].v;
		}
		if (rank == list[ie].rank) {
			assert(b->xv[o + list[ie].pos][dd + 3] * list[i].v < 0.);
			b->xv[o + list[ie].pos][dd + 3] = list[i].v;
			b->xv[o2 + list[ie].pos][dd + 3] = list[i].v;
			P_sum -= list[i].v;
		}
	}

	bcast_beads(b);
}

void rnemd_print(struct beads *b)
{
	if (!inuse) return;
	double sum = 0.;

	MPI_Reduce(&P_sum, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, comm_grid);

	if (ismaster)
		printf("rnemd: P_sum=%lg t=%lg\n", sum, b->time - t0);
}

