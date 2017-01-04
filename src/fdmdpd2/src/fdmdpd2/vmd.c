/*
 * vmd.c - Write current system snapshot as VTF file format to disk
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

/* $Id: vmd.c 296 2011-06-15 12:17:05Z hoemberg $ */

/*
 * The full system configuration is frequently appended to a VTF-file readable
 * by VMD. The file format is borrowed from Espresso. A full description can
 * be found at their website, http://www.espresso.mpg.de/
 *
 * The VTF header contains the particle properties, as well as the
 * connectivity. Unfortunately only linear molecules are supported at the
 * moment. In the part with the coordinates every molecule is shifted, so that
 * its center-of-mass is in the center box.
 */

#include "fdmdpd2.h"

/* ------------------------------------------------------------------------- */

/* file handle for writing the VTF file */
static FILE *FH = NULL;

/* temp array for coordinates */
static VEC *coords = NULL;

/* particle names for the different beads */
static const char names[3][20] = {"A", "B", "C"};

/* ------------------------------------------------------------------------- */

/* TODO: extend this also for branched molecules */
static void vmd_print_connectivity(const struct beads *restrict b,
			PASSPORT pp[restrict])
{
	int i, j, k;
	int idx = 0;
	int chain = 0;

	for (i = 0; i < b->blocks; i++) {
		for (j = 0; j < b->n[i]; j++, chain++) {
			int idx0 = idx;
			for (k = 0; k < b->N[i]; k++, idx++) {
				int t = GET_TYPE(pp[idx]);
				if (!GET_EXISTS(pp[idx])) continue;
				fprintf(FH, "a %llu r %lg n %s resid %d "
				"res %s\n", GET_ID(pp[idx]), .8, names[t],
			       	chain, b->resname[i]);
			}
			fprintf(FH, "b %d::%d\n", idx0, idx - 1);
		}
	}
}

void vmd_alloc(struct beads *restrict b, const char *fn)
{
	int i;
	PASSPORT *passport = b->passport + col_index * b->groupsize;
	PASSPORT *passport_out;

	passport_out = calloc(b->nN, sizeof(*passport_out));
	if (passport_out == NULL) novm("passport_out");

	coords = calloc(b->nN, sizeof(*coords));
	if (coords == NULL) novm("coords");

	for (i = 0; i < b->groupsize; i++) {
		int id = GET_ID(passport[i]);
		passport_out[id] = passport[i];
	}

	if (ismaster) {
		MPI_Reduce(MPI_IN_PLACE, passport_out, b->nN, MPI_LONG_LONG,
						MPI_BOR, 0, comm_grid);
		FH = fopen(fn, "w");
		if (FH == NULL)
			fatal(EIO, "Couldn't open %s for writing", fn);
		vmd_print_connectivity(b, passport_out);
	} else {
		MPI_Reduce(passport_out, NULL, b->nN, MPI_LONG_LONG,
			       			MPI_BOR, 0, comm_grid);
	}

	free(passport_out);
}

void vmd_free(void)
{
	if (ismaster)
		fclose(FH);
	if (coords != NULL)
		free(coords);
}

void vmd_append(struct beads *restrict b)
{
	VEC2 *xv = b->xv + col_index * b->groupsize;
	VEC *n = b->nx + col_index * b->groupsize;
	PASSPORT *passport = b->passport + col_index * b->groupsize;
	int i, j, d, idx;

	memset(coords, 0, b->nN * sizeof(*coords));
	for (i = 0, idx = 0; i < b->local_n; i++) {
		int block = b->local_n_b[i];
		int N = b->N[block]; /* current chain length */
		VEC xv_tmp[N];
		VEC xv_com = {0., 0., 0.};
		int idx0 = idx;

		for (j = 0; j < N; j++, idx++) {
			for (d = 0; d < 3; d++) {
				xv_tmp[j][d] = xv[idx][d] + n[idx][d] * b->l[d];
				xv_com[d] += xv_tmp[j][d];
			}
		}
		for (d = 0; d < 3; d++)
			xv_com[d] = b->l[d] * floor(xv_com[d] / (N * b->l[d]));

		for (j = 0, idx = idx0; j < N; j++, idx++) {
			int id = GET_ID(passport[idx]);
			for (d = 0; d < 3; d++)
				coords[id][d] = xv_tmp[j][d] - xv_com[d];
		}
	}

	if (ismaster) {
		MPI_Reduce(MPI_IN_PLACE, coords, ARRAY_SIZE(*coords) * b->nN,
					MPI_DOUBLE, MPI_SUM, 0, comm_grid);
		fprintf(FH, "\ntimestep\npbc %lg %lg %lg\n",
					b->l[0], b->l[1], b->l[2]);
		for (j = 0; j < b->nN; j++)
			fprintf(FH, "%lg %lg %lg\n",
				coords[j][0], coords[j][1], coords[j][2]);
	} else {
		MPI_Reduce(coords, NULL, ARRAY_SIZE(*coords) * b->nN,
					MPI_DOUBLE, MPI_SUM, 0, comm_grid);
	}
}

