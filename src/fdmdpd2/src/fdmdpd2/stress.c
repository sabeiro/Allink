/*
 * stress.c - Calculate stress tensor and its autocorrelation
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

/* $Id: stress.c 296 2011-06-15 12:17:05Z hoemberg $ */

#include "fdmdpd2.h"


/* --------------------------- common variables ---------------------------- */
static FILE *STRESS = NULL;		/* output file handle */
static char *buf = NULL;		/* output buffer */
static int inuse = 0;			/* is this module used? */

/* ------------------------------------------------------------------------- */

/*
 * The virial theorem is used to calculate the diagonal components of the
 * collective stress tensor.
 */
void get_pressure(const struct beads *restrict b, TENSOR *p)
{
	VEC2 *xv = b->xv + col_index * b->groupsize;
	double V = b->l[0] * b->l[1] * b->l[2];
	int i;

	for (i = 0; i < ARRAY_SIZE(*p); i++)
		(*p)[i] = b->virial[i];

	for (i = 0; i < b->groupsize; i++) {
		(*p)[0] += xv[i][3] * xv[i][3];
		(*p)[1] += xv[i][4] * xv[i][4];
		(*p)[2] += xv[i][5] * xv[i][5];
		/* The ensemble average of the off-diagonal elements for the
		 * ideal gas vanish. For decreasing the fluctuations we do
		 * it without.*/
	}
	for (i = 0; i < ARRAY_SIZE(*p); i++)
		(*p)[i] /= V;
	
	MPI_Allreduce(MPI_IN_PLACE, *p, ARRAY_SIZE(*p), MPI_DOUBLE, MPI_SUM, comm_grid);
	debug("pressure: %lg %lg %lg %lg %lg %lg",
			(*p)[0], (*p)[1], (*p)[2], (*p)[3], (*p)[4], (*p)[5]);
}

/* ------------------------------------------------------------------------- */

void stress_autocorr_alloc(const char *fn, size_t bufsize)
{
	inuse = 1;

	if (ismaster) {
		STRESS = fopen(fn, "a");
		buf = malloc(bufsize);
		if (buf == NULL) novm("stress: buf");
		setbuffer(STRESS, buf, bufsize);
		fprintf(STRESS, "# time step P_xx P_yy P_zz P_xy P_xz P_yz\n");
	}
}

void stress_autocorr_free(void)
{
	if (!inuse) return;

	if (ismaster) {
		fclose(STRESS);
		free(buf);
	}
}

void stress_autocorr_measure(struct beads *b)
{
	TENSOR p;

	if (!inuse) return;
	get_pressure(b, &p);

	if (ismaster)
		fprintf(STRESS, "%lg %lu %lg %lg %lg %lg %lg %lg\n", b->time,
			       	b->step, p[0], p[1], p[2], p[3], p[4], p[5]);
}

