/*
 * mcom.c - MDPD Center-Of-Mass (MCOM) force field
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

/* $Id: mcom.c 305 2011-07-05 14:27:10Z hoemberg $ */

#include "fdmdpd2.h"
#include "rand.h"
#include <gsl/gsl_bspline.h>

#define STEP 16

/* ------------------------------------------------------------------------- */

/* global variables for the mcom_print */
static FILE *FH;			/* file handle (for the master) ... */
static char *buf;			/* ... and its dedicated I/O buffer */
static VEC2 *xvbuf;			/* buffer for positions, velocities, */
static VEC *fbuf;			/* forces */
static VEC *rbuf;			/* orientations */

void mcom_print_init(struct beads *b, const char *fn, size_t bufsize)
{
	int n;

	MPI_Allreduce(&b->local_n, &n, 1, MPI_INT, MPI_SUM, comm_grid);
	
	xvbuf = calloc(n, sizeof(*xvbuf));
	if (xvbuf == NULL) novm("mcom_print_init: xvbuf");
	fbuf = calloc(n, sizeof(*fbuf));
	if (fbuf == NULL) novm("mcom_print_init: fbuf");
	rbuf = calloc(n, sizeof(*rbuf));
	if (rbuf == NULL) novm("mcom_print_init: rbuf");

	if (ismaster) {
		FH = fopen(fn, "w");
		if (FH == NULL)
			fatal(EIO, "Couldn't open %s for writing", fn);

		buf = malloc(bufsize * sizeof(*buf));
		if (buf == NULL) novm("mcom_print_init: I/O buf");
		setbuffer(FH, buf, bufsize);
	}
}

void mcom_print_free(void)
{
	free(xvbuf);
	free(rbuf);
	free(fbuf);
	
	if (ismaster) {
		fclose(FH);
		free(buf);
	}
}

/*
 * calculates the resulting forces on the local molecules' center-of-mass.
 * Both A and B contribute to the forces.
 */
void mcom_calc_resulting_forces(struct beads *b, VEC *b_f)
{
	int i, j, d, idx;
	VEC *f = b->f + col_index * b->groupsize;
		
	memset(b_f, 0, b->local_n * sizeof(*b_f));
	for (i = 0, idx = 0; i < b->local_n; i++) {
		for (j = 0; j < b->N[b->local_n_b[i]]; j++, idx++) {
			for (d = 0; d < 3; d++)
				b_f[i][d] += f[idx][d];
		}
	}
}

/*
 * calculates the end-to-end vector of the local molecules
 */
void mcom_calc_e2e_vector(struct beads *b, VEC *b_r)
{
	int i, j, d, idx;
	VEC *nx = b->nx + col_index * b->groupsize;
	VEC2 *xv = b->xv + col_index * b->groupsize;

	for (i = 0, idx = 0; i < b->local_n; i++) {
		VEC2 rf = {0., 0., 0.};

		for (j = 0; j < b->N[b->local_n_b[i]]; j++, idx++) {
			if (j == 0) {
				for (d = 0; d < 3; d++) {
					rf[d] = -1. * (xv[idx][d] + nx[idx][d] * b->l[d]);
				}
			} else if (j == b->N[b->local_n_b[i]] - 1) {
				for (d = 0; d < 3; d++) {
					rf[d] += xv[idx][d] + nx[idx][d] * b->l[d];
				}
			}
		}
		for (d = 0; d < 3; d++)
			b_r[i][d] = rf[d];
	}
}

/*
 * calculates the local molecules' center-of-mass coordinates and velocities.
 * Both are calculated only from the A beads.
 */
void mcom_calc_com_xv(struct beads *b, VEC2 *b_xv)
{
	int i, j, d, idx, count;
	VEC2 *xv = b->xv + col_index * b->groupsize;
	PASSPORT *p = b->passport + col_index * b->groupsize;
	VEC *nx = b->nx + col_index * b->groupsize;
	
	for (i = 0, idx = 0; i < b->local_n; i++) {
		VEC2 cm = {0., 0., 0., 0., 0., 0.};

		for (j = 0, count = 0; j < b->N[b->local_n_b[i]]; j++, idx++) {
			if (GET_TYPE(p[idx]) == 0) {
				for (d = 0; d < 3; d++) {
					cm[d] += xv[idx][d] + nx[idx][d] * b->l[d];
					cm[d + 3] += xv[idx][d + 3];
				}
				count++;
			}
		}

		for (d = 0; d < ARRAY_SIZE(cm); d++) {
			cm[d] /= count;
			if (d < 3) {
				while (cm[d] > b->l[d]) cm[d] -= b->l[d];
				while (cm[d] < 0.) cm[d] += b->l[d];
			}
			b_xv[i][d] = cm[d];
		}
	}
}

static void get_n_base(struct beads *b, int *n, int *n_base)
{
	int i, rank, size, *n_all;
	MPI_Comm_rank(comm_grid, &rank);
	MPI_Comm_size(comm_grid, &size);

	n_all = alloca(size * sizeof(*n_all));
	n_all[rank] = b->local_n;

	MPI_Allgather(MPI_IN_PLACE, 0, MPI_INT, n_all, 1, MPI_INT, comm_grid);
	for (i = 0, *n_base = 0, *n = 0; i < size; i++) {
		(*n_base) += (i < rank) ? n_all[i] : 0;
		(*n) += n_all[i];
	}
}

void mcom_print(struct beads *b)
{
	int i, n, n_base;

	get_n_base(b, &n, &n_base);
	mcom_calc_resulting_forces(b, fbuf + n_base);
	mcom_calc_e2e_vector(b, rbuf + n_base);
	mcom_calc_com_xv(b, xvbuf + n_base);

	MPI_Reduce(ismaster ? MPI_IN_PLACE : fbuf, ismaster ? fbuf : NULL,
		ARRAY_SIZE(*fbuf) * n, MPI_DOUBLE, MPI_SUM, 0, comm_grid);
	MPI_Reduce(ismaster ? MPI_IN_PLACE : rbuf, ismaster ? rbuf : NULL,
		ARRAY_SIZE(*rbuf) * n, MPI_DOUBLE, MPI_SUM, 0, comm_grid);
	MPI_Reduce(ismaster ? MPI_IN_PLACE : xvbuf, ismaster ? xvbuf : NULL,
		ARRAY_SIZE(*xvbuf) * n, MPI_DOUBLE, MPI_SUM, 0, comm_grid);

	if (ismaster) {
		fprintf(FH, "# L=%lg %lg %lg N=%d t=%lg\n",
					b->l[0], b->l[1], b->l[2], n, b->time);
		for (i = 0; i < n; i++) {
			fprintf(FH, "%lg %lg %lg %lg %lg %lg ",
			       	xvbuf[i][0], xvbuf[i][1], xvbuf[i][2],
			       	xvbuf[i][3], xvbuf[i][4], xvbuf[i][5]);
			fprintf(FH, "%lg %lg %lg %lg %lg %lg\n",
				rbuf[i][0], rbuf[i][1], rbuf[i][2],
				fbuf[i][0], fbuf[i][1], fbuf[i][2]);
		}
	}
}

/* ------------------------------------------------------------------------- */

/*
 * struct mcom - reentrant coefficients for the COM interactions
 */
struct mcom
{
	double *coeffs;				/* interaction coefficients */
	gsl_matrix *dB;				/* B-spline derivatives */
	gsl_bspline_workspace *bw;		/* GSL B-splines */
	gsl_bspline_deriv_workspace *bdw;	/* GSL B-splines derivs. */
	double gamma[SQR(TYPE_MAX)];	/* DPD friction constant */
	double zeta[SQR(TYPE_MAX)];	/* and noise coefficient */
	double gamma_p[SQR(TYPE_MAX)];	/* Transversal DPD friction constant */
	double zeta_p[SQR(TYPE_MAX)];	/* and noise coefficient */
	double rc2;				/* squared cutoff */
	double rc;				/* cutoff */
	double dpdrc;				/* cutoff for thermostat */
	double m;				/* mass of the particles */
	int ncoeffs;				/* number of coefficients */
	int diml;				/* dimensionality: 2 or 3 */
	char fn[NAME_MAX];			/* filename with coefficients*/
	struct nblist nbl;			/* neighbor list */
};

/* standard DPD weighting function */
static double __attribute__((const)) dpd_wf_std(double r, double rc)
{
	return r < rc ? 1. - r / rc : 0.;
}

#define TDPD_WF(r,rc) dpd_wf_std(r,rc)
#define DPD_WF(r,rc) dpd_wf_std(r,rc)

static double mcom_forces2d(struct beads *b, struct mcom *mcom, struct nblist *nbl)
{
	int i, d, k, l;
	const VEC lh = {.5 * b->l[0], .5 * b->l[1], .5 * b->l[2]};
	double energy = 0.;
	size_t istart, iend;

	for (i = 0; i < nbl->count; i += STEP) {
		int j, n, idx[STEP], difleaflet[STEP];
		int p_row[STEP], p_col[STEP], t_row[STEP], t_col[STEP];
		double dr2[STEP], theta[2 * STEP];
		VEC2 drdv[STEP];

		for (j = 0, n = 0; (j < STEP) && (i + j < nbl->count); j++) {
			p_row[j] = GET_POS(nbl->passport[2 * (i + j) + 0]);
			p_col[j] = GET_POS(nbl->passport[2 * (i + j) + 1]);
			t_row[j] = GET_TYPE(nbl->passport[2 * (i + j) + 0]);
			t_col[j] = GET_TYPE(nbl->passport[2 * (i + j) + 1]);
			difleaflet[j] = (t_row[j] == t_col[j]) ? 0 : 1;
			
			dr2[j] = 0.;
			
			for (d = 1; d < 3; d++) {
				drdv[j][d] = b->xv[p_col[j]][d] - b->xv[p_row[j]][d];
				drdv[j][d] += (drdv[j][d] < -lh[d]) ? b->l[d] : 0.;
				drdv[j][d] -= (drdv[j][d] >  lh[d]) ? b->l[d] : 0.;
				dr2[j] += SQR(drdv[j][d]);
			}
			
			if (dr2[j] < mcom->rc2) {
				idx[n++] = j;
				SET_OK(nbl->passport[2 * (i + j)]);
				for (d = 4; d < 6; d++)
					drdv[j][d] = b->xv[p_col[j]][d] -
							b->xv[p_row[j]][d];
			} else {
				UNSET_OK(nbl->passport[2 * (i + j)]);
			}
		}
		
		rng_uniform_vector(rng, 2 * n, theta);

		for (k = 0; k < n; k++) {
			j = idx[k];
			double r = sqrt(dr2[j]);
			double force = 0., omega_r, tdiss, trand, sp, sp2;

			/* conservative force */
			gsl_bspline_deriv_eval_nonzero(r, 1, mcom->dB, &istart,
						&iend, mcom->bw, mcom->bdw);
			int o = istart + mcom->ncoeffs * difleaflet[j];
			for (l = 0; l < mcom->bw->k; l++) {
				force -= mcom->coeffs[o + l] * gsl_matrix_get(mcom->dB, l, 1);
				energy += mcom->coeffs[o + l] * gsl_matrix_get(mcom->dB, l, 0);
			}
			force /= r;
			
			/* Conservative forces are added to the virial */
			b->virial[1] += SQR(drdv[j][1]) * force / mcom->m;
			b->virial[2] += SQR(drdv[j][2]) * force / mcom->m;
			b->virial[5] += drdv[j][1] * drdv[j][2] * force / mcom->m;
			
			theta[2 * k + 0] -= .5;
			theta[2 * k + 1] -= .5;

			sp = drdv[j][1] * drdv[j][4] + drdv[j][2] * drdv[j][5];
			sp2 = (drdv[j][1] * theta[2 * k + 0] +
					drdv[j][2] * theta[2 * k + 1]) / r;
			
			int tmp = t_row[j] * TYPE_MAX + t_col[j];

			/* stochastic forces from standard DPD */
			omega_r = DPD_WF(r, mcom->dpdrc) / r;
			force += (b->zeta[tmp] * sp2 -
			       	b->gamma[tmp] * sp * omega_r) * omega_r;

			/* stochastic forces from transversal DPD */
			omega_r = TDPD_WF(r, mcom->dpdrc) / r;
			force -= (b->zeta_p[tmp] * sp2 -
			       	b->gamma_p[tmp] * sp * omega_r) * omega_r;
			tdiss = SQR(TDPD_WF(r, mcom->dpdrc)) * b->gamma_p[tmp];
			trand = TDPD_WF(r, mcom->dpdrc) * b->zeta_p[tmp];

			tdiss /= mcom->m;
			trand /= mcom->m;
			force /= mcom->m;
			for (d = 1; d < 3; d++) {
				/* longitudinal */
				b->f[p_row[j]][d] -= force * drdv[j][d];
				b->f[p_col[j]][d] += force * drdv[j][d];
				/* dissipative, transversal */
				b->f[p_row[j]][d] += tdiss * drdv[j][d + 3];
				b->f[p_col[j]][d] -= tdiss * drdv[j][d + 3];
				/* random, transversal */
				b->f[p_row[j]][d] -= trand * theta[2 * k + d - 1];
				b->f[p_col[j]][d] += trand * theta[2 * k + d - 1];
			}
		}
	}

	return energy;
}

void mcom_calc_forces(struct beads *b, struct mcom *mcom, enum nblist_flags fl)
{
	if (mcom == NULL) return;

	nblist_update(b, &mcom->nbl, fl);
	b->e[0] = mcom_forces2d(b, mcom, &mcom->nbl);
}

/* ------------------------------------------------------------------------- */

/* open an external file containing the coefficients for the b-splines */
static void mcom_get_coefficients(struct mcom *mcom, const char *fn)
{
	FILE *FH = fopen(fn, "r");
	char buf[LINE_MAX];
	int i = 0, k;

	debug("Trying to open '%s'", fn);	
	if (FH == NULL)
		fatal(ENOENT, "Couldn't open '%s' for reading", fn);

	if (get_data_line(FH, "# k=%d ncoeffs=%d rc=%lg dim=%d",
			&k, &mcom->ncoeffs, &mcom->rc, &mcom->diml) != 4)
		fatal(EINVAL, "Couldn't understand mcom header");
	
	mcom->rc2 = SQR(mcom->rc);
	debug("k=%d ncoeffs=%d rc=%lg dim=%d", k, mcom->ncoeffs, mcom->rc,
								mcom->diml);

	mcom->bw = gsl_bspline_alloc(k, mcom->ncoeffs - k + 2);
	gsl_bspline_knots_uniform(0., mcom->rc, mcom->bw);
	mcom->bdw = gsl_bspline_deriv_alloc(k);

	mcom->dB = gsl_matrix_alloc(k, 2);

	mcom->coeffs = calloc(2 * mcom->ncoeffs, sizeof(*mcom->coeffs));
	if (mcom->coeffs == NULL) novm("mcom->coeffs");

	if (mcom->diml != 2)
		fatal(EINVAL, "only diml=2 supported");

	while (fgets(buf, sizeof(buf), FH) != NULL && i < mcom->ncoeffs) {
		int tmp;

		if (strchr(buf, '#') == buf) continue;
		if (sscanf(buf, "%d %lg %lg", &tmp, mcom->coeffs + i,
					mcom->coeffs + mcom->ncoeffs + i) != 3)
			fatal(EINVAL, "Error reading coeff %d", i);
		if (tmp != i) fatal(EINVAL, "Error reading '%s'", buf);
		i++;
	}
	fclose(FH);

	if (i < mcom->ncoeffs)
		fatal(EINVAL, "Only %d coefficients found, but expected"
				" %d in file '%s'", i, mcom->ncoeffs, fn);
}

static void print_results(struct mcom *mcom)
{
	int i, l;
	double r;

	gsl_matrix *dB = gsl_matrix_alloc(mcom->ncoeffs, mcom->bw->k);

	for (r = 0.; r < mcom->rc; r += .01) {
		printf("%lg", r);
		gsl_bspline_deriv_eval(r, 2, dB, mcom->bw, mcom->bdw);
		gsl_vector_const_view vu = gsl_matrix_const_column(dB, 0);
		gsl_vector_const_view vf = gsl_matrix_const_column(dB, 1);
		for (l = 0; l < 2; l++) {
			double u = 0., du = 0.;
			for (i = 0; i < mcom->ncoeffs; i++) {
				u += mcom->coeffs[l * mcom->ncoeffs + i] * gsl_vector_get(&vu.vector, i);
				du += mcom->coeffs[l * mcom->ncoeffs + i] * gsl_vector_get(&vf.vector, i);
			}
			printf("\t%lg\t%lg", u, du);
		}
		printf("\n");
	}

	printf("## # ncoeffs=%d dim=%d\n## # i\tf_i\tg_i\n", mcom->ncoeffs, mcom->diml);
	for (i = 0; i < mcom->ncoeffs; i++) {
		printf("## %d\t%lg\t%lg\n", i, mcom->coeffs[i], mcom->coeffs[i + mcom->ncoeffs]);
	}
	gsl_matrix_free(dB);
}

/*
 * if a line starting with "# mcom" is found in the system file, then the
 * mcom force field is initialized here.
 */
void mcom_init(FILE *FH, struct beads *b)
{
	int n;
	char buf[LINE_MAX];
	
	struct mcom *mcom = malloc(sizeof(*mcom));
	if (mcom == NULL) novm("mcom");
	b->mcom = mcom;

	/* read everything from another external file */
	fgets(buf, sizeof(buf), FH);
	*strchr(buf, '\n') = 0x0;
	sscanf(buf, "# mcom m=%lg dpdrc=%lg fn=%n", &mcom->m, &mcom->dpdrc, &n);
	if (mcom->m <= 0. || mcom->m > 1000.)
		fatal(EDOM, "mcom mass must be positive");

	strncpy(mcom->fn, buf + n, NAME_MAX);
	printf("fn=%s\n", mcom->fn);
	mcom_get_coefficients(mcom, mcom->fn);

	if (mcom->dpdrc > mcom->rc)
		fatal(EDOM, "dpdrc must be smaller than rc");
	if (mcom->dpdrc <= 0.)
		fatal(EDOM, "dpdrc must be positive");
	
	if (ismaster) print_results(mcom);
}

void mcom_print_header(FILE *FH, struct mcom *mcom)
{
	if (mcom != NULL)
		fprintf(FH, "# mcom m=%lg dpdrc=%lg fn=%s\n", mcom->m, 
							mcom->dpdrc, mcom->fn);
}

void mcom_free(struct mcom *mcom)
{
	if (mcom == NULL) return;
	
	gsl_bspline_free(mcom->bw);
	gsl_bspline_deriv_free(mcom->bdw);
	gsl_matrix_free(mcom->dB);
	free(mcom->coeffs);
	nblist_free(&mcom->nbl);
	free(mcom);
}

/* initialize the neighbor list */
void mcom_init_nblist(struct beads *b, struct mcom *mcom)
{
	if (mcom == NULL) return;

	nblist_alloc(b, &mcom->nbl, groups_per_row * b->groupsize,
			b->passport, b->xv, groups_per_col * b->groupsize,
			b->passport + groups_per_row * b->groupsize,
			b->xv + groups_per_row * b->groupsize,
			mcom->rc + b->rs, NBL_NEWTON3);
}

double mcom_get_rc(struct mcom *mcom)
{
	return mcom->rc;
}

int mcom_get_diml(struct mcom *mcom)
{
	return mcom->diml;
}

/* ------------------------------------------------------------------------- */

static void set_tdpd(struct mcom *mcom, int type1, int type2, double g, double h, double dt)
{
	debug("setting TDPD %d <--> %d to %lg, %lg", type1, type2, g, h);
	mcom->gamma[type1 * TYPE_MAX + type2] = g;
	mcom->zeta[type1 * TYPE_MAX + type2] = sqrt(12. * 2. * g / dt);
	mcom->gamma_p[type1 * TYPE_MAX + type2] = h;
	mcom->zeta_p[type1 * TYPE_MAX + type2] = sqrt(12. * 2. * h / dt);
}

/* Setup transversal DPD thermostat */
static void parse_tdpd(struct mcom *mcom, cfg_t *cfg, double def_gamma,
						double def_gammap, double dt)
{
#if (TYPE_MAX == 2)
	int map[SQR(TYPE_MAX)] = {0, 1, 1, 2};
#else
	int map[SQR(TYPE_MAX)] = {0, 1, 2, 1, 3, 4, 2, 4, 5};
#endif
	int i, j;
	
	if (cfg_size(cfg, "tdpd") > 0) {
		cfg_t *sec = cfg_getsec(cfg, "tdpd");

		if (cfg_size(sec, "gamma") != TYPE_MAX * (TYPE_MAX + 1) / 2)
			fatal(EINVAL, "Incorrect number of TDPD frict"
			"ion coefficients specified. Is: %u. Should be: %d",
			cfg_size(sec, "gamma"), TYPE_MAX * (TYPE_MAX + 1) / 2);
		if (cfg_size(sec, "gamma_p") != TYPE_MAX * (TYPE_MAX + 1) / 2)
			fatal(EINVAL, "Incorrect number of TDPD transve"
			"rsal friction coefficients specified. Is: %u. Should "
			"be: %d", cfg_size(sec, "gamma_p"),
			TYPE_MAX * (TYPE_MAX + 1) / 2);
		for (i = 0; i < TYPE_MAX; i++) {
			for (j = 0; j < TYPE_MAX; j++) {
				set_tdpd(mcom, i, j,
				       	cfg_getnfloat(sec, "gamma",
						map[i * TYPE_MAX + j]),
					cfg_getnfloat(sec, "gamma_p",
					       	map[i * TYPE_MAX + j]), dt);
			}
		}
		if (ismaster) printf("fine-tuned TDPD thermostat.\n");
	} else {
		for (i = 0; i < TYPE_MAX; i++) {
			for (j = 0; j < TYPE_MAX; j++) { 
				set_tdpd(mcom, i, j, def_gamma, def_gammap, dt);
			}
		}
		if (ismaster)
			printf("Transverse DPD thermostat: "
			"gamma||=%lg gamma_|_=%lg\n", def_gamma, def_gammap);
	}
}

/*
 * Setup thermostat. If there are only the options in section "main", then
 * the old-style behavior is used, i.e. there is only one dissipation constant
 * for all interactions. In the new style, each kind of pairwise interaction
 * can have its own dissipation constant. In case of the Langevin thermostat
 * each direction of each species can have its own constant. The new style
 * overrides the old style.
 */
void mcom_parse_thermostat(struct mcom *mcom, cfg_t *cfg)
{
	if (mcom == NULL) return;

	cfg_t *m = cfg_getsec(cfg, "main");
	double gamma = cfg_getfloat(m, "gamma");
	double gamma_p = cfg_getfloat(m, "gamma_p");
	double dt = cfg_getfloat(m, "dt");
	const char *s = cfg_getstr(m, "thermostat");
	
	if (strcasecmp(s, "DPD") == 0 || strcasecmp(s, "TDPD") == 0) {
		parse_tdpd(mcom, cfg, gamma, gamma_p, dt);
		
	} else
		fatal(EINVAL, "Invalid thermostat '%s' specified", s);
}

