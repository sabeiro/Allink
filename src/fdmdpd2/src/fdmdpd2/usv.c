/*
 * usv.c - Umbrella Sampling for Vesicles
 * Copyright (C) 2011 Martin Hoemberg <mhoembe@gwdg.de> and
 *                    Marc Fuhrmans <fuhrmans@theorie.physik.uni-goe.de>
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

/* $Id: usv.c 296 2011-06-15 12:17:05Z hoemberg $ */

#include "fdmdpd2.h"
//#include "dsyevc3.c"

/* ------------------------------------------------------------------------- */

/* options for the tensor constraint */
static cfg_opt_t tensor_opts[] = {
	CFG_FLOAT_LIST("gyr0", "{}", CFGF_NONE),
	CFG_FLOAT("k", 0.0, CFGF_NONE),
	CFG_BOOL("couple_to_com", cfg_false, CFGF_NONE),
	CFG_INT("couple_to_block", 0, CFGF_NONE),
	CFG_END()
}; 

/* collective options for the usv.o module */
cfg_opt_t usv_opts[] = {
	CFG_SEC("tensor", tensor_opts, CFGF_TITLE | CFGF_MULTI),
	CFG_END()
};

/* list of implemented constraint types */
enum cnstrnttype {
	CNSTRNT_TENSOR
};

/*
 * structure describing a single constraint. The fields which are specific
 * to a single constraint are also visible to other types of constraint.
 */
struct cnstrnt
{
	const char *name;	/* name of this constraint */
	double gyr[9];		/* current gyration tensor */
	double gyr0[9];		/* values for a tensor constraint */
	double k;		/* spring constant */
	double sqr_tr0;		/* (tr(gyr0))^2 */
	int couplecom;		/* couple to the centers of mass ? */
	int block;		/* couple to which block ? */
	enum cnstrnttype type;	/* type of this constraint */
};

/* collective structure for all constraints */
struct usv
{
	struct cnstrnt *c;	/* array of the umbrella potentials */
	VEC *mcom;		/* center of mass of the lipids */
	double *ct;		/* cos(theta_i) */
	double com[3];		/* center of mass of the whole system */
	int n;			/* # constraint sections in cfg file */
	int calcmcom;		/* compute molecular centers of mass? */
};

/* ------------------------------------------------------------------------- */

/* print out all entries of a tensor, nicely formatted */
static char buf[LINE_MAX];
static char *print_tensor(const double t[])
{
	snprintf(buf, sizeof(buf), "( %lg %lg %lg / %lg %lg %lg / %lg %lg %lg )",
		t[0], t[1], t[2], t[3], t[4], t[5], t[6], t[7], t[8]);
	return buf;
}

/* print out status message for a tensor constraint */
static void tensor_printstatus(const struct cnstrnt *c)
{
	printf("S=%s\n", print_tensor(c->gyr));
}	

/* initialize tensor constraint */
static void tensor_init(struct cnstrnt *c, cfg_t *t)
{
	int i;

	c->type = CNSTRNT_TENSOR;
	c->name = cfg_title(t);
	c->k = cfg_getfloat(t, "k");
	c->couplecom = (cfg_getbool(t, "couple_to_com") == cfg_true) ? 1 : 0;
	c->block = cfg_getint(t, "couple_to_block");

	/* read in tensor values */
	if (cfg_size(t, "gyr0") == 3) {
		c->gyr0[0] = cfg_getnfloat(t, "gyr0", 0); /* xx */
		c->gyr0[4] = cfg_getnfloat(t, "gyr0", 1); /* yy */
		c->gyr0[8] = cfg_getnfloat(t, "gyr0", 2); /* zz */
	} else if (cfg_size(t, "gyr0") == 9) {
		for (i = 0; i < 9; i++)
			c->gyr0[i] = cfg_getnfloat(t, "gyr0", i);
	} else {
		fatal(EINVAL, "Invalid number of entries in 'gyr0'");
	}
	c->sqr_tr0 = SQR(c->gyr0[0] + c->gyr0[4] + c->gyr0[8]);

	printf("USV: enabled tensor constraint '%s' with k=%lg "
		"S=%s couple_to_com=%d couple_to_block=%d\n",
		c->name, c->k, print_tensor(c->gyr0), c->couplecom, c->block);
}

/* apply forces arising from tensor constraint */
/* TODO: implement coupling to centers of mass -- this option is currently ignored */
static double tensor_force(struct beads *b, struct usv *u, struct cnstrnt *c)
{
	int i, j, k, l, idx = 0, cnt = 0;
	double dgyr[9], e = 0.;

	for (i = 0; i < 9; i++)
		dgyr[i] = c->gyr[i] - c->gyr0[i];

	for (i = 0; i < b->local_n; i++) {
		int N = b->N[b->local_n_b[i]];
		
		if (b->local_n_b[i] == c->block) {
			for (j = idx; j < idx + N; j++) {
				for (k = 0; k < 3; k++) {
					double tmp = 0.;
					for (l = 0; l < 3; l++)
						tmp += dgyr[3 * k + l] * (b->x_intra[idx][l] - u->com[l]);
					b->f_intra[idx][k] -= 2. * c->k * tmp / c->sqr_tr0;
				}
			}
			cnt += N;
		}
		idx += N;
	}

	MPI_Allreduce(MPI_IN_PLACE, &cnt, 1, MPI_INT, MPI_SUM, comm_grid);

	/* return the total energy only in the master process */
	if (ismaster) {
		for (i = 0; i < 9; i++)
			e += SQR(c->gyr[i] - c->gyr0[i]);
		e *= .5 * c->k * cnt / c->sqr_tr0;
	}
	return e;
}

/* ------------------------------------------------------------------------- */

/* update the center of mass for the whole system */
static void getcom(struct beads *b, double com[])
{
	int i, j, k, cnt = 0;
	
	for (i = 0; i < 3; i++)
		com[i] = 0.;

	for (i = 0; i < b->local_n; i++) {
		int N = b->N[b->local_n_b[i]];

		for (j = cnt; j < cnt + N; j++) {
			for (k = 0; k < 3; k++)
				com[k] += b->x_intra[j][k];
		}
		cnt += N;
	}
	
	MPI_Allreduce(MPI_IN_PLACE, com, 3, MPI_DOUBLE, MPI_SUM, comm_grid);
	MPI_Allreduce(MPI_IN_PLACE, &cnt, 1, MPI_INT, MPI_SUM, comm_grid);

	for (i = 0; i < 3; i++)
		com[i] /= cnt;
}

/* calculate gyration tensor of a specific block and return squared radius of gyration */
static double getgyrtensor(struct beads *b, double gyr[], const double com[], const int block)
{
	int i, j, k, l, cnt = 0, idx = 0;
	double t[3];

	for (i = 0; i < 9; i++)
		gyr[i] = 0.;

	for (i = 0; i < b->local_n; i++) {
		int N = b->N[b->local_n_b[i]];

		if (b->local_n_b[i] == block) {
			for (j = idx; j < idx + N; j++) {
				for (k = 0; k < 3; k++)
					t[k] = b->x_intra[j][k] - com[k];
				for (k = 0; k < 3; k++) {
					for (l = 0; l < 3; l++) {
						gyr[3 * k + l] += t[k] * t[l];
					}
				}
			}
			cnt += N;
		}
		idx += N;
	}

	MPI_Allreduce(MPI_IN_PLACE, gyr, 9, MPI_DOUBLE, MPI_SUM, comm_grid);
	MPI_Allreduce(MPI_IN_PLACE, &cnt, 1, MPI_INT, MPI_SUM, comm_grid);

	for (i = 0; i < 9; i++)
		gyr[i] /= cnt;

	return gyr[0] + gyr[4] + gyr[8]; /* trace of gyration tensor */
}

/* calculate the centers of mass of ALL local molecules */
static void calc_mcom(struct beads *b, struct usv *u)
{
	int i, j, k, idx = 0;
	
	for (i = 0; i < b->local_n; i++) {
		int block = b->local_n_b[i];
		int N = b->N[block];
		VEC com = { 0 };

		for (j = idx; j < idx + N; j++) {
			for (k = 0; k < 3; k++) {
				com[k] += b->x_intra[j][k];
			}
		}
		for (j = 0; j < 3; j++)
			u->mcom[i][j] = com[j] / N;
		idx += N;
	}
}

/* ------------------------------------------------------------------------- */

static char strange_error[] = "Unknown type. This should never happen!\n\n";

/* print current radius of gyration, the eigenvalues, and the com */
void usv_printstatus(const struct usv *u)
{
	if (u == NULL || !ismaster) return;
	int i;

	printf("USV: com=%lg %lg %lg\n", u->com[0], u->com[1], u->com[2]);

	for (i = 0; i < u->n; i++) {
		printf("USV %s: ", u->c[i].name);

		switch (u->c[i].type) {
		case CNSTRNT_TENSOR:
			tensor_printstatus(u->c + i);
			break;
		default:
			fprintf(stderr, strange_error);
		}
	}
}

/* calculate all forces due to the Umbrella Sampling potential */
void usv_force(struct beads *b, struct usv *u)
{
	if (u == NULL) return;
	int i;

	getcom(b, u->com);

	/* compute the molecular centers of mass only, if at least a single
	   constraint requires it. */
	if (u->calcmcom)
		calc_mcom(b, u);

	for (i = 0; i < u->n; i++) {
		getgyrtensor(b, u->c[i].gyr, u->com, u->c[i].block);

		if (u->c[i].k == 0.) continue;

		switch (u->c[i].type) {
		case CNSTRNT_TENSOR:
			b->e[3] += tensor_force(b, u, u->c + i);
			break;
		default:
			fprintf(stderr, strange_error);
		}
	}
}

/* init struct */
/* TODO: implement facility for other types of constraint */
void usv_init(struct beads *b, cfg_t *cfg)
{
	if (cfg_size(cfg, "usv") == 0) return;

	cfg_t *t = cfg_getsec(cfg, "usv");
	const int constraints = cfg_size(t, "tensor");
	int i;

	if (constraints > 0) {
		struct usv *u = malloc(sizeof(*b->usv));
		if (u == NULL) novm("b->usv");

		u->n = constraints;
		u->c = calloc(u->n, sizeof(*u->c));
		if (u->c == NULL) novm("b->usv->c");
			
		/* perform general initialization */
		getcom(b, u->com);
		u->mcom = calloc(b->local_n, sizeof(*u->mcom));
		if (u->mcom == NULL) novm("u->mcom");

		/* parse the tensor constraints */
		for (i = 0; i < u->n; i++) {
			cfg_t *tt = cfg_getnsec(t, "tensor", i);

			tensor_init(u->c + i, tt);

			if (u->c[i].couplecom)
				u->calcmcom = 1;
		}
		b->usv = u;
	}
}

/* free struct */
void usv_free(struct beads *b)
{
	if (b->usv != NULL) {
		free(b->usv->c);
		free(b->usv->mcom);
		free(b->usv);
	}
}

