/*
 * main.c - A Parallel Force Decomposition MDPD simulation
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

/* $Id: main.c 366 2013-08-08 14:19:08Z fuhrmans $ */

#include "fdmdpd2.h"
#include "rand.h"
#include <signal.h>
#include <fenv.h>

/* ------------------------------------------------------------------------- */

static cfg_opt_t main_opts[] = {
	CFG_INT("seed", 0, CFGF_NONE),
	CFG_FLOAT("total_time", 1., CFGF_NONE),
	CFG_STR("thermostat", "DPD", CFGF_NONE),
	CFG_FLOAT("gamma", 0.5, CFGF_NONE),
	CFG_FLOAT("gamma_p", 0.0, CFGF_NONE),
	CFG_FLOAT("dt", 0.005, CFGF_NONE),
	CFG_FLOAT("rs", 0., CFGF_NONE),
	CFG_INT("bilayer_normal", 0, CFGF_NONE),
	CFG_END()
};

static cfg_opt_t task_opts[] = {
	CFG_FLOAT("dt", 0., CFGF_NONE),
	CFG_STR("filename", NULL, CFGF_NONE),
	CFG_INT("bufsize", BUFSIZ, CFGF_NONE),
	CFG_INT("snapshots", 1000, CFGF_NONE),
	CFG_END()
};

static cfg_opt_t dpd_opts[] = {
	CFG_FLOAT_LIST("gamma", "{}", CFGF_NONE),
	CFG_END()
};

static cfg_opt_t tdpd_opts[] = {
	CFG_FLOAT_LIST("gamma", "{}", CFGF_NONE),
	CFG_FLOAT_LIST("gamma_p", "{}", CFGF_NONE),
	CFG_END()
};

static cfg_opt_t langevin_opts[] = {
	CFG_FLOAT_LIST("gamma_x", "{}", CFGF_NONE),
	CFG_FLOAT_LIST("gamma_y", "{}", CFGF_NONE),
	CFG_FLOAT_LIST("gamma_z", "{}", CFGF_NONE),
	CFG_END()
};

static cfg_opt_t opts[] = {
	CFG_SEC("main", main_opts, CFGF_NONE),
	CFG_SEC("task", task_opts, CFGF_TITLE | CFGF_MULTI),
	CFG_SEC("nblist", nblist_opts, CFGF_NONE),
	CFG_SEC("integrator", integrator_opts, CFGF_NONE),
	CFG_SEC("Tens", Tens_opts, CFGF_NONE),
	//CFG_SEC("Diff", Diff_opts, CFGF_NONE),
	CFG_SEC("Widom", Widom_opts, CFGF_NONE),
	CFG_SEC("isf", dyn_opts, CFGF_NONE),
	CFG_SEC("rnemd", rnemd_opts, CFGF_NONE),
	CFG_SEC("dpd", dpd_opts, CFGF_NODEFAULT),
	CFG_SEC("tdpd", tdpd_opts, CFGF_NODEFAULT),
	CFG_SEC("langevin", langevin_opts, CFGF_NODEFAULT),
	CFG_SEC("usv", usv_opts, CFGF_NODEFAULT),
	CFG_SEC("mdpd_sg", mdpd_sg_opts, CFGF_NODEFAULT),
	CFG_END()
};

/* ------------------------------------------------------------------------- */

MPI_Comm comm_grid;	/* MPI communicator for all processes */
MPI_Comm comm_row;	/* MPI communicator for the current row */
MPI_Comm comm_col;	/* MPI communicator for the current column */
int groups_per_row;	/* number of groups in one row */
int groups_per_col;	/* number of groups in one column */
int row_index;		/* row index of this process */
int col_index;		/* column index of this process */
int bilayer_orientation;/* bitfield of the axes' indices */
int ismaster;		/* is this the master process? */
struct rng *rng = NULL;	/* random number generator for thermostat, etc. */

/* ------------------------------------------------------------------------- */

/*
 * It is a good idea to terminate the program cleanly, before we are running
 * out of time and the queueing system kills the processes. One way of doing
 * this is to let the queueing system send a signal some minutes before the
 * end of the timelimit. For instance, in a PBS job script a SIGTERM can be
 * sent 300 seconds before the end of the timelimit by including a line like 
 * 
 * #PBS -l signal=15@300
 *
 * Such a signal is normally caught by mpiexec/mpirun and forwarded to all MPI
 * processes, which react correspondingly. Here, a signal handler is installed
 * that catches SIGTERM, SIGINT, SIGUSR1, and SIGURG, so that all processes
 * are able to shutdown cleanly.
 *
 * However, in practice this doesn't always work. Different versions of mpiexec
 * forward only some signals (some of them don't forward at all), different
 * versions of the queueing system send the signals to different processes,
 * and the details of the job script also seem to matter. The following
 * combinations have proven to work:
 *
 * - HLRN UV System, MPT 2.01: signal SIGURG (23)
 * - HLRN XE/ICE1/ICE2, MVAPICH2 1.4.1: signal SIGTERM (15)
 * - JSC JUROPA, Parastation MPI2: signal SIGTERM (15)
 */

static volatile sig_atomic_t __terminated = 0;

static void handle_sigterm(int sig)
{
	__terminated = 1;
}

static void fillmyset(sigset_t *mask)
{
	sigemptyset(mask);
	sigaddset(mask, SIGTERM);
	sigaddset(mask, SIGINT);
	sigaddset(mask, SIGUSR1);
	sigaddset(mask, SIGURG);
}

static void install_sigterm_handler(void)
{
	struct sigaction setup_action;
	
	fillmyset(&setup_action.sa_mask);
	setup_action.sa_flags = 0;
	setup_action.sa_handler = handle_sigterm;

	sigaction(SIGUSR1, &setup_action, NULL);
	sigaction(SIGTERM, &setup_action, NULL);
	sigaction(SIGINT, &setup_action, NULL);
	sigaction(SIGURG, &setup_action, NULL);

	debug("Installed SIGTERM/SIGUSR1/SIGINT/SIGURG handler");
}

static int terminated(void)
{
	sigset_t sigset;

	fillmyset(&sigset);
	sigprocmask(SIG_BLOCK, &sigset, NULL);
	
	int t = __terminated;
	MPI_Allreduce(MPI_IN_PLACE, &t, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	__terminated = t;

	sigprocmask(SIG_UNBLOCK, &sigset, NULL);
	return t;
}

/* ------------------------------------------------------------------------- */

/*
 * Create 2diml MPI process grid topoplogy
 */
static void create_grid_topology(const int rows, const int cols)
{
	int dim[2] = {rows, cols};
	int periods[2] = {0, 0};
	int remain[2];
	int size;
	char name[MPI_MAX_PROCESSOR_NAME];
	
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	if (unlikely(size != rows * cols))
		fatal(EINVAL, "Number of CPUs mismatches %d*%d != %d",
					rows, cols, size);

	MPI_Cart_create(MPI_COMM_WORLD, 2, dim, periods, 1, &comm_grid);
	remain[0] = 0; remain[1] = 1;
	MPI_Cart_sub(comm_grid, remain, &comm_row);
	remain[0] = 1; remain[1] = 0;
	MPI_Cart_sub(comm_grid, remain, &comm_col);

	MPI_Comm_rank(comm_col, &row_index);
	MPI_Comm_rank(comm_row, &col_index);

	groups_per_row = size / rows;
	groups_per_col = size / cols; 
	
	MPI_Get_processor_name(name, &size);
	debug("coords: row=%d/%d, col=%d/%d, name='%s'",
				row_index, rows, col_index, cols, name);

	ismaster = (row_index == 0 && col_index == 0) ? 1 : 0;
}

static void free_grid_topology(void)
{
	MPI_Comm_free(&comm_row);
	MPI_Comm_free(&comm_col);
	MPI_Comm_free(&comm_grid);
}

/* ------------------------------------------------------------------------- */

/*
 * Communicate forces and sum up the intra- and intermolecular contributions.
 * After this function this process has full knowledge of the forces acting
 * on its local monomers (but not about non-local monomers).
 */
static void comm_forces(const struct beads *restrict b, VEC *rbuf,
				MPI_Request sreq[], MPI_Request rreq[])
{
	int i;

	sreq[col_index] = MPI_REQUEST_NULL;
	rreq[col_index] = MPI_REQUEST_NULL;
	sreq[groups_per_row + row_index] = MPI_REQUEST_NULL;
	rreq[groups_per_row + row_index] = MPI_REQUEST_NULL;

	for (i = 0; i < groups_per_row; i++) {
		if (i == col_index) continue;
		MPI_Isend(b->f + i * b->groupsize,
				ARRAY_SIZE(*b->f) * b->groupsize, MPI_DOUBLE,
				i, 0xa, comm_row, sreq + i);
		MPI_Irecv(rbuf + i * b->groupsize,
				ARRAY_SIZE(*b->f) * b->groupsize, MPI_DOUBLE,
				i, 0xa, comm_row, rreq + i);
	}
	
	for (i = 0; i < groups_per_col; i++) {
		if (i == row_index) continue;
		MPI_Isend(b->f + (i + groups_per_row) * b->groupsize,
				ARRAY_SIZE(*b->f) * b->groupsize, MPI_DOUBLE,
				i, 0xb, comm_col, sreq + i + groups_per_row);
		MPI_Irecv(rbuf + (i + groups_per_row) * b->groupsize,
				ARRAY_SIZE(*b->f) * b->groupsize, MPI_DOUBLE,
				i, 0xb, comm_col, rreq + i + groups_per_row);
	}
}

static int wait_request(MPI_Request req[], int *idx)
{
	MPI_Status status;
	MPI_Waitany(groups_per_row + groups_per_col, req, idx, &status);
	return *idx;
}

static void add_forces(const struct beads *restrict b, VEC *rbuf,
					MPI_Request sreq[], MPI_Request rreq[])
{
	int i, d, idx;
	VEC *f_row = b->f + b->groupsize * col_index;
	VEC *f_col = b->f + b->groupsize * (groups_per_row + row_index);
	MPI_Status status[groups_per_row + groups_per_col];

	/* add intramolecular forces */
	for (i = 0; i < b->groupsize; i++) {
		for (d = 0; d < 3; d++) {
		  f_row[i][d] += f_col[i][d] + b->f_intra[i][d];
		}
	}

	while (wait_request(rreq, &idx) != MPI_UNDEFINED) {
		VEC *r = rbuf + idx * b->groupsize;

		for (i = 0; i < b->groupsize; i++) {
			for (d = 0; d < 3; d++) {
				f_row[i][d] += r[i][d];
			}
		}
	}

	/* make sure that the sends have completed */
	MPI_Waitall(groups_per_row + groups_per_col, sreq, status);
}

static void set_xf_intra(struct beads *b)
{
	int i, o = col_index * b->groupsize, d;

	for (i = 0; i < b->groupsize; i++) {
		for (d = 0; d < 3; d++)  {
			b->x_intra[i][d] = b->xv[i + o][d] +
						b->l[d] * b->nx[i + o][d];
			b->f_intra[i][d] = 0.;
		}
	}
}

/*
 * Depending on the block different intramolecular interactions are being
 * calculated. If the name of the block is unknown, then the standard
 * interactions for linear molecules are used.
 */
static void calc_intra(struct beads *b)
{
	int i, bp;

	for (i = 0, bp = 0; i < b->local_n; i++) {
		const char *name = b->resname[b->local_n_b[i]];
		const int N = b->N[b->local_n_b[i]];

		if (strcasestr(name, "PEP") == name && b->NPep > 0) {
			PepCalcForces(b, bp, name, N);

		} else if(strcasestr(name, "TT") == name) {
			mdpd_two_tails(b, b->mdpd, N, bp);

		} else {
			mdpd_linear_chain(b, b->mdpd, N, bp);
		}

		bp += N;
	}
}

/*
 * This is the general force calculation function. No calculations are done
 * here, but the corresponding functions are called. If a different force field
 * or additional interactions are needed, then this is the right place to call
 * the corresponding functions.
 */
void forces_calc(struct beads *b)
{
	MPI_Request sreq[groups_per_row + groups_per_col];
	MPI_Request rreq[groups_per_row + groups_per_col];
	static double maxdispl = 12345.; /* dummy initial value */
	static double volume = 0.;
	enum nblist_flags flags;
	VEC *rbuf;

	/* determine the current state of the neighbor list */
	if (volume != b->l[0] * b->l[1] * b->l[2]) {
		volume = b->l[0] * b->l[1] * b->l[2];
		maxdispl = 0.;
		flags = NBL_VOLUME_CHANGED;
	} else if (maxdispl > .5 * b->rs) {
		maxdispl = 0.;
		flags = NBL_INVALID;
	} else {
		/* TODO: implement a less pessimistic assumption */
		maxdispl += sqrt(b->v2max) * b->dt;
		flags = NBL_VALID;
	}
	
	/* reset the force buffers and calculate new forces */
	memset(b->f, 0, (groups_per_row + groups_per_col) *
					b->groupsize * sizeof(*b->f));
	memset(&b->virial, 0, sizeof(b->virial));
	memset(b->e, 0, ARRAY_SIZE(b->e) * sizeof(*b->e));

	/* TODO: put this calloc() somewhere else, so that this buffer is
	 * not alloc'd / free'd every timestep. */
	rbuf = calloc((groups_per_row + groups_per_col) * b->groupsize,
								sizeof(*rbuf));
	if (unlikely(rbuf == NULL)) novm("rbuf");

	/* non-bonded, i.e. pairwise forces first! */
	mdpd_nonbonded_forces(b, b->mdpd, flags);
#ifdef MCOM
	mcom_calc_forces(b, b->mcom, flags);
#endif
	comm_forces(b, rbuf, sreq, rreq);

	/* while the communications of the non-bonded interactions are
	 * "on the way", the other interactions are calculated */
	set_xf_intra(b);

	calc_intra(b);
	RigidCalcForces(b);
	ExtCalcForces(b);

        surfpot_calc_forces(b, b->surfpot);
        surfbump_calc_forces(b, b->surfbump);
        surfcuff_calc_forces(b, b->surfcuff);
        apply_blockforce(b, b->blockforce);
	usv_force(b, b->usv);
        usp_wrapper(b, b->usp);
        multiblockforce_wrapper(b, b->multiblockforce);

	add_forces(b, rbuf, sreq, rreq);

	free(rbuf);
}
	
/* ------------------------------------------------------------------------- */

/*
 * allocate memory for beads struct
 */
void alloc_beads(struct beads *b)
{
	size_t count = (groups_per_row + groups_per_col) * b->groupsize;
	int i;

	b->f = calloc(count, sizeof(*b->f));
	if (b->f == NULL) novm("beads->f");
	b->xv = calloc(count, sizeof(*b->xv));
	if (b->xv == NULL) novm("beads->xv");
	b->nx = calloc(count, sizeof(*b->nx));
	if (b->nx == NULL) novm("beads->nx");
	b->passport = calloc(count, sizeof(*b->passport));
	if (b->passport == NULL) novm("beads->passport");
	b->v2 = calloc(count, sizeof(*b->v2));
	if (b->v2 == NULL) novm("beads->v2");

	b->x_intra = calloc(b->groupsize, sizeof(*b->x_intra));
	if (b->x_intra == NULL) novm("beads->x_intra");
	b->f_intra = calloc(b->groupsize, sizeof(*b->f_intra));
	if (b->f_intra == NULL) novm("beads->f_intra");

	count = 0;
	for (i = 0; i < b->blocks; i++)
		count += b->n[i];
	
	b->local_n_b = calloc(count, sizeof(*b->local_n_b));
	if (b->local_n_b == NULL) novm("beads->local_n_b");
	b->local_n = 0;
}

/*
 * free allocated memory of beads struct
 */
void free_beads(struct beads *restrict b)
{
	free(b->f);
	free(b->xv);
	free(b->nx);
	free(b->passport);
	free(b->v2);
	free(b->x_intra);
	free(b->f_intra);
	free(b->local_n_b);

	free(b->n);
	free(b->N);
	free(b->resname);

	mdpd_free(b->mdpd);
#ifdef MCOM
	mcom_free(b->mcom);
#endif
        surfpot_free(b->surfpot);
        surfbump_free(b->surfbump);
        surfcuff_free(b->surfcuff);
        blockforce_free(b->blockforce);
        usp_free(b->usp);
        multiblockforce_free(b->multiblockforce);
	RigidFree(b);
	PepFree(b);
	ExtFree(b);
}

/*
 * allocate memory for the different blocks / molecular species
 */
void alloc_blocks(struct beads *restrict b)
{
	b->n = calloc(b->blocks, sizeof(*b->n));
	if (b->n == NULL) novm("beads->n");
	b->N = calloc(b->blocks, sizeof(*b->N));
	if (b->N == NULL) novm("beads->N");
	b->resname = calloc(b->blocks, sizeof(*b->resname));
	if (b->resname == NULL) novm("beads->resname");
}

/*
 * parse config file for thermostat parameters and perform setup
 */
static void parse_thermostat(struct beads *b, cfg_t *cfg)
{
	mdpd_parse_thermostat(b->mdpd, cfg);
#ifdef MCOM	
	mcom_parse_thermostat(b->mcom, cfg);
#endif
}

/*
 * Read in the global configuration file and setup some general things
 */
static void parse_mainoptions(struct beads *restrict b, cfg_t *cfg)
{
	cfg_t *m = cfg_getsec(cfg, "main");

	/* time steps */
	b->dt = cfg_getfloat(m, "dt");
	b->delta_t = cfg_getfloat(m, "total_time");
	debug("delta_t=%lg dt=%lg", b->delta_t, b->dt);

	/* bilayer orientation */
	switch (cfg_getint(m, "bilayer_normal")) {
	case 0: /* x is normal */
		SET_NORMAL(0);
		SET_TANG1(1);
		SET_TANG2(2);
		break;
	case 1: /* y is normal */
		SET_NORMAL(1);
		SET_TANG1(2);
		SET_TANG2(0);
		break;
	case 2: /* z is normal */
		SET_NORMAL(2);
		SET_TANG1(0);
		SET_TANG2(1);
		break;
	}
	debug("Bilayer Orientation: normal=%d tangential=%d %d", NORMAL,
								TANG1, TANG2);
	/* nblist shell */
	b->rs = cfg_getfloat(m, "rs");
	
	/* random number generator */
	unsigned long seed = cfg_getint(m, "seed");
	rng = rng_alloc(seed);
	debug("RNG: seed=%lu", seed);
}

/*
 * read in header of system file and call handlers for the force fields
 */ 
static void parse_header(FILE *FH, struct beads *restrict b)
{
	char buf[LINE_MAX];
	long fpos;
		
	if (get_data_line(FH, "# L=%lg %lg %lg t=%lg blocks=%d", &b->l[0],
			&b->l[1], &b->l[2], &b->time, &b->blocks) != 5)
		fatal(EINVAL, "First line of system file");
	debug("L=%lg %lg %lg t=%.8f blocks=%d", b->l[0], b->l[1], b->l[2],
							b->time, b->blocks);
	b->step = floor(b->time / b->dt + .5);

	do {
		fpos = ftell(FH);
		fgets(buf, sizeof(buf), FH);
		fseek(FH, fpos, SEEK_SET);

		if (strstr(buf, "# v=") == buf) {
			mdpd_init(FH, b);
#ifdef MCOM	
		} else if (strstr(buf, "# mcom") == buf) {
			mcom_init(FH, b);
#endif
		} else if (strstr(buf, "# Rigid") == buf) {
			RigidLoad(b, FH);
		} else if (strstr(buf, "# Pep") == buf) {
			PepLoad(b, FH);
		} else if (strstr(buf, "# Ext") == buf) {
			ExtLoad(b, FH);
		} else if (strstr(buf, "# surfpot") == buf) {
			surfpot_read_header(FH, b);
		} else if (strstr(buf, "# surfbump") == buf) {
			surfbump_read_header(FH, b);
		} else if (strstr(buf, "# surfcuff") == buf) {
			surfcuff_read_header(FH, b);
		} else if (strstr(buf, "# blockforce") == buf) {
			blockforce_read_header(FH, b);
		} else if (strstr(buf, "# multiblockforce") == buf) {
			multiblockforce_read_header(FH, b);
		} else if (strstr(buf, "# usp") == buf) {
			usp_read_header(FH, b);
                }
	} while (fpos < ftell(FH));
}

/*
 * read in blocks definitions / molecular species definitions
 */
static void parse_blocks(FILE *FH, struct beads *restrict b)
{
	int i, j;
	char buf[LINE_MAX];
	long offset = ftell(FH);

	for (i = 0, b->nN = 0; i < b->blocks; i++) {
		if (get_data_line(FH, "# n=%d N=%d name=%s",
					&b->n[i], &b->N[i], buf) != 3)
			fatal(EINVAL, "Invalid header of block %d", i);

		strncpy(b->resname[i], buf, sizeof(b->resname[i]) - 1);
		debug("chains: block=%d n=%d N=%d name='%s'", i, b->n[i],
					b->N[i], b->resname[i]);
		for (j = 0; j < b->n[i] * b->N[i]; j++)
			fgets(buf, LINE_MAX, FH);

		b->nN += b->n[i] * b->N[i];
	}
	
	fseek(FH, offset, SEEK_SET);
}

/* ------------------------------------------------------------------------- */

/*
 * struct twoints - small struct to sort a list of integers carrying some data
 */
struct twoints
{
	int idx;	/* sorting index */
	int data;	/* arbitrary data */	
};

/* sort index descending */
static int cmp_twoints_desc(const void *restrict _a, const void *restrict _b)
{
	const struct twoints *a = _a;
	const struct twoints *b = _b;
	if (a->idx > b->idx) return -1;
	if (a->idx < b->idx) return 1;
	return 0; 
}

/*
 * Determine the size of a group and distribute the groups among the
 * processes. The group size (in beads) is always an integer multiple of the
 * size of one molecule. The number of beads per process ("groupsize") is
 * determined. Approximate load-balancing is achieved with a greedy algorithm.
 */
static void create_groups(struct beads *restrict b, int chain_count[])
{
	struct twoints lbs[b->blocks]; /* blocks sorted descending by N */
	struct twoints *p; /* processes sorted by bead count */
	int i, n, size, rank;
	
	MPI_Comm_size(comm_grid, &size);
	MPI_Comm_rank(comm_grid, &rank);
	p = alloca(size * sizeof(*p));

	for (i = 0; i < b->blocks; i++) {
		lbs[i].data = i; /* block id */
		lbs[i].idx = b->N[i]; /* chain length in this block */
	}
	qsort(lbs, b->blocks, sizeof(*lbs), cmp_twoints_desc);

	for (i = 0; i < size; i++) {
		p[i].idx = 0; /* beads in this group */
		p[i].data = i; /* group id */
	}

	memset(chain_count, 0, b->blocks * sizeof(*chain_count));

	for (i = 0; i < b->blocks; i++) {
		for (n = 0; n < b->n[lbs[i].data]; n++) {
			p[size - 1].idx += b->N[lbs[i].data];
			if (p[size - 1].data == rank)
				chain_count[lbs[i].data]++;
			qsort(p, size, sizeof(*p), cmp_twoints_desc);
		}
	}

	b->groupsize = p[0].idx; /* largest bead count */
	debug("%d group(s) with %d bead(s) per group", size, b->groupsize);
}

static int is_chain_local(int chain_count[], int block)
{
	int rank, transmit;

	/*
	 * we already know how many chains from each type each process will
	 * take care of. The idea is to distribute the molecules in the way
	 * they come, until each process is full. The process with the lowest
	 * rank that is not yet full will get the molecule.
	 */
	MPI_Comm_rank(comm_grid, &rank);
	transmit = (chain_count[block] > 0) ? rank : 100000;
	MPI_Allreduce(MPI_IN_PLACE, &transmit, 1, MPI_INT, MPI_MIN, comm_grid);
	assert(transmit != 100000);

	if (rank == transmit) {
		chain_count[block]--;
		return 1;
	}
	return 0;
}

void bcast_beads(struct beads *b)
{
	int n = b->groupsize;
	int o = b->groupsize * groups_per_row;

	MPI_Allgather(MPI_IN_PLACE, 0, MPI_DOUBLE, b->xv,
		n * ARRAY_SIZE(*b->xv), MPI_DOUBLE, comm_row);
	MPI_Allgather(MPI_IN_PLACE, 0, MPI_DOUBLE, b->nx,
		n * ARRAY_SIZE(*b->nx), MPI_DOUBLE, comm_row);
	MPI_Allgather(MPI_IN_PLACE, 0, MPI_DOUBLE, b->v2,
		n, MPI_DOUBLE, comm_row);
	MPI_Allgather(MPI_IN_PLACE, 0, MPI_UNSIGNED_LONG_LONG, b->passport,
		n, MPI_UNSIGNED_LONG_LONG, comm_row);

	MPI_Allgather(MPI_IN_PLACE, 0, MPI_DOUBLE, b->xv + o,
		n * ARRAY_SIZE(*b->xv), MPI_DOUBLE, comm_col);
	MPI_Allgather(MPI_IN_PLACE, 0, MPI_DOUBLE, b->nx + o,
		n * ARRAY_SIZE(*b->nx), MPI_DOUBLE, comm_col);
	MPI_Allgather(MPI_IN_PLACE, 0, MPI_DOUBLE, b->v2 + o,
		n, MPI_DOUBLE, comm_col);
	MPI_Allgather(MPI_IN_PLACE, 0, MPI_UNSIGNED_LONG_LONG, b->passport + o,
		n, MPI_UNSIGNED_LONG_LONG, comm_col);
}

/*
 * set the passport, position, and velocity of a given (physical) bead. If the
 * bead belongs to the ones owned by this process, it is stored twice.
 */ 
void set_bead(struct beads *restrict b, VEC2 xv, int id, int pos, int type)
{
	int d;

	assert(pos >= 0);
	assert(pos < b->groupsize * (groups_per_row + groups_per_col));
	assert(id >= 0);
	assert(id < b->nN);
	assert(type >= 0);

	SET_ID(b->passport[pos], id);
	SET_TYPE(b->passport[pos], (unsigned long long)type);
	SET_EXISTS(b->passport[pos]);
	SET_POS(b->passport[pos], pos);

	for (d = 0; d < ARRAY_SIZE(*b->nx); d++) {
		b->nx[pos][d] = floor(xv[d] / b->l[d]);
		b->xv[pos][d] = fabs(xv[d] - b->nx[pos][d] * b->l[d]);
		b->xv[pos][d + 3] = xv[d + 3];
	}

	b->v2[pos] = SQR(xv[3]) + SQR(xv[4]) + SQR(xv[5]);
	if (b->v2[pos] > b->v2max)
		b->v2max = b->v2[pos];

	/* store this bead twice? */
	int tmp = pos - b->groupsize * col_index;
	if (tmp >= 0 && tmp < b->groupsize) {
		pos = b->groupsize * (groups_per_row + row_index) + tmp;
		set_bead(b, xv, id, pos, type);
	}
}

/*
 * read in all beads in all blocks
 */
static void parse_beads(FILE *FH, struct beads *restrict b, int chain_count[])
{
	int i, j, k, pos, idx;
	char buf[LINE_MAX];

	idx = 0;
	pos = b->groupsize * col_index;

	for (i = 0; i < b->blocks; i++) { /* blocks */
		fgets(buf, sizeof(buf), FH);
		for (j = 0; j < b->n[i]; j++) { /* molecules */
			if (is_chain_local(chain_count, i)) {
				for (k = 0; k < b->N[i]; k++) { /* beads */
					VEC2 xv;
					int type;

					if (fscanf(FH, "%lf %lf %lf %lf %lf "
					"%lf %d\n", &xv[0], &xv[1], &xv[2],
					&xv[3], &xv[4], &xv[5], &type) != 7)
						fatal(EINVAL,
						"error while reading coordinat"
						"es of particle %d", idx);
					set_bead(b, xv, idx++, pos++, type);
				}
				b->local_n_b[b->local_n++] = i;
			} else { /* ignore this molecule */
				for (k = 0; k < b->N[i]; k++, idx++)
					fgets(buf, sizeof(buf), FH);
			}
		}
	}
	bcast_beads(b);
}

static void parse_systemfile(const char *fn, struct beads *restrict b)
{
	FILE *FH = fopen(fn, "r");

	debug("Starting to read '%s'...", fn);
	if (FH == NULL) fatal(ENOENT, "Couldn't open file %s for reading", fn);
       
	parse_header(FH, b);
	alloc_blocks(b);
	parse_blocks(FH, b);

	int chain_count[b->blocks];
	create_groups(b, chain_count);
	alloc_beads(b);
	parse_beads(FH, b, chain_count);
	mdpd_init_nblist(b, b->mdpd);
#ifdef MCOM
	mcom_init_nblist(b, b->mcom);
#endif
	debug("done!");
	fclose(FH);
}

/* ------------------------------------------------------------------------- */

static PASSPORT *passport_out = NULL;
static VEC2 *xv_out = NULL;
static FILE *mFH;
static char *vbuf = NULL;

static void populate_xvout(struct beads *b)
{
	int i, d;
	VEC2 *xv = b->xv + col_index * b->groupsize;
	VEC *n = b->nx + col_index * b->groupsize;
	PASSPORT *passport = b->passport + col_index * b->groupsize;

	if (unlikely(xv_out == NULL))
		 fatal(EINVAL, "Can't save -- xv_out not alloc'd");

	memset(xv_out, 0, b->nN * sizeof(*xv_out));
	memset(passport_out, 0, b->nN * sizeof(*passport_out));
	for (i = 0; i < b->groupsize; i++) {
		int id = GET_ID(passport[i]);
		if (!GET_EXISTS(passport[i])) continue;

		for (d = 0; d < 3; d++) {
			xv_out[id][d] = xv[i][d] + n[i][d] * b->l[d];
			xv_out[id][d + 3] = xv[i][d + 3];
		}
		passport_out[id] = passport[i];
	}

	MPI_Reduce(ismaster ? MPI_IN_PLACE : xv_out, ismaster ? xv_out : NULL,
		ARRAY_SIZE(*xv) * b->nN, MPI_DOUBLE, MPI_SUM, 0, comm_grid);
	MPI_Reduce(ismaster ? MPI_IN_PLACE : passport_out,
					ismaster ? passport_out : NULL, b->nN,
					MPI_LONG_LONG, MPI_BOR, 0, comm_grid);
}

/*
 * print the whole system to an open file handle (only master process)
 */
static void print_system(FILE *FH, const struct beads *b)
{
	int i, j, o;

	fprintf(FH, "# L=%lg %lg %lg t=%.3f blocks=%d\n",
			b->l[0], b->l[1], b->l[2], b->time, b->blocks);
	mdpd_print_header(FH, b->mdpd); 
#ifdef MCOM
	mcom_print_header(FH, b->mcom);
#endif
	RigidPrintHeader(FH, b);
	PepPrintHeader(FH, b);
	ExtPrintHeader(FH, b);
        surfpot_print_header(FH, b->surfpot);
        surfbump_print_header(FH, b->surfbump);
        surfcuff_print_header(FH, b->surfcuff);
        usp_print_header(FH, b->usp);
        multiblockforce_print_header(FH, b->multiblockforce);

	for (i = 0, o = 0; i < b->blocks; o += b->n[i] * b->N[i], i++) {
		VEC2 *x = xv_out + o;
		PASSPORT *p = passport_out + o;

		fprintf(FH, "# n=%d N=%d name=%s\n",
			b->n[i], b->N[i], b->resname[i]);

		for (j = 0; j < b->n[i] * b->N[i]; j++)
			fprintf(FH, "%lg %lg %lg %lg %lg %lg %llu\n",
					x[j][0], x[j][1], x[j][2], x[j][3],
					x[j][4], x[j][5], GET_TYPE(p[j]));
	}
}

void save_beads(const char *fn, struct beads *b)
{
	populate_xvout(b);

	if (ismaster) {
		FILE *FH = fopen(fn, "w");
		if (FH == NULL)
			fatal(EIO, "Couldn't open %s for writing", fn);

		print_system(FH, b);
		fclose(FH);
		printf("File '%s' successfully written to disk.\n", fn);
	}
}

static void save_beads_frequently(struct beads *b, const int snaphots,
							size_t vbufsize)
{
	populate_xvout(b);
	
	if (ismaster) {
		static int count = 0;

		if (mFH == NULL) {
			char fn[PATH_MAX];
			sprintf(fn, "output%09lu.dat", b->step);
			mFH = fopen(fn, "w");
			if (mFH == NULL)
				fatal(EIO, "Couldn't open '%s' for "
							"writing", fn);
			setbuffer(mFH, vbuf, vbufsize);
		}

		print_system(mFH, b);
		
		if (++count >= snaphots) {
			count = 0;
			fclose(mFH);
			mFH = NULL;
		}
	}
}

static void save_beads_init(struct beads *restrict b, size_t vbufsize)
{
	xv_out = calloc(b->nN, sizeof(*xv_out));
	if (xv_out == NULL) novm("xv_out");

	passport_out = calloc(b->nN, sizeof(*passport_out));
	if (passport_out == NULL) novm("passport_out");

	if (ismaster) {
		vbuf = malloc(vbufsize);
		if (vbuf == NULL) novm("vbuf");
	}
}

static void save_beads_free(void)
{
	free(xv_out);
	free(passport_out);
	if (mFH) fclose(mFH);
	if (vbuf) free(vbuf);
}

/* ------------------------------------------------------------------------- */

/*
 * Prints out non-bonded, bonded, kinetic, and total energy.
 */
static void measure_energy(struct beads *restrict b)
{
	int i, dofs;

	for (i = 0; i < b->groupsize; i++)
		b->e[2] += b->v2[b->groupsize * col_index + i];
	b->e[2] *= .5;
		
	MPI_Reduce(ismaster ? MPI_IN_PLACE : b->e, ismaster ? b->e : NULL,
			ARRAY_SIZE(b->e), MPI_DOUBLE, MPI_SUM, 0, comm_grid);

	if (ismaster) {
		for (i = 0, dofs = 0; i < b->blocks; i++)
			dofs += 3 * b->n[i] * b->N[i];
		printf("time=%lg inter=%lg intra=%lg kin=%lg sum=%lg <kT>=%lg "
			"US=%lg\n", b->time, b->e[0], b->e[1], b->e[2],
			b->e[0] + b->e[1] + b->e[2], (2.*b->e[2])/dofs,
			b->e[3]);
	}
}

/*
 * Print out diagonal entries of the pressure tensor.
 */
static void measure_pressure(struct beads *restrict b)
{
	TENSOR p;
	get_pressure(b, &p);
	if (ismaster)
		printf("l=%lg %lg %lg P=%lg %lg %lg %lg %lg %lg\n",
			       	b->l[0], b->l[1], b->l[2],
			       	p[0], p[1], p[2], p[3], p[4], p[5]);
	double V = b->l[0] * b->l[1] * b->l[2];
}

/* print summed forces of each block */
static void block_forces(struct beads *b)
{
	int i, j, k, idx;
	VEC *f = b->f + b->groupsize * col_index;
	VEC fa[b->blocks];

	for (k = 0; k < b->blocks; k++)
		fa[k][0] = fa[k][1] = fa[k][2] = 0.;
	for (i = 0, idx = 0; i < b->local_n; i++) {
		k = b->local_n_b[i];
		for (j = 0; j < b->N[k]; j++, idx++) {
			fa[k][0] += f[idx][0];
			fa[k][1] += f[idx][1];
			fa[k][2] += f[idx][2];
		}
	}
	MPI_Reduce(ismaster ? MPI_IN_PLACE : fa, ismaster ? fa : NULL,
		b->blocks * ARRAY_SIZE(*f), MPI_DOUBLE, MPI_SUM, 0, comm_grid);

	if (ismaster) {
		printf("block forces:");
		for (k = 0; k < b->blocks; k++)
			printf(" %lg %lg %lg", fa[k][0], fa[k][1], fa[k][2]);
		printf("\n");
	}
}

static void measure(struct beads *restrict b)
{
	measure_energy(b);
	measure_pressure(b);
	rnemd_print(b);
	mdpd_measure_energy(b->mdpd);
	mdpd_sg_measure(b, b->mdpd);
	block_forces(b);
	usv_printstatus(b->usv);
        usp_eval_print(b,b->usp);
}

/* ------------------------------------------------------------------------- */

enum taskname
{
	MEASUREMENT = 1,
	VMD,
	SAVE_BEADS,
	MEASURE_ISF,
	SAVE_ISF,
	STRESS_AUTOCORR,
	RNEMD_SWEEP,
#ifdef MCOM
	PRINT_COMS,
#endif
	TENS_SAVE,
	TENS_MEASURE,
	MDPD_SG_MOVE,
};

struct task
{
	enum taskname taskname; 
	unsigned long nextstep;
	unsigned long dstep;
	char *filename;
	double dt;
	int snapshots;
	size_t bufsize;
};

#define TASKS_MAX 64
static struct task tasks[TASKS_MAX];
static int num_tasks = 0;

static int cmptasks(const void *x, const void *y)
{
	struct task *a = (struct task *)x, *b = (struct task *)y;
	return a->nextstep - b->nextstep;
}

static void task_perform(struct beads *restrict b)
{
	int i = 0;
	char fn[PATH_MAX];

	TensDisable();
	
	while (tasks[i].nextstep == b->step && tasks[i].taskname != 0) {
		switch (tasks[i].taskname) {
		case VMD:
		        vmd_append(b);
			break;
		case SAVE_BEADS:
			save_beads_frequently(b, tasks[i].snapshots, tasks[i].bufsize);
			break;
		case MEASURE_ISF:
			dyn_measure(b);
			break;
		case SAVE_ISF:	
			sprintf(fn, "output%09lu-isf.dat", b->step);
			dyn_save(fn, b);
			break;
		case MEASUREMENT:
			measure(b);
			break;
		case TENS_SAVE:
		        TensWrite(b);
			break;
		case TENS_MEASURE:
			TensEnable(b);
			break;
		case STRESS_AUTOCORR:
			stress_autocorr_measure(b);
			break;
		case RNEMD_SWEEP:
			rnemd_sweep(b);
			break;
#ifdef MCOM
		case PRINT_COMS:
			mcom_print(b);
			break;
#endif
		case MDPD_SG_MOVE:
			mdpd_sg_move(b, b->mdpd);
			break;
		}

		tasks[i].nextstep += tasks[i].dstep;
		i++;
	}

	if (i > 0)
		qsort(tasks, num_tasks, sizeof(*tasks), cmptasks);
}

static void task_alloc(struct beads *restrict b, cfg_t *cfg)
{
	int i;
	memset(tasks, 0, ARRAY_SIZE(tasks) * sizeof(*tasks));
		
	for (i = 0; i < cfg_size(cfg, "task") && i < TASKS_MAX; i++) {
		cfg_t *t = cfg_getnsec(cfg, "task", i);
		const char *title = cfg_title(t);
		double dt = cfg_getfloat(t, "dt");
		tasks[i].dt = dt;
		tasks[i].filename = cfg_getstr(t, "filename");
		const char *fn = tasks[i].filename;
		tasks[i].bufsize = cfg_getint(t, "bufsize");
		tasks[i].snapshots = cfg_getint(t, "snapshots");

		if (strcasecmp(title, "vmd") == 0) {
			tasks[i].taskname = VMD;
			/* TODO: make VMD use a variable output buffer */
			vmd_alloc(b, fn);
			vmd_append(b);
			debug("task VMD: fn=%s dt=%lg", fn, dt);
			
		} else if (strcasecmp(title, "save_beads") == 0) {
			tasks[i].taskname = SAVE_BEADS;
			save_beads_init(b, tasks[i].bufsize);
			debug("task SAVE_BEADS: fn=%s dt=%lg", fn, dt);
			
		} else if (strcasecmp(title, "TensSave") == 0) {
		        tasks[i].taskname = TENS_SAVE;
			debug("task TENS_SAVE dt=%lg", dt);

		} else if (strcasecmp(title, "TensMeasure") == 0) {
			tasks[i].taskname = TENS_MEASURE;
			debug("task TENS_MEASURE dt=%lg", dt);
			TensLoad(b, cfg);

		} else if (strcasecmp(title, "measure") == 0) {
			tasks[i].taskname = MEASUREMENT;
			debug("task MEASUREMENT: dt=%lg", dt);
			
		} else if (strcasecmp(title, "measure_isf") == 0) {
			tasks[i].taskname = MEASURE_ISF;
			dyn_init(tasks[i].dt, tasks[i].bufsize);
			dyn_load(fn, b);
			dyn_save("verify-isf.dat", b);
			debug("task MEASURE_ISF: fn=%s dt=%lg", fn, dt);
			
		} else if (strcasecmp(title, "save_isf") == 0) {
			tasks[i].taskname = SAVE_ISF;
			debug("task SAVE_ISF: fn=%s dt=%lg", fn, dt);

		} else if (strcasecmp(title, "stress_autocorr") == 0) {
			tasks[i].taskname = STRESS_AUTOCORR;
			debug("task STRESS_AUTOCORR: fn=%s", fn);
			stress_autocorr_alloc(fn, tasks[i].bufsize);

		} else if (strcasecmp(title, "rnemd") == 0) {
			tasks[i].taskname = RNEMD_SWEEP;
			debug("task RNEMD_SWEEP dt=%lg", dt);
			rnemd_init(b);
#ifdef MCOM
		} else if (strcasecmp(title, "print_coms") == 0) {
			tasks[i].taskname = PRINT_COMS;
			debug("task PRINT_COMS dt=%lg", dt);
			mcom_print_init(b, fn, tasks[i].bufsize);
#endif
		} else if (strcasecmp(title, "mdpd_sg_move") == 0) {
			tasks[i].taskname = MDPD_SG_MOVE;
			debug("task MDPD_SG_MOVE dt=%lg", dt);
			mdpd_sg_init(b, b->mdpd, cfg);

		} else
			fatal(EINVAL, "task '%s' is unknown", title);
		
		tasks[i].dstep = floor(tasks[i].dt / b->dt + .5);
		tasks[i].nextstep = floor(b->time / tasks[i].dt + .5) *
								tasks[i].dstep;
		if (tasks[i].nextstep < b->step)
			tasks[i].nextstep += tasks[i].dstep;

		num_tasks++;
	}
	qsort(tasks, num_tasks, sizeof(*tasks), cmptasks);
}

static void task_free(struct beads *restrict b)
{
	int i;

	for (i = 0; i < num_tasks; i++) {
		switch (tasks[i].taskname) {
		case VMD:
			vmd_free();
			break;
		case MEASURE_ISF:
			dyn_save("resume-isf.dat", b);
			dyn_free();
			break;
		case SAVE_BEADS:
		  save_beads("resume.dat", b);
			save_beads_free();
			break;
		case STRESS_AUTOCORR:
			stress_autocorr_free();
			break;
		case RNEMD_SWEEP:
			rnemd_print(b);
			rnemd_free();
			break;
#ifdef MCOM
		case PRINT_COMS:
			mcom_print_free();
			break;
#endif
		case TENS_MEASURE:
			TensFree(b);
			break;
		case MDPD_SG_MOVE:
			mdpd_sg_free(b->mdpd);
			break;
		}
	}
}

/* ------------------------------------------------------------------------- */

static void main_loop(struct beads *restrict b)
{
	const double time_end = b->time + b->delta_t;
	double accum = -1. * MPI_Wtime();

	forces_calc(b);
	for (; b->time < time_end && !terminated();
					b->step++, b->time = b->step * b->dt) {
		task_perform(b);
		
		vv_step1(b);
		RigidVv1(b);
		
		forces_calc(b);

		vv_step2(b);
		RigidVv2(b);
		//if(DiffFillProf(b)) break;
	}

	accum += MPI_Wtime();
	printf("Elapsed time for main_loop(): %lf\n", accum);

	if (terminated())
		fprintf(stderr, "Caught terminate signal. Exiting...\n");
}

int main(int argc, char **argv)
{
	const char svn[] = "$Id: main.c 366 2013-08-08 14:19:08Z fuhrmans $";
	int rank;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (rank == 0) {
		printf("FDMDPD2 Copyright 2008-2011 Martin Hoemberg and "
		"Giovanni Marelli\n%s, built %s %s\nThis program comes with "
		"ABSOLUTELY NO WARRANTY. This is free software,\nand you are "
		"welcome to redistribute it under certain conditions;\nsee "
		"COPYING for details.\n\n", svn, __DATE__, __TIME__);
	}

#if defined(DEBUG) && defined(FE_NOMASK_ENV)
	feenableexcept(FE_UNDERFLOW | FE_OVERFLOW | FE_INVALID | FE_DIVBYZERO);
#endif

	if (argc == 5) {
		int rows = atoi(argv[1]);
		int cols = atoi(argv[2]);
		char *configfile = argv[3];
		char *sysfile = argv[4];
		struct beads b = { 0 };

		debug("conf='%s' sys='%s' rows=%d cols=%d",
				configfile, sysfile, rows, cols);

		create_grid_topology(rows, cols);

		cfg_t *cfg = cfg_init(opts, CFGF_NOCASE);
		if (cfg_parse(cfg, configfile) != 0)
			fatal(EINVAL, "cfg_parse() failed");
	
		parse_mainoptions(&b, cfg);
		parse_systemfile(sysfile, &b);
		parse_thermostat(&b, cfg);
		vv_alloc(&b, cfg);
		task_alloc(&b, cfg);
		WidomLoad(&b, cfg);
		//DiffLoad(&b, cfg);
		usv_init(&b, cfg);

		install_sigterm_handler();
		//InitGraphics(argc,argv,&b);
		if(!WidomLoop(&b)) return EXIT_SUCCESS;
		if(!TensLoop(&b))  return EXIT_SUCCESS;
		main_loop(&b);
		task_free(&b);
		free_beads(&b);
		free_grid_topology();
		rng_free(rng);
	} else {
		fprintf(stderr, "Syntax: %s [rows] [cols] [config.txt]"
						" [resume.dat]\n\n", argv[0]);
	}

	MPI_Finalize();
	return EXIT_SUCCESS;
}

