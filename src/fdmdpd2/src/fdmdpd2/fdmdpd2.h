/*
 * fdmdpd2.h - header file for fdmdpd2 main program
 * Copyright (C) 2009-2010 Martin Hoemberg <mhoembe@gwdg.de>
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

/* $Id: fdmdpd2.h 366 2013-08-08 14:19:08Z fuhrmans $ */

#ifndef __FDMDPD2_H__
#define __FDMDPD2_H__

#include "common.h"
#include <mpi.h>
#include <confuse.h>

typedef double VEC[3];
typedef double VEC2[6];

#include "Rigid.h"

typedef unsigned long long PASSPORT;
typedef char RESNAME[9];

/* symmetric 3x3 tensor. The order of the entries is:
 * xx, yy, zz, xy, xz, yz */
typedef double TENSOR[6];

/* number of different bead species, currently only "A" and "B" */
#ifndef TYPE_MAX
#define TYPE_MAX 2
#endif

/*
 *  bit field with the particle properties
 *  ID_MASK: 24 bits for unique particle id
 *  POS_MASK: 24 bits for particle position in memory
 *  TYPE_MASK: 3 bits for particle type (A, B, C, ...)
 *  EXISTS_MASK: 1 bit, which is set unless this is a ghost particle, i.e.
 *               it is used to balance the number of beads per process.
 *  OK_MASK: this bit is only meaningful for the first bead of a pair in
 *           nblist->passport. If it's set, then this pair is within
 *           interaction range, i.e. dr < r_c.
 */ 
#define ID_MASK		(0xffffffLL <<  0)
#define POS_MASK	(0xffffffLL << 24)
#define TYPE_MASK	(0x7LL << 48)
#define EXISTS_MASK	(0x1LL << 51)
#define OK_MASK		(0x1LL << 52)
#define SET_ID(p,i)	((p)=((p)&(~ID_MASK))|((PASSPORT)(i)&(ID_MASK)))
#define GET_ID(p)	((p)&(ID_MASK))
#define SET_POS(p,i)	\
		((p)=((p)&(~POS_MASK))|(((PASSPORT)(i)<<24)&(POS_MASK)))
#define GET_POS(p)	(((p)&(POS_MASK))>>24)
#define SET_TYPE(p,i)	\
		((p)=((p)&(~TYPE_MASK))|(((PASSPORT)(i)<<48)&(TYPE_MASK)))
#define GET_TYPE(p)	(((p)&(TYPE_MASK))>>48)
#define SET_EXISTS(p)	((p)|=(EXISTS_MASK))
#define UNSET_EXISTS(p)	((p)&=(~EXISTS_MASK))
#define GET_EXISTS(p)	(!!((p)&(EXISTS_MASK)))
#define GET_OK(p)	(!!((p)&(OK_MASK)))
#define SET_OK(p)	((p)|=(OK_MASK))
#define UNSET_OK(p)	((p)&=(~OK_MASK))

/*
 * Orientation of the Coordinate System
 * Here a global bilayer normal direction, as well as two tangential directions
 * are explicitly designated. It is assumed in several places, that the
 * system under study possesses this normal vector to the bilayer.
 *
 * Although it is an int, only the 6 lowest bits are used. The lowest bits
 * store the index of the normal direction, the next two bits store the
 * index of the first tangential direction, the last two bits store the index
 * of the second tangential direction. 
 */ 
extern int bilayer_orientation;
#define SET_NORMAL(a)	\
		(bilayer_orientation = (bilayer_orientation & 0x3c) | (a << 0))
#define SET_TANG1(a)	\
		(bilayer_orientation = (bilayer_orientation & 0x33) | (a << 2))
#define SET_TANG2(a)	\
		(bilayer_orientation = (bilayer_orientation & 0x0f) | (a << 4))
#define NORMAL	((bilayer_orientation & 0x03) >> 0)
#define TANG1	((bilayer_orientation & 0x0c) >> 2)
#define TANG2	((bilayer_orientation & 0x30) >> 4)

/* shall we calculate the stress tensor for a given interaction? */
#define STRESS_OMIT 0
#define STRESS_CALC 1

/* forward declarations for the structs of the various force fields */
#ifdef MCOM
struct mcom;
#endif
struct mdpd;
struct surfpot;
struct surfbump;
struct surfcuff;
struct blockforce;
struct multiblockforce;
struct usv;

/*
 * struct beads - particle and simulation box properties
 * The origin is always in the lower left corner, i.e. the simulation box
 * runs from 0. to L in each direction.
 */
struct beads
{
	VEC2 *xv;		/* phase space, position and velocity */
	VEC *nx;		/* image counter for the pbc */
	VEC *f;			/* forces */
	VEC *x_intra;		/* intramolecular coordinates */
	VEC *f_intra;		/* intramolecular and external forces */
	RESNAME *resname;	/* residue names for VMD in each block */
	PASSPORT *passport;	/* particle passport. This is a bitfield! */
	int *local_n_b;		/* block indeces of the local chains */
	int *n;			/* # chains in each block */
	int *N;			/* beads per chain in each block */
	double *v2;		/* squared velocities of the particles */
        struct WIDOM Widom;     /* Widom insertion/deletion loop */
        struct TENS_PROF Tens;  /* Stress profile */
        struct DIFF_PROF Diff;  /* Diffusion profile */
        struct NANO *Nano;      /* Rigid external body */
        struct PEPTIDE *Pep;    /* Infromation about the architecture of the PEP block */
        struct PEPTIDE *Ext;    /* External rules to the system */
	struct mdpd *mdpd;	/* parameters for the MDPD force field */
#ifdef MCOM
	struct mcom *mcom;	/* parameters for the MCOM force field */
#endif
        struct surfpot *surfpot;/* surface potential parameters */
        struct surfbump *surfbump;/* surface potential with bump parameters */
        struct surfcuff *surfcuff;/* cylindrical cuff potential */
        struct blockforce *blockforce; /* parameters for blockforce */
        struct multiblockforce *multiblockforce; /* parameters for multiblockforce */
	struct usv *usv;	/* Umbrella Sampling for vesicles */
	struct usp *usp;	/* Umbrella Sampling for peptides */
	unsigned long step;	/* current intergration step */
	VEC l;			/* lengths of the simulation box */
	TENSOR virial;		/* interaction part of pressure tensor */
	double e[5];		/* inter- intramolecular, kinetic, US energy, external */
	double time;		/* system time */
	double delta_t;		/* time interval for the whole simulation */
	double dt;		/* integration time step */
	double v2max;		/* maximum v^2 in all processses */
	double rs;		/* radius of neighbor shell */
	int local_n;		/* # chains in this process */
	int nN;			/* # all particles in all chains (really!) */
	int blocks;		/* # blocks (different molecular species) */
	int groupsize;		/* number of particles in one group */
        int NNano;              /* Number of nanoparticles */
        int NPep;               /* Number of peptides */
        int NExt;               /* Number of additional interaction fields. */
};

/* ------------------------------------------------------------------------- */

enum nblist_flags
{
	/* flags for nbl_alloc */
	NBL_NOFLAGS = 0x0000,
	NBL_NEWTON3 = 0x0001,		/* use Newton's 3rd law */
	NBL_IGNORE_CC = 0x0002,		/* ignore C-C interactions */
	/* flags for nbl_update */
	NBL_VALID = 0x1000,		/* valid: list is still valid */
	NBL_INVALID = 0x1001,		/* invalid: particles moved too far */
	NBL_VOLUME_CHANGED = 0x1002,	/* invalid: volume has changed */
};

/*
 * struct cpa - cell pointer array
 * This is internally used for the neighbor finding algorithm
 */
struct cpa
{
	VEC2 *xv;
	const PASSPORT *pp;
	PASSPORT *C;
	int *Pc;
	int n;
};

/*
 * struct nblist - list of pairwise interactions
 * This struct facilitates a cache-friendly computation of non-bonded, pairwise
 * interactions. "xv", "f", and "passport" contain 2*"count" elements, so that
 * the corresponding information about one pair is stored sequentially in
 * memory.
 */
struct nblist
{
	VEC2 *drdv;		/* buffer with coord and velocity diffs */
	PASSPORT *passport;	/* passports of listed particles */
	int *Pm;
	void (*loop_col)(struct nblist *nbl, PASSPORT prow, int m1, int m2);
	struct cpa row;		/* cell pointer array for row particles */
	struct cpa col;		/* cell pointer array for column particles */
	double rc2;		/* effective cutoff (rs+RC)^2 */
	int s;			/* number of slices */
	int max;		/* max. number of particles in list */
	int count;		/* current number of particles in list */
	int max_M;		/* maximum number of cells */
};

/* ------------------------------------------------------------------------- */

/* declarations from dynamics.c */
extern void dyn_init(const double _delta_a, const size_t bufsize);
extern void dyn_load(const char *fn, struct beads *restrict b);
extern void dyn_save(const char *fn, struct beads *restrict b);
extern void dyn_measure(struct beads *restrict b);
extern void dyn_free(void);
extern cfg_opt_t dyn_opts[];

/* declarations from error.c */
#ifdef DEBUG
extern void debug(const char *fmt, ...) __attribute__((format(printf, 1, 2)));
#else
#define debug(...) {}
#endif
extern int get_data_line(FILE *restrict FH, const char *restrict fmt, ...)
					__attribute__((format(scanf, 2, 3)));
extern void novm(const char *fmt, ...);
extern void fatal(const int errnum, const char *fmt, ...)
				__attribute__((noreturn,format(printf, 2, 3)));

/* declarations from exforce.c */
extern void blockforce_print_header(FILE *FH, struct blockforce *blockforce);
extern void blockforce_read_header(FILE *FH, struct beads *b);
extern double apply_blockforce(struct beads *b, struct blockforce *blockforce);
extern void blockforce_free(struct blockforce *blockforce);
extern void multiblockforce_print_header(FILE *FH, struct multiblockforce *multiblockforce);
extern void multiblockforce_read_header(FILE *FH, struct beads *b);
extern double apply_multiblockforce(struct beads *b, struct multiblockforce *multiblockforce);
extern void multiblockforce_free(struct multiblockforce *multiblockforce);
extern void multiblockforce_wrapper(struct beads *b, struct multiblockforce *multiblockforce);
extern void multiblockforce_get_coms(struct beads *b, struct multiblockforce *multiblockforce);
extern void multiblockforce_apply_forces(struct beads *b, struct multiblockforce *multiblockforce);

/* declarations from main.c */
extern void alloc_beads(struct beads *restrict b);
extern void free_beads(struct beads *restrict b);
extern void alloc_blocks(struct beads *restrict b);
extern void set_bead(struct beads *restrict b, VEC2 xv, int id, int pos, int type);
extern void bcast_beads(struct beads *restrict b);

extern MPI_Comm comm_grid;	/* MPI communicator for all processes */
extern MPI_Comm comm_row;	/* MPI communicator for the current row */
extern MPI_Comm comm_col;	/* MPI communicator for the current column */
extern int groups_per_row;	/* number of groups in one row */
extern int groups_per_col;	/* number of groups in one column */
extern int row_index;		/* row index of this process */
extern int col_index;		/* column index of this process */
extern int ismaster;		/* is this the master process? */
extern struct rng *rng;		/* random number generator for thermostats */

#ifdef MCOM
/* declarations from mcom.c */
extern void mcom_print_init(struct beads *b, const char *fn, size_t bufsize);
extern void mcom_print_free(void);
extern void mcom_print(struct beads *b);
extern void mcom_calc_forces(struct beads *b, struct mcom *mcom, enum nblist_flags fl);
extern void mcom_init(FILE *FH, struct beads *b);
extern void mcom_free(struct mcom *mcom);
extern void mcom_init_nblist(struct beads *b, struct mcom *mcom);
extern void mcom_print_header(FILE *FH, struct mcom *mcom);
extern double mcom_get_rc(struct mcom *mcom);
extern int mcom_get_diml(struct mcom *mcom);
extern void mcom_calc_resulting_forces(struct beads *b, VEC *b_f);
extern void mcom_calc_e2e_vector(struct beads *b, VEC *b_r);
extern void mcom_calc_com_xv(struct beads *b, VEC2 *b_xv);
extern void mcom_parse_thermostat(struct mcom *mcom, cfg_t *cfg);
#endif

/* declarations from mdpd.c */
extern void mdpd_linear_chain(struct beads *b, struct mdpd *mdpd, const int N, const int bp);
extern void mdpd_two_tails(struct beads *b, struct mdpd *mdpd, const int N, const int bp);
extern void mdpd_nonbonded_forces(struct beads *b, struct mdpd *mdpd, enum nblist_flags fl);
extern void mdpd_free(struct mdpd *mdpd);
extern void mdpd_measure_energy(struct mdpd *mdpd);
extern void mdpd_init(FILE *FH, struct beads *b);
extern void mdpd_init_nblist(struct beads *b, struct mdpd *mdpd);
extern void mdpd_print_header(FILE *FH, struct mdpd *mdpd);
extern void mdpd_parse_thermostat(struct mdpd *mdpd, cfg_t *cfg);
extern cfg_opt_t us_opts[];
extern double mdpd_getrhocoex(struct mdpd *mdpd);
extern double mdpd_getRe(struct mdpd *mdpd);
extern double harmonic_spring(struct beads *b, int bp, int i, int j, double ks, double l0, int calcstress);
extern double bond_angle(struct beads *b, int bp, int i, int j, int k, double kb);
extern void mdpd_sg_init(struct beads *b, struct mdpd *mdpd, cfg_t *cfg);
extern void mdpd_sg_free(struct mdpd *mdpd);
extern void mdpd_sg_move(struct beads *b, struct mdpd *mdpd);
extern void mdpd_sg_measure(struct beads *b, struct mdpd *mdpd);
extern cfg_opt_t mdpd_sg_opts[];

/* declarations from nblist.c */
extern void nblist_alloc(const struct beads *b, struct nblist *nbl,
		int row_n, PASSPORT row_pp[const], VEC2 row_xv[const],
		int col_n, PASSPORT col_pp[const], VEC2 col_xv[const],
		double rc, enum nblist_flags flags);
extern void nblist_update(struct beads *b, struct nblist *nbl, enum nblist_flags f);
extern void nblist_rebuild(struct nblist *restrict nbl);
extern void nblist_free(struct nblist *restrict nbl);
extern cfg_opt_t nblist_opts[];

/* declarations from rnemd.c */
extern void rnemd_init(struct beads *b);
extern void rnemd_free(void);
extern void rnemd_sweep(struct beads *b);
extern void rnemd_print(struct beads *b);
extern cfg_opt_t rnemd_opts[];

/* declarations from stress.c */
extern void get_pressure(const struct beads *restrict b, TENSOR *p);
extern void stress_autocorr_alloc(const char *fn, size_t bufsize);
extern void stress_autocorr_free(void);
extern void stress_autocorr_measure(struct beads *restrict b);
extern void stress_autocorr_print(const char *fn, const double dt);

/* declarations from surfpot.c */
extern void surfpot_print_header(FILE *FH, struct surfpot *surfpot);
extern void surfpot_read_header(FILE *FH, struct beads *b);
extern void surfpot_calc_forces(struct beads *b, struct surfpot *surfpot);
extern void surfpot_free(struct surfpot *surfpot);
extern void surfbump_print_header(FILE *FH, struct surfbump *surfbump);
extern void surfbump_read_header(FILE *FH, struct beads *b);
extern void surfbump_calc_forces(struct beads *b, struct surfbump *surfbump);
extern void surfbump_free(struct surfbump *surfbump);
extern void print_surfbump_forces(struct beads *b, struct surfbump *surfbump);
extern void surfcuff_print_header(FILE *FH, struct surfcuff *surfcuff);
extern void surfcuff_read_header(FILE *FH, struct beads *b);
extern void surfcuff_calc_forces(struct beads *b, struct surfcuff *surfcuff);
extern void surfcuff_free(struct surfcuff *surfcuff);
extern void surfcuff_print_forces(struct beads *b, struct surfcuff *surfcuff);

/* declarations from usv.c */
extern void usv_printstatus(const struct usv *u);
extern void usv_force(struct beads *b, struct usv *u);
extern void usv_init(struct beads *b, cfg_t *cfg);
extern void usv_free(struct beads *b);
extern cfg_opt_t usv_opts[];

/* declarations from usp.c */
extern void usp_print_header(FILE *FH, struct usp *usp);
extern void usp_read_header(FILE *FH, struct beads *b);
extern void usp_free(struct usp *usp);
extern void usp_wrapper(struct beads *b, struct usp *usp);
extern void usp_get_coms(struct beads *b, struct usp *usp);
extern void usp_apply_forces(struct beads *b, struct usp *usp);
extern void usp_eval_print(struct beads *b, struct usp *usp);

/* declarations from vmd.c */
extern void vmd_alloc(struct beads *restrict b, const char *fn);
extern void vmd_append(struct beads *restrict b);
extern void vmd_free(void);

/* declarations from vv.c */
extern void vv_alloc(struct beads *restrict b, cfg_t *cfg);
extern void vv_step1(struct beads *restrict b);
extern void vv_step2(struct beads *restrict b);
extern cfg_opt_t integrator_opts[];

/* declarations from main.c */
extern void forces_calc(struct beads *b);

/* declarations from Tens.c */
extern cfg_opt_t Tens_opts[];
extern void TensLoad(struct beads *restrict b,cfg_t *cfg);
extern void (*TensStress)(struct beads *b,double Force,PASSPORT prow,PASSPORT pcol);
extern void (*TensStressPos)(struct beads *b, double ExtForce,double *Pos1,double *Pos2,int t1,int t2);
extern void (*TensStressPre)(struct beads *b,double *PreExt,PASSPORT prow,PASSPORT pco);
extern void TensWrite(struct beads *b);
extern void TensFree(struct beads *b);
extern void TensDisable(void);
extern void TensEnable(struct beads *b);
extern int TensLoop(struct beads *b);

/* declarations from Diff.c */
extern cfg_opt_t Diff_opts[];
extern void DiffLoad(struct beads *restrict b,cfg_t *cfg);
extern int DiffFillProf(struct beads *restrict b);
extern void DiffFree(struct beads *restrict b);

/* declarations from Rigid.c */
extern void RigidVv1(struct beads *restrict b);
extern void RigidLoad(struct beads *restrict b,FILE *FH);
extern void RigidVv2(struct beads *restrict b);
extern void RigidCalcForces(struct beads *restrict b);
extern void RigidPrintHeader(FILE *FH, const struct beads *b);
extern void RigidFree(struct beads *b);
extern void RigidChNrg(struct beads *restrict b,int Ch);
extern void RigidMonNrg(struct beads *restrict b,int Ch);

/* declaration from Peptide.c */
extern void PepLoad(struct beads *b, FILE *FH);
extern void PepFree(struct beads *b);
extern void PepCalcForces(struct beads *restrict b,int MemPos,const char *name,int NPart);
extern void PepPrintHeader(FILE *FH, const struct beads *b);
extern int __attribute__((format(scanf,3,4))) Fetch(const char *str,const char *mask,const char *fmt, ... );

/* declaration from Extern.c */
extern void ExtLoad(struct beads *b, FILE *FH);
extern void ExtFree(struct beads *b);
extern void ExtCalcForces(struct beads *restrict b);
extern void ExtPrintHeader(FILE *FH, const struct beads *b);

/* in Widom.c */
extern cfg_opt_t Widom_opts[];
//extern void WidomLoopRem(struct beads *restrict b);
//extern void WidomLoad(struct beads *b,cfg_t *cfg);
//extern void WidomNrg(struct beads *b);
extern int WidomLoop(struct beads *b);
extern void WidomFree(struct beads *b);
extern int MCStepIn(struct beads *b);
extern int MCStepOut(struct beads *b);
#endif

