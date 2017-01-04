/*
 * mdpd.c - The MDPD force field
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

/* $Id: mdpd.c 366 2013-08-08 14:19:08Z fuhrmans $ */

#include "fdmdpd2.h"
#include "rand.h"

/* 
 * Unfortunately it is difficult to generalize all expressions in this file
 * to the case of an arbitrary number of coarse-grained species. Therefore
 * only the two most important cases of two (A, B) and three species (A, B, C)
 * are treated here. Generalizations are welcome!
 */
#if ((TYPE_MAX != 2) && (TYPE_MAX != 3))
#error "Only TYPE_MAX=2 or 3 is allowed at the moment"
#endif

#define VIR2(a,b)	((a) * TYPE_MAX + (b))
#define VIR3(a,b,c)	(((a) * TYPE_MAX + (b)) * TYPE_MAX + (c))

/*
 * cache-blocking size used in calculation of non-bonded interactions.
 * According to M. Puetz from IBM is STEP=16 the best choice for POWER6.
 * However, STEP=64 works well on Intel Xeon at the HLRN-II facility.
 */
#define STEP 16

/*
 * weighting functions. They must be selected at compile-time, so that inlining
 * can be performed. It seems that one really needs different weighting
 * functions for the 2nd and the 3rd order interactions, to get meaningful
 * physics. 
 */
#define W_FUNC2(r,a)	w_cubicspline((r),(a))		/* 2nd order */
#define DW_FUNC2(r,a)	dw_cubicspline((r),(a))
#ifdef POLYMER_SYS
#define W_FUNC3(r,b)	w_cubicspline((r),(b))		/* 3rd order */
#define DW_FUNC3(r,b)	dw_cubicspline((r),(b))
#else
#define W_FUNC3(r,b)	w_quadratic((r),(b))		/* 3rd order */
#define DW_FUNC3(r,b)	dw_quadratic((r),(b))
#endif

/* ------------------------------------------------------------------------- */

cfg_opt_t mdpd_sg_opts[] = {
	CFG_FLOAT_LIST("fugacity", "{}", CFGF_NONE),
	CFG_STR_LIST("str", "{}", CFGF_NONE),
	CFG_END()
};

typedef double DENSITY[2 * TYPE_MAX];	/* 2nd order, 3rd order */

/*
 * All interaction coefficients of the MDPD force field are bundled in this
 * struct. For brevity the neighbor list is also contained in this struct.
 */ 
struct mdpd
{
	DENSITY *rho;		/* weighted densities */
	double gamma[SQR(TYPE_MAX)];	/* DPD friction constant */
	double zeta[SQR(TYPE_MAX)];	/* and noise coefficient */
	double gamma_p[SQR(TYPE_MAX)];	/* Transversal DPD friction constant */
	double zeta_p[SQR(TYPE_MAX)];	/* and noise coefficient */
	double lgamma[3 * TYPE_MAX];	/* Langevin friction constant */
	double lzeta[3 * TYPE_MAX];	/* and noise coefficient */
	double a2;		/* junction point for 2nd order wght. func */
	double a3;		/* junction point for 3rd order wght. func */
	double Re;		/* rms chain length in units of the cutoff */
	double ks;		/* spring constant for harmonic potential */
	double kb;		/* spring constant for bond angle potential */
	double l0;		/* minimum of harmonic spring potential */
	double v[SQR(TYPE_MAX)]; /* 2nd order "virial" coefficients */
	double w[CUBE(TYPE_MAX)]; /* 3rd order "virial" coefficients */
	double rho2[SQR(TYPE_MAX)]; /* integrated densities 2nd order */
	double rho3[CUBE(TYPE_MAX)]; /* integrated densities 3rd order */
	int Re_N;		/* # beads in a chain of 1 Re */
	struct nblist nbl;	/* neighbor list */
	int sg_N;		/* # alternative configurations */
	struct rng *sg_rng;	/* rng for semigrand ensemble */
	char **sg_str;		/* alternative configurations */
	double *sg_xi;		/* fugacity fractions, Frenkel Eq. (9.1.15) */
};

static void (*calc_nonbonded_forces)(struct beads *b, struct mdpd *mdpd) = NULL;

/* ------------------------------------------------------------------------- */

/*
 * a single harmonic spring between beads "i" and "j"
 */
extern double harmonic_spring(struct beads *b, int bp, int i, int j, double ks,
					       	double l0, int calcstress)
{
	PASSPORT *p0 = b->passport + col_index * b->groupsize + bp;
	VEC *x = b->x_intra + bp;
	VEC *f = b->f_intra + bp;
	VEC dr;
	double force, dr2 = 0., energy, r;
	int d;
	
	for (d = 0; d < 3; d++) {
		dr[d] = x[j][d] - x[i][d];
		dr2 += SQR(dr[d]);
	}
	r = sqrt(dr2);
	assert(r < 20.); /* this number is somewhat arbitrary */

	energy = .5 * ks * SQR(r - l0);
	force = -ks * (1. - l0 / r);
	if (calcstress == STRESS_CALC) {
               for (d = 0; d < 3; d++) {
                       f[i][d] -= force * dr[d];
                       f[j][d] += force * dr[d];
                       b->virial[d] += SQR(dr[d]) * force;
               }
               b->virial[3] += dr[0] * dr[1] * force;
               b->virial[4] += dr[0] * dr[2] * force;
               b->virial[5] += dr[1] * dr[2] * force;

               TensStress(b, force * r, p0[i], p0[j]);
       } else {
               for (d = 0; d < 3; d++) {
                       f[i][d] -= force * dr[d];
                       f[j][d] += force * dr[d];
               }
	}
	return energy;
}

/*
 * a single bond angle potential (3 particles needed: "i", "j", "k")
 */
extern double bond_angle(struct beads *b, int bp, int i, int j, int k, double kb)
{
	PASSPORT *p0 = b->passport + col_index * b->groupsize + bp;
	VEC *x = b->x_intra + bp;
	VEC *f = b->f_intra + bp;
	VEC d0, d1;
	double n0, n1;
	double lim_n0, lim_n1;
	double ct = 0.;
	double energy;
	const double max_f = 500.; /* FIXME: this shouldn't be hard-coded */
	int d;
	double Pre[12];

	for (d = 0; d < 3; d++) {
		d0[d] = x[j][d] - x[i][d];
		d1[d] = x[k][d] - x[j][d]; 
	}
	n0 = sqrt(SQR(d0[0]) + SQR(d0[1]) + SQR(d0[2]));
	n1 = sqrt(SQR(d1[0]) + SQR(d1[1]) + SQR(d1[2]));
	for (d = 0; d < 3; d++) {
		d0[d] /= n0;
		d1[d] /= n1;
		ct += d0[d] * d1[d];
	}
	energy = kb * (1. - ct);

	/* we limit the forces to max_f */
	lim_n0 = MAX(n0, 2. * kb / max_f);
	lim_n1 = MAX(n1, 2. * kb / max_f);
	for (d = 0; d < 3; d++) {
		f[i][d] -= kb * (d1[d] - d0[d] * ct) / lim_n0;
		f[j][d] += kb * (d1[d] - d0[d] * ct) / lim_n0;
		Pre[d]   = d0[d] * kb * (d1[d] - d0[d] * ct);
		b->virial[d] += Pre[d];

		f[j][d] -= kb * (d0[d] - d1[d] * ct) / lim_n1;
		f[k][d] += kb * (d0[d] - d1[d] * ct) / lim_n1;
		Pre[d + 6] = d1[d] * kb * (d0[d] - d1[d] * ct);
		b->virial[d] +=  Pre[d + 6];
	}

	Pre[3] = d0[0] * kb * (d1[1] - d0[1] * ct);
	b->virial[3] += Pre[3];
	Pre[9] = d1[0] * kb * (d0[1] - d1[1] * ct);
	b->virial[3] += Pre[9];
	Pre[4] = d0[0] * kb * (d1[2] - d0[2] * ct);
	b->virial[4] += Pre[4];
	Pre[10] = d1[0] * kb * (d0[2] - d1[2] * ct);
	b->virial[4] += Pre[10];
	Pre[5] = d0[1] * kb * (d1[2] - d0[2] * ct);
	b->virial[5] += Pre[5];
	Pre[11] = d1[1] * kb * (d0[2] - d1[2] * ct);
	b->virial[5] += Pre[11];

	TensStressPre(b, Pre, p0[i], p0[j]);
	TensStressPre(b, Pre + 6, p0[j], p0[k]);

	return energy;
}

/*
 * Linear lipid molecule (default architecture)
 */
void mdpd_linear_chain(struct beads *b, struct mdpd *mdpd, int N, int bp)
{
	if (mdpd == NULL) return;
	int j;
	for (j = 0; j < N - 1; j++)
		b->e[1] += harmonic_spring(b, bp, j, j + 1, mdpd->ks, mdpd->l0, STRESS_CALC);
	for (j = 1; j < N - 1; j++)
		b->e[1] += bond_angle(b, bp, j - 1, j, j + 1, mdpd->kb);
}

/*
 * Two-tailed lipids with the following architecture (here: N_A = 8, N = 11):
 * 0A-1A-2A-3A-9B-10B  (IfAA = 1)
 *             |
 * 4A-5A-6A-7A-8B 
 * 0A-1A-2A-3A         (IfAA = -1)
 *          |
 * 4A-5A-6A-7A-8B-9B-10B
 */
void mdpd_two_tails(struct beads *b, struct mdpd *mdpd, int N, int bp)
{
  PASSPORT *p0 = b->passport + col_index * b->groupsize + bp;
  int j, lim, IfAA=-1;

  if (mdpd == NULL) return;

  for (lim = 1; lim < N; lim++)
    if (GET_TYPE(p0[lim - 1]) == 0 && GET_TYPE(p0[lim]) == 1)
      break;
  // short backbone
  for (j = 0; j < lim/2; j++)
    b->e[1] += harmonic_spring(b, bp, j, j+1, mdpd->ks, mdpd->l0, STRESS_CALC);
  // long backbone
  for (j = lim/2; j < N - 1; j++)
    b->e[1] += harmonic_spring(b, bp, j, j+1, mdpd->ks, mdpd->l0, STRESS_CALC);
  // join the tails
  b->e[1] += harmonic_spring(b, bp,lim/2-1,lim+IfAA, mdpd->ks, mdpd->l0, STRESS_CALC);
  for (j = 1; j < lim / 2 - 1; j++)
    b->e[1] += bond_angle(b, bp, j - 1, j, j + 1, mdpd->kb);
  for (j = lim / 2 + 1; j < N - 1 ; j++)
    b->e[1] += bond_angle(b, bp, j - 1, j, j + 1, mdpd->kb);
}
/* ------------------------------------------------------------------------- */

/*
 * weighting function: we use a constant part for 0<r<a and and a cubic spline
 * in the domain a<r<1.
 */
static double w_cubicspline(const double r, const double a)
{
	const double b = 1.;
	double Ai = .5 * 15. / M_PI / (CUBE(a) + CUBE(a + b) + CUBE(b));
	return (r < a) ? Ai : (Ai / CUBE(a - b)) *
		((((-2. * r) + 3. * (a + b)) * r - 6. * a * b) * r +
				 3. * a * b * b - b * b * b);
}

/*
 * derivative of the weighting function: dw/dr
 */
static double dw_cubicspline(const double r, const double a)
{
	const double b = 1.;
	double Ai = 15. / (M_PI * (((((4. * a) + 6.) * a) + 6.) * a + 4.));
	return (r > a && r < b) ? (6. * Ai / CUBE(a - b) * (-r * r + (a + b) *
			       			r - a * b)) : 0.; 
}

/*
 * quadratic weighting function. The Fourier transform of this weighting
 * function is stricly positive, i.e. we are avoiding cluster crystals.
 */
static double w_quadratic(const double r, const double b)
{
	return (r < b) ? 15. / (2. * M_PI * CUBE(b)) * SQR(1. - r / b) : 0.;
}

/*
 * derivative of the quadratic weighting function
 */
static double dw_quadratic(const double r, const double b)
{
	return (r < b) ? 15. / (M_PI * SQR(SQR(b))) * (r / b - 1.) : 0.;
}

/*
 * calculate weighted density by iterating over all interacting pairs in this
 * process.
 */
static void weight_density(const struct beads *b, struct mdpd *mdpd)
{
	const VEC lh = {.5 * b->l[0], .5 * b->l[1], .5 * b->l[2]};
	struct nblist *nbl = &mdpd->nbl;
	int i, d;
	
	for (i = 0; i < nbl->count; i += STEP) {
		int j, n, idx[STEP];
		int p_row[STEP], p_col[STEP];
		double dr2[STEP];
		VEC2 *drdv = nbl->drdv + i;

		for (j = 0, n = 0; (j < STEP) && (i + j < nbl->count); j++) {
			VEC tmp;
			
			p_row[j] = GET_POS(nbl->passport[2 * (i + j) + 0]);
			p_col[j] = GET_POS(nbl->passport[2 * (i + j) + 1]);
			
			/* compute dr vector and minimum image convention */
			dr2[j] = 0.;
			for (d = 0; d < 3; d++) {
				tmp[d] = b->xv[p_col[j]][d]-b->xv[p_row[j]][d];
				tmp[d] += (tmp[d] < -lh[d]) ? b->l[d] : 0.;
				tmp[d] -= (tmp[d] >  lh[d]) ? b->l[d] : 0.;
				dr2[j] += SQR(tmp[d]);
			}

			/* if this is within interaction range, write back
			 * dr and dv vectors
			 */
			if (dr2[j] < 1.) {
				idx[n++] = j;
				SET_OK(nbl->passport[2 * (i + j)]);
				for (d = 0; d < 3; d++)
					drdv[j][d] = tmp[d];
				for (d = 3; d < 6; d++)
					drdv[j][d] = b->xv[p_col[j]][d] -
							b->xv[p_row[j]][d];
			} else {
				UNSET_OK(nbl->passport[2 * (i + j)]);
			}
		}

		/*
		 * calculate weighted densities
		 */
		for (j = 0; j < n; j++) {
			int k = idx[j];
			int t_row = GET_TYPE(nbl->passport[2 * (i + k) + 0]);
			int t_col = GET_TYPE(nbl->passport[2 * (i + k) + 1]);
			double w2 = W_FUNC2(sqrt(dr2[k]), mdpd->a2);
			double w3 = W_FUNC3(sqrt(dr2[k]), mdpd->a3);
			mdpd->rho[p_row[k]][t_col] += w2;
			mdpd->rho[p_col[k]][t_row] += w2;
			mdpd->rho[p_row[k]][t_col + TYPE_MAX] += w3;
			mdpd->rho[p_col[k]][t_row + TYPE_MAX] += w3;
			//printf("%d-%d %lf %f %lf %lf \n",p_row[k],p_col[k],mdpd->rho[p_row[k]][t_col],mdpd->rho[p_col[k]][t_row],mdpd->rho[p_row[k]][t_col + TYPE_MAX],mdpd->rho[p_col[k]][t_row + TYPE_MAX]);
		}
	}
}

/* ------------------------------------------------------------------------- */

/* DPD weighting function. this is w^R(r). */
static double __attribute__((const)) dpd_wf(double r, double rc)
{
//	return 1.;		/* constant */
	return 1. - r / rc;	/* standard */ 
}

/*
 * This function calculates all non-bonded interactions.
 */
static void calc_nonbonded_forces_dpd(struct beads *b, struct mdpd *mdpd)
{
	struct nblist *nbl = &mdpd->nbl;
	int i, d;

	for (i = 0; i < nbl->count; i += STEP) {
		int k, j, n, idx[STEP];
		double dr2[STEP], sp[STEP], theta[STEP];
		int p_row[STEP], p_col[STEP], t_row[STEP], t_col[STEP];
		double rho[TYPE_MAX * STEP];
		VEC2 *drdv = nbl->drdv + i;

		for (j = 0, n = 0; (j < STEP) && (i + j < nbl->count); j++) {

			if (!GET_OK(nbl->passport[2 * (i + j) + 0]))
				continue;
			
			p_row[j] = GET_POS(nbl->passport[2 * (i + j) + 0]);
			t_row[j] = GET_TYPE(nbl->passport[2 * (i + j) + 0]);
			p_col[j] = GET_POS(nbl->passport[2 * (i + j) + 1]);
			t_col[j] = GET_TYPE(nbl->passport[2 * (i + j) + 1]);

			/* only the third order densities are needed here */
			for (d = 0; d < TYPE_MAX; d++)
				rho[TYPE_MAX * j + d] = /* repulsive only */
					mdpd->rho[p_row[j]][TYPE_MAX + d] +
					mdpd->rho[p_col[j]][TYPE_MAX + d];
			/* printf("%d %d %d %lf %lf\n",j,p_row[j],p_col[j],mdpd->rho[p_row[j]][TYPE_MAX + 0],mdpd->rho[p_col[j]][TYPE_MAX + 0]); */
			dr2[j] = SQR(drdv[j][0]) + SQR(drdv[j][1]) +
							SQR(drdv[j][2]);
			sp[j] = drdv[j][0] * drdv[j][3] +
				drdv[j][1] * drdv[j][4] +
				drdv[j][2] * drdv[j][5];
			idx[n++]=j;
		} 

		rng_uniform_vector(rng, n, theta);
		
		for (k = 0; k < n; k++) {
			j = idx[k];
			double r = sqrt(dr2[j]);
			double force = 0., omega_r;
			
			/* conservative force */
			for (d = 0; d < TYPE_MAX; d++)
				force += mdpd->w[VIR3(d, t_row[j], t_col[j])] *
							rho[TYPE_MAX * j + d];
			force *= DW_FUNC3(r, mdpd->a3) * 2. / 3.;
			force += DW_FUNC2(r, mdpd->a2) * mdpd->v[VIR2(t_row[j], t_col[j])];
			force *= -1. / r;

			/* printf("%d-%d (%lf-%lf) %lf %lf %lf  %lf_%lf\n",p_row[j],p_col[j]-b->nN,force,r,force * drdv[j][0],force * drdv[j][1],force * drdv[j][2],rho[TYPE_MAX * j + 0],rho[TYPE_MAX * j + 1]); */
			
			/* fixme: TensStress() lacks the stochastic forces */
			TensStress(b, force * r, nbl->passport[2 * (i + j) + 0],
     						nbl->passport[2 * (i + j) + 1]);
			/* dissipative and random force */
			omega_r = dpd_wf(r, 1.) / r;
			theta[k] -= .5;
			int tmp = t_row[j] * TYPE_MAX + t_col[j];
			force += (mdpd->zeta[tmp] * theta[k] -
			       	mdpd->gamma[tmp] * sp[j] * omega_r) * omega_r;

			/* 
			 * The random and the dissipative forces have average
			 * values of zero, so that they do not contribute to
			 * the virial pressure. However, they do decorrelate in
			 * time and contribute to the viscosity. Hence, we take
			 * the conservative as well as the dissipative forces
			 * into account.
			 *
			 * This contrasts with former versions of this program.
			 */
			for (d = 0; d < 3; d++)
				b->virial[d] += SQR(drdv[j][d]) * force;
			b->virial[3] += drdv[j][0] * drdv[j][1] * force;
			b->virial[4] += drdv[j][0] * drdv[j][2] * force;
			b->virial[5] += drdv[j][1] * drdv[j][2] * force;

			for (d = 0; d < ARRAY_SIZE(*b->f); d++) {
				b->f[p_row[j]][d] -= force * drdv[j][d];
				b->f[p_col[j]][d] += force * drdv[j][d];
			}
		}
	}
}

#define VIR_TDPD_R(alpha,beta) (omega_r * ((mdpd->zeta[tmp] - mdpd->zeta_p[tmp]) * \
			sp2 * drdv[j][alpha] + mdpd->zeta_p[tmp] * \
			theta[3 * k + alpha] * r) * drdv[j][beta])
#define VIR_TDPD_D(alpha,beta) (-1. * SQR(omega_r) * ((mdpd->gamma[tmp] - \
			mdpd->gamma_p[tmp]) * drdv[j][alpha] * drdv[j][beta] * \
			sp[j] + dr2[j] * mdpd->gamma_p[tmp] * drdv[j][3 + alpha] \
			* drdv[j][beta]))
static void calc_nonbonded_forces_tdpd(struct beads *b, struct mdpd *mdpd)
{
	struct nblist *nbl = &mdpd->nbl;
	int i, d;

	for (i = 0; i < nbl->count; i += STEP) {
		int k, j, n, idx[STEP];
		double dr2[STEP], sp[STEP], theta[3 * STEP];
		int p_row[STEP], p_col[STEP], t_row[STEP], t_col[STEP];
		double rho[TYPE_MAX * STEP];
		VEC2 *drdv = nbl->drdv + i;

		for (j = 0, n = 0; (j < STEP) && (i + j < nbl->count); j++) {

			if (!GET_OK(nbl->passport[2 * (i + j) + 0]))
				continue;
			
			p_row[j] = GET_POS(nbl->passport[2 * (i + j) + 0]);
			t_row[j] = GET_TYPE(nbl->passport[2 * (i + j) + 0]);
			p_col[j] = GET_POS(nbl->passport[2 * (i + j) + 1]);
			t_col[j] = GET_TYPE(nbl->passport[2 * (i + j) + 1]);

			/* only the third order densities are needed here */
			for (d = 0; d < TYPE_MAX; d++)
				rho[TYPE_MAX * j + d] = /* repulsive only */
					mdpd->rho[p_row[j]][TYPE_MAX + d] +
					mdpd->rho[p_col[j]][TYPE_MAX + d];
	
			dr2[j] = SQR(drdv[j][0]) + SQR(drdv[j][1]) +
							SQR(drdv[j][2]);
			sp[j] = drdv[j][0] * drdv[j][3] +
				drdv[j][1] * drdv[j][4] +
				drdv[j][2] * drdv[j][5];
			idx[n++]=j;
		} 

		rng_uniform_vector(rng, 3 * n, theta);

		for (k = 0; k < n; k++) {
			j = idx[k];
			double r = sqrt(dr2[j]);
			double force = 0., omega_r;
			
			/* conservative force */
			for (d = 0; d < TYPE_MAX; d++)
				force += mdpd->w[VIR3(d, t_row[j], t_col[j])] *
							rho[TYPE_MAX * j + d];
			force *= DW_FUNC3(r, mdpd->a3) * 2. / 3.;
			force += DW_FUNC2(r, mdpd->a2) * mdpd->v[VIR2(t_row[j], t_col[j])];
			force *= -1. / r;

			for (d = 0; d < 3; d++)
				b->virial[d] += SQR(drdv[j][d]) * force;
			b->virial[3] += drdv[j][0] * drdv[j][1] * force;
			b->virial[4] += drdv[j][0] * drdv[j][2] * force;
			b->virial[5] += drdv[j][1] * drdv[j][2] * force;

			/* FIXME: TensStress() lacks the stochastic forces */
			TensStress(b, force * r, nbl->passport[2 * (i + j) + 0],
					nbl->passport[2 * (i + j) + 1]);
			omega_r = dpd_wf(r, 1.) / r;
			theta[3 * k + 0] -= .5;
			theta[3 * k + 1] -= .5;
			theta[3 * k + 2] -= .5;
			double sp2 = (drdv[j][0] * theta[3 * k + 0] +
					drdv[j][1] * theta[3 * k + 1] +
					drdv[j][2] * theta[3 * k + 2]) / r;
			int tmp = t_row[j] * TYPE_MAX + t_col[j];
			force += ((mdpd->zeta[tmp] - mdpd->zeta_p[tmp]) * sp2 -
			       	(mdpd->gamma[tmp] - mdpd->gamma_p[tmp]) * sp[j] * omega_r) * omega_r;

			/* The transverse forces would render the stress tensor
			 * non-symmetric. To avoid this unphysical behavior,
			 * we symmetrize the tensor manually. */
			b->virial[0] += VIR_TDPD_R(0, 0) + VIR_TDPD_D(0, 0);
			b->virial[1] += VIR_TDPD_R(1, 1) + VIR_TDPD_D(1, 1);
			b->virial[2] += VIR_TDPD_R(2, 2) + VIR_TDPD_D(2, 2);
			b->virial[3] += .5 * (VIR_TDPD_R(0, 1) + VIR_TDPD_D(0, 1) + VIR_TDPD_R(1, 0) + VIR_TDPD_D(1, 0));
			b->virial[4] += .5 * (VIR_TDPD_R(0, 2) + VIR_TDPD_D(0, 2) + VIR_TDPD_R(2, 0) + VIR_TDPD_D(2, 0));
			b->virial[5] += .5 * (VIR_TDPD_R(1, 2) + VIR_TDPD_D(1, 2) + VIR_TDPD_R(2, 1) + VIR_TDPD_D(2, 1));
			
			double tdiss = SQR(dpd_wf(r, 1.)) * mdpd->gamma_p[tmp];
			double trand = dpd_wf(r, 1.) * mdpd->zeta_p[tmp];
			for (d = 0; d < ARRAY_SIZE(*b->f); d++) {
				/* longitudinal */
				b->f[p_row[j]][d] -= force * drdv[j][d];
				b->f[p_col[j]][d] += force * drdv[j][d];
				/* dissipative, transverse */
				b->f[p_row[j]][d] += tdiss * drdv[j][d + 3];
				b->f[p_col[j]][d] -= tdiss * drdv[j][d + 3];
				/* random, transverse */
				b->f[p_row[j]][d] -= trand * theta[3 * k + d];
				b->f[p_col[j]][d] += trand * theta[3 * k + d];
			}
		}
	}
}

static void langevin_thermostat(struct beads *b, struct mdpd *mdpd)
{
	int i, d, o = b->groupsize * col_index;
	VEC r, fsum = { 0 };

	for (i = 0; i < b->groupsize; i++) {
		int t = GET_TYPE(b->passport[o + i]);

		rng_uniform_vector(rng, 3, r);
		for (d = 0; d < 3; d++) {
			double tmp = mdpd->lzeta[3 * t + d] * (r[d] - .5) - 
				mdpd->lgamma[3 * t + d] * b->xv[o + i][d + 3];
			b->f[o + i][d] += tmp;
			fsum[d] += tmp;
		}
	}

	/*
	 * sum up all thermostat forces and subtract the mean
	 * to avoid a motion of the system's center of mass
	 */
	MPI_Allreduce(MPI_IN_PLACE, fsum, 3, MPI_DOUBLE, MPI_SUM, comm_grid);
	for (d = 0; d < 3; d++) fsum[d] /= b->nN;

	for (i = 0; i < b->groupsize; i++)
		for (d = 0; d < 3; d++)
			b->f[o + i][d] -= fsum[d];
}

static void calc_nonbonded_forces_langevin(struct beads *b, struct mdpd *mdpd)
{
	struct nblist *nbl = &mdpd->nbl;
	int i, j, k, d, n;
	int idx[STEP], p_row[STEP], p_col[STEP], t_row[STEP], t_col[STEP];
	double r, rho[TYPE_MAX * STEP];

	for (i = 0; i < nbl->count; i += STEP) {
		VEC2 *drdv = nbl->drdv + i;

		for (j = 0, n = 0; (j < STEP) && (i + j < nbl->count); j++) {
			if (!GET_OK(nbl->passport[2 * (i + j) + 0])) continue;
			
			p_row[j] = GET_POS(nbl->passport[2 * (i + j) + 0]);
			t_row[j] = GET_TYPE(nbl->passport[2 * (i + j) + 0]);
			p_col[j] = GET_POS(nbl->passport[2 * (i + j) + 1]);
			t_col[j] = GET_TYPE(nbl->passport[2 * (i + j) + 1]);

			/* only the third order densities are needed here */
			for (d = 0; d < TYPE_MAX; d++)
				rho[TYPE_MAX * j + d] = /* repulsive only */
					mdpd->rho[p_row[j]][TYPE_MAX + d] +
					mdpd->rho[p_col[j]][TYPE_MAX + d];
			idx[n++]=j;
		} 

		for (k = 0; k < n; k++) {
			double force = 0.;
			j = idx[k];
			r = sqrt(SQR(drdv[j][0]) + SQR(drdv[j][1]) + SQR(drdv[j][2]));

			/* conservative force */
			for (d = 0; d < TYPE_MAX; d++)
				force += mdpd->w[VIR3(d, t_row[j], t_col[j])] *
							rho[TYPE_MAX * j + d];
			force *= DW_FUNC3(r, mdpd->a3) * 2. / 3.;
			force += DW_FUNC2(r, mdpd->a2) * mdpd->v[VIR2(t_row[j], t_col[j])];
			force *= -1. / r;
			
			/* The stochastic forces do not contribute to the
			 * virial, and we use the conservative part only. */
			b->virial[3] += drdv[j][0] * drdv[j][1] * force;
			b->virial[4] += drdv[j][0] * drdv[j][2] * force;
			b->virial[5] += drdv[j][1] * drdv[j][2] * force;
			TensStress(b, force * r, nbl->passport[2 * (i + j) + 0],
					nbl->passport[2 * (i + j) + 1]);
			
			for (d = 0; d < ARRAY_SIZE(*b->f); d++) {
				b->virial[d] += SQR(drdv[j][d]) * force;
				b->f[p_row[j]][d] -= force * drdv[j][d];
				b->f[p_col[j]][d] += force * drdv[j][d];
			}
		}
	}
	langevin_thermostat(b, mdpd);
}

static void calc_nonbonded_forces_dpdlangevin(struct beads *b, struct mdpd *mdpd)
{
	calc_nonbonded_forces_dpd(b, mdpd);
	langevin_thermostat(b, mdpd);
}

/*
 * the non-bonded densities are summed up and the energy is calculated
 */
static double add_local_densities(struct beads *b, struct mdpd *mdpd)
{
	DENSITY *rho_row = mdpd->rho + b->groupsize * col_index;
	DENSITY *rho_col = mdpd->rho + b->groupsize * (groups_per_row + row_index);
	double energy = 0.;
	int i, j, k;

	memset(mdpd->rho2, 0, ARRAY_SIZE(mdpd->rho2) * sizeof(*mdpd->rho2));
	memset(mdpd->rho3, 0, ARRAY_SIZE(mdpd->rho3) * sizeof(*mdpd->rho3));
	
	for (i = 0; i < b->groupsize; i++) {
		for (j = 0; j < ARRAY_SIZE(*mdpd->rho); j++) {
		        rho_row[i][j] += rho_col[i][j];
			rho_col[i][j] = rho_row[i][j];
		}
		
		int t = GET_TYPE(b->passport[col_index * b->groupsize + i]);
		for (j = 0; j < TYPE_MAX; j++) {
			mdpd->rho2[VIR2(t, j)] += rho_row[i][j];
			for (k = 0; k < TYPE_MAX; k++) {
				mdpd->rho3[VIR3(t, j, k)] +=
						rho_row[i][TYPE_MAX + j] *
						rho_row[i][TYPE_MAX + k];
			}
		}
	}

	for (i = 0; i < SQR(TYPE_MAX); i++)
	  energy += .5 * mdpd->v[i] * mdpd->rho2[i];
	for (i = 0; i < CUBE(TYPE_MAX); i++)
	  energy += 1. / 3. * mdpd->w[i] * mdpd->rho3[i];
	return energy;
}

static double compute_density(struct beads *b, struct mdpd *mdpd)
{
	DENSITY *rho_row = mdpd->rho;
	DENSITY *rho_col = mdpd->rho + groups_per_row * b->groupsize;
	int i;

	memset(mdpd->rho, 0, (groups_per_row + groups_per_col) *
                                        b->groupsize * sizeof(*mdpd->rho));
	weight_density(b, mdpd);

	for (i = 0; i < groups_per_row; i++) {
		if (i == col_index) {
			MPI_Reduce(MPI_IN_PLACE, rho_row + col_index * b->groupsize,
			       	ARRAY_SIZE(*rho_row) * b->groupsize,
				MPI_DOUBLE, MPI_SUM, i, comm_row);
		} else {
			MPI_Reduce(rho_row + i * b->groupsize, NULL,
				ARRAY_SIZE(*rho_row) * b->groupsize,
				MPI_DOUBLE, MPI_SUM, i, comm_row);
		}
	}
	for (i = 0; i < groups_per_col; i++) {
		MPI_Reduce(i == row_index ? MPI_IN_PLACE : rho_col +
			i * b->groupsize, i == row_index ? rho_col +
			row_index * b->groupsize : NULL, ARRAY_SIZE(*rho_col) *
			b->groupsize, MPI_DOUBLE, MPI_SUM, i, comm_col);
	}
	double U = add_local_densities(b, mdpd);

	MPI_Allgather(MPI_IN_PLACE, 0, MPI_DOUBLE, rho_row,
		ARRAY_SIZE(*rho_row) * b->groupsize, MPI_DOUBLE, comm_row);
	MPI_Allgather(MPI_IN_PLACE, 0, MPI_DOUBLE, rho_col,
		ARRAY_SIZE(*rho_col) * b->groupsize, MPI_DOUBLE, comm_col);
	return U;
}

/*
 * calculate all forces and potential energies. This is the point where the
 * interactions come into play. No communication of the forces is done here!
 */
void mdpd_nonbonded_forces(struct beads *b, struct mdpd *mdpd, enum nblist_flags fl)
{
	if (mdpd == NULL) return;
	nblist_update(b, &mdpd->nbl, fl);
	b->e[0] += compute_density(b, mdpd);
	calc_nonbonded_forces(b, mdpd);
}

/* ------------------------------------------------------------------------- */

/*
 * Print accumulated densities that can be used for histogram reweighting.
 */
void mdpd_measure_energy(struct mdpd *mdpd)
{
	int i;

	if (mdpd == NULL) return;

	MPI_Reduce(ismaster ? MPI_IN_PLACE : mdpd->rho2,
			ismaster ? mdpd->rho2 : NULL, ARRAY_SIZE(mdpd->rho2),
			MPI_DOUBLE, MPI_SUM, 0, comm_grid);
	MPI_Reduce(ismaster ? MPI_IN_PLACE : mdpd->rho3,
			ismaster ? mdpd->rho3 : NULL, ARRAY_SIZE(mdpd->rho3),
			MPI_DOUBLE, MPI_SUM, 0, comm_grid);

	if (ismaster) {
		double r2[TYPE_MAX * (TYPE_MAX + 1) / 2];
		double r3[TYPE_MAX * (TYPE_MAX + 1) * (TYPE_MAX + 2) / 6];

		/* print out accumulated densities */
#if (TYPE_MAX == 2)
		r2[0] = mdpd->rho2[VIR2(0, 0)];
		r2[1] = mdpd->rho2[VIR2(0, 1)] + mdpd->rho2[VIR2(1, 0)];
		r2[2] = mdpd->rho2[VIR2(1, 1)];
		r3[0] = mdpd->rho3[VIR3(0, 0, 0)];
		r3[1] = mdpd->rho3[VIR3(0, 0, 1)] +
			mdpd->rho3[VIR3(0, 1, 0)] + mdpd->rho3[VIR3(1, 0, 0)];
		r3[2] = mdpd->rho3[VIR3(1, 1, 0)] +
                        mdpd->rho3[VIR3(1, 0, 1)] + mdpd->rho3[VIR3(0, 1, 1)];
		r3[3] = mdpd->rho3[VIR3(1, 1, 1)];
#elif (TYPE_MAX == 3)
		r2[0] = mdpd->rho2[VIR2(0, 0)];
		r2[1] = mdpd->rho2[VIR2(0, 1)] + mdpd->rho2[VIR2(1, 0)];
		r2[2] = mdpd->rho2[VIR2(0, 2)] + mdpd->rho2[VIR2(2, 0)];
		r2[3] = mdpd->rho2[VIR2(1, 1)];
		r2[4] = mdpd->rho2[VIR2(1, 2)] + mdpd->rho2[VIR2(2, 1)];
		r2[5] = mdpd->rho2[VIR2(2, 2)];

		r3[0] = mdpd->rho3[VIR3(0, 0, 0)];
		r3[1] = mdpd->rho3[VIR3(0, 0, 1)] +
			mdpd->rho3[VIR3(0, 1, 0)] + mdpd->rho3[VIR3(1, 0, 0)];
		r3[2] = mdpd->rho3[VIR3(0, 0, 2)] +
			mdpd->rho3[VIR3(0, 2, 0)] + mdpd->rho3[VIR3(2, 0, 0)];
		r3[3] = mdpd->rho3[VIR3(0, 1, 1)] +
			mdpd->rho3[VIR3(1, 0, 1)] + mdpd->rho3[VIR3(1, 1, 0)];
		r3[4] = mdpd->rho3[VIR3(0, 1, 2)] + mdpd->rho3[VIR3(2, 0, 1)] +
			mdpd->rho3[VIR3(1, 2, 0)] + mdpd->rho3[VIR3(0, 2, 1)] +
			mdpd->rho3[VIR3(2, 1, 0)] + mdpd->rho3[VIR3(1, 0, 2)];
		r3[5] = mdpd->rho3[VIR3(0, 2, 2)] +
			mdpd->rho3[VIR3(2, 0, 2)] + mdpd->rho3[VIR3(2, 2, 0)];
		r3[6] = mdpd->rho3[VIR3(1, 1, 1)];
		r3[7] = mdpd->rho3[VIR3(1, 1, 2)] +
			mdpd->rho3[VIR3(2, 1, 1)] + mdpd->rho3[VIR3(1, 2, 1)];
		r3[8] = mdpd->rho3[VIR3(1, 2, 2)] +
			mdpd->rho3[VIR3(2, 1, 2)] + mdpd->rho3[VIR3(2, 2, 1)];
		r3[9] = mdpd->rho3[VIR3(2, 2, 2)];
#endif
		/* the densities are weighted by R_e^3/N^2 or R_e^6/N^3 */
		for (i = 0; i < ARRAY_SIZE(r2); i++) 
			r2[i] *= CUBE(mdpd->Re) / SQR(mdpd->Re_N);
		for (i = 0; i < ARRAY_SIZE(r3); i++)
			r3[i] *= CUBE(SQR(mdpd->Re) / mdpd->Re_N);
#if (TYPE_MAX == 2)
		printf("sum rho: AA=%lg AB=%lg BB=%lg AAA=%lg AAB=%lg ABB=%lg "
		"BBB=%lg\n", r2[0], r2[1], r2[2], r3[0], r3[1], r3[2], r3[3]);
#elif (TYPE_MAX == 3)
		printf("sum rho: AA=%lg AB=%lg AC=%lg BB=%lg BC=%lg CC=%lg ",
			r2[0], r2[1], r2[2], r2[3], r2[4], r2[5]);
		printf("AAA=%lg AAB=%lg AAC=%lg ABB=%lg ABC=%lg ",
			r3[0], r3[1], r3[2], r3[3], r3[4]);
		printf("ACC=%lg BBB=%lg BBC=%lg BCC=%lg CCC=%lg\n",
			r3[5], r3[6], r3[7], r3[8], r3[9]);	
#endif
	}
}

/*
 * this function reads the two lines in the header of the system file that
 * specify the interaction parameters. This is the first hint, that this
 * force field will be used, so that this function must be called first.
 * The mdpd struct is allocated here.
 */
void mdpd_init(FILE *FH, struct beads *b)
{
	/* double v[TYPE_MAX * (TYPE_MAX + 1) / 2]; */
	/* double w[TYPE_MAX * (TYPE_MAX + 1) * (TYPE_MAX + 2) / 6]; */
	double v[6];
	double w[10];
	struct mdpd *mdpd = malloc(sizeof(*mdpd));
	int i;

	if (mdpd == NULL) novm("mdpd");
	memset(mdpd, 0, sizeof(*mdpd));
	b->mdpd = mdpd;
	int IfTwoType = 1;
	int IfThreeType = 1;
	int NSecond = 6;
	fscanf(FH,"# v=");
	for(int i=0;i<6;i++){
	  if(fscanf(FH,"%lf ",v+i) != 1){
	    //fsetpos(FileToRead,&PosTemp);
	    if(i == 2) IfThreeType = 0;
	    NSecond = 4;
	    break;
	  }
	  if(i >  2) IfTwoType = 0;
	}
	int NThird = 10;
	if(IfTwoType) NThird = 4;
	fscanf(FH,"w=");
	for(int i=0;i<NThird;i++){
	  if(fscanf(FH,"%lf ",w+i) != 1){
	    break;
	  }
	}
	/* fscanf(FH, "# v="); */
	/* for (i = 0; i < ARRAY_SIZE(v); i++) { */
	/* 	if (fscanf(FH, "%lg ", &v[i]) != 1) */
	/* 		fatal(EINVAL, "MDPD virial coefficients (2)"); */
	/* } */
	/* fscanf(FH, "w="); */
	/* for (i = 0; i < ARRAY_SIZE(w); i++) { */
	/* 	if (fscanf(FH, "%lg ", &w[i]) != 1) */
	/* 		fatal(EINVAL, "MDPD virial coefficients (3)"); */
	/* } */
#if (TYPE_MAX == 2)
	debug("MDPD: v_AA=%lg v_AB=%lg v_BB=%lg", v[0], v[1], v[2]);
	debug("MDPD: w_AAA=%lg w_AAB=%lg w_ABB=%lg w_BBB=%lg", w[0], w[1],
								w[2], w[3]);
#elif (TYPE_MAX == 3)
	debug("MDPD: v_AA=%lg v_AB=%lg v_AC=%lg v_BB=%lg v_BC=%lg v_CC=%lg",
					v[0], v[1], v[2], v[3], v[4], v[5]);
	debug("MDPD: w_AAA=%lg w_AAB=%lg w_AAC=%lg w_ABB=%lg w_ABC=%lg",
						w[0], w[1], w[2], w[3], w[4]);
	debug("MDPD: w_ACC=%lg w_BBB=%lg w_BBC=%lg w_BCC=%lg w_CCC=%lg",
						w[5], w[6], w[7], w[8], w[9]);
#endif
	if (ismaster) {
		printf("COEFFICIENTS: v=");
		for (i = 0; i < NSecond; i++)
			printf("%lg ", v[i]);
		printf("w=");
		for (i = 0; i < NThird; i++)
			printf("%lg ", w[i]);
		printf("\n");
	}

	/*
	 * The parameters Re, N, ks, kb, and l0 must be chosen consistently,
	 * so that the rms end-to-end distance is indeed given by Re. However,
	 * there is currently no analytical expression available, so that these
	 * coefficients must be evaluated numerically.
	 *
	 * The program /src/ideal-chain/ic can be used for this purpose.
	 */ 
	if (get_data_line(FH, "# a2=%lg a3=%lg Re=%lg N=%d ks=%lg kb=%lg"
		" l0=%lg", &mdpd->a2, &mdpd->a3, &mdpd->Re, &mdpd->Re_N,
	       	&mdpd->ks, &mdpd->kb, &mdpd->l0) != 7)
		fatal(EINVAL, "MDPD interaction coefficients");
	if (mdpd->ks == 0.) mdpd->ks = 3. * (mdpd->Re_N - 1.) / SQR(mdpd->Re);
	debug("MDPD: a2=%lg a3=%lg Re=%lg N=%d ks=%lg kb=%lg l0=%lg",
				mdpd->a2, mdpd->a3, mdpd->Re, mdpd->Re_N,
			       	mdpd->ks, mdpd->kb, mdpd->l0);

	double sv = CUBE(mdpd->Re) / SQR(mdpd->Re_N);
	double sw = CUBE(SQR(mdpd->Re) / mdpd->Re_N);
	/* TODO: it should be possible to find a general analytical expression
	 * relating the indices of the virial coefficients to array indices. */
	if(IfTwoType){
	  mdpd->v[VIR2(0, 0)] = v[0] * sv;
	  mdpd->v[VIR2(0, 1)] = v[1] * sv;
	  mdpd->v[VIR2(1, 0)] = v[1] * sv;
	  mdpd->v[VIR2(1, 1)] = v[2] * sv;
	  mdpd->w[VIR3(0, 0, 0)] = w[0] * sw;
	  mdpd->w[VIR3(0, 0, 1)] = w[1] * sw;
	  mdpd->w[VIR3(0, 1, 0)] = w[1] * sw;
	  mdpd->w[VIR3(0, 1, 1)] = w[2] * sw;
	  mdpd->w[VIR3(1, 0, 0)] = w[1] * sw;
	  mdpd->w[VIR3(1, 0, 1)] = w[2] * sw;
	  mdpd->w[VIR3(1, 1, 0)] = w[2] * sw;
	  mdpd->w[VIR3(1, 1, 1)] = w[3] * sw;
	}
	else{
	  mdpd->v[VIR2(0, 0)] = v[0] * sv;
	  mdpd->v[VIR2(0, 1)] = v[1] * sv;
	  mdpd->v[VIR2(0, 2)] = v[2] * sv;
	  mdpd->v[VIR2(1, 0)] = v[1] * sv;
	  mdpd->v[VIR2(1, 1)] = v[3] * sv;
	  mdpd->v[VIR2(1, 2)] = v[4] * sv;
	  mdpd->v[VIR2(2, 0)] = v[2] * sv;
	  mdpd->v[VIR2(2, 1)] = v[4] * sv;
	  mdpd->v[VIR2(2, 2)] = v[5] * sv;
	  mdpd->w[VIR3(0, 0, 0)] = w[0] * sw;
	  mdpd->w[VIR3(0, 0, 1)] = w[1] * sw;
	  mdpd->w[VIR3(0, 0, 2)] = w[2] * sw;
	  mdpd->w[VIR3(0, 1, 0)] = w[1] * sw;
	  mdpd->w[VIR3(0, 1, 1)] = w[3] * sw;
	  mdpd->w[VIR3(0, 1, 2)] = w[4] * sw;
	  mdpd->w[VIR3(0, 2, 0)] = w[2] * sw;
	  mdpd->w[VIR3(0, 2, 1)] = w[4] * sw;
	  mdpd->w[VIR3(0, 2, 2)] = w[5] * sw;
	  mdpd->w[VIR3(1, 0, 0)] = w[1] * sw;
	  mdpd->w[VIR3(1, 0, 1)] = w[3] * sw;
	  mdpd->w[VIR3(1, 0, 2)] = w[4] * sw;
	  mdpd->w[VIR3(1, 1, 0)] = w[3] * sw;
	  mdpd->w[VIR3(1, 1, 1)] = w[6] * sw;
	  mdpd->w[VIR3(1, 1, 2)] = w[7] * sw;
	  mdpd->w[VIR3(1, 2, 0)] = w[4] * sw;
	  mdpd->w[VIR3(1, 2, 1)] = w[7] * sw;
	  mdpd->w[VIR3(1, 2, 2)] = w[8] * sw;
	  mdpd->w[VIR3(2, 0, 0)] = w[2] * sw;
	  mdpd->w[VIR3(2, 0, 1)] = w[4] * sw;
	  mdpd->w[VIR3(2, 0, 2)] = w[5] * sw;
	  mdpd->w[VIR3(2, 1, 0)] = w[4] * sw;
	  mdpd->w[VIR3(2, 1, 1)] = w[7] * sw;
	  mdpd->w[VIR3(2, 1, 2)] = w[8] * sw;
	  mdpd->w[VIR3(2, 2, 0)] = w[5] * sw;
	  mdpd->w[VIR3(2, 2, 1)] = w[8] * sw;
	  mdpd->w[VIR3(2, 2, 2)] = w[9] * sw;
	}
}

/*
 * print the two header lines for the system file that specify the interaction
 * parameters. This is somehow the counterpart of mdpd_init().
 */
void mdpd_print_header(FILE *FH, struct mdpd *mdpd)
{
	if (mdpd == NULL) return;
	double v[TYPE_MAX * (TYPE_MAX + 1) / 2];
	double w[TYPE_MAX * (TYPE_MAX + 1) * (TYPE_MAX + 2) / 6];
	int i;

	double sv = SQR(mdpd->Re_N) / CUBE(mdpd->Re);
	double sw = CUBE(mdpd->Re_N) / SQR(CUBE(mdpd->Re));

	/* TODO: it should be possible to find a general analytical expression
	 * relating the indices of the virial coefficients to array indices. */
#if (TYPE_MAX == 2)
	v[0] = mdpd->v[VIR2(0, 0)] * sv; /* v_AA */
	v[1] = mdpd->v[VIR2(0, 1)] * sv; /* v_AB */
	v[2] = mdpd->v[VIR2(1, 1)] * sv; /* v_BB */
	w[0] = mdpd->w[VIR3(0, 0, 0)] * sw; /* w_AAA */
	w[1] = mdpd->w[VIR3(0, 0, 1)] * sw; /* w_AAB */
	w[2] = mdpd->w[VIR3(0, 1, 1)] * sw; /* w_ABB */
	w[3] = mdpd->w[VIR3(1, 1, 1)] * sw; /* w_BBB */
#elif (TYPE_MAX == 3)
	v[0] = mdpd->v[VIR2(0, 0)] * sv; /* v_AA */
	v[1] = mdpd->v[VIR2(0, 1)] * sv; /* v_AB */
	v[2] = mdpd->v[VIR2(0, 2)] * sv; /* v_AC */
	v[3] = mdpd->v[VIR2(1, 1)] * sv; /* v_BB */
	v[4] = mdpd->v[VIR2(1, 2)] * sv; /* v_BC */
	v[5] = mdpd->v[VIR2(2, 2)] * sv; /* v_CC */
	w[0] = mdpd->w[VIR3(0, 0, 0)] * sw; /* w_AAA */
	w[1] = mdpd->w[VIR3(0, 0, 1)] * sw; /* w_AAB */
	w[2] = mdpd->w[VIR3(0, 0, 2)] * sw; /* w_AAC */
	w[3] = mdpd->w[VIR3(0, 1, 1)] * sw; /* w_ABB */
	w[4] = mdpd->w[VIR3(0, 1, 2)] * sw; /* w_ABC */
	w[5] = mdpd->w[VIR3(0, 2, 2)] * sw; /* w_ACC */
	w[6] = mdpd->w[VIR3(1, 1, 1)] * sw; /* w_BBB */
	w[7] = mdpd->w[VIR3(1, 1, 2)] * sw; /* w_BBC */
	w[8] = mdpd->w[VIR3(1, 2, 2)] * sw; /* w_BCC */
	w[9] = mdpd->w[VIR3(2, 2, 2)] * sw; /* w_CCC */
#endif
	fprintf(FH, "# v=");
	for (i = 0; i < ARRAY_SIZE(v); i++) 
		fprintf(FH, "%lg ", v[i]);
	fprintf(FH, "w=");
	for (i = 0; i < ARRAY_SIZE(w); i++) 
		fprintf(FH, "%lg ", w[i]);
	fprintf(FH, "\n# a2=%lg a3=%lg Re=%lg N=%d ks=%lg kb=%lg l0=%lg\n",
				mdpd->a2, mdpd->a3, mdpd->Re, mdpd->Re_N,
				mdpd->ks, mdpd->kb, mdpd->l0);
}

/*
 * frees the memory of an mdpd struct
 */
void mdpd_free(struct mdpd *mdpd)
{
	if (mdpd == NULL) return;

	nblist_free(&mdpd->nbl);
	free(mdpd->rho);
	free(mdpd);
}

/* ------------------------------------------------------------------------- */

static void set_dpd(struct mdpd *mdpd, int type1, int type2, double g, double dt)
{
	debug("setting DPD %d <--> %d to %lg", type1, type2, g);
	mdpd->gamma[type1 * TYPE_MAX + type2] = g;
	mdpd->zeta[type1 * TYPE_MAX + type2] = sqrt(12. * 2. * g / dt);
}

/* Setup standard DPD thermostat */
static void parse_dpd(struct mdpd *mdpd, cfg_t *cfg, double def_gamma, double dt)
{
#if (TYPE_MAX == 2)
	int map[SQR(TYPE_MAX)] = {0, 1, 1, 2};
#else
	int map[SQR(TYPE_MAX)] = {0, 1, 2, 1, 3, 4, 2, 4, 5};
#endif
	int i, j;
	
	if (cfg_size(cfg, "dpd") > 0) {
		cfg_t *sec = cfg_getsec(cfg, "dpd");

		if (cfg_size(sec, "gamma") != TYPE_MAX * (TYPE_MAX + 1) / 2)
			fatal(EINVAL, "Incorrect number of DPD friction"
			" coefficients specified. Is: %u. Should be: %d",
			cfg_size(sec, "gamma"), TYPE_MAX * (TYPE_MAX + 1) / 2);
		for (i = 0; i < TYPE_MAX; i++)
			for (j = 0; j < TYPE_MAX; j++)
				set_dpd(mdpd, i, j, cfg_getnfloat(sec, "gamma",
						map[i * TYPE_MAX + j]), dt);
		if (ismaster) printf("fine-tuned DPD thermostat.\n");
	} else {
		for (i = 0; i < TYPE_MAX; i++) {
			for (j = 0; j < TYPE_MAX; j++) { 
				set_dpd(mdpd, i, j, def_gamma, dt);
			}
		}
		if (ismaster)
			printf("DPD thermostat: gamma=%lg\n", def_gamma);
	}
}

static void set_tdpd(struct mdpd *mdpd, int type1, int type2, double g, double h, double dt)
{
	debug("setting TDPD %d <--> %d to %lg, %lg", type1, type2, g, h);
	mdpd->gamma[type1 * TYPE_MAX + type2] = g;
	mdpd->zeta[type1 * TYPE_MAX + type2] = sqrt(12. * 2. * g / dt);
	mdpd->gamma_p[type1 * TYPE_MAX + type2] = h;
	mdpd->zeta_p[type1 * TYPE_MAX + type2] = sqrt(12. * 2. * h / dt);
}

/* Setup transversal DPD thermostat */
static void parse_tdpd(struct mdpd *mdpd, cfg_t *cfg, double def_gamma,
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
				set_tdpd(mdpd, i, j,
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
				set_tdpd(mdpd, i, j, def_gamma, def_gammap, dt);
			}
		}
		if (ismaster)
			printf("Transverse DPD thermostat: "
			"gamma||=%lg gamma_|_=%lg\n", def_gamma, def_gammap);
	}
}

static void set_langevin(struct mdpd *mdpd, int dir, int type, double g, double dt)
{
	debug("setting Langevin dir %d, type %d to %lg", dir, type, g);
	mdpd->lgamma[type * 3 + dir] = g;
	mdpd->lzeta[type * 3 + dir] = sqrt(12. * 2. * g / dt);
}

/* Setup Langevin thermostat */
static void parse_langevin(struct mdpd *mdpd, cfg_t *cfg, double def_gamma, double dt)
{
	const char *p[3] = {"gamma_x", "gamma_y", "gamma_z"};
	int i, j;
	
	if (cfg_size(cfg, "langevin") > 0) {
		cfg_t *sec = cfg_getsec(cfg, "langevin");

		for (i = 0; i < 3; i++) {
			if (cfg_size(sec, p[i]) != TYPE_MAX &&
						cfg_size(sec, p[i]) != 0)
				fatal(EINVAL, "Incorrect number of Langevin "
				"%s friction coefficients specified. Is "
				"%u. Should be %d or 0", p[i], 
				cfg_size(sec, p[i]), TYPE_MAX);
			for (j = 0; j < TYPE_MAX; j++) {
				set_langevin(mdpd, i, j, cfg_getnfloat(sec, p[i],
						       			j), dt);
			}
		}
		if (ismaster) printf("fine-tuned LANGEVIN thermostat.\n");
	} else {
		for (i = 0; i < 3; i++) {
			for (j = 0; j < TYPE_MAX; j++) {
				set_langevin(mdpd, i, j, def_gamma, dt);
			}
		}
		if (ismaster)
			printf("LANGEVIN thermostat: gamma=%lg\n", def_gamma);
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
void mdpd_parse_thermostat(struct mdpd *mdpd, cfg_t *cfg)
{
	if (mdpd == NULL) return;

	cfg_t *m = cfg_getsec(cfg, "main");
	double gamma = cfg_getfloat(m, "gamma");
	double gamma_p = cfg_getfloat(m, "gamma_p");
	double dt = cfg_getfloat(m, "dt");
	const char *s = cfg_getstr(m, "thermostat");
	
	if (strcasecmp(s, "DPD") == 0) {
		parse_dpd(mdpd, cfg, gamma, dt);
		calc_nonbonded_forces = calc_nonbonded_forces_dpd;

	} else if (strcasecmp(s, "TDPD") == 0) {
		parse_tdpd(mdpd, cfg, gamma, gamma_p, dt);
		calc_nonbonded_forces = calc_nonbonded_forces_tdpd;
		
	} else if (strcasecmp(s, "Langevin") == 0) {
		parse_langevin(mdpd, cfg, gamma, dt);
		calc_nonbonded_forces = calc_nonbonded_forces_langevin;

	} else if (strcasecmp(s, "DPDLangevin") == 0) {
		parse_langevin(mdpd, cfg, gamma, dt);
		parse_dpd(mdpd, cfg, gamma, dt);
		calc_nonbonded_forces = calc_nonbonded_forces_dpdlangevin;

	} else
		fatal(EINVAL, "Unknown thermostat '%s' specified", s);
}

/* ------------------------------------------------------------------------- */

void mdpd_sg_init(struct beads *b, struct mdpd *mdpd, cfg_t *cfg)
{
	if (mdpd == NULL) return;
	int i, j;
		
	if (cfg_size(cfg, "mdpd_sg") > 0) {
		cfg_t *sec = cfg_getsec(cfg, "mdpd_sg");
		double denom;
	
		mdpd->sg_rng = rng_alloc(b->step);
		if (mdpd->sg_rng == NULL) fatal(ENOMEM, "sg_rng");

		if (cfg_size(sec, "str") != cfg_size(sec, "fugacity")) 
			fatal(EINVAL, "Incosistent number of strings (%d) and"
			" fugacities (%d) specified", (int)cfg_size(sec, "str"),
			(int)cfg_size(sec, "fugacity"));
		mdpd->sg_N = cfg_size(sec, "str");
		if (mdpd->sg_N < 2)
			fatal(EINVAL, "You need at least 2 different strings"
			" to run semigrand ensemble (is: %d)", mdpd->sg_N);

		mdpd->sg_str = calloc(mdpd->sg_N, sizeof(*mdpd->sg_str));
		if (mdpd->sg_str == NULL) novm("mdpd->sg_str");
		mdpd->sg_xi = calloc(mdpd->sg_N, sizeof(*mdpd->sg_xi));
		if (mdpd->sg_xi == NULL) novm("mdpd->sg_xi");

		for (i = 0, denom = 0.; i < mdpd->sg_N; i++) {
			mdpd->sg_xi[i] = cfg_getnfloat(sec, "fugacity", i);
			denom += mdpd->sg_xi[i];

			mdpd->sg_str[i] = (i == 0) ?
				calloc(mdpd->sg_N * b->N[0], sizeof(**mdpd->sg_str)) :
				mdpd->sg_str[i - 1] + b->N[0];
			if (mdpd->sg_str[i] == NULL) novm("mdpd->sg_str[0]");

			char *s = cfg_getnstr(sec, "str", i);
			if (strlen(s) != b->N[0])
				fatal(EINVAL, "String %d has invalid length: %d",
					i, (int)strlen(s));
			for (j = 0; j < strlen(s); j++) {
				if (s[j] == '0' || s[j] == 'A')
					mdpd->sg_str[i][j] = 0;
				else if (s[j] == '1' || s[j] == 'B')
					mdpd->sg_str[i][j] = 1;
#if (TYPE_MAX >= 3)
				else if (s[j] == '2' || s[j] == 'C')
					mdpd->sg_str[i][j] = 2;
#endif
				else
					fatal(EINVAL, "Invalid character '%c' in string '%s'",
							s[j], s);
			}
		}
		if (mdpd->sg_xi[0] != 1.) 
			fatal(EINVAL, "1st fugacity must be unity");
		for (i = 0; i < mdpd->sg_N; i++)
			mdpd->sg_xi[i] /= denom;
	} else 
		fatal(EINVAL, "No mdpd_sg section found");
}

void mdpd_sg_free(struct mdpd *mdpd)
{
	if (mdpd == NULL) return;

	if (mdpd->sg_N > 0) {
		free(mdpd->sg_str[0]);
		free(mdpd->sg_str);
		free(mdpd->sg_xi);
		rng_free(mdpd->sg_rng);
	}
}

struct chainprop
{
	PASSPORT *prow;	/* passports in the row for the selected molecule */
	PASSPORT *pcol;	/* passports in the column for the selected molecule */
	int row_index;	/* row index of the process that owns this chain */
	int col_index;	/* column index */
	int offset;	/* offset of the first bead within its group */
};

/* select randomly a molecule and decide about the root process */
static void find_molecule(struct beads *b, struct mdpd *mdpd, struct chainprop *cp)
{
	int i, idx;
	int molid = (int)(rng_uniform(mdpd->sg_rng) * b->n[0]);
	PASSPORT *p = b->passport + b->groupsize * col_index;

	cp->row_index = -1;
	cp->col_index = -1;
	cp->offset = -1;

	for (i = 0, idx = 0; i < b->local_n; i++) {
		int block = b->local_n_b[i];
		int N = b->N[block];

		if (block == 0) {
			if (GET_ID(p[idx]) == molid * b->N[0]) {
				cp->row_index = row_index;
				cp->col_index = col_index;
				cp->offset = idx;
				break;
			}
		}
		idx += N;
	}

	MPI_Allreduce(MPI_IN_PLACE, &cp->row_index, 1, MPI_INT, MPI_MAX, comm_grid);
	MPI_Allreduce(MPI_IN_PLACE, &cp->col_index, 1, MPI_INT, MPI_MAX, comm_grid);
	MPI_Allreduce(MPI_IN_PLACE, &cp->offset, 1, MPI_INT, MPI_MAX, comm_grid);
	assert(cp->row_index >= 0);
	assert(cp->col_index >= 0);
	assert(cp->offset >= 0 && cp->offset < b->groupsize);

	cp->prow = (row_index == cp->row_index) ? b->passport + b->groupsize *
					cp->col_index + cp->offset : NULL;
	cp->pcol = (col_index == cp->col_index) ? b->passport + b->groupsize *
			(groups_per_row + cp->row_index) + cp->offset : NULL;
}

/* send passports to all other participating processes */
static void bcast_passports(struct beads *b, struct chainprop *cp)
{
	if (cp->prow)
		MPI_Bcast(cp->prow, b->N[0], MPI_UNSIGNED_LONG_LONG, cp->col_index, comm_row);
	if (cp->pcol)
		MPI_Bcast(cp->pcol, b->N[0], MPI_UNSIGNED_LONG_LONG, cp->row_index, comm_col);
}

/* switch the bead types to a new, random configuration */
static void switch_passports(struct beads *b, struct mdpd *mdpd, struct chainprop *cp, int *rsn, int *rso)
{
	int i, cnt = 0, so, sn;
	char oldstr[b->N[0]];

	for (i = 0; i < b->N[0]; i++)
		oldstr[i] = GET_TYPE(cp->prow[i]);

	for (i = 0; i < mdpd->sg_N; i++) {
		if (memcmp(oldstr, mdpd->sg_str[i], sizeof(oldstr)) == 0)
			so = i;
	}	

	do {
		sn = (int)(mdpd->sg_N * rng_uniform(rng));
		assert(++cnt < 100);
		assert(sn >= 0 && sn < mdpd->sg_N);
	} while (so == sn);

	for (i = 0; i < b->N[0]; i++) {
		SET_TYPE(cp->prow[i], mdpd->sg_str[sn][i]);
		SET_TYPE(cp->pcol[i], mdpd->sg_str[sn][i]);
	}
	(*rsn) = sn;
	(*rso) = so;
}

/* only the molecules in the **first** (0) block are affected */
void mdpd_sg_move(struct beads *b, struct mdpd *mdpd)
{
	if (mdpd == NULL) return;

	struct chainprop cp;
	PASSPORT prow_backup[b->N[0]], pcol_backup[b->N[0]];
	int acc = 0;

	find_molecule(b, mdpd, &cp);

	if (cp.prow)
		memcpy(prow_backup, cp.prow, b->N[0] * sizeof(*cp.prow));
	if (cp.pcol)
		memcpy(pcol_backup, cp.pcol, b->N[0] * sizeof(*cp.pcol));

	if (cp.prow && cp.pcol) {
		int sn, so;
		switch_passports(b, mdpd, &cp, &sn, &so);
		double r = mdpd->sg_xi[sn] / mdpd->sg_xi[so];
		bcast_passports(b, &cp);
		double dU = b->e[0] - compute_density(b, mdpd);
		
		double z = r * exp(dU);
		if (z >= 1. || z > rng_uniform(rng)) {
			acc = 1;
			printf("SGCMC: move accepted (%d-->%d, U=%lg, r=%d, c=%d)\n",
				so, sn, -1. * dU, cp.row_index, cp.col_index);
		} else {
			acc = 0;
			printf("SGCMC: move rejected (%d-->%d, U=%lg, r=%d, c=%d)\n",
				so, sn, -1. * dU, cp.row_index, cp.col_index);
		}
	} else {
		bcast_passports(b, &cp);
		compute_density(b, mdpd);
	}

	if (cp.prow) {
		MPI_Bcast(&acc, 1, MPI_INT, cp.col_index, comm_row);
		if (!acc) memcpy(cp.prow, prow_backup, b->N[0] * sizeof(*cp.prow));
	}
	if (cp.pcol) {
		MPI_Bcast(&acc, 1, MPI_INT, cp.row_index, comm_col);
		if (!acc) memcpy(cp.pcol, pcol_backup, b->N[0] * sizeof(*cp.pcol));
	}
}

void mdpd_sg_measure(struct beads *b, struct mdpd *mdpd)
{
	if (mdpd == NULL || mdpd->sg_N < 2) return;

	int i, j, idx, cnt[b->N[0]];
	char str[b->N[0]];
	PASSPORT *p = b->passport + col_index * b->groupsize;

	for (j = 0; j < mdpd->sg_N; j++)
		cnt[j] = 0;
	for (i = 0, idx = 0; i < b->local_n; i++) {
		int block = b->local_n_b[i];
		if (block > 0) continue;

		for (j = 0; j < b->N[block]; j++)
			str[j] = GET_TYPE(p[idx + j]);
		for (j = 0; j < mdpd->sg_N; j++) {
			if (memcmp(str, mdpd->sg_str[j], sizeof(str)) == 0)
				cnt[j]++;
		}
		idx += b->N[block];
	}
	MPI_Reduce(ismaster ? MPI_IN_PLACE : cnt, ismaster ? cnt : NULL,
				b->N[0], MPI_INT, MPI_SUM, 0, comm_grid);
	if (ismaster) {
		printf("SGCMC: composition");
		for (i = 0, j = 0; i < mdpd->sg_N; i++) {
			printf(" %d", cnt[i]);
			j += cnt[i];
		}
		printf("\n");
		assert(j == b->n[0]);
	}
}

/* ------------------------------------------------------------------------- */

/* initialize the neighbor list */
void mdpd_init_nblist(struct beads *b, struct mdpd *mdpd)
{
	enum nblist_flags flags = NBL_NEWTON3;
	size_t count = (groups_per_row + groups_per_col) * b->groupsize;
	if (mdpd == NULL) return;
	
	mdpd->rho = calloc(count, sizeof(*mdpd->rho));
	if (mdpd->rho == NULL) novm("mdpd->rho");

#if TYPE_MAX == 3
	/* in case of the phantom solvent, the C-C interactions should
	 * not be included in the neighbor list. */
	if (mdpd->w[VIR3(2, 2, 2)] == 0. && mdpd->v[VIR2(2, 2)] == 0.) {
		flags |= NBL_IGNORE_CC;
		printf("MDPD: Disabling C-C interactions!\n");
	}
#endif
	nblist_alloc(b, &mdpd->nbl, groups_per_row * b->groupsize,
			b->passport, b->xv, groups_per_col * b->groupsize,
			b->passport + groups_per_row * b->groupsize,
			b->xv + groups_per_row * b->groupsize,
			1. + b->rs, flags);
}

/* return the coexistence density */
double mdpd_getrhocoex(struct mdpd *mdpd)
{
	double sv = SQR(mdpd->Re_N) / CUBE(mdpd->Re);
	double sw = CUBE(mdpd->Re_N) / SQR(CUBE(mdpd->Re));
	double v = mdpd->v[VIR2(0, 0)]*sv;
	double w = mdpd->w[VIR3(0, 0, 0)]*sw;
	return 3.*(- .5*v + sqrt(SQR(v)*.25-8.*w/3.))/(4.*w);
}
double mdpd_getRe(struct mdpd *mdpd)
{
	return mdpd->Re;
}
double mdpd_getVar(struct mdpd *mdpd,int Nppc)
{
  return sqrt(SQR(mdpd->Re)/(double)(Nppc-1)/3.)/2.;
}
static int ReallocPairlist(struct nblist *nbl){
  if (unlikely(nbl->count >= nbl->max)) {
    nbl->max += 0x200000;
    nbl->passport = realloc(nbl->passport,
			    2 * nbl->max * sizeof(*nbl->passport));
    if (unlikely(nbl->passport == NULL))
      novm("nbl->passport");
    nbl->drdv = realloc(nbl->drdv, nbl->max * sizeof(*nbl->drdv));
    if (unlikely(nbl->drdv == NULL)) novm("nbl->drdv");
  }
  return 0;
}
//doesn't consider the interaction with the other molecules
//how to touch the pairlist?
double mdpd_EnergySingleCh(struct beads *restrict b,int Ch){
  //int Nppc = b->N[b->local_n_b[0]];
  int Nppc = b->N[0];
  int pFirst = col_index * b->groupsize;
  int pLast = (col_index + 1) * b->groupsize;
  int pFirstCh = pFirst+Ch*Nppc;
  struct mdpd *mdpd = b->mdpd;
  struct nblist *nbl = &mdpd->nbl;
  double dr[3];
  double Nrg = 0.;
  nbl->count = 0;
  // filling the pairlist
  for(int p1=pFirstCh;p1<pFirstCh+Nppc;p1++){
    //printf("%lf %d\n",b->xv[p1][0],Ch);
    for(int p2=pFirst;p2<pLast;p2++){
      //exclude the double interactions
      if(p2 >= pFirst+Ch*Nppc && p2 < pFirst+(Ch+1)*Nppc)
	if(p1 >= p2) 
	  continue;
      for(int d=0;d<3;d++)
	dr[d] = remainder(b->xv[p1][d] - b->xv[p2][d],b->l[d]);
      double Dist = SQR(dr[0])+SQR(dr[1])+SQR(dr[2]);
      if(Dist > SQR(1.)) continue;
      int c = nbl->count++;
      if(ReallocPairlist(nbl)) return 1;
      nbl->passport[2 * c + 0] = b->passport[p1];
      nbl->passport[2 * c + 1] = b->passport[p2];
    }
  }
  //  for(int i=0;i<nbl->count;i++) printf("%d-%d) %d %d %d\n",Ch,i,GET_ID(nbl->passport[i*2+0]),GET_ID(nbl->passport[i*2+1]));
  memset(mdpd->rho, 0, (groups_per_row + groups_per_col) *b->groupsize * sizeof(*mdpd->rho));
  weight_density(b, mdpd);
  DENSITY *rho_row = mdpd->rho + b->groupsize * col_index;
  DENSITY *rho_col = mdpd->rho + b->groupsize * (groups_per_row + row_index);
  double energy = 0.;
  memset(mdpd->rho2,0,ARRAY_SIZE(mdpd->rho2)*sizeof(*mdpd->rho2));
  memset(mdpd->rho3,0,ARRAY_SIZE(mdpd->rho3)*sizeof(*mdpd->rho3));
  for (int i = pFirstCh; i < pFirstCh+Nppc; i++) {
    for (int j = 0; j < ARRAY_SIZE(*mdpd->rho); j++) {
      rho_row[i][j] += rho_col[i][j];
      rho_col[i][j] = rho_row[i][j];
    }
    int t = GET_TYPE(b->passport[col_index * b->groupsize + i]);
    for (int j = 0; j < TYPE_MAX; j++) {
      mdpd->rho2[VIR2(t, j)] += rho_row[i][j];
      for (int k = 0; k < TYPE_MAX; k++) {
  	mdpd->rho3[VIR3(t, j, k)] +=
  	  rho_row[i][TYPE_MAX + j] *
  	  rho_row[i][TYPE_MAX + k];
      }
    }
  }
  for (int i = 0; i < SQR(TYPE_MAX); i++)
    energy += .5 * mdpd->v[i] * mdpd->rho2[i];
  for (int i = 0; i < CUBE(TYPE_MAX); i++)
    energy += 1. / 3. * mdpd->w[i] * mdpd->rho3[i];
  b->e[0] += energy;
  Nrg += energy;
  for (int j = 0; j < Nppc - 1; j++){
    b->e[1] += harmonic_spring(b, Ch*Nppc, j, j + 1, mdpd->ks, mdpd->l0,STRESS_OMIT);
  }
  for (int j = 1; j < Nppc - 1; j++)
    b->e[1] += bond_angle(b, Ch, j - 1, j, j + 1, mdpd->kb);
  return Nrg;
}
double mdpd_EnergySingleMon(struct beads *restrict b,int p1){
  //int Nppc = b->N[b->local_n_b[0]];
  int pFirst = col_index * b->groupsize;
  int pLast = (col_index + 1) * b->groupsize;
  struct mdpd *mdpd = b->mdpd;
  struct nblist *nbl = &mdpd->nbl;
  double dr[3];
  nbl->count = 0;
  int NCutOff = 0;
  for(int p2=pFirst;p2<pLast;p2++){
    //exclude the double interactions
    if(p1 == p2) continue;
    for(int d=0;d<3;d++)
      dr[d] = remainder(b->xv[p1][d] - b->xv[p2][d],b->l[d]);
    double Dist = SQR(dr[0])+SQR(dr[1])+SQR(dr[2]);
    if(Dist > SQR(1.)) continue;
    int c = nbl->count++;
    if(ReallocPairlist(nbl)) return 0.;
    nbl->passport[2 * c + 0] = b->passport[p1];
    nbl->passport[2 * c + 1] = b->passport[p2];
    NCutOff++;
  }
  //for(int i=0;i<nbl->count;i++) printf("%d-%d) %d %d %d\n",Ch,i,GET_ID(nbl->passport[i*2+0]),GET_ID(nbl->passport[i*2+1]));
  memset(mdpd->rho, 0, (groups_per_row + groups_per_col) *
  	 b->groupsize * sizeof(*mdpd->rho));
  weight_density(b, mdpd);
  DENSITY *rho_row = mdpd->rho + b->groupsize * col_index;
  DENSITY *rho_col = mdpd->rho + b->groupsize * (groups_per_row + row_index);
  double energy = 0.;
  memset(mdpd->rho2, 0, ARRAY_SIZE(mdpd->rho2) * sizeof(*mdpd->rho2));
  memset(mdpd->rho3, 0, ARRAY_SIZE(mdpd->rho3) * sizeof(*mdpd->rho3));
  for (int j = 0; j < ARRAY_SIZE(*mdpd->rho); j++) {
    rho_row[p1][j] += rho_col[p1][j];
    rho_col[p1][j] = rho_row[p1][j];
  }
  int t = GET_TYPE(b->passport[col_index * b->groupsize + p1]);
  for (int j = 0; j < TYPE_MAX; j++) {
    mdpd->rho2[VIR2(t, j)] += rho_row[p1][j];
    for (int k = 0; k < TYPE_MAX; k++) {
      mdpd->rho3[VIR3(t, j, k)] +=
	rho_row[p1][TYPE_MAX + j] *
	rho_row[p1][TYPE_MAX + k];
    }
  }
  for (int i = 0; i < SQR(TYPE_MAX); i++)
    energy += .5 * mdpd->v[i] * mdpd->rho2[i];
  for (int i = 0; i < CUBE(TYPE_MAX); i++)
    energy += 1. / 3. * mdpd->w[i] * mdpd->rho3[i];
  //b->e[0] += energy;
  printf("%d) %d %lf\n",p1,NCutOff,energy);
  return energy;
}
