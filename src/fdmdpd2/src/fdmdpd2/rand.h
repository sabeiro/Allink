/*
 * rand.h - Random number generator include file
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

#ifndef __RAND_H__
#define __RAND_H__

struct rng;

extern struct rng *rng_alloc(unsigned long seed);
extern void rng_free(struct rng *restrict rng);
extern double rng_uniform(struct rng *restrict rng);
extern int rng_uniform_vector(struct rng *restrict rng, int n,
							double *restrict r);
extern double rng_gaussian(struct rng *r, double m, double s);

#endif

