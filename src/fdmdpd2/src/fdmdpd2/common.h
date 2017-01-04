/*
 * common.h - Common include file for FDMDPD2
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

#ifndef __COMMON_H__
#define __COMMON_H__

#define _GNU_SOURCE

//#ifndef DEBUG
//#define NDEBUG
//#endif

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <complex.h>
#include <assert.h>
#include <math.h>
#include <errno.h>
#include <limits.h>

/* some useful macros */
#define SQR(X)		((X) * (X))
#define CUBE(X)		((X) * (X) * (X))
#define MAX(X,Y)	(((X) > (Y)) ? (X) : (Y))
#define MIN(X,Y)	(((X) < (Y)) ? (X) : (Y))
#define ARRAY_SIZE(arr) (sizeof(arr) / sizeof((arr)[0]))

/* branch prediction macros, uses a GCC extension. */
#define	likely(exp)	__builtin_expect(!!(exp), 1)
#define	unlikely(exp)	__builtin_expect(!!(exp), 0)


#endif

