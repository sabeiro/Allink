/*
 * error.c - Error handling
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

/* $Id: error.c 126 2010-05-05 07:40:12Z hoemberg $ */

#include "common.h"

/*
 * print debugging message to stderr
 */
#ifdef DEBUG
void debug(const char *fmt, ...)
{
	va_list args;
	va_start(args, fmt);
#ifdef MPI_INCLUDED
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	fprintf(stderr, "%05d: debug: ", rank);
#endif
	vfprintf(stderr, fmt, args);
	fprintf(stderr, "\n");
	va_end(args);
}
#endif

/*
 * print error message to stderr and return error code
 */
int error(const int errnum, const char *fmt, ...)
{
	va_list args;
	va_start(args, fmt);
#ifdef MPI_INCLUDED
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	fprintf(stderr, "%05d: ",rank);
#endif
	vfprintf(stderr, fmt, args);
	fprintf(stderr, ": %s\n", strerror(errnum));
	va_end(args);

	return -errnum;
}

/*
 * print out of memory message and return error code
 */ 
int novm(const char *text)
{
	return error(ENOMEM, text);
}

/*
 * This function reads in a single line from the given file handle and
 * parses the information.
 */
int get_data_line(FILE *FH, const char *restrict fmt, ...)
{
	va_list ap;
	int ret;
	char *s, line[LINE_MAX];

	va_start(ap, fmt);
	s = fgets(line, LINE_MAX, FH);
	if (s) {
		ret = vsscanf(line, fmt, ap);
	} else {
		ret = 0;
	}
	va_end(ap);
	return ret;
}

