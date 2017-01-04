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

/* $Id: error.c 296 2011-06-15 12:17:05Z hoemberg $ */

#include "fdmdpd2.h"

/*
 * Terminate the MPI program with a SIGABRT of the calling process. If 
 * 'ulimit -c' is set appropriately, one gets a single core dump of the calling
 * process; all others processes receive a SIGTERM and later a SIGKILL.
 */
static void __attribute__((noreturn)) crash(void)
{
	fprintf(stderr, "\nAborting.\n\n");
	abort(); /* in contrast to MPI_Abort() this does a core dump */
}

/*
 * Print nicely formatted message to STDERR
 */
static void printmsg(int errnum, const char *premsg, const char *fmt, va_list args)
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	fprintf(stderr, "%05d: ", rank);
	if (premsg != NULL)
		fprintf(stderr, "%s ", premsg);
	vfprintf(stderr, fmt, args);
	fprintf(stderr, "\n");
}

/*
 * print debugging message to stderr
 */
#ifdef DEBUG
void debug(const char *fmt, ...)
{
	va_list args;
	va_start(args, fmt);
	printmsg(0, "debug", fmt, args);
	va_end(args);
}
#endif

/*
 * print out of memory message and return error code
 */ 
void novm(const char *fmt, ...)
{
	va_list args;
	va_start(args, fmt);
	printmsg(ENOMEM, "out of vm at ", fmt, args);
	va_end(args);
	crash();
}

/*
 * Report a critical error and die horribly.
 */
void fatal(const int errnum, const char *fmt, ...)
{
	va_list args;
	va_start(args, fmt);
	printmsg(errnum, NULL, fmt, args);
	va_end(args);
	crash();
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

