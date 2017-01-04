/*
 * corr.c - Calculate Correlation, Autocorrelation and Covariance by FFT
 * (C) Copyright 2010 by Martin Hömberg <mhoembe@gwdg.de>
 *
 * Syntax: ./corr [-cov] [-sm] maxlag re1{,im1} re2{,im2} < in > out
 *
 * cov:    calculate covariance (instead of correlation)
 * maxlag: maximum lag
 * col1:   index of first column (counting starts by one)
 * col2:   index of second column
 * sm:     subtract mean from input quantities
 *
 */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <fftw3.h>
#include <math.h>
#include <errno.h>
#include <limits.h>
#include <string.h>
#include <argp.h>
#include <stdarg.h>
#include <assert.h>

/*
 * print debugging message to stderr
 */
static void debug(const char *fmt, ...)
{
#ifdef DEBUG
	va_list args;
	va_start(args, fmt);
	vfprintf(stderr, fmt, args);
	fprintf(stderr, "\n");
	va_end(args);
#endif
}

/*
 * print error message to stderr and return error code
 */
static int error(const int errnum, const char *fmt, ...)
{
	va_list args;
	va_start(args, fmt);
	vfprintf(stderr, fmt, args);
	fprintf(stderr, ": %s\n", strerror(errnum));
	va_end(args);

	return -errnum;
}

/*
 * print out of memory message and return error code
 */ 
static int novm(const char *text)
{
	return error(ENOMEM, text);
}

static struct argp_option argp_option[] = {
	{"cov", 'C', NULL, 0, "Calculate covariance instead of correlation"},
	{"mean", 'M', NULL, 0, "Subtract mean from input quantities"},
	{"lag", 'L', "time", 0, "Maximum lag"},
	{NULL, 'A', "col1{,col2}", 0, "Column for first quantity \"A\""},
	{NULL, 'B', "col1{,col2}", 0, "Column for second quantity \"B\""},
	{0}
};

/* static variables for command line options */
static int do_cov = 0;			/* calculate covariance? */
static int do_subtract_mean = 0;	/* subtract mean from input? */
static int maxlag = 10000;		/* max. lag */
static int cre1 = 0;
static int cim1 = -1;
static int cre2 = 0;
static int cim2 = -1;


static error_t parse_opt(int key, char *arg, struct argp_state *state)
{
	char *str;
	
	switch (key) {
	case 'C':
		do_cov = 1;
		break;
	case 'M':
		do_subtract_mean = 1;
		break;
	case 'L':
		maxlag = atoi(arg);
		debug("Using maxlag=%d", maxlag);
		break;
	case 'A':
		cre1 = atoi(arg) - 1;
		if (str = strpbrk(arg, ";:,|")) {
			cim1 = atoi(str + 1) - 1;
			debug("A is complex. Real/imag parts from columns "
						"%d/%d", cre1 + 1, cim1 + 1);
		} else {
			debug("A is real. Data from column %d", cre1 + 1);
		}
		break;
	case 'B':
		cre2 = atoi(arg) - 1;
		if (str = strpbrk(arg, ";:,|")) {
			cim2 = atoi(str + 1) - 1;
			debug("B is complex. Real/imag parts from columns "
						"%d/%d", cre2 + 1, cim2 + 1);
		} else {
			debug("B is real. Data from column %d", cre2 + 1);
		}
		break;
	default:
		return ARGP_ERR_UNKNOWN;
	}
	return 0;
}

const char *argp_program_version = "corr 1.0";
const char *argp_program_bug_address = "<mhoembe@gwdg.de>";
static const char doc[] = "Calculate (Auto-)Correlation and Covariance";
static struct argp argp = {argp_option, parse_opt, 0, doc};


static complex double *bufA = NULL;	/* buffer for first quantity */
static complex double *bufB = NULL;	/* buffer for second quantity */
static int maxcount = 0;		/* max. size for buffer */
static int count = 0;			/* current number of entries */
static double cov = 0.;			/* covariance */

static int parse_input(FILE *FH)
{
	int omc = 0;

	while (!feof(FH)) {
		double re1 = 0., im1 = 0., re2 = 0., im2 = 0.;
		char buf[LINE_MAX], *str;
		int n = 0;

		str = fgets(buf, sizeof(buf), FH);
		if (str == NULL) continue;

		if (strchr(buf, '#') == buf) {
			omc++;
			continue;
		}
		
		str = strtok(buf, " \n\t");
		if (str == NULL)
			return error(EINVAL, "first strtok() returned NULL");
		do {
			if (n == cre1) sscanf(str, "%lf", &re1);
			if (n == cim1) sscanf(str, "%lf", &im1);
			if (n == cre2) sscanf(str, "%lf", &re2);
			if (n == cim2) sscanf(str, "%lf", &im2);

			n++;
			str = strtok(NULL, " \n\t");
		} while (str != NULL);

		if (count >= maxcount) {
			maxcount += 65536;

			bufA = realloc(bufA, maxcount * sizeof(*bufA));
			if (bufA == NULL) return novm("bufA");
			bufB = realloc(bufB, maxcount * sizeof(*bufB));
			if (bufB == NULL) return novm("bufB");
		}

		bufA[count] = re1 + I * im1;
		bufB[count] = re2 + I * im2;
		count++;
	}
	debug("Found %d data lines; %d lines omitted.", count, omc);

	if (count < maxlag)
		return error(EINVAL, "not enough datapoints");
	return 0;
}

static int subtract_mean(void)
{
	int i;
	complex double m1 = 0., m2 = 0.;

	for (i = 0; i < count; i++) {
		m1 += bufA[i];
		m2 += bufB[i];
	}
	m1 /= count;
	m2 /= count;
	for (i = 0; i < count; i++) {
		bufA[i] -= m1;
		bufB[i] -= m2;
	}
	debug("Subtracted means. m1=%lg+I*%lg, m2=%lg+I*%lg.",
				creal(m1), cimag(m1), creal(m2), cimag(m2));
	return 0;
}

static int covariance(void)
{
	int i;
	for (i = 0; i < count; i++)
		cov += creal(bufA[i] * conj(bufB[i]));
	cov /= count;
	return 0;
}

static int correlate(void)
{
	assert(maxlag <= count);
	int fftsize = count + maxlag;
	int i;
	
	complex double *in = fftw_malloc(fftsize * sizeof(*in));
	complex double *out = fftw_malloc(fftsize * sizeof(*out));
	complex double *tmp = fftw_malloc(fftsize * sizeof(*out));
	if (in == NULL) return novm("in");
	if (out == NULL) return novm("out");
	if (tmp == NULL) return novm("tmp");
	
	fftw_plan p1, p2;
       
	p1 = fftw_plan_dft_1d(fftsize, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	p2 = fftw_plan_dft_1d(fftsize, out, in, FFTW_BACKWARD, FFTW_ESTIMATE);

	/* FFT first quantity */
	for (i = 0; i < count; i++)
		in[i] = bufA[i];
	for (; i < fftsize; i++)
		in[i] = 0.;
	fftw_execute(p1);
	memcpy(tmp, out, fftsize * sizeof(*out));

	/* FFT second quantity */
	for (i = 0; i < count; i++)
		in[i] = bufB[i];
	for (; i < fftsize; i++)
		in[i] = 0.;
	fftw_execute(p1);

	/* multiply both */
	for (i = 0; i < fftsize; i++) 
		out[i] = conj(out[i]) * tmp[i];

	/* backwards transform */
	fftw_execute(p2);

	/* copy back results */
	for (i = 0; i < maxlag; i++) {
		bufA[i] = in[i] / fftsize / count;
		bufB[i] = 0.;
	}
	for (; i < count; i++) {
		bufA[i] = 0.;
		bufB[i] = 0.;
	}

	fftw_destroy_plan(p1);
	fftw_destroy_plan(p2);
	fftw_free(in);
	fftw_free(out);
	fftw_free(tmp);

	return 0;
}

static int print_corr(void)
{
	int i;

	for (i = 0; i < maxlag; i++) {
		double re = creal(bufA[i]);
		double im = cimag(bufA[i]);
		if (!do_cov) {
			re /= cov;
			im /= cov;
		}
		if (cim1 < 0 && cim2 < 0) {
			printf("%d %lg\n", i, re);
		} else {
			printf("%d %lg %lg\n", i, re, im);
		}
	}
	return 0;
}


int main(int argc, char **argv)
{
	int ret;
	
	argp_parse(&argp, argc, argv, 0, 0, 0);
	
	ret = parse_input(stdin);
	if (ret) goto error;

	if (do_subtract_mean)
		subtract_mean();
	covariance();

	ret = correlate();
	if (ret) goto error;

	print_corr();


	return EXIT_SUCCESS;
error:
	return EXIT_FAILURE;
}

