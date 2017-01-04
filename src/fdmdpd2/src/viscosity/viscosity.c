/*
 * viscosity.c - Calculate stress tensor autocorrelation function
 * (C) Copyright 2010 by Martin Hoemberg
 *
 * Syntax: ./viscosity [tw] [ts] [col] {files}
 * avgbins specifies over many bins a running average is performed.
 * The first column has index 1!
 *
 */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <fftw3.h>
#include <string.h>
#include <limits.h>
#include <errno.h>
#include <fenv.h>
#include <assert.h>
#include "common.h"


static int novm(const char *text)
{
	fprintf(stderr, "out of vm: %s\n", text);
	return -ENOMEM;
}

struct data
{
	complex double *fft;
	double *x;
	double *y;
	double *ac;
	double *ac2;
	double *mx;
	double *vx;
	int count;
	int acc;
	int max;
	int n;
	fftw_plan p1;
	fftw_plan p2;
};

static int alloc_fft(struct data *data, int tw)
{
	int n = tw, k = 0;
	while (n > 0) {
		n >>= 1;
		k++;
	}
	n = 2 << k;
	data->n = n;

	data->y = fftw_malloc(n * sizeof(*data->y));
	if (data->y == NULL) return novm("data->y");
	data->fft = fftw_malloc((n / 2 + 1) * sizeof(*data->fft));
	if (data->fft == NULL) return novm("data->fft");

	data->p1 = fftw_plan_dft_r2c_1d(n, data->y, data->fft, FFTW_ESTIMATE);
	data->p2 = fftw_plan_dft_c2r_1d(n, data->fft, data->y, FFTW_ESTIMATE);

	data->ac = calloc(tw, sizeof(*data->ac));
	if (data->ac == NULL) return novm("data->ac");
	data->ac2 = calloc(tw, sizeof(*data->ac2));
	if (data->ac2 == NULL) return novm("data->ac2");
	
	data->mx = calloc(tw, sizeof(*data->mx));
	if (data->mx == NULL) return novm("data->mx");
	data->vx = calloc(tw, sizeof(*data->vx));
	if (data->vx == NULL) return novm("data->xv");

	return 0;
}

static int read_file(struct data *data, const char *fn, int col)
{
	int k, line;
	char buf[LINE_MAX], *s;

	FILE *FH = fopen(fn, "r");
	if (FH == NULL) {
		fprintf(stderr, "No such file '%s'.\n", fn);
		return -ENOENT;
	}

	data->count = 0;
	line = 0;
	while (!feof(FH)) {
		line++;
		if (data->count >= data->max) {
			data->max += 0x10000;
			data->x = realloc(data->x, data->max * sizeof(*data->x));
			if (data->x == NULL) return novm("x");
		}
		if(!fgets(buf, LINE_MAX, FH)) break;
		s = strchr(buf, '#');
		if (s) *s = 0x00;

		s = strtok(buf, " \t");
		k = 0;
		while (s) {
			k++;
			double a = atof(s);
			if (k == col)
				data->x[data->count] = a;
			s = strtok(NULL, " \t");
		}
		if (k >= col) {
			data->count++;
		} else {
			fprintf(stderr, "line %d has only %d columns. ignored.\n", line, k);
		}
	}
	fclose(FH);
	return 0;
}

static int autocorr(struct data *data, int tw, int ts)
{
	int i, j;

	for (i = 0; i + tw <= data->count; i += ts) {
		fprintf(stderr, ".");
		
		memset(data->y, 0, data->n * sizeof(*data->y));
		memset(data->fft, 0, (data->n / 2 + 1) * sizeof(*data->fft));

		for (j = 0; j < tw; j++)
			data->y[j] = data->x[i + j];
		fftw_execute(data->p1);

		for (j = 0; j < data->n / 2 + 1; j++)
			data->fft[j] = data->fft[j] * conj(data->fft[j]);
		fftw_execute(data->p2);
	
		for (j = 0; j < tw; j++) {
			double tmp = data->y[j] / (tw * (double)data->n);
			data->ac[j] += tmp;
			data->ac2[j] += SQR(tmp);
		}
		data->acc++;
	}
	fprintf(stderr, "\n");
	return 0;
}

static void cleanup(struct data *data)
{
	free(data->x);
	fftw_destroy_plan(data->p1);
	fftw_destroy_plan(data->p2);
	fftw_free(data->y);
	fftw_free(data->fft);
	free(data->ac);
	free(data->ac2);
	free(data->mx);
	free(data->vx);
}

int main(int argc, char **argv)
{
	int ret, c = 1, i;
	int col[10], cols = 0, tw, ts;
	struct data data;
	
	feenableexcept(FE_UNDERFLOW | FE_OVERFLOW | FE_INVALID | FE_DIVBYZERO);

	memset(&data, 0, sizeof(data));

	if (argc < 5) {
		fprintf(stderr, "Syntax: %s [tw] [ts] [col1,col2,...] "
							"{files}\n", argv[0]);
		return EXIT_FAILURE;
	}

	tw = atoi(argv[c++]);
	ts = atoi(argv[c++]);
	fprintf(stderr, "tw=%d ts=%d col=", tw, ts);
	char *s = strtok(argv[c++], ",;");
	while (s && atoi(s) > 0) {
		col[cols] = atoi(s);
		fprintf(stderr, "%d ", col[cols++]);
		s = strtok(NULL, ",;");
	}
	fprintf(stderr, "\n");

	if (cols == 0) {
		fprintf(stderr, "No valid column(s) specified!\n");
		return EXIT_FAILURE;
	}
	
	ret = alloc_fft(&data, tw);
	if (ret) return EXIT_FAILURE;

	for (; c < argc; c++) {
		fprintf(stderr, "Reading File %s\n", argv[c]);
		for (i = 0; i < cols; i++) {
			if (read_file(&data, argv[c], col[i]))
				return EXIT_FAILURE;
			if (autocorr(&data, tw, ts))
				return EXIT_FAILURE;
		}
	}
	
	for (c = 0; c < tw; c++) {
		double m = data.ac[c] / data.acc;
		double v = data.ac2[c] / data.acc - m * m;
		v = (v > 0. && data.acc > 1) ? sqrt(v / (data.acc - 1.)) : 0.;
		data.mx[c] = m;
		data.vx[c] = v;
	}
		
	printf("# t C(t) mean-error \\bar C(t) avg.dev.\n# n=%d\n", data.acc);
	for (c = 0; c < tw; c++) {
		int lower = (int)floor(.9 * c);
		int upper = (int)floor(1.1 * c) + 1;
		double m = 0., v = 0.;

		if (upper <= tw) {
			for (i = lower; i < upper; i++)
				m += data.mx[i];
			m /= upper - lower;
			for (i = lower; i < upper; i++)
				v += fabs(data.mx[i] - m);
			v /= upper - lower;
		}

		printf("%d %lg %lg %lg %lg\n", c, data.mx[c], data.vx[c], m, v);
	} 

	cleanup(&data);

	fprintf(stderr, "Done!\n");
	return EXIT_SUCCESS;
}

