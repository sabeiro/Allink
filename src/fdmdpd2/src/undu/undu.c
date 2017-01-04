/*
 * undu.c - calculate height and thickness fluctuation spectrum
 * (C) Copyright 2009,2011 Martin Hoemberg <mhoembe@gwdg.de>
 *
 * Syntax: ./undu [ny] [nz] [output.dat] > spectrum.dat
 */

#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <complex.h>
#include <math.h>
#include <limits.h>
#include <fenv.h>
#include <errno.h>
#include <fftw3.h>
#include <assert.h>

#define SQR(x)		((x)*(x))

#define NORMAL		0
#define TANG1		1
#define TANG2		2
#define BLOCKS_MAX	2

struct data
{
	complex double *out;
	double *hm;
	double *h;
	double *cnt;
	double *in;
	double *tq;
	double *hq;
	double *norm;		/* normalization factor (sinc(x)) */
	FILE *fhtq;
	FILE *fhhq;
	double A;		/* mean area */
	double A2;		/* second moment of area */
	double L[3];
	double time;
	double Lavg[3];
	double b1L;
	double b2L;
	int filecount;		/* number of files */
	int n1, n2;		/* maximum quantum number */
	int rsize, qsize;	/* realspace / fourier space size */
	fftw_plan plan;
};

/* ------------------------------------------------------------------------- */

/*
 * print error message to stderr and return error code
 */
static void __attribute__((format(printf, 2, 3),noreturn)) error(const int errnum, const char *fmt, ...)
{
	va_list args;

	va_start(args, fmt);
	vfprintf(stderr, fmt, args);
	fprintf(stderr, ": %s\n", strerror(errnum));
	va_end(args);

	exit(EXIT_FAILURE);
}

/* ------------------------------------------------------------------------- */

static void allocate(int n1, int n2, struct data *data)
{
	data->n1 = n1;
	data->n2 = n2;
	data->rsize = n1 * n2;
	data->qsize = n1 * (n2 / 2 + 1);
	
	data->in = fftw_malloc(n1 * n2 * sizeof(*data->in));
	if (data->in == NULL) error(ENOMEM, "data->in");
	data->out = fftw_malloc(n1 * (n2 / 2 + 1) * sizeof(*data->out));
	if (data->out == NULL) error(ENOMEM, "data->out");

	data->plan = fftw_plan_dft_r2c_2d(n1, n2, data->in, data->out, FFTW_MEASURE);

	data->cnt = calloc(n1 * n2, sizeof(*data->cnt));
	if (data->cnt == NULL) error(ENOMEM, "data->cnt");
	data->hm = calloc(n1 * n2, sizeof(*data->hm));
	if (data->hm == NULL) error(ENOMEM, "data->hm");
	data->h = calloc(BLOCKS_MAX * n1 * n2, sizeof(*data->h));
	if (data->h == NULL) error(ENOMEM, "data->h");

	data->tq = calloc(data->qsize, sizeof(*data->tq));
	if (data->tq == NULL) error(ENOMEM, "data->tq");
	data->hq = calloc(data->qsize, sizeof(*data->hq));
	if (data->hq == NULL) error(ENOMEM, "data->hq");

	data->fhtq = fopen("dyndump_tq.dat", "a");
	if (data->fhtq == NULL) error(EIO, "dyndump_tq.dat");
	data->fhhq = fopen("dyndump_hq.dat", "a");
	if (data->fhhq == NULL) error(EIO, "dyndump_hq.dat");

	data->norm = calloc(data->qsize, sizeof(*data->norm));
	if (data->norm == NULL) error(ENOMEM, "data->norm");
}


static void cleanup(struct data *data)
{
	fftw_destroy_plan(data->plan);
	fftw_free(data->in);
	fftw_free(data->out);
	fftw_cleanup();
	free(data->cnt);
	free(data->hm);
	free(data->h);
	free(data->tq);
	free(data->hq);
	free(data->norm);
	fclose(data->fhtq);
	fclose(data->fhhq);
}

static double sinc(const double x) {
	return (x == 0.) ? 1. : sin(M_PI * x) / (M_PI * x);
}

static void calculate_normalization(struct data *data)
{
	int i, n1, n2;

	for (i = 0; i < data->n1; i++) {
		n1 = (i < data->n1 / 2) ? i : (i - data->n1);
		for (n2 = 0; n2 < data->n2 / 2 + 1; n2++)
			data->norm[i * (data->n2 / 2 + 1) + n2] = 1. /
					(data->n1 * sinc(n1 / data->n1) *
					data->n2 * sinc(n2 / data->n2));
	}
}

/* ------------------------------------------------------------------------- */

static void print_results(struct data *data)
{
	int i, n2;
	double q1 = 2. * M_PI / (data->Lavg[TANG1] / data->filecount);
	double q2 = 2. * M_PI / (data->Lavg[TANG2] / data->filecount);
	data->A /= data->filecount;
	data->A2 /= data->filecount;
	printf("# q1 q2 |q| hq tq\n# <A>=%lg <A2>=%lg nx=%d ny=%d\n",
				data->A, data->A2, data->n1, data->n2);
	for (i = 0; i < data->n1; i++) {
		int n1 = (i < data->n1 / 2) ? i : (i - data->n1);
		for (n2 = 0; n2 < data->n2 / 2 + 1; n2++) {
			double hq = data->hq[i * (data->n2 / 2 + 1) + n2] / data->filecount;
			double tq = data->tq[i * (data->n2 / 2 + 1) + n2] / data->filecount;
			printf("%lg %lg %lg %lg %lg\n", q1 * n1, q2 * n2, hypot(q1 * n1, q2 * n2), hq, tq);
		}
	}
}

/* ------------------------------------------------------------------------- */

static inline int getbin(double x, double L, double bL)
{
	while (x >= L) x -= L;
	while (x < 0.) x += L;
	return (int)(bL * x);
}

static void parsecoord(const char *line, struct data *data)
{
	double x[6];
	int t;
	
	if (sscanf(line, "%lg %lg %lg %lg %lg %lg %d\n",
		&x[0], &x[1], &x[2], &x[3], &x[4], &x[5], &t) != 7)
		error(EINVAL, "Invalid data line '%s'", line);
	if (t == 1) return; /* only hydrophobic beads */

	int bin1 = getbin(x[TANG1], data->L[TANG1], data->b1L);
	int bin2 = getbin(x[TANG2], data->L[TANG2], data->b2L);
	assert(bin1 >= 0 && bin1 < data->n1);
	assert(bin2 >= 0 && bin2 < data->n2);

	data->hm[bin1 * data->n2 + bin2] += x[NORMAL];
	data->cnt[bin1 * data->n2 + bin2]++;
}

static int isnewsnapshot(const char *line, struct data *data, int *blocks)
{
	if (strstr(line, "# L=") != line) return 0;
	if (sscanf(line, "# L=%lg %lg %lg t=%lg blocks=%d", &data->L[0],
		&data->L[1], &data->L[2], &data->time, blocks) == 5) {
		assert(*blocks <= BLOCKS_MAX);
		data->b1L = data->n1 / data->L[TANG1];
		data->b2L = data->n2 / data->L[TANG2];
		fprintf(stderr, "%s", line);
		return 1;
	}
	error(EINVAL, "Line(1) '%s' not understood", line);
}

static int isnewblock(const char *line, struct data *data, int *n, int *N, char *name)
{
	if (strstr(line, "# n=") != line) return 0;
	if (sscanf(line, "# n=%d N=%d name=%s\n", n, N, name) == 3) {
		memset(data->hm, 0, data->rsize * sizeof(*data->hm));
		memset(data->cnt, 0, data->rsize * sizeof(*data->cnt));
		return 1;
	}
	error(EINVAL, "Line(2) '%s' not understood", line);
}

static void finishblock(int block, char *s, struct data *data)
{
	int i;
	for (i = 0; i < data->rsize; i++) {
		double r = (data->cnt[i] > 0.) ? (data->hm[i] / data->cnt[i]) : 0.;
		data->h[block * data->rsize + i] = r;
	}
}

static void dyndump(FILE *FH, struct data *data)
{
	int i;

	fprintf(FH, "# time=%lg\n", data->time);
	for (i = 0; i < data->qsize; i++) {
		complex double z = data->out[i];
		fprintf(FH, "%lg %lg\n", creal(z), cimag(z));
	}
}

static void finishsnapshot(struct data *data)
{
	int i;

	/* thickness spectrum */
	for (i = 0; i < data->rsize; i++)
		data->in[i] = .5 * (data->h[i] - data->h[data->rsize + i]);
	fftw_execute(data->plan);
	for (i = 0; i < data->qsize; i++) {
		data->out[i] *= data->norm[i];
		data->tq[i] += creal(data->out[i] * conj(data->out[i]));
	}
	dyndump(data->fhtq, data);

	/* height spectrum */
	for (i = 0; i < data->rsize; i++)
		data->in[i] = .5 * (data->h[i] + data->h[data->rsize + i]);
	fftw_execute(data->plan);
	for (i = 0; i < data->qsize; i++) {
		data->out[i] *= data->norm[i];
		data->hq[i] += creal(data->out[i] * conj(data->out[i]));
	}
	dyndump(data->fhhq, data);

	data->A += data->L[TANG1] * data->L[TANG2];
	data->A2 += SQR(data->L[TANG1] * data->L[TANG2]);

	for (i = 0; i < 3; i++)
		data->Lavg[i] += data->L[i];

	data->filecount++;
}

#define crBegin static int state=0; switch(state) { case 0:
#define crReturn(x) do { state=__LINE__; return x; \
	                         case __LINE__:; } while (0)
#define crFinish } do { state=0; } while(0)

static void parseall(const char *line, struct data *data)
{
	static int block, blocks, i, n, N;
	static char name[40];

	crBegin;

	if (isnewsnapshot(line, data, &blocks)) {
		for (block = 0; block < blocks; block++) {
			crReturn();
			if (isnewblock(line, data, &n, &N, name)) {
				for (i = 0; i < n * N; i++) {
					crReturn();
					parsecoord(line, data);
				}
				finishblock(block, name, data);
			} else
				return; /* ignore other header lines */
		}
		finishsnapshot(data);
	}
	crFinish;
}

static void parse(const char *fn, struct data *data)
{
	char line[LINE_MAX];
	FILE *FH = (strcmp(fn, "-") == 0) ? stdin : fopen(fn, "r");
	if (FH == NULL)
		error(ENOENT, "Couldn't open %s for reading", fn);
			
	while (fgets(line, LINE_MAX, FH) != NULL)
		parseall(line, data);
	parseall("Premature end", data);
	if (FH != stdin) fclose(FH);
}

/* ------------------------------------------------------------------------- */

int main(int argc, char **argv)
{
#if defined(DEBUG) && defined(FE_NOMASK_ENV)
	feenableexcept(FE_UNDERFLOW | FE_OVERFLOW | FE_INVALID | FE_DIVBYZERO);
#endif
	if (argc == 4) {
		struct data data = { 0 };

		allocate(atoi(argv[1]), atoi(argv[2]), &data);
		calculate_normalization(&data);
		parse(argv[3], &data);
		print_results(&data);
		cleanup(&data);
	} else
		fprintf(stderr, "Syntax: %s [n1] [n2] [output.dat|-]\n\n", argv[0]);
	return EXIT_SUCCESS;
}

