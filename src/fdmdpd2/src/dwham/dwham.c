/*
 * dwham.c - Distributed Weighted Histogram Analysis Method
 * (C) Copyright 2009 Martin Hoemberg <mhoembe@gwdg.de>
 *
 * Syntax: wham [config.dat] [list]
 *
 */

#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include <errno.h>
#include <fenv.h>
#include <mpfr.h>
#include <mpi.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_histogram2d.h>

static mp_rnd_t rnd = GMP_RNDN;

#define CUBE(x)		((x)*(x)*(x))
#define SQR(x)		((x)*(x))

#define COEFF_MAX	20

struct expansion
{
	double aa;
	double ab;
	double bb;
	double aaa;
	double aab;
	double abb;
	double bbb;
};

struct dp
{
	struct expansion r;		/* integrated densities */
	double lambda;			/* order parameter */
	double U;			/* total energy */
	double A;			/* surface area */
	int run;			/* index of the run */
};

struct run
{
	struct expansion vir;		/* virial coefficients */
	double a[COEFF_MAX];		/* US potential coefficients */
	double l0;			/* US potential center */
	double f;			/* free energy */
	double u;			/* average total energy */
	int n;				/* number of data points in this run */
};

struct data
{
	struct dp *dp;			/* process local data points */
	struct run *run;		/* pointer to struct with runs */
	mpfr_t *G;			/* G buffer */
	mpfr_t *F;			/* F buffer */
	struct expansion vir0;		/* reference virial coefficients */
	double lambda_min;		/* lambda histogram starts here */
	double lambda_max;		/* lambda histogram ends here */
	double lambda_bw;		/* lambda histogram bin width */
	double rho_min;			/* rho_coex histogram starts here */
	double rho_max;			/* rho_coex histogram ends here */
	double rho_bw;			/* rho_coex histogram bin width */
	double U_min;			/* internal energy histogram */
	double U_max;			/* start, end, ... */
	double U_bw;			/* ... and binwidth */
	double A_min;			/* start of area histogram ... */
	double A_max;			/* ... end ... */
	double A_bw;			/* ... and binwidth */
	double rho_1;			/* reweighting target rho_coex */
	double kN_1;			/* reweighting target kN */
	double chiN_1;			/* reweighting target chiN */
	int nN;				/* extensive number of bonds */
	int lambda_bins;		/* number of bins in lambda histo */
	int rho_bins;			/* number of bins in rho_coex histo */
	int U_bins;			/* number of bins in U */
	int A_bins;			/* number of bins in A */
	int runs_n;			/* number of runs */
	int runs_max;			/* max number of runs */
	int dpc;			/* total data point counter */
	int dp_count;			/* process local data point count */
	int dp_max;			/* max number of data points */
};

/* ------------------------------------------------------------------------- */

static int error(int errnum, const char *fmt, ...)
					__attribute__((format(printf, 2, 3)));
static int get_data_line(FILE *restrict FH, const char *restrict fmt,
		       const int n, ...) __attribute__((format(scanf, 2, 4)));

/*
 * print debugging message to stderr
 */
#ifdef DEBUG
static void debug(const char *fmt, ...) __attribute__((format(printf, 1, 2)));
static void debug(const char *fmt, ...)
{
	va_list args;
	int rank;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	va_start(args, fmt);
	fprintf(stderr, "%05d: debug: ", rank);
	vfprintf(stderr, fmt, args);
	fprintf(stderr, "\n");
	va_end(args);
}
#else
#define debug(...) {}
#endif

/*
 * print error message to stderr and return error code
 */
static int error(const int errnum, const char *fmt, ...)
{
	va_list args;
	int rank;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	va_start(args, fmt);
	fprintf(stderr, "%05d: ",rank);
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

/*
 * This function reads in a single line from the given file handle and
 * parses the information.
 */
static int get_data_line(FILE *FH, const char *restrict fmt, const int n, ...)
{
	va_list ap;
	int ret;
	char line[LINE_MAX];

	va_start(ap, n);
	fgets(line, LINE_MAX, FH);
	ret = vsscanf(line, fmt, ap);
	va_end(ap);

	if (ret == EOF) {
		debug("EOF encountered.");
	} else if (ret != n) {
		return error(EINVAL, "Error in '%s'. Expected %d; Found %d",
								line, n, ret);
	}
	return 0;
}

/* ------------------------------------------------------------------------- */

static int mpfr_calloc_array(size_t n, mpfr_t **restrict rc)
{
	mpfr_t *c = calloc(n, sizeof(*c));
	size_t i;

	if (c == NULL) return novm("mpfr_calloc_array");

	for (i = 0; i < n; i++) {
		mpfr_init(c[i]);
		mpfr_set_d(c[i], 0., rnd);
	}

	(*rc) = c;

	return 0; 
}

static void mpfr_free_array(size_t n, mpfr_t *rc)
{
	size_t i;
	
	for (i = 0; i < n; i++)
		mpfr_clear(rc[i]);
	free(rc);
}

static void mpfr_exp_d(mpfr_t out, const double x)
{
	mpfr_set_d(out, x, rnd);
	mpfr_exp(out, out, rnd);
}

/* ------------------------------------------------------------------------- */

static void calc_virials(struct expansion *vir, double coex, double kN,
								double chiN)
{
	vir->aa = -2. * (kN + 3.) / coex;
	vir->bb = 0.1;
	vir->ab = (chiN / coex) + .5 * (vir->aa + vir->bb);
	vir->aaa = 1.5 * (kN + 2.) / SQR(coex);
	vir->aab = vir->aaa;
	vir->abb = vir->abb;
	vir->bbb = 0.;
}
	
static int load_config_file(const char *fn, struct data *data)
{
	int ret;
	FILE *FH = fopen(fn, "r");
	if (FH == NULL)
		return error(ENOENT, "Can't open '%s' for reading", fn);

	/* 1. extensive number of bond angles used to calc lambda */ 
	ret = get_data_line(FH, "%d", 1, &data->nN);
	if (ret) return ret;
	debug("extensive nN=%d", data->nN);

	/* 2. binning of lambda */
	ret = get_data_line(FH, "%lg %lg %lg", 3, &data->lambda_min,
					&data->lambda_max, &data->lambda_bw);
	if (ret) return ret;
	data->lambda_bins = (int)((data->lambda_max - data->lambda_min) /
							       data->lambda_bw);
	debug("lambda binning: min=%lg max=%lg bw=%lg", data->lambda_min,
				       data->lambda_max, data->lambda_bw);

	/* 3. binning of rho_coex */
	ret = get_data_line(FH, "%lg %lg %lg", 3, &data->rho_min,
					&data->rho_max, &data->rho_bw);
	if (ret) return ret;
	data->rho_bins = (int)((data->rho_max - data->rho_min) / data->rho_bw);
	debug("reweighting target: rho=%lg..%lg in steps of %lg",
				data->rho_min, data->rho_max, data->rho_bw);

	/* 4. reference virial coefficients */
	double rho0, kN0, chiN0;
	ret = get_data_line(FH, "%lg %lg %lg", 3, &rho0, &kN0, &chiN0);
	if (ret) return ret;
	calc_virials(&data->vir0, rho0, kN0, chiN0);
	debug("reference values: rho=%lg kN=%lg chiN=%lg", rho0, kN0, chiN0);
	
	/* 5. reweighting target for rho_coex, kN, chiN */
	ret = get_data_line(FH, "%lg %lg %lg", 3, &data->rho_1, &data->kN_1,
								&data->chiN_1);
	if (ret) return ret;
	debug("reweighting target: rho_coex=%lg kN=%lg chiN=%lg",
				data->rho_1, data->kN_1, data->chiN_1);

	/* 6. binning for U */
	ret = get_data_line(FH, "%lg %lg %lg", 3, &data->U_min,
					&data->U_max, &data->U_bw);
	if (ret) return ret;
	data->U_bins = (int)((data->U_max - data->U_min) / data->U_bw);
	debug("U binning: U=%lg..%lg in steps of %lg",
				data->U_min, data->U_max, data->U_bw);

	/* 7. binning for A */
	ret = get_data_line(FH, "%lg %lg %lg", 3, &data->A_min,
					&data->A_max, &data->A_bw);
	if (ret) return ret;
	data->A_bins = (int)((data->A_max - data->A_min) / data->A_bw);
	debug("A binning: A=%lg..%lg in steps of %lg",
				data->A_min, data->A_max, data->A_bw);

	fclose(FH);
	return 0;
}

static int realloc_runs(struct data *data)
{
	data->runs_max += 0x1000;
	data->run = realloc(data->run, data->runs_max * sizeof(*data->run));
	if (data->run == NULL)
		return novm("data->run");

	return 0;
}

static int check_for_new_run(const char *fn, struct data *data,
				       struct expansion *vir, double pot[],
				       double l0, int *rid)
{
	int found = 0;
	int i, n, ret;

	/* check if at least one coefficient is nonzero */
	for (i = 0; i < COEFF_MAX; i++) {
		if (pot[i] != 0.)
			found = 1;
	}
	if (!found)
		return error(EINVAL, "It seems that '%s' does not contain "
				"an US potential.\n", fn);

	/* make sure that this is not a warm-up run */
//#ifdef OLD
	if (pot[1] > 4.) {
		(*rid) = -1;
		return 0;
	}
//#endif

	if (vir->aa == 0. && vir->ab == 0.)
		memcpy(vir, &data->vir0, sizeof(*vir));
		//calc_virials(vir, 17., 100., 30.);

	/* Now let's see if there is already such a histogram present */
	for (n = 0; n < data->runs_n; n++) {
		if (data->run[n].l0 == l0 &&
		memcmp(data->run[n].a, pot, COEFF_MAX * sizeof(*pot)) == 0 &&
		memcmp(&data->run[n].vir, vir, sizeof(*vir)) == 0)
			break;
	}

	if (n >= data->runs_n) {
		if (data->runs_n >= data->runs_max) {
			ret = realloc_runs(data);
			if (ret) return ret;
		}

		memcpy(data->run[n].a, pot, COEFF_MAX * sizeof(*pot));
		data->run[n].l0 = l0;
		memcpy(&data->run[n].vir, vir, sizeof(*vir));
		data->runs_n++;
		debug("New run (%d): v=%lg %lg %lg w=%lg %lg %lg %lg l0=%lg",
				n, vir->aa, vir->ab, vir->bb,
			       	vir->aaa, vir->aab, vir->abb, vir->bbb, l0);
		debug("              a_0=%lg a_1=%lg a_2=%lg a_3=%lg",
				pot[0], pot[1], pot[2], pot[3]);
		data->run[n].n = 0;
		data->run[n].f = 0.;
	}

	(*rid) = n;
	return 0;
}

static void parse_harmonic_potential(char *s, double *pot)
{
	double k;
	sscanf(strcasestr(s, "lk="), "lk=%lg", &k);
	pot[1] = .5 * k;
}

static void parse_full_potential(char *s, double *pot)
{
	char *s1 = strcasestr(s, "a_1=");
	sscanf(s1, "a_1=%lg a_2=%lg", &pot[0], &pot[1]); /* FIXME */
}

static int check_dp(struct dp *dp)
{
	struct expansion virnull;
	memset(&virnull, 0, sizeof(virnull));
	
	if (!isnormal(dp->lambda) || dp->lambda < .1 || dp->lambda > .9)
		return error(EINVAL, "dp->lambda=%lg", dp->lambda);
	if (!isnormal(dp->U) || dp->U < -1e6 || dp->U > 1e6)
		return error(EINVAL, "dp->U=%lg", dp->U);
	if (dp->run < 0 || dp->run > 1000)
		return error(EINVAL, "dp->run=%d", dp->run);
	if (memcmp(&virnull, &dp->r, sizeof(virnull)) == 0)
		return error(EINVAL, "dp->r is null");
	if (!isnormal(dp->r.aa) || dp->r.aa < -1e7 || dp->r.aa > 1e7)
		return error(EINVAL, "dp->r.aa=%lg", dp->r.aa);
	if (!isnormal(dp->r.ab) || dp->r.ab < -1e7 || dp->r.ab > 1e7)
		return error(EINVAL, "dp->r.ab=%lg", dp->r.ab);
	if (!isnormal(dp->r.bb) || dp->r.bb < -1e7 || dp->r.bb > 1e7)
		return error(EINVAL, "dp->r.bb=%lg", dp->r.bb);
	if (!isnormal(dp->r.aaa) || dp->r.aaa < -1e7 || dp->r.aaa > 1e7)
		return error(EINVAL, "dp->r.aaa=%lg", dp->r.aaa);
	if (!isnormal(dp->r.aab) || dp->r.aab < -1e7 || dp->r.aab > 1e7)
		return error(EINVAL, "dp->r.aab=%lg", dp->r.aab);
	if (!isnormal(dp->r.abb) || dp->r.abb < -1e7 || dp->r.abb > 1e7)
		return error(EINVAL, "dp->r.abb=%lg", dp->r.abb);
	if (!isnormal(dp->r.bbb) || dp->r.bbb < -1e7 || dp->r.bbb > 1e7)
		return error(EINVAL, "dp->r.bbb=%lg", dp->r.bbb);
	return 0;
}

static int add_dp(struct data *data, struct dp *dp)
{
	if (data->dp_count >= data->dp_max) {
		data->dp_max += 0x100000;
		data->dp = realloc(data->dp, data->dp_max * sizeof(*data->dp));
		if (data->dp == NULL) return novm("data->dp");
	}

	int ret = check_dp(dp);
	if (ret) return ret;

	memcpy(data->dp + data->dp_count++, dp, sizeof(*dp));
	return 0;
}

static int check_histofile(const char *fn, struct data *data)
{
	int rank, size, ret;
	struct expansion vir;
	struct dp dp;
	double pot[COEFF_MAX];
	char buf[LINE_MAX];
	FILE *FH = fopen(fn, "r");
	if (FH == NULL)
		return error(ENOENT, "Can't open '%s' for reading", fn);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	if (rank == 0)	
		printf("Accessing file '%s'...\n", fn);

	memset(&vir, 0, sizeof(vir));
	memset(&dp, 0, sizeof(dp));

	/* Try to get Umbrella Sampling potential coefficients */
	while (!feof(FH)) {
		char *s;

		fgets(buf, sizeof(buf), FH);
		if (feof(FH)) continue;

		s = strcasestr(buf, "COEFFICIENTS: ");
		if (s != NULL) {
			s = strcasestr(s, "v=");
			if (sscanf(s, "v=%lg %lg %lg w=%lg %lg %lg %lg",
				&vir.aa, &vir.ab, &vir.bb, &vir.aaa,
				&vir.aab, &vir.abb, &vir.bbb) != 7)
				return error(EINVAL, "coeff line wrong");
		}

		s = strcasestr(buf, "UMBRELLA SAMPLING: ");
		if (s != NULL) {
			double l0;
			memset(pot, 0, sizeof(*pot) * COEFF_MAX);
			s = strcasestr(s, "l0=");
			if (sscanf(s, "l0=%lg", &l0) != 1)
				return error(EINVAL, "No l0 found");

			/* power series for potential */
			if (strstr(s, "a_") != NULL)
				parse_full_potential(s, pot);
			else if (strstr(s, "lk=") != NULL)
				parse_harmonic_potential(s, pot);
			else
				return error(EINVAL, "Premature end of line");

			ret = check_for_new_run(fn, data, &vir, pot,
								l0, &dp.run);
			if (ret) return ret;
		}

		s = strcasestr(buf, " sum=");
		if (s != NULL) {
			sscanf(s, " sum=%lg", &dp.U);
		}
		
		s = strcasestr(buf, "sum rho");
		if (s != NULL) {
			sscanf(s, "sum rho: AA=%lg AB=%lg BB=%lg AAA=%lg "
				"AAB=%lg ABB=%lg BBB=%lg",
				&dp.r.aa, &dp.r.ab, &dp.r.bb,
				&dp.r.aaa, &dp.r.aab, &dp.r.abb, &dp.r.bbb);
		}

		if (strstr(buf, "l=") != NULL && strstr(buf, " P=") != NULL) {
			double lx, ly, lz;
			s = strstr(buf, "l=");
			sscanf(s, "l=%lg %lg %lg", &lx, &ly, &lz);
			dp.A = ly * lz;
		}

		s = strcasestr(buf, "Sxx=");
		if (s != NULL && dp.run != -1) {
			s = strchr(s, '=') + 1;
			dp.lambda = atof(s);

			if (data->dpc++ % size == rank) {
				ret = add_dp(data, &dp);
				if (ret) return ret;
			}
			data->run[dp.run].n++;
		}
	}

	fclose(FH);
	return 0;
}

static int load_list_file(const char *fn, struct data *data)
{
	int ret;
	int filecount = 0;
	FILE *FH = fopen(fn, "r");
	if (FH == NULL)
		return error(ENOENT, "Can't open '%s' for reading", fn);

	while (!feof(FH)) {
		char buf[LINE_MAX];
		*buf=0x0;
		get_data_line(FH, "%s", 1, buf);
		if (!feof(FH)) {
			ret = check_histofile(buf, data);
			if (ret) return ret;
			filecount++;
		}
	}
	debug("%d files, %d total data points, %.3f%% in this process",
		filecount, data->dpc, 100. * data->dp_count / data->dpc);
	fclose(FH);
	return 0;
}

static int alloc_work_buffers(struct data *data)
{
	int ret;

	ret = mpfr_calloc_array(SQR(data->runs_n), &data->G);
	if (ret) return ret;
	
	ret = mpfr_calloc_array(data->runs_n, &data->F);
	if (ret) return ret;

	return 0;
}

static int setup(int argc, char **argv, struct data *data)
{
	int ret;

	memset(data, 0, sizeof(*data));

	if (argc < 3) {
		fprintf(stderr, "Syntax: dwham [config.dat] [list]\n\n");
		return 1;
	}

	debug("config file='%s' list file='%s'", argv[1], argv[2]);
	
	ret = load_config_file(argv[1], data);
	if (ret) return ret;

	ret = load_list_file(argv[2], data);
	if (ret) return ret;

	ret = alloc_work_buffers(data);
	if (ret) return ret;

	return 0;
}
	
static void free_data(struct data *data)
{
	free(data->run);
	free(data->dp);
	mpfr_free_array(data->runs_n, data->G);
	mpfr_free_array(data->runs_n, data->F);
}

/* ------------------------------------------------------------------------- */

static double hamilton(const struct expansion *vir, const struct expansion *ir)
{
	double e = 0.;

	e += .5 * vir->aa * ir->aa;
	e += .5 * vir->ab * ir->ab;
	e += .5 * vir->bb * ir->bb;
	e += 1. / 3. * vir->aaa * ir->aaa;
	e += 1. / 3. * vir->aab * ir->aab;
	e += 1. / 3. * vir->abb * ir->abb;
	e += 1. / 3. * vir->bbb * ir->bbb;

	return e;
}

static double dhamilton(struct expansion *vir_1, struct expansion *vir_0,
							const struct dp *dp)
{
	return hamilton(vir_1, &dp->r) - hamilton(vir_0, &dp->r);
}

/* the potential is gauged, so that the constant term is zero */
static double pot(struct data *data, const int i, const double l)
{
	double ret = 0.;
	int j;

	for (j = COEFF_MAX; j > 0; j--) {
		ret += data->run[i].a[j - 1];
		ret *= l - data->run[i].l0;
	}
	return ret * data->nN;
}

static double get_de(struct data *data, const struct dp *dp, const int run)
{
	return dhamilton(&data->run[run].vir, &data->vir0, dp)
						 + pot(data, run, dp->lambda);
}

static int reduce_cache(const int n, mpfr_t *c)
{
	int rank, size, i, j;
	mpfr_t tmp;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	mpfr_init(tmp);

	char *buf = calloc(size * LINE_MAX, sizeof(*buf));
	if (buf == NULL) return novm("reduce_cache");

	for (i = 0; i < n; i++) {
		mpfr_snprintf(buf + rank * LINE_MAX, LINE_MAX, "%Re", c[i]);
		if (!mpfr_number_p(c[i]))
			return error(EDOM, "c[%d]=%s is not a number",
					       i, buf + rank * LINE_MAX);
		
		MPI_Allgather(MPI_IN_PLACE, 0, MPI_CHAR, buf, LINE_MAX,
			       			MPI_CHAR, MPI_COMM_WORLD);
		mpfr_set_d(c[i], 0., rnd);
		for (j = 0; j < size; j++) {
			mpfr_set_str(tmp, buf + j * LINE_MAX, 0, rnd);
			mpfr_add(c[i], c[i], tmp, rnd);
		}
		if (!mpfr_number_p(c[i]))
			return error(EDOM, "c[%d] is not a number", i);
	}

	mpfr_clear(tmp);
	free(buf);
	return 0;
}

enum action
{
	CALC_F,
	CALC_FG
};

static int calc_F_G(struct data *data, const gsl_vector *f,
							const enum action act)
{
	int i, k, l, ret;
	mpfr_t z, tmp, *x, *y;
	
	data->run[0].f = 0.;
	for (i = 1; i < data->runs_n; i++)
		data->run[i].f = gsl_vector_get(f, i - 1);

	mpfr_inits(z, tmp, NULL);
	ret = mpfr_calloc_array(data->runs_n, &x);
	if (ret) return ret;
	ret = mpfr_calloc_array(data->runs_n, &y);
	if (ret) return ret;

	for (k = 0; k < data->runs_n; k++)
		mpfr_set_d(data->F[k], 0., rnd);
	if (act == CALC_FG) {
		for (k = 0; k < SQR(data->runs_n); k++)
			mpfr_set_d(data->G[k], 0., rnd);
	}

	for (i = 0; i < data->dp_count; i++) {
		mpfr_set_d(z, 0., rnd);
		for (k = 0; k < data->runs_n; k++) {
			double x_exp = get_de(data, data->dp + i, k);
			double y_exp = x_exp - data->run[k].f;
			mpfr_exp_d(x[k], -1. * x_exp);
			mpfr_exp_d(y[k], -1. * y_exp);
			mpfr_mul_si(tmp, y[k], data->run[k].n, rnd);
			mpfr_add(z, z, tmp, rnd);
		}
		for (k = 0; k < data->runs_n; k++) {
			mpfr_div(tmp, x[k], z, rnd);
			mpfr_add(data->F[k], data->F[k], tmp, rnd);
			for (l = 0; l < data->runs_n && act == CALC_FG; l++) {
				mpfr_mul(tmp, x[k], y[l], rnd);
				mpfr_div(tmp, tmp, z, rnd);
				mpfr_div(tmp, tmp, z, rnd);
				mpfr_add(data->G[k * data->runs_n + l],
						data->G[k * data->runs_n + l],
						tmp, rnd);
			}
		}
		if (i % 5000 == 0)
			debug("%d/%d (%lg%%)", i, data->dp_count,
					100. * i / data->dp_count);
	}
	
	ret = reduce_cache(data->runs_n, data->F);
	if (ret) return ret;

	if (act == CALC_FG) {
		ret = reduce_cache(SQR(data->runs_n), data->G);
		if (ret) return ret;
	}

	mpfr_clears(z, tmp, NULL);
	mpfr_free_array(data->runs_n, x);
	mpfr_free_array(data->runs_n, y);

	return 0;
}

static void update_f(struct data *data, gsl_vector *g)
{
	int k;
	mpfr_t tmp;

	mpfr_init(tmp);
	for (k = 1; k < data->runs_n; k++) {
		double *ge = gsl_vector_ptr(g, k - 1);
		(*ge) = data->run[k].f;
		mpfr_div(tmp, data->F[k], data->F[0], rnd);
		mpfr_log(tmp, tmp, rnd);
		(*ge) += mpfr_get_d(tmp, rnd);
	}
	mpfr_clear(tmp);
	MPI_Bcast(g->data, g->size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

static void update_df(struct data *data, gsl_matrix *J)
{
	int k, l;
	mpfr_t tmp, res;

	gsl_matrix_set_identity(J);

	mpfr_inits(tmp, res, NULL);
	for (l = 1; l < data->runs_n; l++) {
		mpfr_div(tmp, data->G[0 + l], data->F[0], rnd);
		for (k = 1; k < data->runs_n; k++) {
			double *Je = gsl_matrix_ptr(J, k - 1, l - 1);

			mpfr_div(res, data->G[k * data->runs_n + l],
							data->F[k], rnd);
			mpfr_sub(res, tmp, res, rnd);
			mpfr_mul_si(res, res, data->run[l].n, rnd);
			(*Je) += mpfr_get_d(res, rnd);
		}
	}
	mpfr_clears(tmp, res, NULL);
	MPI_Bcast(J->data, J->size1 * J->size2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

static int free_energy_f(const gsl_vector *f, void *params, gsl_vector *g)
{
	int rank;
	struct data *data = params;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0) printf("Entering f...\n");

	int ret = calc_F_G(data, f, CALC_F);
	if (ret) return GSL_ENOMEM;

	update_f(data, g);
	return GSL_SUCCESS;
}

static int free_energy_df(const gsl_vector *f, void *params, gsl_matrix *J)
{
	int rank;
	struct data *data = params;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0) printf("Entering df...\n");
	
	int ret = calc_F_G(data, f, CALC_FG);
	if (ret) return GSL_ENOMEM;

	update_df(data, J);
	return GSL_SUCCESS;
}

static int free_energy_fdf(const gsl_vector *f, void *params, gsl_vector *g,
								gsl_matrix *J)
{
	int rank;
	struct data *data = params;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0) printf("Entering fdf...\n");
	
	int ret = calc_F_G(data, f, CALC_FG);
	if (ret) return GSL_ENOMEM;

	update_df(data, J);
	update_f(data, g);
	return GSL_SUCCESS;
}

/* ------------------------------------------------------------------------- */

static void print_state(int iter, gsl_multiroot_fdfsolver *s)
{
	int i, rank;
	int n = s->x->size;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (rank != 0) return;
	
	printf("iter = %d\n", iter);
	for (i = 0; i < n; i++) {
		printf("f_%04d = %lg, g_%04d = %lg\n", 
			i, gsl_vector_get(s->x, i),
			i, gsl_vector_get(s->f, i));
	}
}

static void load_initial(const char *fn, gsl_vector *x)
{
	FILE *FH = fopen(fn, "rb");

	if (FH == NULL) {
		gsl_vector_set_all(x, 0.);
		return;
	}

	if (gsl_vector_fread(FH, x) != GSL_SUCCESS)
		gsl_vector_set_all(x, 0.);
		
	fclose(FH);
}

static void save_initial(const char *fn, gsl_vector *x)
{
	FILE *FH = fopen(fn, "wb");
	gsl_vector_fwrite(FH, x);
	fclose(FH);
}

static int produce_final_rho(const char *fn, struct data *data)
{
	int i, j, rank;
	FILE *FH;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank != 0) return 0;

	FH = fopen(fn, "w");
	if (FH == NULL)
		return error(EIO, "Couldn't write to '%s'", fn);

	fprintf(FH, "# no dpn f u vAA vAB vBB wAAA wAAB wABB wBBB l0 a1 a2..\n");

	for (i = 0; i < data->runs_n; i++) {
		fprintf(FH, "%d %d %lg %lg ", i, data->run[i].n, data->run[i].f, data->run[i].u);
		fprintf(FH, "%lg %lg %lg %lg %lg %lg %lg ",
				data->run[i].vir.aa, data->run[i].vir.ab,
				data->run[i].vir.bb, data->run[i].vir.aaa,
				data->run[i].vir.aab, data->run[i].vir.abb,
				data->run[i].vir.bbb);
		fprintf(FH, "%lg ", data->run[i].l0);
		for (j = 0; j < COEFF_MAX; j++)
			fprintf(FH, " %lg", data->run[i].a[j]);
		fprintf(FH, "\n");
	}

	fclose(FH);
	return 0;
}

#ifndef NEW_VERSION
static int free_energies(struct data *data)
{
	const gsl_multiroot_fdfsolver_type *T = gsl_multiroot_fdfsolver_hybridj;
	gsl_multiroot_fdfsolver *s =
			gsl_multiroot_fdfsolver_alloc(T, data->runs_n - 1);
	gsl_multiroot_function_fdf fdf;
	int iter = 0;
	int status;
	const char f_fn[] = ".f.dat";
	gsl_vector *initial = gsl_vector_alloc(data->runs_n - 1);

	fdf.f = free_energy_f;
	fdf.df = free_energy_df;
	fdf.fdf = free_energy_fdf;
	fdf.n = data->runs_n - 1;
	fdf.params = data;

	load_initial(f_fn, initial);

	gsl_multiroot_fdfsolver_set(s, &fdf, initial);

	do {
		iter++;
		status = gsl_multiroot_fdfsolver_iterate(s);

		print_state(iter, s);
		save_initial(f_fn, s->x);
	
		if (status)
			break;

		status = gsl_multiroot_test_residual(s->f, 1e-7);
	} while (status == GSL_CONTINUE && iter < 500);
	
	save_initial(f_fn, s->x);

	data->run[0].f = 0.;
	for (iter = 1; iter < data->runs_n; iter++)
		data->run[iter].f = gsl_vector_get(s->x, iter - 1);

	printf("status = %s\n", gsl_strerror(status));
	produce_final_rho("results.dat", data);

	gsl_multiroot_fdfsolver_free(s);
	gsl_vector_free(initial);

	return (status == GSL_SUCCESS) ? 0 : 1;
}
#else

static int reduce_mpfr(mpfr_t x)
{
	int rank, size, j;
	char *s;
	mpfr_t tmp;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	mpfr_init(tmp);

	s = calloc(size * LINE_MAX, sizeof(*s));
	if (s == NULL) return novm("s");

	mpfr_snprintf(s + rank * LINE_MAX, LINE_MAX, "%Re", x);
	if (!mpfr_number_p(x))
		return error(EDOM, "x=%s is not a number",
					       s + rank * LINE_MAX);
		
	MPI_Allgather(MPI_IN_PLACE, 0, MPI_CHAR, s, LINE_MAX,
			       			MPI_CHAR, MPI_COMM_WORLD);
	mpfr_set_d(x, 0., rnd);
	for (j = 0; j < size; j++) {
		mpfr_set_str(tmp, s + j * LINE_MAX, 0, rnd);
		mpfr_add(x, x, tmp, rnd);
	}
	if (!mpfr_number_p(x))
		return error(EDOM, "x is not a number");

	mpfr_clear(tmp);
	free(s);

	return 0;
}

static void print_status(int iter, mpfr_t *ebf, struct data *data)
{
	int i, rank;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	mpfr_t tmp;

	if (rank != 0) return;

	mpfr_init(tmp);	

	
	printf("iter = %d\n", iter);
	for (i = 0; i < data->runs_n; i++) {
		mpfr_log(tmp, ebf[i], rnd);
		
		printf("f_%04d = %lg\n", 
			i, mpfr_get_d(tmp, rnd));
	}

	mpfr_clear(tmp);
}

static int free_energies(struct data *data)
{
	int i, j, k, l, ret;
	int conv;
	int iter = 0, maxiter = 1000;
	double tol = 1e-5;
	mpfr_t ebfk, bottom, tmp, *ebw, *ebf, *ebf2, *fact;
	
	printf("1\n");

	/* pre computing all exponentials */
	ret = mpfr_calloc_array(data->runs_n * data->dp_count, &ebw);
	if (ret) {
		printf("Sorry. Not enough memory available. Try to increase"
		" the number of MPI processes.\n");
		printf("data->dp_count=%d data->runs_n=%d\n",
						data->dp_count, data->runs_n);
		return ret;
	}

	ret = mpfr_calloc_array(data->runs_n, &fact);
	if (ret) return ret;
	ret = mpfr_calloc_array(data->runs_n, &ebf);
	if (ret) return ret;
	ret = mpfr_calloc_array(data->runs_n, &ebf2);
	if (ret) return ret;

	for (l = 0; l < data->dp_count; l++) {
		for (i = 0; i < data->runs_n; i++) {
			double bw = -1. * get_de(data, data->dp + l, i);
			mpfr_exp_d(ebw[l * data->runs_n + i], bw);
		}
	}

	mpfr_inits(ebfk, bottom, tmp, NULL);

	for (k = 0; k < data->runs_n; k++) {
		mpfr_set_d(ebf[k], 1., rnd);
		mpfr_mul_si(fact[k], ebf[k], data->run[k].n, rnd);
	}

	printf("2\n");

	do {
		iter++;
		/* calculate free energies for each run */
		for (k = 0; k < data->runs_n; k++) {
			mpfr_set_d(ebfk, 0., rnd);

			for (i = 0; i < data->runs_n; i++) {
				for (l = 0; l < data->dp_count; l++) {
					if (data->dp[l].run != i) continue;
					mpfr_set_d(bottom, 0., rnd);

					for (j = 0.; j < data->runs_n; j++) {
						mpfr_mul(tmp, ebw[l * data->runs_n + j], fact[j], rnd);
						mpfr_add(bottom, bottom, tmp, rnd);
					}
					mpfr_div(tmp, ebw[l * data->runs_n + k], bottom, rnd);
					mpfr_add(ebfk, ebfk, tmp, rnd);
				}
			}
			reduce_mpfr(ebfk);

			mpfr_set(ebf2[k], ebfk, rnd);
			mpfr_mul(tmp, ebf[0], ebfk, rnd);
			mpfr_d_div(ebf[k], 1., tmp, rnd);
			mpfr_mul_si(fact[k], ebf[k], data->run[k].n, rnd);
		}

		/* test if converged */
		conv = 1;
		double max = 0.;
		for (k = 0; k < data->runs_n; k++) {
			mpfr_mul(tmp, ebf[k], ebf2[k], rnd);
			mpfr_log(tmp, tmp, rnd);
			double delta = fabs(mpfr_get_d(tmp, rnd));
			if (delta > max) max = delta;
			if (delta >= tol) conv = 0;
			mpfr_div(ebf[k], ebf2[0], ebf2[k], rnd);
		}
		//printf("max=%lg\n", max);
		print_status(iter, ebf, data);
	} while (iter < maxiter && !conv);
	print_status(iter, ebf, data);

	for (i = 0; i < data->runs_n; i++) {
		mpfr_log(tmp, ebf[i], rnd);
		data->run[i].f = mpfr_get_d(tmp, rnd);
	}
	produce_final_rho("results.dat", data);

	mpfr_clears(ebfk, bottom, tmp, NULL);
	mpfr_free_array(data->runs_n * data->dp_count, ebw);
	mpfr_free_array(data->runs_n, fact);
	mpfr_free_array(data->runs_n, ebf);
	mpfr_free_array(data->runs_n, ebf2);

	return 0;
}
#endif

static int average_energies(struct data *data)
{
	int i, count[data->runs_n];
	double u[data->runs_n];
	
	for (i = 0; i < data->runs_n; i++) {
		u[i] = 0.;
		count[i] = 0;
	}

	for (i = 0; i < data->dp_count; i++) {
		struct dp *dp = data->dp + i;
		u[dp->run] += dp->U;
		count[dp->run]++;
	}

	MPI_Allreduce(MPI_IN_PLACE, count, data->runs_n, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, u, data->runs_n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	for (i = 0; i < data->runs_n; i++)
		data->run[i].u = u[i] / count[i];

	return 0;
}

/* ------------------------------------------------------------------------- */

static int print_histo_1d(const char *fn, const double l, const double dx,
							const int n, mpfr_t *y)
{
	FILE *FH = fopen(fn, "w");
	mpfr_t tmp;
	int i;

	if (FH == NULL)
		return error(EIO, "Can't open '%s' for writing", fn);

	mpfr_init(tmp);

	for (i = 0; i < n; i++) {
		mpfr_log(tmp, y[i], rnd);
		double x = l + (i + .5) * dx;
		double log_y = -1. * mpfr_get_d(tmp, rnd);
		mpfr_fprintf(FH, "%lg %lg %Re\n", x, log_y, y[i]);
	}

	mpfr_clear(tmp);
	fclose(FH);
	return 0;
}

static int print_histo_2d(const char *fn, const double lx, const double dx,
		const int nx, const double ly, const double dy, const int ny,
		mpfr_t *z)
{
	FILE *FH = fopen(fn, "w");
	mpfr_t tmp;
	int i, j;

	if (FH == NULL)
		return error(EIO, "Can't open '%s' for writing", fn);

	mpfr_init(tmp);

	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			double x = lx + (i + .5) * dx;
			double y = ly + (j + .5) * dy;
			mpfr_log(tmp, z[i * ny + j], rnd);
			double log_z = -1. * mpfr_get_d(tmp, rnd);
			mpfr_fprintf(FH, "%lg %lg %lg %Re\n", x, y, log_z,
							       z[i * ny + j]);
		}
		fprintf(FH, "\n");
	}

	mpfr_clear(tmp);
	fclose(FH);
	return 0;
}

/* Free Energy as a function of rho_coex */
static int rw_f_rho(struct data *data, mpfr_t *rho, mpfr_t *wlr, mpfr_t *wur,
								mpfr_t *wla)
{
	struct expansion vir;
	mpfr_t tmp, z;
	int i, k, b, ret;
	
	mpfr_inits(tmp, z, NULL);

	/* target bin for reweighting at a specific rho_coex */
	const int trb = (int)((data->rho_1 - data->rho_min) / data->rho_bw);

	/* the normalization constant C is still missing */
	for (i = 0; i < data->dp_count; i++) {
		const struct dp *dp = data->dp + i;

		mpfr_set_d(z, 0., rnd);
		for (k = 0; k < data->runs_n; k++) {
			double y_exp = get_de(data, dp, k) - data->run[k].f;
			mpfr_exp_d(tmp, -1. * y_exp);
			mpfr_mul_si(tmp, tmp, data->run[k].n, rnd);
			mpfr_add(z, z, tmp, rnd);
		}


		for (k = 0; k < data->rho_bins; k++) {
			double coex = data->rho_min + (k + .5) * data->rho_bw;
			calc_virials(&vir, coex, data->kN_1, data->chiN_1);
			double x_exp = dhamilton(&vir, &data->vir0, dp);
			mpfr_exp_d(tmp, -1. * x_exp);
			mpfr_div(tmp, tmp, z, rnd);

			/* Basic Free Energy as a function of rho_coex */
			mpfr_add(rho[k], rho[k], tmp, rnd);

			/* Contribution to w(lambda, rho) */
			b = (int)((dp->lambda - data->lambda_min) / data->lambda_bw);
			if (b >= 0 && b < data->lambda_bins) {
				b = b * data->rho_bins + k; 
				mpfr_add(wlr[b], wlr[b], tmp, rnd);
			}
			
			/* Contribution to w(U, rho) */
			b = (int)((dp->U - data->U_min) / data->U_bw);
			if (b >= 0 && b < data->U_bins) {
				b = b * data->rho_bins + k;
				mpfr_add(wur[b], wur[b], tmp, rnd);
			}

			/* Contribution to w(lambda, A, rho = rho_1) */
			b = (int)((dp->A - data->A_min) / data->A_bw);
			int b2 = (int)((dp->lambda - data->lambda_min) / data->lambda_bw);
			if (b >= 0 && b < data->A_bins && b2 >= 0 &&
					b2 < data->lambda_bins && k == trb)
			{
				b = b2 * data->A_bins + b;
				mpfr_add(wla[b], wla[b], tmp, rnd);
			}
					
		}
		if (i % 2000 == 0)
			debug("RW %d/%d (%lg%%)", i, data->dp_count,
					100. * i / data->dp_count);
	}
	
	ret = reduce_cache(data->rho_bins, rho);
	if (ret) return ret;
	ret = reduce_cache(data->rho_bins * data->lambda_bins, wlr);
	if (ret) return ret;
	ret = reduce_cache(data->rho_bins * data->U_bins, wur);
	if (ret) return ret;
	ret = reduce_cache(data->lambda_bins * data->A_bins, wla);
	if (ret) return ret;

	for (i = 0; i < data->rho_bins; i++) {
		mpfr_d_div(tmp, 1., rho[i], rnd);
		if (mpfr_number_p(tmp)) {
			for (k = 0; k < data->lambda_bins; k++) {
				b = k * data->rho_bins + i;
				mpfr_div(wlr[b], wlr[b], rho[i], rnd);
			}
			for (k = 0; k < data->U_bins; k++) {
				b = k * data->rho_bins + i;
				mpfr_div(wur[b], wur[b], rho[i], rnd);
			}
			if (i == trb) { /* for reweighting at a specific rho */
				for (k = 0; k < data->lambda_bins * data->A_bins; k++) {
					mpfr_div(wla[k], wla[k], rho[i], rnd);
				}
			}
		}
	}
	
	mpfr_clears(tmp, z, NULL);
	return 0;
}

/* Heat Capacity */
static int rw_heat_capacity(const char *fn, struct data *data, mpfr_t *wur)
{
	mpfr_t tmp, z, mean_U, mean_U2;
	int i, k;
	FILE *FH;

	mpfr_inits(tmp, z, mean_U, mean_U2, NULL);

	FH = fopen(fn, "w");
	if (FH == NULL)
		return error(EIO, "Cannot open '%s' for writing", fn);
	fprintf(FH, "# rho <U>/kT C/k\n");
	
	for (i = 0; i < data->rho_bins; i++) {
		mpfr_set_d(mean_U, 0., rnd);
		mpfr_set_d(mean_U2, 0., rnd);
		mpfr_set_d(z, 0., rnd);

		for (k = 0; k < data->U_bins; k++) {
			const double U = data->U_min + (k + .5) * data->U_bw;
			mpfr_set(tmp, wur[k * data->rho_bins + i], rnd);
			mpfr_add(z, z, tmp, rnd); /* normalisation */
			mpfr_mul_d(tmp, tmp, U, rnd);
			mpfr_add(mean_U, mean_U, tmp, rnd); /* 1st moment */
			mpfr_mul_d(tmp, tmp, U, rnd);
			mpfr_add(mean_U2, mean_U2, tmp, rnd); /* 2nd moment */
		}
		mpfr_d_div(tmp, 1., z, rnd);
		if (mpfr_number_p(tmp)) {
			mpfr_mul(mean_U, mean_U, tmp, rnd);
			mpfr_mul(mean_U2, mean_U2, tmp, rnd);
		} else {
			mpfr_set_d(mean_U, 0., rnd);
			mpfr_set_d(mean_U2, 0., rnd);
		}

		double c1 = mpfr_get_d(mean_U, rnd);
		double c2 = mpfr_get_d(mean_U2, rnd);

		fprintf(FH, "%lg %lg %lg\n",
				data->rho_min + (i + .5) * data->rho_bw, 
				c1, c2 - c1 * c1);
	}
	mpfr_clears(tmp, z, mean_U, mean_U2, NULL);
	fclose(FH);
	return 0;
}

/* Heat Capacity as a function of lambda at a specific rho */
static int rw_heat_capacity2(const char *fn, struct data *data, mpfr_t *wlu)
{
	mpfr_t tmp, z, mean_U, mean_U2;
	int i, k;
	FILE *FH;

	mpfr_inits(tmp, z, mean_U, mean_U2, NULL);

	FH = fopen(fn, "w");
	if (FH == NULL)
		return error(EIO, "Cannot open '%s' for writing", fn);
	fprintf(FH, "# lambda <U>/kT C/k\n# rho_1 = %lg\n", data->rho_1);
	
	for (i = 0; i < data->lambda_bins; i++) {
		mpfr_set_d(mean_U, 0., rnd);
		mpfr_set_d(mean_U2, 0., rnd);
		mpfr_set_d(z, 0., rnd);

		for (k = 0; k < data->U_bins; k++) {
			const double U = data->U_min + (k + .5) * data->U_bw;
			mpfr_set(tmp, wlu[i * data->U_bins + k], rnd);
			mpfr_add(z, z, tmp, rnd); /* normalisation */
			mpfr_mul_d(tmp, tmp, U, rnd);
			mpfr_add(mean_U, mean_U, tmp, rnd); /* 1st moment */
			mpfr_mul_d(tmp, tmp, U, rnd);
			mpfr_add(mean_U2, mean_U2, tmp, rnd); /* 2nd moment */
		}
		mpfr_d_div(tmp, 1., z, rnd);
		if (mpfr_number_p(tmp)) {
			mpfr_mul(mean_U, mean_U, tmp, rnd);
			mpfr_mul(mean_U2, mean_U2, tmp, rnd);
		} else {
			mpfr_set_d(mean_U, 0., rnd);
			mpfr_set_d(mean_U2, 0., rnd);
		}

		double c1 = mpfr_get_d(mean_U, rnd);
		double c2 = mpfr_get_d(mean_U2, rnd);

		fprintf(FH, "%lg %lg %lg\n",
				data->lambda_min + (i + .5) * data->lambda_bw,
				c1, c2 - c1 * c1);
	}
	mpfr_clears(tmp, z, mean_U, mean_U2, NULL);
	fclose(FH);
	return 0;
}

/* Area fluctuations as a function of lambda at specific rho_coex */
static int rw_area(const char *fn, struct data *data, mpfr_t *wla)
{
	mpfr_t tmp, z, mean_A, mean_A2;
	int i, k;
	FILE *FH;

	mpfr_inits(tmp, z, mean_A, mean_A2, NULL);

	FH = fopen(fn, "w");
	if (FH == NULL)
		return error(EIO, "Cannot open '%s' for writing", fn);
	fprintf(FH, "# lambda <A>/s^2 <A^2>/s^4\n# rho_1 = %lg\n", data->rho_1);
	
	for (i = 0; i < data->lambda_bins; i++) {
		mpfr_set_d(mean_A, 0., rnd);
		mpfr_set_d(mean_A2, 0., rnd);
		mpfr_set_d(z, 0., rnd);

		for (k = 0; k < data->A_bins; k++) {
			const double A = data->A_min + (k + .5) * data->A_bw;
			mpfr_set(tmp, wla[i * data->A_bins + k], rnd);
			mpfr_add(z, z, tmp, rnd); /* normalisation */
			mpfr_mul_d(tmp, tmp, A, rnd);
			mpfr_add(mean_A, mean_A, tmp, rnd); /* 1st moment */
			mpfr_mul_d(tmp, tmp, A, rnd);
			mpfr_add(mean_A2, mean_A2, tmp, rnd); /* 2nd moment */
		}
		mpfr_d_div(tmp, 1., z, rnd);
		if (mpfr_number_p(tmp)) {
			mpfr_mul(mean_A, mean_A, tmp, rnd);
			mpfr_mul(mean_A2, mean_A2, tmp, rnd);
		} else {
			mpfr_set_d(mean_A, 0., rnd);
			mpfr_set_d(mean_A2, 0., rnd);
		}

		double c1 = mpfr_get_d(mean_A, rnd);
		double c2 = mpfr_get_d(mean_A2, rnd);

		fprintf(FH, "%lg %lg %lg\n",
				data->lambda_min + (i + .5) * data->lambda_bw,
				c1, c2 - c1 * c1);
	}
	mpfr_clears(tmp, z, mean_A, mean_A2, NULL);
	fclose(FH);
	return 0;
}

/* calculate order parameter S as a function of rho_coex */
static int rw_order_param(const char *fn, struct data *restrict data,
							mpfr_t *restrict wlr)
{
	mpfr_t tmp, z, mean_S, mean_S2;
	int i, k;
	FILE *FH;

	mpfr_inits(tmp, z, mean_S, mean_S2, NULL);

	FH = fopen(fn, "w");
	if (FH == NULL)
		return error(EIO, "Cannot open '%s' for writing", fn);
	fprintf(FH, "# rho <S> <DS^2>\n");
	
	for (i = 0; i < data->rho_bins; i++) {
		mpfr_set_d(mean_S, 0., rnd);
		mpfr_set_d(mean_S2, 0., rnd);
		mpfr_set_d(z, 0., rnd);

		for (k = 0; k < data->lambda_bins; k++) {
			const double S = 0 + (k + .5) * data->lambda_bw;
			mpfr_set(tmp, wlr[k * data->lambda_bins + i], rnd);
			mpfr_add(z, z, tmp, rnd); /* normalisation */
			mpfr_mul_d(tmp, tmp, S, rnd);
			mpfr_add(mean_S, mean_S, tmp, rnd); /* 1st moment */
			mpfr_mul_d(tmp, tmp, S, rnd);
			mpfr_add(mean_S2, mean_S2, tmp, rnd); /* 2nd moment */
		}
		mpfr_d_div(tmp, 1., z, rnd);
		if (mpfr_number_p(tmp)) {
			mpfr_mul(mean_S, mean_S, tmp, rnd);
			mpfr_mul(mean_S2, mean_S2, tmp, rnd);
		} else {
			mpfr_set_d(mean_S, 0., rnd);
			mpfr_set_d(mean_S2, 0., rnd);
		}

		double c1 = mpfr_get_d(mean_S, rnd);
		double c2 = mpfr_get_d(mean_S2, rnd);

		fprintf(FH, "%lg %lg %lg\n",
				data->rho_min + (i + .5) * data->rho_bw, 
				c1, sqrt(c2 - c1 * c1));
	}
	mpfr_clears(tmp, z, mean_S, mean_S2, NULL);
	fclose(FH);
	return 0;
}

static int reweighting(struct data *data)
{
	int ret, rank;
	mpfr_t *rho, *wlr, *wur, *wla;

	ret = mpfr_calloc_array(data->rho_bins, &rho);
	if (ret) return ret;
	ret = mpfr_calloc_array(data->rho_bins * data->lambda_bins, &wlr);
	if (ret) return ret;
	ret = mpfr_calloc_array(data->rho_bins * data->U_bins, &wur);
	if (ret) return ret;
	ret = mpfr_calloc_array(data->lambda_bins * data->A_bins, &wla);
	if (ret) return ret;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (rank == 0) {
		/* Free Energy as a function of rho_coex */
		printf("Calculating F(rho_coex)...\n");
		ret = rw_f_rho(data, rho, wlr, wur, wla);
		if (ret) return ret;
		ret = print_histo_1d("f_rho.dat", data->rho_min, data->rho_bw,
						       data->rho_bins, rho);
		if (ret) return ret;

		/* Free Energy as a function of rho_coex and lambda */
		printf("Calculating w(lambda, rho_coex)...\n");
		ret = print_histo_2d("f_lambda_rho.dat", data->lambda_min,
				data->lambda_bw, data->lambda_bins,
				data->rho_min, data->rho_bw, data->rho_bins,
				wlr);
		if (ret) return ret;

		/* heat capacity as a function of rho_coex */
		printf("Calculating C(rho_coex)...\n");
		ret = rw_heat_capacity("c_gamma.dat", data, wur);
		if (ret) return ret;

		/* Area fluctuations as function of lambda at specific rho */
		printf("Calculating A(lambda,rho_coex=rho_1)...\n");
		ret = rw_area("a_lambda.dat", data, wla);
		if (ret) return ret;

		/* Order Parameter S as a function of rho_coex */
		printf("Calculating S(rho_coex)...\n");
		ret = rw_order_param("s_rho.dat", data, wlr);
		if (ret) return ret;
		
		printf("Reweighting done!\n");
	} else {
		ret = rw_f_rho(data, rho, wlr, wur, wla);
		if (ret) return ret;
	}

	mpfr_free_array(data->rho_bins, rho);
	mpfr_free_array(data->rho_bins * data->lambda_bins, wlr);
	mpfr_free_array(data->rho_bins * data->U_bins, wur);
	mpfr_free_array(data->lambda_bins * data->A_bins, wla);

	return 0;
}

/* ------------------------------------------------------------------------- */

extern int main(int argc, char **argv)
{
	int ret;
	struct data data;


#if defined(DEBUG) && defined(FE_NOMASK_ENV)
	//feenableexcept(FE_UNDERFLOW | FE_OVERFLOW | FE_INVALID | FE_DIVBYZERO);
	feenableexcept(FE_OVERFLOW | FE_INVALID | FE_DIVBYZERO);
#endif

	MPI_Init(&argc, &argv);

	ret = setup(argc, argv, &data);
	if (ret) goto error;

	ret = average_energies(&data);
	if (ret) goto error;

	ret = free_energies(&data);
	if (ret) goto error;

	ret = reweighting(&data);
	if (ret) goto error;

	free_data(&data);

	MPI_Finalize();
	mpfr_free_cache();

	return EXIT_SUCCESS;

error:
	MPI_Abort(MPI_COMM_WORLD, ret);
	return EXIT_FAILURE;
}

