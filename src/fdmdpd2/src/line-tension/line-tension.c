/*
 * line-tension.c - Calculate line tension of curved membrane
 * Copyright (C) 2009 Martin Hoemberg <mhoembe@gwdg.de>
 *
 * Syntax: line-tension [config.txt] {output1.dat} {output2.dat} ...
 */

#include <stdio.h>
#include <stdlib.h>
#include <fenv.h>
#include <math.h>
#include <string.h>
#include <stdarg.h> 
#include <limits.h>
#include <errno.h>
#include <fftw3.h>

/* number of different bead species, currently only "A" and "B" */
#define TYPE_MAX 2

typedef unsigned long long PASSPORT;
typedef double VEC[3];
typedef double VEC2[6];
typedef double DENSITY[2 * TYPE_MAX];	/* 2nd order, 3rd order */
typedef char RESNAME[9];

/* some useful macros */
#define SQR(X)		((X) * (X))
#define CUBE(X)		((X) * (X) * (X))
#define MAX(X,Y)	(((X) > (Y)) ? (X) : (Y))
#define MIN(X,Y)	(((X) < (Y)) ? (X) : (Y))
#define ARRAY_SIZE(arr) (sizeof(arr) / sizeof((arr)[0]))
					      
/* branch prediction macros, uses a GCC extension. */
#define	likely(exp)	__builtin_expect(!!(exp), 1)
#define	unlikely(exp)	__builtin_expect(!!(exp), 0)

/* offsets of neighboring hexagons in a unit cell */
static const int offy[8] = {1, 2, 1, 2, 1, 0, 1, 0};
static const int offz[8] = {0, 1, 0, 0, 0, 1, 0, 0};

/* index differences of neighbor hexagons, needed for triangulation */
static const int nb_py[3][8] = {
	{ 1, 1, 1, 1, 0, 0, 0, 0},
	{ 2, 2, 2, 2, 1, 1, 1, 1},
	{ 1, 2, 2, 1, 1, 0, 0, 1}
};
static const int nb_pz[3][8] = {
	{ 0, 0, 0,-1, 1, 0, 0, 0},
	{ 1, 0, 0, 0, 0, 0, 0,-1},
	{ 1, 1, 1, 0, 1, 1, 1, 0}
};

/* struct needed for the hexagonal triangulation */
struct hexgrid
{
	double *x;			/* local height */
	double *xn;			/* counter */
	double a;			/* lattice constant in y direction */
	double b1;			/* b * cos(theta) */
	double b2;			/* b * sin(theta) */
	double w0;			/* .5 * a + b1 */
	double cy;			/* unit cell, 2 * a + 2 * b1 */
	double cz;			/* unit cell, 2 * b2 */
	int my;				/* # hexagons in x, even number! */
	int mz;				/* # hegagons in z direction */
	int maxbins;			/* max. number of bins */
};

enum leaflet
{
	UPPER_LEAFLET,
	LOWER_LEAFLET,
	BOTH_LEAFLETS
};

/* struct containing system information */
struct beads
{
	VEC2 *xv;			/* coordinates and velocities */
	VEC *nx;			/* pbc counters */
	double *lambdas;		/* order parameter of each molecule */
	VEC *normal;			/* average normal of each molecule */
	VEC *com;			/* center-of-mass */
	int *type;			/* A or B */
	enum leaflet *leaflet;		/* leaflet of each molecule */
	VEC l;				/* box size */
	double time;			/* current system time */
	int n;				/* number of molecules */
	int N;				/* number of beads per molecule */
};

/* main struct containing everything */
struct data
{
	char **filename;		/* array with filenames */
	double *fft_in;			/* FFT input buffer */
	double *S;			/* FFT power spectrum */
	double *S2;			/* FFT power spectrum */
	fftw_complex *fft_out;		/* FFT output buffer */
	double dy0;			/* rough bin width */
	double q1;			/* average 2 * PI / L_z */
	double Lz;			/* average L_z */
	fftw_plan fftplan;		/* FFT plan */
	int fftn;			/* FFT size, number of z bins */
	int lines;			/* number of lines per system */
	int files;			/* # data files */
	int dump;			/* dump flatted profile? */
	struct hexgrid hg;		/* triangulation struct */
};

/* ------------------------------------------------------------------------- */

static int error(int errnum, const char *fmt, ...)
					__attribute__((format(printf, 2, 3)));
static int get_data_line(FILE *restrict FH, const char *restrict fmt, ...)
					__attribute__((format(scanf, 2, 3)));

/*
 * print debugging message to stderr
 */
#ifdef DEBUG
static void debug(const char *fmt, ...) __attribute__((format(printf, 1, 2)));
static void debug(const char *fmt, ...)
{
	va_list args;

	va_start(args, fmt);
	fprintf(stderr, "DEBUG: ");
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

/*
 * This function reads in a single line from the given file handle and
 * parses the information.
 */
static int get_data_line(FILE *FH, const char *restrict fmt, ...)
{
	va_list ap;
	int ret;
	char line[LINE_MAX];

	va_start(ap, fmt);
	fgets(line, LINE_MAX, FH);
	ret = vsscanf(line, fmt, ap);
	va_end(ap);

	return ret;
}

/* ------------------------------------------------------------------------- */

static int init_hexgrid(struct hexgrid *restrict hg, double a)
{
	memset(hg, 0, sizeof(*hg));

	hg->a = a;
	debug("hexgrid: a=%lg", a);

	return 0;
}

static void free_hexgrid(struct hexgrid *restrict hg)
{
	if (hg->x)
		free(hg->x);
	if (hg->xn)
		free(hg->xn);
}

static void normalize(VEC n)
{
	int k;
	double norm = sqrt(SQR(n[0]) + SQR(n[1]) + SQR(n[2]));
	for (k = 0; k < 3; k++) n[k] /= norm;
}

/*
 * return the correct index of a hexagon from the given lateral coordinates
 */ 
static int get_hexagon(const struct hexgrid *restrict hg, const VEC2 xv)
{
	int ny = (int)floor(xv[1] / hg->cy);
	int nz = (int)floor(xv[2] / hg->cz);
	double dy = xv[1] - ny * hg->cy;
	double dz = xv[2] - nz * hg->cz;
	int idx = 0;

	ny *= 2;

	idx |= (dy < hg->cy - dy) ? 0x4 : 0x0;
	dy = (dy < hg->cy - dy) ? dy : hg->cy - dy;
	idx |= (dz < hg->cz - dz) ? 0x2 : 0x0;
	dz = (dz < hg->cz - dz) ? dz : hg->cz - dz;
	idx |= (dz < (hg->w0 - dy) * hg->b2 / hg->b1) ? 0x1 : 0x0;

	ny = (ny + offy[idx]) % hg->my;
	nz = (nz + offz[idx]) % hg->mz;

	return ny * hg->mz + nz;
}

/*
 * return the local normal vector from the given lateral coordinates
 */
static int get_normal(const struct hexgrid *restrict hg, const VEC x, VEC n)
{
	int ny = (int)floor(x[1] / hg->cy);
	int nz = (int)floor(x[2] / hg->cz);
	double dy = x[1] - ny * hg->cy;
	double dz = x[2] - nz * hg->cz;
	VEC a, b;
	int idx = 0, k;
	int p_idx[3];

	ny *= 2;

	idx |= (dy < hg->cy - dy) ? 0x4 : 0x0;
	idx |= (dz < hg->cz / hg->cy * dy) ? 0x2 : 0x0;
	idx |= (dz < hg->cz * (1. - dy / hg->cy)) ? 0x1 : 0x0;

	/* all three points for this triangle */
	for (k = 0; k < 3; k++) { 
		int ty = ny + nb_py[k][idx] + hg->my;
		int tz = nz + nb_pz[k][idx] + hg->mz;
		ty %= hg->my;
		tz %= (ty % 2 == 0) ? hg->mz : hg->mz - 1;
		p_idx[k] = ty * hg->mz + tz;
	}

	a[0] = hg->x[p_idx[1]] - hg->x[p_idx[0]];
	a[1] = (nb_py[1][idx] - nb_py[0][idx]) * (hg->a + hg->b1);
	a[2] = (nb_pz[1][idx] - nb_pz[0][idx]) * 2. * hg->b2;
	
	b[0] = hg->x[p_idx[2]] - hg->x[p_idx[0]];
	b[1] = (nb_py[2][idx] - nb_py[0][idx]) * (hg->a + hg->b1);
	b[2] = (nb_pz[2][idx] - nb_pz[0][idx]) * 2. * hg->b2;

	n[0] = a[1] * b[2] - a[2] * b[1];
	n[1] = a[2] * b[0] - a[0] * b[2];
	n[2] = a[0] * b[1] - a[1] * b[0];

	normalize(n);
	return 0;
}

/*
 * setup hexagonal grid and setup triangulation
 */
static int refresh_hexgrid(struct hexgrid *restrict hg, const VEC l)
{
	int my = 2 * (int)(l[1] / hg->a / 3.732050808 + .5);
	int mz, mz_best;
	double best_ratio = 1000.;
	double ty, tz, b, ct;
	int start = MAX(1, hg->mz - 10);
	int end = (hg->mz == 0) ? 1000 : hg->mz + 10;

	ty = l[1] / my - hg->a;
	for (mz = start; mz < end; mz++) {
		tz = .5 * l[2] / mz;
		b = sqrt(ty * ty + tz * tz);
		ct = ty / b;
	
		if (fabs(ct - sqrt(.5)) < best_ratio) {
			best_ratio = fabs(ct - sqrt(.5));
			mz_best = mz;
		}
	}

	tz = .5 * l[2] / mz_best;
	b = sqrt(ty * ty + tz * tz);
	ct = ty / b;

	hg->my = my;
	hg->mz = mz_best;
	hg->b1 = b * ct;
	hg->b2 = b * sqrt(1. - SQR(ct));
	hg->w0 = .5 * hg->a + hg->b1;
	hg->cy = 2. * (hg->a + hg->b1);
	hg->cz = 2. * hg->b2;

	if (hg->maxbins < hg->my * hg->mz) {
		hg->maxbins = hg->my * hg->mz;
		
		hg->x = calloc(hg->maxbins, sizeof(*hg->x));
		if (hg->x == NULL) return novm("hg->x");
		hg->xn = calloc(hg->maxbins, sizeof(*hg->xn));
		if (hg->xn == NULL) return novm("hg->xn");
	}

	debug("Triangulation: my=%d mz=%d a=%lg b=%lg, theta=%lg",
			       	my, mz_best, hg->a, b, acos(ct) * 180. / M_PI);
	return 0;
}

/*
 * calculate surface by mapping the beads to the next grid points
 */
static int set_surface(const struct hexgrid *restrict hg,
						const struct beads *restrict b)
{
	int i;

	memset(hg->x, 0, hg->maxbins * sizeof(*hg->x));
	memset(hg->xn, 0, hg->maxbins * sizeof(*hg->xn));

	for (i = 0; i < b->n * b->N; i++) {
		if (b->type[i] == 0) {
			int n = get_hexagon(hg, b->xv[i]);
			hg->xn[n]++;
			hg->x[n] += b->xv[i][0];
		}
	}

	for (i = 0; i < hg->my * hg->mz; i++)
		hg->x[i] /= hg->xn[i];
	return 0;
}

static int print_surface(FILE *FH, const struct hexgrid *restrict hg)
{
	int i, j;

	for (i = 0; i< hg->my / 2; i++) {
		for (j = 0; j < hg->mz; j++) {
			fprintf(FH, "%lg %lg %lg\n",
				hg->x[2 * i * hg->mz + j],
				i * 2. * (hg->a + hg->b1),
				j * 2. * hg->b2);
		}
		fprintf(FH, "\n");

		for (j = 0; j < hg->mz; j++) {
			fprintf(FH, "%lg %lg %lg\n",
				hg->x[(2 * i + 1) * hg->mz + j],
				(i * 2. + 1) * (hg->a + hg->b1),
				(j * 2. + 1) * hg->b2);
		}
		fprintf(FH, "\n");
	}

	return 0;
}

/* ------------------------------------------------------------------------- */

/*
 * Read in the global configuration file and perform a lot of initial setup
 */
static int parse_configfile(const char *fn, struct data *restrict data)
{
	FILE *FH = fopen(fn, "r");
	double a;
	int ret;

	if (FH == NULL)
		return error(ENOENT, "Couldn't open file %s for reading", fn);

	/* Triangulation with hexagonal lattice */
	/* only lattice constant in y direction is expected */
	if (get_data_line(FH, "%lg", &a) != 1) goto error;
	ret = init_hexgrid(&data->hg, a);
	if (ret) return ret;

	/* dump flatted profiles? */
	if (get_data_line(FH, "%d", &data->dump) != 1) goto error;

	/* read rough estimate of bin widths */
	if (get_data_line(FH, "%lg", &data->dy0) != 1) goto error;
	
	/* FFT size, number of z bins */
	if (get_data_line(FH, "%d", &data->fftn) != 1) goto error;
	data->fft_in = calloc(data->fftn, sizeof(*data->fft_in));
	if (data->fft_in == NULL) return novm("data->fft_in");
	data->fft_out = calloc(data->fftn / 2 + 1, sizeof(*data->fft_out));
	if (data->fft_out == NULL) return novm("data->fft_out");
	data->fftplan = fftw_plan_dft_r2c_1d(data->fftn, data->fft_in,
			data->fft_out, FFTW_MEASURE);
	data->S = calloc(data->fftn / 2 + 1, sizeof(*data->S));
	if (data->S == NULL) return novm("data->S");
	data->S2 = calloc(data->fftn / 2 + 1, sizeof(*data->S2));
	if (data->S2 == NULL) return novm("data->S2");

	/* Number of lines per system, can bei either 2 or 4 */
	if (get_data_line(FH, "%d", &data->lines) != 1) goto error;
	if (data->lines != 2 && data->lines != 4)
		return error(EINVAL, "Only lines=2 or lines=4 allowed!");

	fclose(FH);
	return 0;
error:
	fclose(FH);
	return error(EINVAL, "Error while reading from '%s'", fn);
}

static int load_config(int argc, char **argv, struct data *restrict data)
{
	int ret;

	if (argc < 3) {
		fprintf(stderr, "Syntax: %s [config.txt] {file}...\n\n",
							       argv[0]);
		return -EINVAL;
	}
	data->files = argc - 2;
	data->filename = argv + 2;
	
	debug("conf='%s' files=%d", argv[1], data->files);
	
	ret = parse_configfile(argv[1], data);
	if (ret) return ret;

	return 0;
}

/* ------------------------------------------------------------------------- */

static int read_beads_file(const char *fn, struct beads *restrict b)
{
	FILE *FH = fopen(fn, "r");
	int blocks, i, d;
	char buf[LINE_MAX];

	if (FH == NULL)
		return error(ENOENT, "Couldn't open '%s' for reading!", fn);

	get_data_line(FH, "# L=%lg %lg %lg t=%lg blocks=%d",
			&b->l[0], &b->l[1], &b->l[2], &b->time, &blocks);
	debug("L=%lg %lg %lg time=%lg", b->l[0], b->l[1], b->l[2], b->time);

	if (blocks != 1)
		return error(EINVAL, "Only blocks=1 supported (is: %d)",
		     						blocks);

	fgets(buf, sizeof(buf), FH);
	fgets(buf, sizeof(buf), FH);

	get_data_line(FH, "# n=%d N=%d name=%s", &b->n, &b->N, buf);
	debug("Found n=%d N=%d molecules called '%s'", b->n, b->N, buf);

	b->xv = calloc(b->n * b->N, sizeof(*b->xv));
	if (b->xv == NULL) return novm("b->xv");
	
	b->nx = calloc(b->n * b->N, sizeof(*b->nx));
	if (b->nx == NULL) return novm("b->nx");
	
	b->type = calloc(b->n * b->N, sizeof(*b->type));
	if (b->type == NULL) return novm("b->type");

	b->lambdas = calloc(b->n, sizeof(*b->lambdas));
	if (b->lambdas == NULL) return novm("b->lambdas");

	b->normal = calloc(b->n, sizeof(*b->normal));
	if (b->normal == NULL) return novm("b->normal");

	b->com = calloc(b->n, sizeof(*b->com));
	if (b->com == NULL) return novm("b->com");
	
	b->leaflet = calloc(b->n, sizeof(*b->leaflet));
	if (b->leaflet == NULL) return novm("b->leaflet");

	for (i = 0; i < b->n * b->N; i++) {
		VEC x;
		int ret = get_data_line(FH, "%lg %lg %lg %lg %lg %lg %d",
				&x[0], &x[1], &x[2],
				&b->xv[i][3], &b->xv[i][4], &b->xv[i][5],
				&b->type[i]);
		if (ret != 7)
			return error(EINVAL, "Parsing failed. Expected 7. "
						       "Found %d.", ret);

		for (d = 0; d < 3; d++) {
			double n = 0.;
			while (x[d] >= b->l[d]) {
				x[d] -= b->l[d];
				n += 1.;
			}
			while (x[d] < 0.) {
				x[d] += b->l[d];
				n -= 1.;
			}
			b->xv[i][d] = x[d];
			b->nx[i][d] = n;
		}
	}

	fclose(FH);
	return 0;
}

static void free_beads(struct beads *restrict b)
{
	free(b->xv);
	free(b->nx);
	free(b->type);
	free(b->lambdas);
	free(b->normal);
	free(b->com);
	free(b->leaflet);
}

/* ------------------------------------------------------------------------- */

/* calculate S_{xx} order parameter */
static double get_sxx(VEC2 x[const restrict], const int i, const int j,
	       							const VEC n)
{
	int d;
	double sp = 0., norm2 = 0.;
	
	/* the local normal vector must be normalized. */
	for (d = 0; d < 3; d++) {
		double tmp = x[j][d] - x[i][d];
		norm2 += SQR(tmp);
		sp += tmp * n[d];
	}
	return .5 * (3. * SQR(sp) / norm2 - 1.);
}

static int calc_order_parameters(struct hexgrid *restrict hg,
						struct beads *restrict b)
{
	int i, j, d;
	VEC n, sum_n, com;
	int upper = 0, lower = 0;

	for (i = 0; i < b->n; i++) {
		int count = 0;
		double s = 0.;

		for (d = 0; d < 3; d++) {
			sum_n[d] = 0.;
			com[d] = 0.;
		}

		for (j = i * b->N; j < (i + 1) * b->N; j++) {
			for (d = 0; d < 3; d++)
				com[d] += b->xv[j][d] + b->nx[j][d] * b->l[d];
				
			/* only A-A bonds are taken into account */
			if (b->type[j] != 0 || b->type[j + 1] != 0)
				continue;
			
			get_normal(hg, b->xv[j], n);
			s += get_sxx(b->xv, j, j + 1, n);
			count++;

			for (d = 0; d < 3; d++)
				sum_n[d] += n[d];
		}

		b->lambdas[i] = s / count;
		normalize(sum_n);

		for (d = 0; d < 3; d++) {
			b->normal[i][d] = sum_n[d];
			b->com[i][d] = com[d] / b->N;
		}

		/* Count how many lipids are on each leaflet */
		double cc = 0.;
		count = 0;
		for (j = i * b->N; j < (i + 1) * b->N; j++) {
			if (b->type[j] == 1) {
				cc += b->xv[j][0];
				count++;
			}
		}
		cc /= count;

		for (d = 0; d < 3; d++) {
			com[d] = b->com[i][d];
			while (com[d] >= b->l[d]) com[d] -= b->l[d];
			while (com[d] < 0.) com[d] += b->l[d];
		}

		if (cc > hg->x[get_hexagon(hg, com)]) {
			b->leaflet[i] = UPPER_LEAFLET;
			upper++;
		} else {
			b->leaflet[i] = LOWER_LEAFLET;
			lower++;
		}
	}

	printf("Upper=%d Lower=%d Ratio=%lg\n", upper, lower, (double)upper / lower);

	FILE *surface = fopen("surface.dat", "w");
	print_surface(surface, hg);
	fclose(surface);

	return 0; 
}

static int rotate_molecules(struct hexgrid *restrict hg, struct beads *restrict b)
{
	int i, d;
	VEC nges = {0., 0., 0.};
	const double threshold = 0.60;

	/* calculate collective rotation angle */
	for (i = 0; i < b->n; i++) {
		if (b->lambdas[i] < threshold)
		       continue;
		for (d = 0; d < 3; d++) 
			nges[d] += b->normal[i][d];
	}
	
	normalize(nges);
	printf("Total normal: (%lg|%lg|%lg)\n", nges[0], nges[1], nges[2]);
	double theta = atan2(nges[1], nges[0]);
	double st, ct;
	sincos(theta, &st, &ct);
	
	for (i = 0; i < b->n; i++) {
		double x = b->com[i][0];
		double y = b->com[i][1];
		b->com[i][0] = x * ct + y * st;
		b->com[i][1] = y * ct - x * st;
	}

	return 0;
}


#ifdef BLBUB

	int i, j, d;
	VEC com;

	
	for (i = 0; i < b->n; i++) {
		double st, ct;
		sincos(theta, &st, &ct);

		for (d = 0; d < 3; d++) {
			com[d] = b->com[i][d];
			while (com[d] >= b->l[d]) com[d] -= b->l[d];
			while (com[d] < 0.) com[d] += b->l[d];
		}
		
		int bin = get_hexagon(hg, com);
		double dx = hg->x[bin];

		for (j = i * b->N; j < (i + 1) * b->N; j++) {
			double x = b->xv[j][0] + b->nx[j][0] * b->l[0] - b->com[i][0];
			double y = b->xv[j][1] + b->nx[j][1] * b->l[1] - b->com[i][1];
			b->xv[j][0] = x * ct + y * st + b->com[i][0] - dx;
			b->xv[j][1] = y * ct - x * st + b->com[i][1];
			b->xv[j][0] += .5 * b->l[0];
		
			for (d = 0; d < 3; d++) {
				while (b->xv[j][d] >= b->l[d]) {
					b->xv[j][d] -= b->l[d];
					b->nx[j][d] += 1.;
				}
				while (b->xv[j][d] < 0.) {
					b->xv[j][d] += b->l[d];
					b->nx[j][d] -= 1.;
				}
			}
		}
	}

	return 0;
}
#endif

static int dump_molecules(const char *fn, struct beads *restrict b)
{
	int i;
	FILE *FH = fopen(fn, "w");
	static char *type[2] = {"H", "C"};

	if (FH == NULL)
		return error(EIO, "Couldn't open '%s' for writing", fn);

	fprintf(FH, "%d\n\n", b->n * b->N);

	for (i = 0; i < b->n * b->N; i++)
		fprintf(FH, "%s %lg %lg %lg\n", type[b->type[i]], b->xv[i][0],
						b->xv[i][1], b->xv[i][2]);
	fprintf(FH, "\n");

	fclose(FH);

	return 0;
}

/* ------------------------------------------------------------------------- */

static double integral(const int lower, const int upper, const double m0,
       const double m1, const double g, double profile[], const double dy)
{
	int i;
	double sum = 0.;
	
	for (i = lower; i < upper; i++) {
		if ((i + 1) * dy < g) {
			sum += (m0 - profile[i]) * dy;
		} else if (i * dy < g && (i + 1) * dy > g) {
			double f = (g - i * dy) / dy;
			sum += f * (m0 - profile[i]) * dy;
			sum += (1. - f) * (m1 - profile[i]) * dy;
		} else if (i * dy > g) {
			sum += (m1 - profile[i]) * dy;
		}
	}

	return sum;
}

/* Gibbs dividing surface */
static double gds(int lower, int upper, double m0, double m1,
					double dy, double profile[])
{
	double sum, g;
	double lb = lower * dy;
	double ub = upper * dy;
	const double eps = .0001;
	int step = 0;

	do {
		step++;
		g = .5 * (lb + ub);
		
		sum = integral(lower, upper, m0, m1, g, profile, dy);
		debug("step=%d %lg<%lg<%lg sum=%lg", step, lb, g, ub, sum);

		if (m0 > m1) {
			if (sum > eps) ub = g;
			if (sum < -1. * eps) lb = g;
		} else {
			if (sum > eps) lb = g;
			if (sum < -1. * eps) ub = g;
		}

	} while (fabs(sum) > eps && step < 20);
	return g;
}

static int makehisto(struct beads *restrict b, double zl, double zu,
				int bins, double histo[],
			       	enum leaflet leaflet)
{
	double histoc[bins];
	double dy = b->l[2] / bins; /* XXX */
	int i;
	
	memset(histo, 0, bins * sizeof(*histo));
	memset(histoc, 0, bins * sizeof(*histoc));

	for (i = 0; i < b->n; i++) {
		double com_z = b->com[i][1]; /* XXX */

		while (com_z >= b->l[1]) com_z -= b->l[1]; /* XXX */
		while (com_z < 0.) com_z += b->l[1]; /* XXX */

		if (com_z < zl || com_z >= zu)
			continue;
		if (leaflet != BOTH_LEAFLETS && leaflet != b->leaflet[i])
			continue;
		
		double com = b->com[i][2]; /* XXX */
		int bin;

		while (com >= b->l[2]) com -= b->l[2]; /* XXX */
		while (com < 0.) com += b->l[2]; /* XXX */
		bin = com / dy;

		if (bin < 0 || bin >= bins)
			return error(EINVAL, "This should never happen, "
			"bin=%d/%d com=%lg (%lg)", bin, bins, com, b->com[i][1]);

		histo[bin] += b->lambdas[i];
		histoc[bin]++;
	}

	for (i = 0; i < bins; i++) {
		if (histoc[i] > 0) {
			histo[i] /= histoc[i];
		}
	}

	return 0;
}

/*
 * calculate the two bulk average quantities
 */
static int calc_bulk_values(struct beads *restrict b, const int max1,
		       const int max2, const double dy, double r_lambda[])
{	
	int i, d;
	double lower[3], upper[3];
	double c[3], lambda[3];

	lower[0] = 0.;
	upper[0] = (MIN(max1, max2) - 3) * dy;
	lower[1] = (MIN(max1, max2) + 3) * dy;
	upper[1] = (MAX(max1, max2) - 3) * dy;
	lower[2] = (MAX(max1, max2) + 3) * dy;
	upper[2] = b->l[2]; /* XXX */
	for (d = 0; d < ARRAY_SIZE(lower); d++) {
		c[d] = 0.;
		lambda[d] = 0.;
	}

	for (i = 0; i < b->n; i++) {
		double com = b->com[i][2]; /* XXX */

		while (com >= b->l[2]) com -= b->l[2]; /* XXX */
		while (com < 0.) com += b->l[2]; /* XXX */

		for (d = 0; d < ARRAY_SIZE(lower); d++) {
			if (com >= lower[d] && com < upper[d]) {
				lambda[d] += b->lambdas[i];
				c[d]++;
			}
		}
	}

	if (c[0] == 0. || c[1] == 0. || c[2] == 0.) goto error;
	c[0] += c[2];
	lambda[0] += lambda[2];
	
	for (d = 0; d < 2; d++) {
		r_lambda[d] = lambda[d] / c[d];
		printf("lambda_%d = %lg\n", d, r_lambda[d]);
		if (r_lambda[d] < 0. || r_lambda[d] > 1.) goto error;
	}
	if (r_lambda[0] < 0.6 && r_lambda[1] < 0.6) goto error;

	return 0;
error:
	printf("Something went wrong in %s()!\n", __func__);
	printf("* max1=%d max2=%d dy=%lg\n", max1, max2, dy);
	printf("* c=(%lg|%lg|%lg) lambda=(%lg|%lg|%lg)\n",
		       c[0], c[1], c[2], lambda[0], lambda[1], lambda[2]);
	printf("* r_lambda=");
	for (d = 0; d < 2; d++)
		printf("%lg ", r_lambda[d]);
	printf("\n\n");
	return -EINVAL;
}

static int find_gds_from_histo(struct data *restrict data,
		struct beads *restrict b, double lambda[const restrict],
		double prof[restrict], enum leaflet leaflet, const int bins,
		const int boundary[static 4])
{
	int i;
	double histo[bins];
	double dy = b->l[2] / bins; /* XXX */
	double dz = b->l[1] / data->fftn; /* XXX */

	for (i = 0; i < data->fftn; i++) {
		makehisto(b, i * dz, (i + 1) * dz, bins, histo, leaflet);

		prof[i] = gds(boundary[0], boundary[1], lambda[0], lambda[1], dy, histo);
		prof[i + data->fftn] = gds(boundary[2], boundary[3], lambda[1], lambda[0], dy, histo);

		debug("%d\t%lg\t%lg", i, prof[i], prof[i + data->fftn]);
	}

	return 0;
}

static int find_profile(struct data *restrict data, struct beads *restrict b,
		double prof[])
{
	int i, ret;
	int bins = (int)(b->l[2] / data->dy0 + .5); /* XXX */
	double histo[bins], dy;
	int max1 = 0, max2 = 0;

	dy = b->l[2] / bins; /* XXX */

	do {
		double sd1 = 0., sd2 = 0.;
		max1 = 0;
		max2 = 0;

		makehisto(b, 0., b->l[1], bins, histo, BOTH_LEAFLETS); /* XXX */
	
		for (i = 0; i < bins; i++) {
			printf("%i %lg %lg\n", i , i * dy, histo[i]);
		}
	
		const int md = 5; /* minimal distance between minima */
		for (i = 1; i < bins - 1; i++) {
			double diff = (histo[i + 1] - histo[i - 1]) / (2. * dy); 
	
			if (fabs(diff) > sd1) {
				sd1 = fabs(diff);
				max1 = i;
			}
		}
		for (i = 1; i < bins - 1; i++) {
			double diff = (histo[i + 1] - histo[i - 1]) / (2. * dy);
			
			if (fabs(diff) > sd2 && abs(i - max1) > md) {
				sd2 = fabs(diff);
				max2 = i;
			}
		}
		debug("MAX1: i=%d sd=%lg MAX2: i=%d sd=%lg", max1, sd1, max2, sd2);

		if (max1 < 5 || max2 < 5) {
			for (i = 0; i < b->n; i++) {
				b->com[i][2] += 5. * dy; /* XXX */
			}
		}
	} while (max1 < 5 || max2 < 5);

	double lambda[2];
	ret = calc_bulk_values(b, max1, max2, dy, lambda);
	if (ret) {
		for (i = 0; i < bins; i++)
			printf("%d %lg %lg\n", i, i * dy, histo[i]);
		return ret;
	}

	int boundary[4] = {
		MAX(0, max1 - (int)(3. / dy)),
		MIN(bins, max1 + (int)(3. / dy)),
		MAX(0, max2 - (int)(3. / dy)),
		MIN(bins, max2 + (int)(3. / dy))
	};

	memset(prof, 0, sizeof(prof));
	if (data->lines == 2) {
		find_gds_from_histo(data, b, lambda, prof, BOTH_LEAFLETS,
							bins, boundary);
	} else if (data->lines == 4) {
		find_gds_from_histo(data, b, lambda, prof, UPPER_LEAFLET,
							bins, boundary);
		find_gds_from_histo(data, b, lambda, prof + 2 * data->fftn,
					       LOWER_LEAFLET, bins, boundary);
	}

	data->q1 += 2. * M_PI / b->l[1]; /* XXX */
	data->Lz += b->l[1]; /* XXX */

	return 0;
}

/* ------------------------------------------------------------------------- */

static double sinc(double x) __attribute__((const));
static double sinc(const double x)
{
	if (x == 0.)
		return 1.;
	else
		return sin(x * M_PI) / (x * M_PI);
}

/*
 * Calculate FFT of the profile. Either 2 or 4 lines can be used
 */
static int fft_prof(struct data *restrict data, double prof[const restrict])
{
	int i, o;

	for (o = 0; o < data->lines * data->fftn; o += data->fftn) {
		double m = 0.;
		for (i = 0; i < data->fftn; i++)
			m += prof[i + o];
		m /= data->fftn;
		for (i = 0; i < data->fftn; i++)
			data->fft_in[i] = prof[i + o] - m;
		
		fftw_execute(data->fftplan);

		for (i = 0; i < data->fftn / 2 + 1; i++) {
			double re = data->fft_out[i][0];
			double im = data->fft_out[i][1];
			double S = (re * re + im * im) / SQR(data->fftn);
			double damp = SQR(sinc((double)i / data->fftn));
			data->S[i] += S; // / damp;
			data->S2[i] += SQR(S); // / damp);
		}
	}

	return 0;
}

static int fft_print_spectrum(const char *fn, struct data *restrict data)
{
	int i;
	FILE *FH = fopen(fn, "w");

	for (i = 0; i < data->fftn / 2 + 1; i++) {
		data->S[i] /= data->lines * data->files;
		data->S2[i] /= data->lines * data->files;
	}
		
	data->q1 /= data->files;
	data->Lz /= data->files;

	fprintf(FH, "# q S dS n\n");
	fprintf(FH, "# Lz=%lg\n", data->Lz);

	for (i = 0; i < data->fftn / 2 + 1; i++) {
		double q = data->q1 * i;
		double dS = sqrt(data->S2[i] - SQR(data->S[i])) /
		       			sqrt(data->lines * data->files - 1);
		fprintf(FH, "%lg %lg %lg %d\n", q, data->S[i], dS, i);
	}

	fclose(FH);
	return 0;
}

/* ------------------------------------------------------------------------- */

static int main_loop(struct data *restrict data)
{
	int i, ret;
	struct beads b;

	for (i = 0; i < data->files; i++) {
		const char *fn = data->filename[i];
		double prof[data->lines * data->fftn];
		fprintf(stderr, "Opening file %d/%d '%s'...\n", i + 1, data->files, fn);

		FILE *FH = fopen(fn, "r");
		if (FH == NULL) {
			printf("Couldn't open '%s' for reading!\n", fn);
			continue;
		}
		fclose(FH);

		ret = read_beads_file(fn, &b);
		if (ret) return ret;

		ret = refresh_hexgrid(&data->hg, b.l);
		if (ret) return ret;
		
		ret = set_surface(&data->hg, &b);
		if (ret) return ret;

		ret = calc_order_parameters(&data->hg, &b);
		if (ret) return ret;

//		ret = rotate_molecules(&data->hg, &b);
//		if (ret) return ret;
		
		if (data->dump) {
			char dfn[LINE_MAX];
			snprintf(dfn, sizeof(dfn), "%s.dump.xyz", fn);

			ret = dump_molecules(dfn, &b);
			if (ret) return ret;
		}

		ret = find_profile(data, &b, prof);
		if (ret) return ret;

		/* the format in prof is y0, y1, y2, y3 (left side),
		 * y0, y1, y2, y3 (right side) */
		ret = fft_prof(data, prof);
		if (ret) return ret;
		
		free_beads(&b);
	}

	ret = fft_print_spectrum("spectrum.dat", data);
	if (ret) return ret;

	return 0;
}


/* ------------------------------------------------------------------------- */

static void free_data(struct data *restrict data)
{
	free_hexgrid(&data->hg);
	free(data->fft_in);
	free(data->fft_out);
	fftw_destroy_plan(data->fftplan);
	free(data->S);
	free(data->S2);
}

int main(int argc, char **argv)
{
	int ret;
	struct data data;

#if defined(DEBUG) && defined(FE_NOMASK_ENV)
	feenableexcept(FE_UNDERFLOW | FE_OVERFLOW | FE_INVALID | FE_DIVBYZERO);
#endif
	ret = load_config(argc, argv, &data);
	if (ret) goto error;

	ret = main_loop(&data);
	if (ret) goto error;

	free_data(&data);
	fftw_cleanup();

	return EXIT_SUCCESS;
error:
	return EXIT_FAILURE;
}

