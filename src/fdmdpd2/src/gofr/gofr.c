/*
 * gofr.c - computation of radial distribution function (RDF)
 * (C) Copyright 2007-2010 Martin Hoemberg <mhoembe@gwdg.de>
 *
 * This program reads an FDMDPD2 configuration file, and calculates g(r) for
 * a specified combination of blocks and types of particles. This program
 * uses OpenMP, so that the number of concurrent threads can influenced by
 * setting the environment variable OMP_NUM_THREADS prior to execution.
 *
 * The data are always read from STDIN.
 *
 * Parameters:
 *    -B1 / -B2: 1st / 2nd block 
 *    -T1 / -T2: type of 1st/2nd particles (in block B1 / B2)
 *     If one of B1/B2/T1/T2 is left unspecified, then all matching
 *     particles are taken into consideration.
 * 
 *    -D: dimensions, e.g. -D=yz
 *     if "xyz" is specified, the 3diml g(r) is calculated.
 *     Default value is: "yz", i.e. the 2diml in-plane g(r)
 *
 *    -W: radial binwidth, default 0.01
 *
 *    -R: cutoff, max "r"; default 20.0
 *
 */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include <assert.h>
#include <argp.h>


struct list
{
	struct list *next;
	int val;
};

struct filter
{
	struct list *type;
	struct list *block;
	int filter_blocks;
	int filter_types;
};

struct cfg
{
	double *x;
	double *g;
	short *t;
	double L[3];
	int n[20];
	int N[20];
	int max;
};

static double dr = 0.01; 		/* radial binwidth */
static double rmax = 20.0;		/* max. r */
static int bins = 0;			/* # bins */
static char names[20][20];		/* names of the blocks */
static int blocks = 0;			/* # blocks */
static int files = 0;			/* # files */
static struct filter filter1 = { 0 };	/* filter for 1st particles */
static struct filter filter2 = { 0 };	/* filter for 2nd particles */
static char *fn = NULL;			/* filename */

/* ------------------------------------------------------------------------- */

/* 
 * find a specific entry in the list; if it is found, this entry is moved to
 * the front, and '1' is returned. If nothing is found, '0' is returned.
 */
static int isinlist(struct list *head, int val)
{
	struct list *l = head;
	while (l != NULL) {
		if (l->val == val) {
			if (l != head) {
				l->val = head->val;
				head->val = val;
			}
			return 1;
		}
		l = l->next;
	}
	return 0;
}

/* ------------------------------------------------------------------------- */

/* is the block 'b' included in this filter? */
static int useblock(struct filter *f, int b)
{
	return (f->filter_blocks) ? isinlist(f->block, b) : 1;
}

/* is the type 't' included in this filter? */
static int usetype(struct filter *f, int t)
{
	return (f->filter_types) ? isinlist(f->type, t) : 1;
}

/* ------------------------------------------------------------------------- */

static struct argp_option argp_option[] = {
	{"dim", 'D', "xyz", 0, "dimensions in which g(r) shall be calculated"},
	{"binwidth", 'W', "value", 0, "radial binwdith, default is 0.01"},
	{"cutoff", 'R', "value", 0, "maximum r, default 20.0"},
	{"block1", 'A', "block1,block2,...", 0, "possible blocks for 1st particles"},
	{"block2", 'B', "block1,block2,...", 0, "possible blocks for 2nd particles"},
	{"type1", 'S', "type1,type2,...", 0, "possible types for 1st particles"},
	{"type2", 'T', "type1,type2,...", 0, "possible types for 2nd particles"},
	{0}
};

static int create_list(char *arg, struct list **l)
{
	struct list *list = NULL;
	const char delim[] = ",;";
	char *tok, s[LINE_MAX];
	int i, cnt = 0;

	strncpy(s, arg, sizeof(s));
	tok = strtok(s, delim);

	while (tok != NULL) {
		cnt++;
		tok = strtok(NULL, delim);
	}

	fprintf(stderr, "'%s' has %d entries.\n", arg, cnt);

	list = calloc(cnt, sizeof(*list));
	if (list == NULL) return ENOMEM;

	for (i = 0; i < cnt; i++)
		list[i].next = (i < cnt - 1) ? list + i + 1 : NULL;

	strncpy(s, arg, sizeof(s));
	tok = strtok(s, delim);

	i = 0;
	while (tok != NULL) {
		list[i++].val = atoi(tok);
		tok = strtok(NULL, delim);
	}
	*l = list;
	return 0;
}

static error_t parse_opt(int key, char *arg, struct argp_state *state)
{
	switch (key) {
	case 'D':
		/* not implemented yet. */
		break;
	case 'W':
		dr = atof(arg);
		fprintf(stderr, "Using binwidth dr=%lg.\n", dr);
		break;
	case 'R':
		rmax = atof(arg);
		fprintf(stderr, "Using rmax=%lg.\n", rmax);
		break;
	case 'A':
		filter1.filter_blocks = 1;
		return create_list(arg, &filter1.block);
	case 'B':
		filter2.filter_blocks = 1;
		return create_list(arg, &filter2.block);
	case 'S':
		filter1.filter_types = 1;
		return create_list(arg, &filter1.type);
	case 'T':
		filter2.filter_types = 1;
		return create_list(arg, &filter2.type);
	case ARGP_KEY_ARG:
		fn = arg;
		break;
	case ARGP_KEY_END:
		if (state->arg_num < 1)
			argp_usage(state);
		break;
	default:
		return ARGP_ERR_UNKNOWN;
	}
	return 0;
}

const char *argp_program_version = "gofr 1.0";
const char *argp_program_bug_address = "<mhoembe@gwdg.de>";
static const char doc[] = "Calculate radial pair correlation function";
static struct argp argp = {argp_option, parse_opt, "[filename|-]", doc};


/* ------------------------------------------------------------------------- */


static int readconfig(FILE *FH, struct cfg *cfg)
{
	char buf[LINE_MAX];
	double tmp;
	int i, b, max;

	cfg->n[0] = 0;
	if (fgets(buf, sizeof(buf), FH) == NULL)
		return 0;

	if (sscanf(buf, "# L=%lg %lg %lg t=%lg blocks=%d",
		       	&cfg->L[0], &cfg->L[1], &cfg->L[2], &tmp, &b) != 5) {
		fprintf(stderr, "Error0 reading '%s'\n", buf);
		return 1;
	}

	if (blocks == 0)
		blocks = b;

	b = 0;
	max = 0;
	do {
		fgets(buf, sizeof(buf), FH);
		if (strstr(buf, "# n=") == buf) {
			if (sscanf(buf, "# n=%d N=%d name=%s", &cfg->n[b], &cfg->N[b], names[b]) != 3) {
				fprintf(stderr, "Error1 reading '%s'\n", buf);
				return 1;
			}
			fprintf(stderr, "%s", buf);

			if (cfg->max < max + cfg->n[b] * cfg->N[b]) {
				cfg->max += cfg->n[b] * cfg->N[b];
				cfg->x = realloc(cfg->x, 3 * cfg->max * sizeof(*cfg->x));
				if (cfg->x == NULL) {
					fprintf(stderr, "novm: cfg->x");
					return 1;
				}
				cfg->t = realloc(cfg->t, cfg->max * sizeof(*cfg->t));
				if (cfg->t == NULL) {
					fprintf(stderr, "novm: cfg->t");
					return 1;
				}
			}

			for (i = 0; i < cfg->N[b] * cfg->n[b]; i++) {
				fgets(buf, sizeof(buf), FH);
				if (sscanf(buf, "%lg %lg %lg %lg %lg %lg %hd",
					       	&cfg->x[3 * (max + i) + 0],
						&cfg->x[3 * (max + i) + 1],
						&cfg->x[3 * (max + i) + 2],
						&tmp, &tmp, &tmp,
						&cfg->t[max + i]) != 7) {
					fprintf(stderr, "Error2 reading '%s'",
									buf);
					return 1;
				}
			}

			max += cfg->n[b] * cfg->N[b];
			b++;
		}
	} while (b < blocks);

	if (cfg->g == NULL) {
		cfg->g = calloc(bins, sizeof(*cfg->g));
		if (cfg->g == NULL) {
			fprintf(stderr, "novm: cfg->g");
			return 1;
		}
	} else {
		memset(cfg->g, 0, bins * sizeof(*cfg->g));
	}

	return 0;
}

static void pair(struct cfg *cfg, int i, int j, double pref)
{
	double dy, dz;

//	dx = cfg->x[3 * i + 0] - cfg->x[3 * j + 0];
	dy = cfg->x[3 * i + 1] - cfg->x[3 * j + 1];
	dz = cfg->x[3 * i + 2] - cfg->x[3 * j + 2];

//	while (dx > .5 * cfg->L[0]) dx -= L[0];
//	while (dx < -.5 * cfg->L[0]) dx += L[0];
	while (dy > .5 * cfg->L[1]) dy -= cfg->L[1];
	while (dy < -.5 * cfg->L[1]) dy += cfg->L[1];
	while (dz > .5 * cfg->L[2]) dz -= cfg->L[2];
	while (dz < -.5 * cfg->L[2]) dz += cfg->L[2];

	double r = hypot(dy, dz);
	if (r < rmax) {
		int bin = (int)floor(r / dr);
		assert(bin >= 0 && bin < bins);
		cfg->g[bin] += pref;
	}
}

static void calc_gofr(struct cfg *cfg)
{
	int b1, b2, i, j, offset1 = 0, offset2 = 0, n1 = 0, n2 = 0;
	
	for (b1 = 0; b1 < blocks; b1++) {
		for (i = 0; i < cfg->n[b1] * cfg->N[b1]; i++) {
			if (useblock(&filter1, b1) && usetype(&filter1, cfg->t[offset1 + i]))
				n1++;
			if (useblock(&filter2, b1) && usetype(&filter2, cfg->t[offset1 + i]))
				n2++;
		}

		if (useblock(&filter1, b1) || useblock(&filter2, b1)) {
			offset2 = offset1;
			for (b2 = b1; b2 < blocks; b2++) {
				if ((useblock(&filter1, b2) && useblock(&filter2, b1)) || (useblock(&filter1, b1)) && useblock(&filter2, b2)) {
					for (i = 0; i < cfg->n[b1] * cfg->N[b1]; i++) {
						int t1 = cfg->t[offset1 + i];
						if (!usetype(&filter1, t1) && !usetype(&filter2, t1)) continue;
						j = (b1 == b2) ? i + 1 : 0;
						for (; j < cfg->n[b2] * cfg->N[b2]; j++) {
							int t2 = cfg->t[offset2 + j];
							double x = 0.;
							if (usetype(&filter1, t1) && usetype(&filter2, t2))
								x += 1.;
							if (usetype(&filter1, t2) && usetype(&filter2, t1))
								x += 1.;
							if (x > 0.)
								pair(cfg, offset1 + i, offset2 + j, x);
						}
					}
				}
				offset2 += cfg->n[b2] * cfg->N[b2];
			}
		}
		offset1 += cfg->n[b1] * cfg->N[b1];
	}

	for (i = 0; i < bins; i++)
		cfg->g[i] *= cfg->L[1] * cfg->L[2] / (n1 * n2);
}

int main(int argc, char **argv)
{
	int i;

	argp_parse(&argp, argc, argv, 0, 0, 0);
	bins = (int)ceil(rmax / dr) + 1;

	double *g = calloc(bins, sizeof(*g));
	if (g == NULL) {
		fprintf(stderr, "novm: g\n");
		return EXIT_FAILURE;
	}

	FILE *FH = (strcmp(fn, "-") == 0) ? stdin : fopen(fn, "r");
	if (FH == NULL) {
		fprintf(stderr, "Couldn't open '%s' for reading", fn);
		return EXIT_FAILURE;
	}

	#pragma omp parallel private(i)
	{
		struct cfg cfg = { 0 };

		do {
			#pragma omp critical
			{
				if (readconfig(FH, &cfg)) abort();
			}

			if (cfg.n[0] > 0) {
				calc_gofr(&cfg);
				#pragma omp atomic
				files++;
		
				for (i = 0; i < bins; i++) {
					# pragma omp atomic
					g[i] += cfg.g[i];
				}
			}
		} while (cfg.n[0] > 0);

		free(cfg.x);
		free(cfg.t);
		free(cfg.g);
	}

	printf("# r g(r)\n");
	for (i = 1; i < bins; i++) {
		double l = i * dr;
//		double r = pow(l*l*l + 1.5*l*l*dr + 1.5*l*dr*dr + .5*dr*dr*dr, 1./3.);
		double r = sqrt(l * l + l * dr + .5 * dr * dr);
		double dv = 2. * M_PI * r * dr;
		printf("%lg %lg\n", r, g[i] / (dv * files)); 
	}
	free(g);
	fclose(FH);
	if (filter1.filter_blocks) free(filter1.block);
	if (filter1.filter_types) free(filter1.type);
	if (filter2.filter_blocks) free(filter2.block);
	if (filter2.filter_types) free(filter2.type);

	return EXIT_SUCCESS;
}

