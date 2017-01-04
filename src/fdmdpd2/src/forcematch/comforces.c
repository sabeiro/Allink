/*
 * comforces.c - program to merge fdmdpd's output generated with "mcom_print2"
 * Copyright 2010 Martin Hoemberg
 *
 * The output of this program can be used either in forcematch.c and friction.c
 *
 * Syntax:  ./comforces resume.dat forces.dat.0?-0?
 *
 */

#include "common.h"

typedef double VEC[3];
typedef double VEC2[6];

static int blocks = 0;			/* # blocks */
static int global_n[50];		/* # molecules in each block */
static int global_N[50];		/* # beads/molecule in each block */
static int bound[51];			/* # of beads per block */
static int files = 0;			/* # data files */
static FILE **fhh;			/* file handles to data files */
static int nmols = -1;			/* # molecules (total) */
static VEC2 *com;			/* com coordinates and velocities */
static VEC *re;				/* end-to-end vectors */
static int *comcnt;			/* # contributing beads to com */
static VEC *force;			/* resulting forces on each com */
static VEC *pforce;			/* pairwise forces (huge!) */

/*
 * read system file to determine the number of blocks, the number of molecules,
 * and their sizes.
 */
static int read_blocks(const char *fn)
{
	char buf[LINE_MAX];
	int n, N;
	
	FILE *FH = fopen(fn, "r");
	if (FH == NULL) return error(ENOENT, "block file");

	bound[0] = 0;
	while (fgets(buf, sizeof(buf), FH)) {
		if (strncmp(buf, "# n=", 4) == 0) {
			if (sscanf(buf, "# n=%d N=%d ", &n, &N) != 2)
			       return error(EINVAL, "block header");
			global_n[blocks] = n;
			global_N[blocks] = N;
			fprintf(stderr, "block=%d n=%d N=%d\n", blocks, n, N);
			bound[blocks + 1] = bound[blocks] + n * N;
			blocks++;
		}
	}
	if (blocks == 0)
		return error(EINVAL, "no blocks found");

	fclose(FH);
	return 0;
}

/*
 * open the force data files
 */
static int open_datafiles(int n, char **fn)
{
	int i;
	fhh = malloc(n * sizeof(*fhh));
	files = n;
	for (i = 0; i < n; i++) {
		fhh[i] = fopen(fn[i], "r");
		if (fhh[i] == NULL)
			return error(ENOENT, "file '%s'", fn[i]);
	}
	return 0;
}

/*
 * calculate the id of the molecule from the id of the bead
 */ 
static int id2mol(int id)
{
	int b, found = 0, mol = 0, k = id, bl;
	
	for (b = 0; b < blocks && !found; b++) {
		if (id >= bound[b] && id < bound[b + 1]) {
			found = 1;
			bl = b;
		} else {
			mol += global_n[b];
			k -= global_n[b] * global_N[b];
		}
	}
	assert(found);
	return mol + k / global_N[bl];
}

/*
 * is this bead the first bead in a molecule?
 */
static int is_first_bead(int id)
{
	int b, k = id;
	
	for (b = 0; b < blocks; b++) {
		if (id >= bound[b] && id < bound[b + 1]) {
			return (k % global_N[b] == 0);
		} else {
			k -= global_n[b] * global_N[b];
		}
	}
	fprintf(stderr, "This should never happen.\n");
	return -1;
}

/*
 * is this bead the last bead in a molecule?
 */
static int is_last_bead(int id)
{
	int b, k = id;
	
	for (b = 0; b < blocks; b++) {
		if (id >= bound[b] && id < bound[b + 1]) {
			return (k % global_N[b] == global_N[b] - 1);
		} else {
			k -= global_n[b] * global_N[b];
		}
	}
	fprintf(stderr, "This should never happen2.\n");
	return -1;
}

/*
 * parse the header of a new time step in the data files
 */
static int parse_header(void)
{
	int i;
	char buf[LINE_MAX];

	for (i = 0; i < files; i++) {
		double L[3], t;
		int n;

		if (fgets(buf, sizeof(buf), fhh[i]) == NULL)
			return 1;
		if (sscanf(buf, "# L=%lg %lg %lg N=%d t=%lg", &L[0], &L[1], &L[2], &n, &t) != 5)
			return error(EINVAL, "header of file %d", i);
		if (i == 0) {
			printf("%s", buf);
			fprintf(stderr, "%s", buf);
		}

		if (nmols < 0) {
			nmols = n;

			com = calloc(nmols, sizeof(*com));
			if (com == NULL) return novm("parse_header com");
			comcnt = calloc(nmols, sizeof(*comcnt));
			if (comcnt == NULL) return novm("parse_header comcnt");
			re = calloc(nmols, sizeof(*re));
			if (re == NULL) return novm("parse_header re");

			force = calloc(nmols, sizeof(*force));
			if (force == NULL) return novm("force");
			pforce = calloc(SQR(nmols), sizeof(*pforce));
			if (pforce == NULL) return novm("pforce");

		} else if (nmols != n) {
			return error(EINVAL, "nmols == %d, n = %d", nmols, n);
		}
	}

	return 0;
}

/*
 * calc positions and velocities, as well as the end-to-end vectors
 * of the molecules.
 */
static int calc_xv(void)
{
	int i, id, molid, type;
	char buf[LINE_MAX];
	VEC2 xv;

	memset(com, 0, nmols * sizeof(*com));
	memset(comcnt, 0, nmols * sizeof(*comcnt));
	memset(re, 0, nmols * sizeof(*re));

	for (i = 0; i < files; i++) {
		debug("calc_xv() at file %d.\n", i);

		do {
			if (fgets(buf, sizeof(buf), fhh[i]) == NULL)
				return error(EINVAL, "calc_xv premature eof");
			if (strncmp(buf, "# --START--", 11) == 0)
				break;
			if (sscanf(buf, "%d %lg %lg %lg %lg %lg %lg %d\n",
				&id, &xv[0], &xv[1], &xv[2],
				&xv[3], &xv[4], &xv[5], &type) != 8)
				return error(EINVAL, "calc_xv line '%s'", buf);
			molid = id2mol(id);

			/* center of mass coordinate and velocity */
			if (type == 0) { /* only A beads contribute to com */
				com[molid][0] += xv[0];
				com[molid][1] += xv[1];
				com[molid][2] += xv[2];
				com[molid][3] += xv[3];
				com[molid][4] += xv[4];
				com[molid][5] += xv[5];
				comcnt[molid]++;
			}

			/* end-to-end vector */
			if (is_first_bead(id)) {
				re[molid][0] -= xv[0];
				re[molid][1] -= xv[1];
				re[molid][2] -= xv[2];
			} else if (is_last_bead(id)) {
				re[molid][0] += xv[0];
				re[molid][1] += xv[1];
				re[molid][2] += xv[2];
			}
		} while (1);
	}
	return 0;
}

/*
 * calculate total pairwise forces between coms
 */
static int calc_f(void)
{
	int i, id1, id2, mol1, mol2;
	char buf[LINE_MAX];
	VEC f;

	memset(force, 0, nmols * sizeof(*force));
	memset(pforce, 0, SQR(nmols) * sizeof(*pforce));

	for (i = 0; i < files; i++) {
		debug("calc_f() at file %d.\n", i);

		do {
			if (fgets(buf, sizeof(buf), fhh[i]) == NULL)
				return error(EINVAL, "calc_f premature eof");
			if (strncmp(buf, "# --END--", 9) == 0)
				break;
			if (sscanf(buf, "%d %d %lg %lg %lg\n", &id1, &id2,
						&f[0], &f[1], &f[2]) != 5)
				return error(EINVAL, "calc_f line '%s'", buf);
			mol1 = id2mol(id1);
			mol2 = id2mol(id2);
			assert(mol1 != mol2);

			force[mol1][0] -= f[0];
			force[mol1][1] -= f[1];
			force[mol1][2] -= f[2];
			force[mol2][0] += f[0];
			force[mol2][1] += f[1];
			force[mol2][2] += f[2];

			pforce[mol1 * nmols + mol2][0] += f[0];
			pforce[mol1 * nmols + mol2][1] += f[1];
			pforce[mol1 * nmols + mol2][2] += f[2];
		} while (1);
	}
	return 0;
}

static void print_coms(void)
{
	int i, j;

	for (i = 0; i < nmols; i++) {
		for (j = 0; j < 6; j++)
			com[i][j] /= comcnt[i];
		printf("%lg %lg %lg %lg %lg %lg ",
			com[i][0], com[i][1], com[i][2],
			com[i][3], com[i][4], com[i][5]);
		printf("%lg %lg %lg %lg %lg %lg\n",
			re[i][0], re[i][1], re[i][2],
			force[i][0], force[i][1], force[i][2]);
	}
}

static void print_forces(void)
{
	int i, j;

	for (i = 0; i < nmols; i++) {
		for (j = i + 1; j < nmols; j++) {
			VEC f;

			f[0] = pforce[j * nmols + i][0] - pforce[i * nmols + j][0];
			f[1] = pforce[j * nmols + i][1] - pforce[i * nmols + j][1];
			f[2] = pforce[j * nmols + i][2] - pforce[i * nmols + j][2];

			if (f[0] != 0. || f[1] != 0. || f[2] != 0.) {
				printf("%d %d %lg %lg %lg\n", i, j, f[0], f[1], f[2]);
			}
		}
	}
}

/* 
 * final clean up
 */ 
static void cleanup(void)
{
	int i;

	for (i = 0; i < files; i++)
		fclose(fhh[i]);
	free(fhh);
	free(com);
	free(comcnt);
	free(re);
	free(force);
	free(pforce);
}

int main(int argc, char **argv)
{
	int ret;

	if (argc < 3) {
		fprintf(stderr, "Syntax: %s [resume.dat] [forces.dat.00-00 ...]\n\n",
								argv[0]);
		return EXIT_FAILURE;
	}

	ret = read_blocks(argv[1]);
	if (ret) return EXIT_FAILURE;

	ret = open_datafiles(argc - 2, argv + 2);
	if (ret) return EXIT_FAILURE;
	
	while (parse_header() == 0) {
		ret = calc_xv();
		if (ret) return ret;

		ret = calc_f();
		if (ret) return ret;

		print_coms();
		print_forces();
		printf("# --END--\n");
	}

	cleanup();
	return EXIT_SUCCESS;
}

