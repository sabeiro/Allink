/*
 * blockavg.c -- Block averaging to determine error of correlated data
 * (C) Copyright 2010 by Martin Hoemberg
 *
 * To understand what this program does, please read the following article:
 * H. Flyvbjerg and H.G. Petersen, J.Chem.Phys. 91, 461--466 (1989)
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#define SQR(x) ((x)*(x))

static int rg_transform(double x[], int n)
{
	int i;

	for (i = 0; i < n / 2; i++) {
		x[i] = .5 * (x[2 * i] + x[2 * i + 1]);
	}
	return n / 2;
}

static double estimate(double x[], double m, int n)
{
	int i;
	double c0;

	for (i = 0, c0 = 0.; i < n; i++)
		c0 += SQR(x[i] - m);
	return c0 / (n * (n - 1.));
}

static double mean(double x[], int n)
{
	int i;
	double m;
	
	for (i = 0, m = 0.; i < n; i++) 
		m += x[i];
	return m / n;
}

static int readdata(int argc, char **argv, double **rx, int *rn)
{
	int cnt = 0, max = 0;
	double tmp, *x = NULL;

	if (argc < 2) {
		fprintf(stderr, "Syntax: %s [-|filename]\n\n", argv[0]);
		return -1;
	}
	FILE *FH = (strcmp(argv[1], "-") == 0) ? stdin : fopen(argv[1], "r");
	if (FH == NULL) {
		fprintf(stderr, "Couldn't open '%s' for reading\n", argv[1]);
		return -1;
	}

	while (!feof(FH)) {
		if (fscanf(FH, " %lg ", &tmp) == 1) {
			if (cnt >= max) {
				max += 1048576;
				x = realloc(x, max * sizeof(*x));
				if (x == NULL) {
					fprintf(stderr, "novm for a\n");
					return -1;
				}
			}
			x[cnt++] = tmp;
		}
	}
	if (FH != stdin) fclose(FH);
	fprintf(stderr, "Read %d datapoints\n", cnt);

	(*rx) = x;
	(*rn) = cnt;
	return 0;
}

int main(int argc, char **argv)
{
	double *x = NULL, m;
	int n = 0, t = 0;

	if (readdata(argc, argv, &x, &n))
		return EXIT_FAILURE;

	m = mean(x, n);
	fprintf(stderr, "m=%lg\n", m);

	printf("# transform sigma error n\n");
	do {
		double sigma = sqrt(estimate(x, m, n));
		double err = sigma / sqrt(2. * n - 2.);

		printf("%d %lg %lg %d\n", t++, sigma, err, n);

		n = rg_transform(x, n);
	} while (n >= 2);

	free(x);
	return EXIT_SUCCESS;
}

