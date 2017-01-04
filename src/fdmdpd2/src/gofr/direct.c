#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <errno.h>


void cylindrical_ft(int N, double in[restrict], double out[restrict], double h)
{
	int n, i;
	double buf[N];

	for (n = 0; n < N; n++) {
		double k =  2. * M_PI * n / (2. * h * N);
		for (i = 0; i < N; i++) {
			double r = h * i;
			buf[i] = r * j0(k * r) * in[i];
		}
		out[n] = 3./8. * (buf[0] + buf[N - 1]);
		out[n] += 7./6. * (buf[1] + buf[N - 2]);
		out[n] += 23./24. * (buf[2] + buf[N - 3]);
		for (i = 3; i < N - 3; i++) out[n] += buf[i];
		out[n] *= h;
	}
}

#ifdef _TEST_
const double dr = 0.1;

/*
 * gnuplot> plot "a" u 4:5 w lp, exp(-x*x/2./10.)/10.
 */
static void test(void)
{
	const int n = 1000;
	double a[n], b[n];
	int i;

	for (i = 0; i < n; i++)
		a[i] = exp(-10. * (i*dr)*(i*dr)*.5);

	cylindrical_ft(n, a, b, dr);

	for (i = 0; i < n; i++) {
		double x = i * dr;
		printf("%d %lg %lg ", i, x, a[i]);
		double k = 2. * M_PI * i / (2. * n * dr); 
		printf("%lg %lg\n", k, b[i]);
	}

	return 0;
}
#endif

int main(int argc, char **argv)
{
#ifdef _TEST_
	test();
#else
	if (argc < 4) {
		printf("Syntax: %s [filename] [dr] [rho_A]\n", argv[0]);
		return 1;
	}
	
	FILE *FH = fopen(argv[1], "r");
	double *hr = NULL;
	int count = 0, max = 0;

	double dr = atof(argv[2]);
	double rho = atof(argv[3]);

	if (FH == NULL) {
		fprintf(stderr, "Couldn't open '%s' for reading\n", argv[1]);
		return ENOENT;
	}

	while(!feof(FH)) {
		char line[512];
		double x, y;
		
		fgets(line, 512, FH);
		sscanf(line, "%lg %lg", &x, &y);

		if (count >= max) {
			max += 0x10000;
			hr = realloc(hr, max * sizeof(*hr));
			if (hr == NULL) {
				fprintf(stderr, "novm: hr\n");
				return ENOMEM;
			}
		}
		hr[count++] = y - 1.;
	}

	double *ck = calloc(count, sizeof(*ck));
	if (ck == NULL) {
		fprintf(stderr, "novm: ck\n");
		return ENOMEM;
	}
	double *cr = calloc(count, sizeof(*cr));
	if (cr == NULL) {
		fprintf(stderr, "novm: cr\n");
		return ENOMEM;
	}

	int i;
	cylindrical_ft(count, hr, ck, dr);
	for (i = 0; i < count; i++)
		ck[i] /= (1. + 2. * M_PI * rho * ck[i]);

	cylindrical_ft(count, ck, cr, M_PI / (count * dr));
	for (i = 0; i < count; i++)
		printf("%lg %lg\n", i * dr, cr[i]);

	free(hr);
	free(cr);
	free(ck);
	fclose(FH);
#endif
	return 0;
}

