#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <fftw3.h>
#include <string.h>
#include <math.h>


int main(int argc, char **argv)
{
	fftw_plan plan;
	double *tmp = NULL, *in, *out;
	int cnt = 0, max = 0;
	int N, i;
	double dx = (argc == 2) ? atof(argv[1]) : 1.;

	fprintf(stderr, "dx=%lg\n", dx);

	while (!feof(stdin)) {
		if (cnt >= max) {
			max += 65536;
			tmp = realloc(tmp, max * sizeof(*tmp));
		}
		if (fscanf(stdin, "%lf ", tmp + cnt) == 1) {
			cnt++;
		} else {
			fprintf(stderr, "fscanf failed.");
		}
	}
	fprintf(stderr, "cnt=%d\n", cnt);

	in = fftw_malloc(cnt * sizeof(*in));
	out = fftw_malloc(cnt * sizeof(*out));
	N = 2 * (cnt - 1);
		
	plan = fftw_plan_r2r_1d(cnt, in, out, FFTW_REDFT00, FFTW_ESTIMATE);

	memcpy(in, tmp, cnt * sizeof(*tmp));

	fftw_execute(plan);

	for (i = 0; i < cnt; i++) {
		double k = 2. * M_PI * i / (N * dx);
		printf("%lg %lg\n", k, out[i]);
	}

	fftw_destroy_plan(plan);
	fftw_free(in);
	fftw_free(out);
	free(tmp);
	fftw_cleanup();

	return 0;
}

