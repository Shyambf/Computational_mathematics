#include <iostream>
#include <omp.h>
#include <string.h>


int eiler(
	double (*f[])(double*, double),
	double y[],
	double star_time,
	double end_time,
	double tau,
	int n,
	FILE* file
);


int runge_kutta_2nd(
	double (*f[])(double*, double),
	double y[],
	double star_time,
	double end_time,
	double tau,
	int n,
	FILE* file
);


int predictor_corrector(
	double (*f[])(double*, double),
	double y[],
	double star_time,
	double end_time,
	double tau,
	int n,
	FILE* file
);

int runge_kutta_4nd(
	double (*f[])(double*, double),
	double y[],
	double star_time,
	double end_time,
	double tau,
	int n,
	FILE* file
);


void ImplictEulerMethod(
	double (*f[])(double*, double),
	double y[],
	double star_time,
	double end_time,
	double tau,
	int n,
	FILE* file
);
