#include "t3.h"
#include <utility>
#include <limits.h>
#include <Windows.h>
#include <math.h>
#include <omp.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void CoordinatedDescent(double (*Func)(double, double), double xmin, double xmax,
	double ymin, double ymax,
	double x, double y, double eps, double h)
{
	int countIteration = 0;
	double temp, resx, resy;
	double tempxmin = xmin, tempxmax = xmax, tempymin = ymin, tempymax = ymax;

	double tn = omp_get_wtime();

	do
	{
		temp = Func(x, y);
		double tempx1 = xmax - (xmax - xmin) / 1.618, tempx2 = xmin + (xmax - xmin) / 1.618;
		double FuncX1 = Func(tempx1, y), FuncX2 = Func(tempx2, y);
		while (fabs(tempx1 - tempx2) >= eps)
		{
			FuncX1 >= FuncX2 ? xmin = tempx1 : xmax = tempx2;
			tempx1 = xmax - (xmax - xmin) / 1.618;
			tempx2 = xmin + (xmax - xmin) / 1.618;
			FuncX1 = Func(tempx1, y);
			FuncX2 = Func(tempx2, y);
		}

		resx = (xmax + xmin) / 2.0;

		double tempy1 = ymax - (ymax - ymin) / 1.618, tempy2 = ymin + (ymax - ymin) / 1.618;
		double FuncY1 = Func(x, tempy1), FuncY2 = Func(x, tempy2);

		while (fabs(tempy1 - tempy2) >= eps)
		{
			FuncY1 >= FuncY2 ? ymin = tempy1 : ymax = tempy2;
			tempy1 = ymax - (ymax - ymin) / 1.618;
			tempy2 = ymin + (ymax - ymin) / 1.618;
			FuncY1 = Func(x, tempy1);
			FuncY2 = Func(x, tempy2);
		}

		resy = (ymax + ymin) / 2.0;
		countIteration++;

		x = resx;
		y = resy;

		xmin = tempxmin;
		xmax = tempxmax;
		ymin = tempymin;
		ymax = tempymax;

	} while (fabs(Func(x, y) - temp) >= h);

	double tk = omp_get_wtime();

	printf("Time = %lf\t CountIteration = %d\n", tk - tn, countIteration);
	printf("X = %lf\tY = %lf\tMin -> Func(X, Y) = %lf\n", resx, resy, Func(resx, resx - resy));
}

void Scanirovanie(double (*Func)(double, double), double xmin, double xmax, double ymin, double ymax, double x, double y, double h1, double h2)
{
	double rez, w, scan_min, x_scan = xmin, y_scan = ymin;
	int kolvo_sections = (int)((xmax - xmin) / h1) + 1;

	scan_min = 99999999.0;
	x = xmin; y = ymin;
	double tn = omp_get_wtime();
	for (int i = 0; i < kolvo_sections; i++)
	{
		y = ymin;
		for (int j = 0; j < kolvo_sections; j++)
		{
			w = Func(x, y);
			if (w < scan_min) {
				scan_min = w;
				x_scan = x;
				y_scan = y;
				rez = w;
			}
			y += h1;
		}
		x += h1;
	}
	printf("h1 = %1.3f rez = %f xsc = %1.3f ysc = % 1.3f\n", h1, rez, x_scan, y_scan);

	double xminh2 = x_scan - h1, xmaxh2 = x_scan + h1;
	if (xminh2 < xmin) xminh2 = xmin;
	if (xmaxh2 > xmax) xmaxh2 = xmax;
	double yminh2 = y_scan - h1, ymaxh2 = y_scan + h1;
	if (yminh2 < ymin) yminh2 = ymin;
	if (ymaxh2 > ymax) ymaxh2 = ymax;
	kolvo_sections = (int)((xmaxh2 - xminh2) / h2) + 1;
	int kolvo_strok = (int)((ymaxh2 - xminh2) / h2) + 1;
	scan_min = 99999999.0;
	x = xminh2; y = yminh2;
	for (int i = 0; i < kolvo_sections; i++) {
		y = yminh2;
		for (int j = 0; j < kolvo_strok; j++) {
			w = Func(x, y);
			if (w < scan_min) {
				scan_min = w;
				x_scan = x;
				y_scan = y;
				rez = w;
			}
			y += h2;
		}
		x += h2;
	}
	printf("h2 = % 1.3f rez = %f xsc = % 1.3f ysc = % 1.3f\n", h2, rez, x_scan, y_scan);
	double tk = omp_get_wtime();
	std::cout << "Time execute: " << tk - tn << std::endl;

}

void Gradient(double (*Func)(double, double), double X, double Y, double eps, double h)
{
	double x1 = X, x2 = Y, x1zapomnit, x2zapomnit, p1, p2;
	double sum = 100;
	double tngrad = omp_get_wtime();
	double kolvo_iteracii = 0;
	while (sum > eps && abs(x1) <= X && abs(x2) <= Y)
	{
		kolvo_iteracii++;
		p1 = (Func(x1 + eps, x2) - Func(x1, x2)) / eps;
		p2 = (Func(x1, x2 + eps) - Func(x1, x2)) / eps;
		x1zapomnit = x1;
		x1 -= h * p1;
		x2zapomnit = x2;
		x2 -= h * p2;
		sum = p1 * p1 + p2 * p2;
	}
	double tkgrad = omp_get_wtime();
	printf("x1 = % 1.9f x2 = %1.9f F = % 1.9f\n", x1zapomnit, x2zapomnit, Func(x1zapomnit, x2zapomnit));
	std::cout << "Time execute: " << tkgrad - tngrad << std::endl;
	std::cout << "Kolvo iteracii: " << kolvo_iteracii << std::endl;

}
