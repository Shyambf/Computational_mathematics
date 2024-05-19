#pragma once

#include<complex>

void CoordinatedDescent(double (*Func)(double, double), double xmin, double xmax, double ymin, double ymax, double x, double y, double eps, double h);

void Scanirovanie(double (*Func)(double, double), double xmin, double xmax, double ymin, double ymax, double x, double y, double h, double h2);

void Gradient(double (*Func)(double, double), double X, double Y, double eps, double h);
