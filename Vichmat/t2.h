#include <fstream>
#include <iostream>
#include <cstring>
#include <cmath>
#include <omp.h>
#include <complex>

void CreateRealFiles(int matrixSizes[], int n);
void CreateComplexFiles(int matrixSizes[], int n);
void YakobiR(int n);
void YakobiC(int n);
void GausZeyd(int n);
void GausZeydC(int n);
void Gaus(int n);
void GausC(int n);
