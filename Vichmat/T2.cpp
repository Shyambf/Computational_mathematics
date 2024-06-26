﻿#include "t2.h"


using namespace std;

void CreateRealFiles(int matrixSizes[], int n)
{
	srand(time(0));

	char filenameA[256];
	char filenameB[256];

	//int matrixSizes[] = { 100, 250, 500 };
	for (int k = 0; k < n; k++)
	{
		int N = matrixSizes[k];
		sprintf_s(filenameA, "txt/matrixA%d.txt", N);
		sprintf_s(filenameB, "txt/matrixB%d.txt", N);
		ofstream matrixA(filenameA);
		ofstream matrixB(filenameB);

		if (!matrixA.is_open() || !matrixB.is_open())
		{
			cout << "Error!" << endl;
			return;
		}

		double* matrA = new double[N];
		double* matrB = new double[N];
		for (int i = 0; i < N; i++)
		{
			double sum = 0.0;
			for (int j = 0; j < N; j++)
			{
				matrA[j] = rand() % 100;
				sum += matrA[j];
			}
			matrA[i] = sum + 1;
			for (int j = 0; j < N; j++)
			{
				matrixA << matrA[j] << "\t";
			}
			matrixA << endl;
			matrB[i] = rand() % 100;
			matrixB << matrB[i] << "\t";
		}
		matrixB.close();
		matrixA.close();
		delete[] matrB;
		delete[] matrA;
	}
}

void CreateComplexFiles(int matrixSizes[], int n)
{
	srand(time(0));

	char filenameA[256];
	char filenameB[256];

	//int matrixSizes[] ={ 100, 250, 500 };

	for (int k = 0; k < n; k++)
	{
		int N = matrixSizes[k];
		sprintf_s(filenameA, "txt/matrixAC%d.txt", N);
		sprintf_s(filenameB, "txt/matrixBC%d.txt", N);
		ofstream matrixA(filenameA);
		ofstream matrixB(filenameB);

		if (!matrixA.is_open() || !matrixB.is_open())
		{
			cout << "Error!" << endl;
			return;
		}

		complex<double>* matrA = new complex<double>[N];
		complex<double>* matrB = new complex<double>[N];
		for (int i = 0; i < N; i++)
		{
			complex<double> sum = 0.0;
			for (int j = 0; j < N; j++)
			{
				matrA[j].real(rand() % 100);
				matrA[j].imag(rand() % 100);
				sum += matrA[j];
			}
			matrA[i] = sum + 1.0;
			for (int j = 0; j < N; j++)
			{
				matrixA << matrA[j] << "\t";
			}
			matrixA << endl;
			matrB[i].real(rand() % 100);
			matrB[i].imag(rand() % 100);
			matrixB << matrB[i] << "\t";
		}
		matrixB.close();
		matrixA.close();
		delete[] matrB;
		delete[] matrA;
	}
}

void YakobiR(int n)
{
	const int N = n;
	const double eps = 0.0001;
	double sum1 = 0.0;
	double sum2 = 0.0;
	double xmax = -1.0;
	double tmp = 0.0;
	double tn, tk;


	double** matrixA = new double* [N];
	double* matrixB = new double[N];
	double* matrixX = new double[N];
	double* matrixXX = new double[N];
	char filenameA[256];
	char filenameB[256];
	char filenameX[256];

	sprintf_s(filenameA, "txt/matrixA%d.txt", N);
	sprintf_s(filenameB, "txt/matrixB%d.txt", N);
	sprintf_s(filenameX, "txt/matrixYakobiXR%d.txt", N);
	ifstream matrA(filenameA);
	ifstream matrB(filenameB);
	ofstream matrX(filenameX);

	if (!matrA.is_open() || !matrB.is_open() || !matrX.is_open())
	{
		cout << "error" << endl;
		return;
	}

	for (int i = 0; i < N; i++)
	{
		matrixA[i] = new double[N];
	}

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			matrA >> matrixA[i][j];
		}
		matrB >> matrixB[i];
		matrixX[i] = matrixXX[i] = 1.0;
	}

	tn = omp_get_wtime();
	do {
		xmax = -1.0;
		for (int i = 0; i < N; i++)
		{
			sum1 = 0.0;
			sum2 = 0.0;
			for (int j = 0; j < n; j++)
			{
				if (j < i)
					sum1 += matrixA[i][j] * matrixX[j];
				if (j > i)
					sum2 += matrixA[i][j] * matrixX[j];
			}
			matrixXX[i] = (1.0 / matrixA[i][i]) * (matrixB[i] - sum1 - sum2);
			tmp = abs(matrixX[i] - matrixXX[i]);
			if (tmp > xmax) xmax = tmp;
		}
		for (int i = 0; i < N; i++)
			matrixX[i] = matrixXX[i];
	} while (xmax > eps);

	tk = omp_get_wtime();
	double deltat = tk - tn;

	for (int i = 0; i < N; i++)
	{
		matrX << "x[" << i << "] = " << matrixX[i] << endl;
	}

	matrA.close();
	matrB.close();
	matrX.close();

	for (int i = 0; i < N; i++)
	{
		delete[] matrixA[i];
	}
	delete[] matrixA;
	delete[] matrixB;
	delete[] matrixX;

	cout << deltat << endl;
}

void YakobiC(int n)
{
	const int N = n;
	const double eps = 0.0001;
	complex<double> sum1 = 0.0;
	complex<double> sum2 = 0.0;
	double xmax = -1.0;
	double tmp = 0.0;
	double tn, tk;


	complex<double>** matrixA = new complex<double>*[N];
	complex<double>* matrixB = new complex<double>[N];
	complex<double>* matrixX = new complex<double>[N];
	complex<double>* matrixXX = new complex<double>[N];
	char filenameA[256];
	char filenameB[256];
	char filenameX[256];

	sprintf_s(filenameA, "txt/matrixAC%d.txt", N);
	sprintf_s(filenameB, "txt/matrixBC%d.txt", N);
	sprintf_s(filenameX, "txt/matrixYakobiXC%d.txt", N);
	ifstream matrA(filenameA);
	ifstream matrB(filenameB);
	ofstream matrX(filenameX);

	if (!matrA.is_open() || !matrB.is_open() || !matrX.is_open())
	{
		cout << "Error" << endl;
		return;
	}

	for (int i = 0; i < N; i++)
	{
		matrixA[i] = new complex<double>[N];
	}

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			matrA >> matrixA[i][j];
		}
		matrB >> matrixB[i];
		matrixX[i] = matrixXX[i] = 1.0;
	}
	tn = omp_get_wtime();

	do {
		xmax = -1.0;
		for (int i = 0; i < N; i++)
		{
			sum1 = 0.0;
			sum2 = 0.0;
			for (int j = 0; j < n; j++)
			{
				if (j < i)
					sum1 += matrixA[i][j] * matrixX[j];
				if (j > i)
					sum2 += matrixA[i][j] * matrixX[j];
			}
			matrixXX[i] = (1.0 / matrixA[i][i]) * (matrixB[i] - sum1 - sum2);
			tmp = abs(matrixX[i] - matrixXX[i]);
			if (tmp > xmax) xmax = tmp;
		}
		for (int i = 0; i < N; i++)
			matrixX[i] = matrixXX[i];
	} while (xmax > eps); //взять real

	tk = omp_get_wtime();
	double deltat = tk - tn;

	for (int i = 0; i < N; i++)
	{
		matrX << "x[" << i << "] = " << matrixX[i] << endl;
	}

	matrA.close();
	matrB.close();
	matrX.close();

	for (int i = 0; i < N; i++)
	{
		delete[] matrixA[i];

	}
	delete[] matrixA;
	delete[] matrixB;
	delete[] matrixX;
	delete[] matrixXX;

	cout << deltat << endl;
}

void GausZeyd(int n)
{
	const int N = n;
	const double eps = 0.0001;
	double sum1 = 0.0;
	double sum2 = 0.0;
	double xmax = -1.0;
	double tmp = 0.0;
	double tn, tk;
	double xv;


	double** matrixA = new double* [N] { 0 };
	double* matrixB = new double[N];
	double* matrixX = new double[N];
	char filenameA[256];
	char filenameB[256];
	char filenameX[256];

	sprintf_s(filenameA, "txt/matrixA%d.txt", N);
	sprintf_s(filenameB, "txt/matrixB%d.txt", N);
	sprintf_s(filenameX, "txt/matrixGausZeydXR%d.txt", N);
	ifstream matrA(filenameA);
	ifstream matrB(filenameB);
	ofstream matrX(filenameX);

	if (!matrA.is_open() || !matrB.is_open() || !matrX.is_open())
	{
		cout << "Ошибка открытия файла" << endl;
		return;
	}

	for (int i = 0; i < N; i++)
	{
		matrixA[i] = new double[N] { 0 };
	}

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			matrA >> matrixA[i][j];
		}
		matrB >> matrixB[i];
		matrixX[i] = 1.0;
	}

	tn = omp_get_wtime();
	do {
		xmax = -1.0;
		for (int i = 0; i < N; i++)
		{
			xv = matrixX[i];
			sum1 = 0.0;
			sum2 = 0.0;
			for (int j = 0; j < n; j++)
			{
				if (j < i)
					sum1 += matrixA[i][j] * matrixX[j];
				if (j > i)
					sum2 += matrixA[i][j] * matrixX[j];
			}
			matrixX[i] = (1.0 / matrixA[i][i]) * (matrixB[i] - sum1 - sum2);
			tmp = abs(matrixX[i] - xv);
			if (tmp > xmax) xmax = tmp;
		}
	} while (xmax > eps);

	tk = omp_get_wtime();
	double deltat = tk - tn;

	for (int i = 0; i < N; i++)
	{
		matrX << "x[" << i << "] = " << matrixX[i] << endl;
	}

	matrA.close();
	matrB.close();
	matrX.close();

	for (int i = 0; i < N; i++)
	{
		delete[] matrixA[i];

	}
	delete[] matrixA;
	delete[] matrixB;
	delete[] matrixX;


	cout << deltat << endl;
}

void GausZeydC(int n)
{
	const int N = n;
	const double eps = 0.0001;
	complex<double> sum1 = 0.0;
	complex<double> sum2 = 0.0;
	double xmax = -1.0;
	double tmp = 0.0;
	double tn, tk;
	complex<double> xv;


	complex<double>** matrixA = new complex<double>*[N] { 0 };
	complex<double>* matrixB = new complex<double>[N];
	complex<double>* matrixX = new complex<double>[N];
	char filenameA[256];
	char filenameB[256];
	char filenameX[256];

	sprintf_s(filenameA, "txt/matrixAC%d.txt", N);
	sprintf_s(filenameB, "txt/matrixBC%d.txt", N);
	sprintf_s(filenameX, "txt/matrixGausZeydXC%d.txt", N);
	ifstream matrA(filenameA);
	ifstream matrB(filenameB);
	ofstream matrX(filenameX);

	if (!matrA.is_open() || !matrB.is_open() || !matrX.is_open())
	{
		cout << "Ошибка открытия файла" << endl;
		return;
	}

	for (int i = 0; i < N; i++)
	{
		matrixA[i] = new complex<double>[N] { 0 };
	}

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			matrA >> matrixA[i][j];
		}
		matrB >> matrixB[i];
		matrixX[i] = 1.0;
	}

	tn = omp_get_wtime();
	do {
		xmax = -1.0;
		for (int i = 0; i < N; i++)
		{
			xv = matrixX[i];
			sum1 = 0.0;
			sum2 = 0.0;
			for (int j = 0; j < n; j++)
			{
				if (j < i)
					sum1 += matrixA[i][j] * matrixX[j];
				if (j > i)
					sum2 += matrixA[i][j] * matrixX[j];
			}
			matrixX[i] = (1.0 / matrixA[i][i]) * (matrixB[i] - sum1 - sum2);
			tmp = abs(matrixX[i] - xv);
			if (tmp > xmax) xmax = tmp;
		}
	} while (xmax > eps);

	tk = omp_get_wtime();
	double deltat = tk - tn;

	for (int i = 0; i < N; i++)
	{
		matrX << "x[" << i << "] = " << matrixX[i] << endl;
	}

	matrA.close();
	matrB.close();
	matrX.close();

	for (int i = 0; i < N; i++)
	{
		delete[] matrixA[i];

	}
	delete[] matrixA;
	delete[] matrixB;
	delete[] matrixX;


	cout << deltat << endl;
}

void Gaus(int n)
{
	const int N = n;
	double tn, tk;
	double w = 0.0, sum = 0.0;
	double deltat;


	double** matrixA = new double* [N] { 0 };
	double* matrixB = new double[N];
	double* matrixX = new double[N] { 0 };
	char filenameA[256];
	char filenameB[256];
	char filenameX[256];

	sprintf_s(filenameA, "txt/matrixA%d.txt", N);
	sprintf_s(filenameB, "txt/matrixB%d.txt", N);
	sprintf_s(filenameX, "txt/matrixGausXR%d.txt", N);
	ifstream matrA(filenameA);
	ifstream matrB(filenameB);
	ofstream matrX(filenameX);

	if (!matrA.is_open() || !matrB.is_open() || !matrX.is_open())
	{
		cout << "Ошибка открытия файла" << endl;
		return;
	}

	for (int i = 0; i < N; i++)
	{
		matrixA[i] = new double[N] { 0 };
	}

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			matrA >> matrixA[i][j];
		}
		matrB >> matrixB[i];
		//matrixX[i] = 1.0;
	}

	tn = omp_get_wtime();
	for (int i = 0; i <= n - 1; i++)
	{
		for (int j = i + 1; j <= n - 1; j++)
		{
			w = matrixA[j][i] / matrixA[i][i];
			for (int k = i; k <= N - 1; k++)
			{
				matrixA[j][k] = matrixA[j][k] - w * matrixA[i][k];
			}
			matrixB[j] = matrixB[j] - w * matrixB[i];
		}
	}

	for (int i = N - 1; i >= 0; i--)
	{
		sum = 0.0;
		for (int j = i + 1; j < N; j++)
		{
			//if (i != NULL)
				sum = sum + matrixA[i][j] * matrixX[j];
		}
		matrixX[i] = (matrixB[i] - sum) / matrixA[i][i];
	}

	tk = omp_get_wtime();

	for (int i = 0; i < N; i++)
	{
		matrX << "x[" << i << "] = " << matrixX[i] << endl;
	}

	matrA.close();
	matrB.close();
	matrX.close();

	for (int i = 0; i < N; i++)
	{
		delete[] matrixA[i];

	}
	delete[] matrixA;
	delete[] matrixB;
	delete[] matrixX;


	deltat = tk - tn;

	cout << deltat << endl;
}

void GausC(int n)
{
	const int N = n;
	double tn, tk;
	complex<double> w = 0.0, sum = 0.0;
	double deltat;


	complex<double>** matrixA = new complex<double>*[N] { 0 };
	complex<double>* matrixB = new complex<double>[N];
	complex<double>* matrixX = new complex<double>[N] { 0 };
	char filenameA[256];
	char filenameB[256];
	char filenameX[256];

	sprintf_s(filenameA, "txt/matrixAC%d.txt", N);
	sprintf_s(filenameB, "txt/matrixBC%d.txt", N);
	sprintf_s(filenameX, "txt/matrixGausXC%d.txt", N);
	ifstream matrA(filenameA);
	ifstream matrB(filenameB);
	ofstream matrX(filenameX);

	if (!matrA.is_open() || !matrB.is_open() || !matrX.is_open())
	{
		cout << "Ошибка открытия файла" << endl;
		return;
	}

	for (int i = 0; i < N; i++)
	{
		matrixA[i] = new complex<double>[N] { 0 };
	}

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			matrA >> matrixA[i][j];
		}
		matrB >> matrixB[i];
		matrixX[i] = 1.0;
	}
	tn = omp_get_wtime();

	for (int i = 0; i <= n - 1; i++)
	{
		for (int j = i + 1; j <= n - 1; j++)
		{
			w = matrixA[j][i] / matrixA[i][i];
			for (int k = i; k <= N - 1; k++)
			{
				matrixA[j][k] = matrixA[j][k] - w * matrixA[i][k];
			}
			matrixB[j] = matrixB[j] - w * matrixB[i];
		}
	}

	for (int i = N - 1; i >= 0; i--)
	{
		sum = 0.0;
		for (int j = i; j <= N - 1; j++)
		{
			if (i != NULL)
				sum = sum + matrixA[i][j + 1] * matrixX[j + 1];
		}
		matrixX[i] = (matrixB[i] - sum) / matrixA[i][i];
	}

	tk = omp_get_wtime();

	for (int i = 0; i < N; i++)
	{
		matrX << "x[" << i << "] = " << matrixX[i] << endl;
	}

	matrA.close();
	matrB.close();
	matrX.close();

	for (int i = 0; i < N; i++)
	{
		delete[] matrixA[i];

	}
	delete[] matrixA;
	delete[] matrixB;
	delete[] matrixX;



	deltat = tk - tn;

	cout << deltat << endl;
}
