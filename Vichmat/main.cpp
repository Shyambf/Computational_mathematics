#include "t1.h"
#include "t2.h"
#include <cstdlib>
#include <iostream>

using namespace std;

double F11(double* y, double t)
{
	return 2.0 * (y[0] - y[0] * y[1]);
}

double F21(double* y, double t)
{
	return -(y[1] - y[0] * y[1]);
}

int main() {
	if (std::system("mkdir txt") == 0) {
		std::cout << "Папка успешно создана." << std::endl;
	}
	else {
		std::cout << "Ошибка при создании папки." << std::endl;
	}
	const int n = 2;
	int count_func = 5;
	double y[n] = { 1.0, 3.0 };
	double t0 = 0.0, tmax = 10.0, tau = 0.001;
	double y_temp[n];
	double (*f[n])(double*, double) = { F11, F21 };

	FILE* file = nullptr;



	fopen_s(&file, "txt/output0.txt", "w");
	memcpy(&y_temp, &y, sizeof(y));
	eiler(f, y_temp, t0, tmax, tau, n, file);
	fclose(file);

	fopen_s(&file, "txt/output1.txt", "w");
	memcpy(&y_temp, &y, sizeof(y));
	runge_kutta_2nd(f, y_temp, t0, tmax, tau, n, file);
	fclose(file);


	fopen_s(&file, "txt/output2.txt", "w");
	memcpy(&y_temp, &y, sizeof(y));
	predictor_corrector(f, y_temp, t0, tmax, tau, n, file);
	fclose(file);

	fopen_s(&file, "txt/output3.txt", "w");
	memcpy(&y_temp, &y, sizeof(y));
	runge_kutta_4nd(f, y_temp, t0, tmax, tau, n, file);
	fclose(file);

	fopen_s(&file, "txt/output4.txt", "w");
	memcpy(&y_temp, &y, sizeof(y));
	ImplictEulerMethod(f, y_temp, t0, tmax, tau, n, file);
	fclose(file);
	rewind(stdin);
	printf("Show graph? [Y/N] - ");
	char Ans_User;
	Ans_User = getchar();
	if (Ans_User == 'Y') {
		char command[50];
		sprintf_s(command, "python plot.py %d %f %lf", count_func, tmax, tau);
		system(command);
	}

	// Theme 2

	int matrR[] = { 100, 250, 500 };
	int matrC[] = { 50, 125, 250 };
	CreateRealFiles(matrR, 3);
	YakobiR(100);
	YakobiR(250);
	YakobiR(500);
	return 1;
	return 1;
}