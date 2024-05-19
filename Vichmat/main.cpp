#include "t1.h"
#include "t2.h"
#include "t3.h"
#include <cstdlib>
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <iomanip>

using namespace std;

double F11(double* y, double t)
{
	return 2.0 * (y[0] - y[0] * y[1]);
}

double F21(double* y, double t)
{
	return -(y[1] - y[0] * y[1]);
}

double Celevaya_func(double x1, double x2)
{
	return 2.2 * x1 + -1.3 * x2 + exp(0.04 * x1 * x1 + 0.12 * x2 * x2);
}

void readLinesFromFile(const string& filename, vector<string>& lines, int start, int size) {
	ifstream file(filename);
	if (!file.is_open()) {
		cerr << "Unable to open file: " << filename << endl;
		return;
	}

	string line;
	int count = 0;
	while (getline(file, line) && count < size + start + 1) {
		if (count > start) {
			lines.push_back(line);
		}
		count++;

	}
	file.close();
}

void printTable(const vector<vector<string>>& data, int size) {
	int otst = 40;

	for (int i = 0; i < size; ++i) {
		for (const auto& fileLines : data) {
			if (i < fileLines.size()) {
				cout << left << setw(otst) << fileLines[i]; // Форматированный вывод
			}
			else {
				cout << left << setw(otst) << ""; // Если нет значения для этой строки, выводим пустую строку
			}
		}
		cout << endl;
	}
	
}
void theme1() {
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
}

void print_data() {
	vector<vector<string>> fileData;
	int start, N;
	char n;
	int size;
	rewind(stdin);
	std::cout << "Ведите тип:\n1 - R\n2 - C\n";
	n = getchar();
	if (n == 'R') {
		std::cout << "Ведите размер:\n1 - 100\n2 - 250\n3 - 500\n";
		scanf_s("%d", &N);
	}
	else {
		std::cout << "Ведите размер:\n1 - 50\n2 - 125\n3 - 250\n";
		scanf_s("%d", &N);
	}
	std::cout << "Ведите с какого элемента начать:\n";
	scanf_s("%d", &start);

	std::cout << "Ведите колличество элементов для сравнения:\n";
	scanf_s("%d", &size);

	system("pause");
	system("cls");
	char filenameA[256];
	char filenameB[256];
	char filenameX[256];

	std::cout << "type: " << n << N << "[" << start << ":" << size + start << "]" << std::endl;
	sprintf_s(filenameA, "txt/matrixYakobiX%c%d.txt", n, N);
	sprintf_s(filenameB, "txt/matrixGausZeydX%c%d.txt", n, N);
	sprintf_s(filenameX, "txt/matrixGausX%c%d.txt", n, N);
	vector<string> filenames = { filenameA, filenameB, filenameX };
	for (const auto& filename : filenames) {
		vector<string> lines;
		readLinesFromFile(filename, lines, start-1, size);
		fileData.push_back(lines);
	}

	printTable(fileData, size);
}

void metods() {

	std::cout << "якоби" << std::endl;
	std::cout << "Время выполнения 100 ";
	YakobiR(100);
	std::cout << std::endl;
	std::cout << "Время выполнения 250 ";
	YakobiR(250);
	std::cout << std::endl;
	std::cout << "Время выполнения 500 ";
	YakobiR(500);
	std::cout << std::endl;
	std::cout << std::endl;


	std::cout << "якоби коплекс" << std::endl;
	std::cout << "Время выполнения 50 ";
	YakobiC(50);
	std::cout << std::endl;
	std::cout << "Время выполнения 125 ";
	YakobiC(125);
	std::cout << std::endl;
	std::cout << "Время выполнения 250 ";
	YakobiC(250);
	std::cout << std::endl;
	std::cout << std::endl;


	std::cout << "гаус зейд" << std::endl;
	std::cout << "Время выполнения 100 ";
	GausZeyd(100);
	std::cout << std::endl;
	std::cout << "Время выполнения 250 ";
	GausZeyd(250);
	std::cout << std::endl;
	std::cout << "Время выполнения 500 ";
	GausZeyd(500);
	std::cout << std::endl;
	std::cout << std::endl;


	std::cout << "гаус зейд комплекс" << std::endl;
	std::cout << "Время выполнения 50 ";
	GausZeydC(50);
	std::cout << std::endl;
	std::cout << "Время выполнения 125 ";
	GausZeydC(125);
	std::cout << std::endl;
	std::cout << "Время выполнения 250 ";
	GausZeydC(250);
	std::cout << std::endl;
	std::cout << std::endl;


	std::cout << "гаус" << std::endl;
	std::cout << "Время выполнения 100 ";
	Gaus(100);
	std::cout << std::endl;
	std::cout << "Время выполнения 250 ";
	Gaus(250);
	std::cout << std::endl;
	std::cout << "Время выполнения 500 ";
	Gaus(500);
	std::cout << std::endl;
	std::cout << std::endl;


	std::cout << "гаус комплекс" << std::endl;
	std::cout << "Время выполнения 50 ";
	GausC(50);
	std::cout << std::endl;
	std::cout << "Время выполнения 125 ";
	GausC(125);
	std::cout << std::endl;
	std::cout << "Время выполнения 250 ";
	GausC(250);
	std::cout << std::endl;
	std::cout << std::endl;
}

int main() {
	setlocale(LC_ALL, "RU");
	if (std::system("mkdir txt") == 0) {
		std::cout << "Папка успешно создана." << std::endl;
	}

	// Theme 2


	int mode;
	int matrR[] = { 100, 250, 500 };
	int matrC[] = { 50, 125, 250 };
	do
	{
		system("cls");
		printf("1 - Создать обычные матрицы\n");
		printf("2 - Создать комплексные матрицы\n");
		printf("3 - Запустить методы\n");
		printf("4 - Сравнить\n");
		printf("5 - Тема 1\n");
		printf("6 - Тема 3\n");
		printf("7 - Выход\n");
		mode = getchar();
		switch (mode)
		{
		case '1':
			system("cls");
			CreateRealFiles(matrR, 3);
			system("pause");
			break;
		case '2':
			system("cls");
			CreateComplexFiles(matrC, 3);
			system("pause");
			break;
		case '3':
			system("cls");
			metods();
			system("pause");
			break;
		case '4':
			system("cls");
			print_data();
			system("pause");
			break;
		case '6':
			std::cout << "Answer CoordinatedDescent: " << "\n";
			CoordinatedDescent(Celevaya_func, -5.0, 5.0, -5.0, 5.0, 0, 0, 0.000001, 0.0001);
			std::cout << "\n\n" << "Answer Scanirovanie: " << "\n";
			Scanirovanie(Celevaya_func, -5.0, 5.0, -5.0, 5.0, 0, 0, 0.1, 0.001);
			std::cout << "\n\n" << "Answer Gradient: " << "\n";
			Gradient(Celevaya_func, 5.0, 5.0, 0.00001, 0.000001);
			system("pause");

			break;
		case '7':
			break;
		}

	} while (mode != '7');
	return 0;

}