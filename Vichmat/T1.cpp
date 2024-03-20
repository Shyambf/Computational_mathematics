#include <iostream>
#include <omp.h>
#include <string.h>

using namespace std;

double F1(double* y, double t)
{
    return 2.0 * (y[0] - y[0] * y[1]);
}

double F2(double* y, double t)
{
    return -(y[1] - y[0] * y[1]);
}


void print_array(double array[], int n, FILE *file) {
    for (int i = 0; i < n; i++) {
        fprintf(file, "%lf ", array[i]);
    }
    fprintf(file, "\n");
}


int eiler(double (*f[])(double*, double), double y[], double star_time, double end_time, double tau, int n, FILE *file) {
    double* yy = new double[n] {0.0};
    double start = omp_get_wtime();
    for (; star_time <= end_time; star_time += tau) {
        print_array(y, n, file);
        for (int i = 0; i < n; i++) {
            yy[i] = y[i] + tau * f[i](y, i);
        }
        for (int i = 0; i < n; i++) {
            y[i] = yy[i];
        }
    }

    printf_s("%f\n", omp_get_wtime() - start);
    delete[] yy;
    return 0;
}


int runge_kutta_2nd(double (*f[])(double*, double), double y[], double star_time, double end_time, double tau, int n, FILE *file) {
    double* yy = new double[n] {0.0};
    double* fk = new double[n] {0.0};
    double* ff = new double[n] {0.0};
    double start = omp_get_wtime();

    for (; star_time <= end_time; star_time += tau) {
        print_array(y, n, file);
        for (int i = 0; i < n; i++) {
            fk[i] = f[i](y, i);
        }

        for (int i = 0; i < n; i++) {
            yy[i] = y[i] + 0.5 * tau * fk[i];
        }

        for (int i = 0; i < n; i++) {
            ff[i] = f[i](yy, i);
        }

        for (int i = 0; i < n; i++) {
            y[i] = y[i] + tau * ff[i];
        }

    }
    printf_s("%f\n", omp_get_wtime() - start);
    delete[] yy;
    delete[] fk;
    delete[] ff;
    return 1;
}

int predictor_corrector(double (*f[])(double*, double), double y[], double star_time, double end_time, double tau, int n, FILE *file) {
    double* yy = new double[n] {0.0};
    double* fk = new double[n] {0.0};
    double* ff = new double[n] {0.0};
    double start = omp_get_wtime();

    for (; star_time <= end_time; star_time += tau) {
        print_array(y, n, file);
        for (int i = 0; i < n; i++) {
            fk[i] = f[i](y, i);
        }

        for (int i = 0; i < n; i++) {
            yy[i] = y[i] + tau * fk[i];
        }

        for (int i = 0; i < n; i++) {
            ff[i] = f[i](yy, i);
        }

        for (int i = 0; i < n; i++) {
            y[i] = y[i] + 0.5 * tau * (fk[i] + ff[i]);
        }

    }
    printf_s("%f\n", omp_get_wtime() - start);
    delete[] yy;
    delete[] fk;
    delete[] ff;
    return 1;
}


int runge_kutta_4nd(double (*f[])(double*, double), double y[], double star_time, double end_time, double tau, int n, FILE* file) {
    const int m = 4;
    double t = star_time;
    double* yy = new double [n] { 0.0 };
    double* ff = new double [n] { 0.0 };
    double** R = new double* [m] {};
    for (int i = 0; i < m; i++)
    {
        R[i] = new double[n] {0.0};
    }

    double start = omp_get_wtime();
    for (double t = star_time; t < end_time; t += tau)
    {
        print_array(y, n, file);
        for (int i = 0; i < n; i++)
        {
            R[0][i] = tau * f[i](y, t);
        }

        for (int i = 0; i < n; i++)
        {
            yy[i] = y[i] + 0.5 * R[0][i];
        }

        for (int i = 0; i < n; i++)
        {
            R[1][i] = tau * f[i](yy, t + 0.5 * tau);

        }

        for (int i = 0; i < n; i++)
        {
            yy[i] = y[i] + 0.5 * R[1][i];
        }

        for (int i = 0; i < n; i++)
        {
            R[2][i] = tau * f[i](yy, t + 0.5 * tau);
        }

        for (int i = 0; i < n; i++)
        {
            yy[i] = y[i] + R[2][i];
        }

        for (int i = 0; i < n; i++)
        {
            R[3][i] = tau * f[i](yy, t + tau);
        }

        for (int i = 0; i < n; i++)
        {
            y[i] += (R[0][i] + 2.0 * R[1][i] + 2.0 * R[2][i] + R[3][i]) / 6.0;
        }
    }
    printf_s("%f\n", omp_get_wtime() - start);
    for (int i = 0; i < m; i++)
    {
        delete[] R[i];
    }
    delete[] R;
    delete[] yy;
    delete[] ff;
    return 1;

}


void ImplictEulerMethod(double (*f[])(double*, double), double y[], double start_time, double end_time, double tau, int n, FILE* file)
{
    double* yy = new double[n] {0.0};
    double* b = new double[n] {0.0};
    double* p = new double[n] {0.0};
    double* opredn = new double[n] { 0.0 };
    double** a = new double* [n];
    for (int i = 0; i < n; i++)
    {
        a[i] = new double[n] {0.0};
    }

    double tn, tk, deltat;
    double deltah = 0.0001, opred;


    double start = omp_get_wtime();
    for (double t = start_time; t < end_time; t += tau)
    {
        print_array(y, n, file);
        for (int i = 0; i < n; i++)
        {
            b[i] = -f[i](y, t);
        }

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                memcpy(yy, y, sizeof(y) * n);
                yy[j] += deltah;
                a[i][j] = (f[i](yy, t) + b[i]) / deltah;
            }
        }

        a[0][0] -= (double)(1.0 / tau);
        a[1][1] -= (double)(1.0 / tau);

        opred = a[0][0] * a[1][1] - a[0][1] * a[1][0];
        if (!opred)
        {
            exit(1);
        }

        opredn[0] = b[0] * a[1][1] - b[1] * a[0][1];
        opredn[1] = b[1] * a[0][0] - b[0] * a[1][0];

        for (int i = 0; i < n; i++)
        {
            p[i] = opredn[i] / opred;
        }

        for (int i = 0; i < n; i++)
        {
            y[i] += p[i];
        }
    }
    printf_s("%f\n", omp_get_wtime() - start);
    for (int i = 0; i < n; i++)
    {
        delete[] a[i];
    }
    delete[] a;
    delete[] yy;
    delete[] opredn;
    delete[] b;
    delete[] p;
}


int main() {

    const int n = 2;
    int count_func = 5;
    double y[n] = { 1.0, 3.0 };
    double t0 = 0.0, tmax = 1000.0, tau = 0.001;
    double y_temp[n];
    double (*f[n])(double*, double) = { F1, F2 };

    FILE* file = nullptr;


    fopen_s(&file, "output0.txt", "w");
    memcpy(&y_temp, &y, sizeof(y)); 
    eiler(f, y_temp, t0, tmax, tau, n, file);
    fclose(file);


    fopen_s(&file, "output1.txt", "w");
    memcpy(&y_temp, &y, sizeof(y));
    runge_kutta_2nd(f, y_temp, t0, tmax, tau, n, file);
    fclose(file);


    fopen_s(&file, "output2.txt", "w");
    memcpy(&y_temp, &y, sizeof(y));
    predictor_corrector(f, y_temp, t0, tmax, tau, n, file);
    fclose(file);

    fopen_s(&file, "output3.txt", "w");
    memcpy(&y_temp, &y, sizeof(y));
    runge_kutta_4nd(f, y_temp, t0, tmax, tau, n, file);
    fclose(file);

    fopen_s(&file, "output4.txt", "w");
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
    return 0;
}
