#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define XN 128
#define YN 128

double** matrix_alloc(int x, int y)
{
    double** wyn = calloc(x, sizeof(double*));
    for (int i = 0; i < x; i++)
    {
        wyn[i] = calloc(y, sizeof(double));
    }
    return wyn;
}

void matrix_free(double** A, int x)
{
    for (int i = 0; i < x; i++)
    {
        free(A[i]);
    }
    free(A);
}

double sum_cal(double** matrix_V, double del, int k)
{
    double sum = 0.0;
    for (int i = 0; i <= XN - k; i+=k) {
        for (int j = 0; j <= YN - k; j+=k) {
            double temp1 = (matrix_V[i+k][j] - matrix_V[i][j]) / (2.0 * k * del) + (matrix_V[i+k][j+k] - matrix_V[i][j+k]) / (2.0 * k * del);
            double temp2 = (matrix_V[i][j+k] - matrix_V[i][j]) / (2.0 * k * del) + (matrix_V[i+k][j+k] - matrix_V[i+k][j]) / (2.0 * k * del);
            sum += (k * k * del * del) / 2.0 * (temp1 * temp1 + temp2 * temp2);
        }
    }
    return sum;
}

double V_cal(int i, int j, int k, double** matrix_V)
{
    return 0.25 * (matrix_V[i+k][j] + matrix_V[i-k][j] + matrix_V[i][j+k] + matrix_V[i][j-k]);
}

void thicken(double** matrix_V, double del, int old_k)
{
    int new_k = old_k / 2;
    for (int i = 0; i < XN; i += old_k) {
        for (int j = 0; j < YN; j += old_k) {
            if (i + old_k < XN) {
                matrix_V[i+old_k][j+new_k] = 0.5 * (matrix_V[i+old_k][j] + matrix_V[i+old_k][j+old_k]);
            }
            if (j + old_k < YN) {
                matrix_V[i+new_k][j+old_k] = 0.5 * (matrix_V[i][j+old_k] + matrix_V[i+old_k][j+old_k]);
            }
            matrix_V[i+new_k][j+new_k] = 0.25 * (matrix_V[i][j] + matrix_V[i+old_k][j] + matrix_V[i][j+old_k] + matrix_V[i+old_k][j+old_k]);
        }
    }
}

void fill(double** matrix_V, double del, int k, double TOL, FILE* file)
{
    static int iter_count = 0;
    double prev_S, curr_S;
    curr_S = sum_cal(matrix_V, del, k);

    do {
        prev_S = curr_S;
        for (int i = k; i <= XN - k; i += k) {
            for (int j = k; j <= YN - k; j += k) {
                matrix_V[i][j] = V_cal(i, j, k, matrix_V);
            }
        }
        curr_S = sum_cal(matrix_V, del, k);
        fprintf(file, "%d %f\n", ++iter_count, curr_S);
    } while(fabs((curr_S - prev_S) / prev_S) >= TOL);
    fprintf(file, "\n\n");
}

void file_write(double** matrix_V, double del, int k, FILE* file)
{
    for (int i = 0; i < XN + 1; i += k) {
        for (int j = 0; j < YN + 1; j += k) {
            fprintf(file, "%f %f %f\n", i * del, j * del, matrix_V[i][j]);
        }
    }
}


int main(int argc, const char* argv[])
{
    const double TOL = 1e-8;
    const int k_max = 16;
    const double del = 0.2;
    const double x_max = del * XN;
    const double y_max = del * YN;
    

    double** matrix_V = matrix_alloc(XN + 1, YN + 1);
    for (int i = 0; i < XN + 1; ++i) {
        matrix_V[i][YN] = -1.0 * sin(2.0 * M_PI * del * i / y_max);
        matrix_V[i][0] = sin(2.0 * M_PI * del * i / y_max);
    }

    for (int i = 0; i < YN + 1; ++i) {
        matrix_V[0][i] = sin(M_PI * del * i / x_max);
        matrix_V[XN][i] = sin(M_PI * del * i / x_max);
    }

    FILE* f_sum = fopen("sum.dat", "w");
    FILE* fk1 = fopen("k1.dat", "w");
    FILE* fk2 = fopen("k2.dat", "w");
    FILE* fk4 = fopen("k4.dat", "w");
    FILE* fk8 = fopen("k8.dat", "w");
    FILE* fk16 = fopen("k16.dat", "w");

    FILE* files[] = {fk16, fk8, fk4, fk2, fk1};

    int k = k_max;
    int i = 0;
    while (k > 0) {
        fill(matrix_V, del, k, TOL, f_sum);
        file_write(matrix_V, del, k, files[i++]);
        thicken(matrix_V, del, k);
        k /= 2;
    }

    fclose(fk1);
    fclose(fk2);
    fclose(fk4);
    fclose(fk8);
    fclose(fk16);
    fclose(f_sum);
    matrix_free(matrix_V, XN + 1);
}

