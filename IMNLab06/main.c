#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mgmres.h"

double** alloc_matrix(int x, int y)
{
    double** result = calloc(x, sizeof(double*));
    for (int i = 0; i < x; i++)
    {
        result[i] = calloc(y, sizeof(double));
    }
    return result;
}

void free_matrix(double** matrix, int x)
{
    for (int i = 0; i < x; i++)
    {
        free(matrix[i]);
    }
    free(matrix);
}

int l_to_j(int l, int nx) {
    return l / (nx + 1);
}

int l_to_i(int l, int j, int nx) {
    return l - j * (nx + 1);
}

double eps_l(int l, int nx, double eps1, double eps2) {
    int j = l_to_j(l, nx);
    int i = l_to_i(l, j, nx);
    return i <= nx / 2 ? eps1:eps2;
}

double rho1(double x, double y, double max_x, double max_y, double sigma) {
    double temp1 = -pow((x - 0.25 * max_x) / sigma, 2);
    double temp2 = -pow((y - 0.5 * max_y) / sigma, 2);
    return exp(temp1 + temp2);
}

double rho2(double x, double y, double max_x, double max_y, double sigma) {
    double temp1 = -pow((x - 0.75 * max_x) / sigma, 2);
    double temp2 = -pow((y - 0.5 * max_y) / sigma, 2);
    return -1.0 * exp(temp1 + temp2);
}

void poisson_calculate(int nx, int ny, double V1, double V2, double V3, double V4, double* V,
               double* a, double* b, double delta, double eps1, double eps2, double sigma, int has_rho);

int main(int argc, const char* argv[]) {
    const double delta = 0.1;
    double eps1 = 1.0;
    double eps2 = 1.0;
    int V1 = 10;
    int V2 = -10;
    int V3 = 10;
    int V4 = -10;

    // Czesc pierwsza (nx, ny = 4, wypisac wektor i macierz)
    int nx = 4;
    int ny = 4;
    double max_x = nx * delta;
    double sigma = max_x / 10.0;
    int N = (nx + 1) * (ny + 1);

    double* a = calloc(5 * N, sizeof(double));
    double* b = calloc(N, sizeof(double));
    double* V = calloc(N, sizeof(double));

    poisson_calculate(nx, ny, V1, V2, V3, V4, V, a, b, delta, eps1, eps2, sigma, 0);

    FILE* file_nz = fopen("nonzero_el.dat", "w");
    fprintf(file_nz, "# l\ti_l\tj_l\tb[l]\n");
    for (int l = 0; l < N; l++) {
        int j = l_to_j(l, nx);
        int i = l_to_i(l, j, nx);
        fprintf(file_nz, "%d\t%d\t%d\t%f\n", l, i, j, b[l]);
    }
    fprintf(file_nz, "\n\n# k\ta[k]\n");
    int count = 0;
    for (int l = 0; l < 5 * N; l++) {
        if (a[l] != 0) {
            fprintf(file_nz, "%d\t%f\n", count++, a[l]);
        }
    }
    fclose(file_nz);
    free(a);
    free(b);
    free(V);


    // Czesc druga - mapy potencjalu, rho = 0
    FILE* file_n50 = fopen("n50.dat", "w");
    FILE* file_n100 = fopen("n100.dat", "w");
    FILE* file_n200 = fopen("n200.dat", "w");

    // nx=50
    nx = 50;
    ny = 50;
    max_x = nx * delta;
    sigma = max_x / 10.0;
    N = (nx + 1) * (ny + 1);

    a = calloc(5 * N, sizeof(double));
    b = calloc(N, sizeof(double));
    V = calloc(N, sizeof(double));

    poisson_calculate(nx, ny, V1, V2, V3, V4, V, a, b, delta, eps1, eps2, sigma, 0);
    for (int l = 0; l < N; l++) {
        int j = l_to_j(l, nx);
        int i = l_to_i(l, j, nx);
        fprintf(file_n50, "%f %f %f\n", i * delta, j * delta, V[l]);
    }
    free(a);
    free(b);
    free(V);

    // nx=100
    nx = 100;
    ny = 100;
    max_x = nx * delta;
    sigma = max_x / 10.0;
    N = (nx + 1) * (ny + 1);

    a = calloc(5 * N, sizeof(double));
    b = calloc(N, sizeof(double));
    V = calloc(N, sizeof(double));

    poisson_calculate(nx, ny, V1, V2, V3, V4, V, a, b, delta, eps1, eps2, sigma, 0);
    for (int l = 0; l < N; l++) {
        int j = l_to_j(l, nx);
        int i = l_to_i(l, j, nx);
        fprintf(file_n100, "%f %f %f\n", i * delta, j * delta, V[l]);
    }
    free(a);
    free(b);
    free(V);

    // nx=200
    nx = 200;
    ny = 200;
    max_x = nx * delta;
    sigma = max_x / 10.0;
    N = (nx + 1) * (ny + 1);

    a = calloc(5 * N, sizeof(double));
    b = calloc(N, sizeof(double));
    V = calloc(N, sizeof(double));

    poisson_calculate(nx, ny, V1, V2, V3, V4, V, a, b, delta, eps1, eps2, sigma, 0);
    for (int l = 0; l < N; l++) {
        int j = l_to_j(l, nx);
        int i = l_to_i(l, j, nx);
        fprintf(file_n200, "%f %f %f\n", i * delta, j * delta, V[l]);
    }
    free(a);
    free(b);
    free(V);

    fclose(file_n50);
    fclose(file_n100);
    fclose(file_n200);


    // Czesc trzecia - mapy potencjalu, rho != 0
    FILE* file_eps1 = fopen("eps1.dat", "w");
    FILE* file_eps2 = fopen("eps2.dat", "w");
    FILE* file_eps10 = fopen("eps10.dat", "w");

    nx = 100;
    ny = 100;
    max_x = nx * delta;
    sigma = max_x / 10.0;
    N = (nx + 1) * (ny + 1);
    V1 = V2 = V3 = V4 = 0;

    eps1 = eps2 = 1.0;

    a = calloc(5 * N, sizeof(double));
    b = calloc(N, sizeof(double));
    V = calloc(N, sizeof(double));

    poisson_calculate(nx, ny, V1, V2, V3, V4, V, a, b, delta, eps1, eps2, sigma, 1);
    for (int l = 0; l < N; l++) {
        int j = l_to_j(l, nx);
        int i = l_to_i(l, j, nx);
        fprintf(file_eps1, "%f %f %f\n", i * delta, j * delta, V[l]);
    }

    eps2 = 2.0;
    poisson_calculate(nx, ny, V1, V2, V3, V4, V, a, b, delta, eps1, eps2, sigma, 1);
    for (int l = 0; l < N; l++) {
        int j = l_to_j(l, nx);
        int i = l_to_i(l, j, nx);
        fprintf(file_eps2, "%f %f %f\n", i * delta, j * delta, V[l]);
    }

    eps2 = 10.0;
    poisson_calculate(nx, ny, V1, V2, V3, V4, V, a, b, delta, eps1, eps2, sigma, 1);
    for (int l = 0; l < N; l++) {
        int j = l_to_j(l, nx);
        int i = l_to_i(l, j, nx);
        fprintf(file_eps10, "%f %f %f\n", i * delta, j * delta, V[l]);
    }
    free(a);
    free(b);
    free(V);

    fclose(file_eps1);
    fclose(file_eps2);
    fclose(file_eps10);
}

void poisson_calculate(int nx, int ny, double V1, double V2, double V3, double V4, double* V, double* a, double* b,
                    double delta, double eps1, double eps2, double sigma, int has_rho) {
    double max_x = nx * delta;
    double max_y = ny * delta;
    int N = (nx + 1) * (ny + 1);

    int* ja = calloc(5 * N, sizeof(int));
    int* ia = calloc(N + 1, sizeof(int));
    for (int i = 0; i < N + 1; i++) {
        ia[i] = -1;
    }

    const int itr_max = 500;
    const int mr = 500;
    const double tol_abs = 1e-8;
    const double tol_rel = 1e-8;

    double* eps = calloc(N + nx + 1, sizeof(double));
    for (int l = 0; l < N + nx + 1; l++) {
        eps[l] = eps_l(l, nx, eps1, eps2);
    }

    double* rho = calloc(N, sizeof(double));
    for (int l = 0; l < N; l++) {
        if (has_rho) {
            int j = l_to_j(l, nx);
            int i = l_to_i(l, j, nx);
            rho[l] = -1.0 * (rho1(i * delta, j * delta, max_x, max_y, sigma)
                    + rho2(i * delta, j * delta, max_x, max_y, sigma));
        }
        else {
            rho[l] = 0.0;
        }
    }

    int k = -1;
    for (int l = 0; l < N; l++) {
        int edge = 0;
        double vb = 0;
        int j = l_to_j(l, nx);
        int i = l_to_i(l, j, nx);

        if (i == 0) {
            edge = 1;
            vb = V1;
        }
        if (j == ny) {
            edge = 1;
            vb = V2;
        }
        if (i == nx) {
            edge = 1;
            vb = V3;
        }
        if (j == 0) {
            edge = 1;
            vb = V4;
        }

        b[l] = rho[l];

        if (edge == 1) {
            b[l] = vb;
        }

        ia[l] = -1;

        // Left outer
        if (l - nx - 1 >= 0 && edge == 0) {
            k++;
            if (ia[l] < 0) {
                ia[l] = k;
            }
            a[k] = eps[l] / (delta * delta);
            ja[k] = l - nx - 1;
        }

        // Under diagonal
        if (l - 1 >= 0 && edge == 0) {
            k++;
            if (ia[l] < 0) {
                ia[l] = k;
            }
            a[k] = eps[l] / (delta * delta);
            ja[k] = l - 1;
        }

        // Diagonal
        k++;
        if (ia[l] < 0) {
            ia[l] = k;
        }
        if (edge == 0) {
            a[k] = -1.0 * (2.0 * eps[l] + eps[l+1] + eps[l+nx+1]) / (delta * delta);
        } else {
            a[k] = 1.0;
        }
        ja[k] = l;

        // Over diagonal
        if (l < N && edge == 0) {
            k++;
            a[k] = eps[l+1] / (delta * delta);
            ja[k] = l + 1;
        }

        // Right outer
        if (l < N - nx - 1 && edge == 0) {
            k++;
            a[k] = eps[l+nx+1] / (delta * delta);
            ja[k] = l + nx + 1;
        }
    }
    int nz_count = k + 1;
    ia[N] = nz_count;
    pmgmres_ilu_cr(N, nz_count, ia, ja, a, V, b, itr_max, mr, tol_abs, tol_rel);

    free(rho);
    free(eps);
    free(ia);
    free(ja);
}
