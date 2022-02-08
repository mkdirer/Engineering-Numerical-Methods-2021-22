#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define SIG_X(max_x) (0.1 * max_x)
#define SIG_Y(max_y) (0.1 * max_y)

typedef struct
{
    double del;
    double ome;
    double eps;
    double TOL;
    int xn;
    int yn;
    double max_x;
    double max_y;
    double V1;
    double V2;
} paramet;

double** matrix_alloc(int x, int y);
void matrix_free(double** matrix, int x);
double den(double x, double y, double max_x, double max_y);
double sum_calc(double** V, int max_x, int max_y, double del, double** den_values);
void global_calc(paramet params, FILE* file_sum, FILE* file_result, FILE* file_sigma);
void local_calc(paramet params, FILE* file);


int main(int argc, const char* argv[])
{
    const double del = 0.1;
    const int xn = 150.0;
    const int yn = 100.0;
    const double V1 = 10.0;
    const double V2 = 0.0;
    const double max_x = del * xn;
    const double max_y = del * yn;
    const double TOL = 1e-8;
    const double eps = 1.0;
    
    
    double omes_l[] = {1.0, 1.4, 1.8, 1.9};
    double omes_g[] = {0.6, 1.0};

    ///////// Relaksacja globalna /////////////
    FILE* f_glob_sig_1 = fopen("glob_sigma_ome1.dat", "w");
    FILE* f_glob_sig_2 = fopen("glob_sigma_ome2.dat", "w");
    FILE* f_glob_sol_1 = fopen("glob_solution_ome1.dat", "w");
    FILE* f_glob_sol_2 = fopen("glob_solution_ome2.dat", "w");
    FILE* f_glob_sum_1 = fopen("glob_sum_ome1.dat", "w");
    FILE* f_glob_sum_2 = fopen("glob_sum_ome2.dat", "w");
   
    

    paramet params = {xn, yn, max_x, max_y, V1, V2, del, omes_g[0], eps, TOL};
    global_calc(params, f_glob_sum_1, f_glob_sol_1, f_glob_sig_1);
    params.ome = omes_g[1];
    global_calc(params, f_glob_sum_2, f_glob_sol_2, f_glob_sig_2);

    fclose(f_glob_sum_1);
    fclose(f_glob_sum_2);

    ////////// Relaksacja lokalna //////////////
    FILE* f_loc_sum_1 = fopen("loc_sum_ome1.dat", "w");
    FILE* f_loc_sum_2 = fopen("loc_sum_ome2.dat", "w");
    FILE* f_loc_sum_3 = fopen("loc_sum_ome3.dat", "w");
    FILE* f_loc_sum_4 = fopen("loc_sum_ome4.dat", "w");

    params.ome = omes_l[0];
    local_calc(params, f_loc_sum_1);
    params.ome = omes_l[1];
    local_calc(params, f_loc_sum_2);
    params.ome = omes_l[2];
    local_calc(params, f_loc_sum_3);
    params.ome = omes_l[3];
    local_calc(params, f_loc_sum_4);

    fclose(f_loc_sum_1);
    fclose(f_loc_sum_2);
    fclose(f_loc_sum_3);
    fclose(f_loc_sum_4);
}

double** matrix_alloc(int x, int y)
{
    double** result = calloc(x, sizeof(double*));
    for (int i = 0; i < x; i++)
    {
        result[i] = calloc(y, sizeof(double));
    }
    return result;
}

void matrix_free(double** matrix, int x)
{
    for (int i = 0; i < x; i++)
    {
        free(matrix[i]);
    }
    free(matrix);
}

double den(double x, double y, double max_x, double max_y)
{
    double den1 = exp(-pow((x - 0.35 * max_x) / (SIG_X(max_x)), 2)
            - pow((y - 0.5 * max_y) / (SIG_Y(max_y)), 2));
    double den2 = -exp(-pow((x - 0.65 * max_x) / (SIG_X(max_x)), 2)
            - pow((y - 0.5 * max_y) / (SIG_Y(max_y)), 2));
    return den1 + den2;
}

double sum_calc(double** V, int max_x, int max_y, double del, double** den_values)
{
    double sum = 0.0;
    for (int i = 0; i < max_x; i++)
    {
        for (int j = 0; j < max_y; j++)
        {
            double temp1 = (V[i+1][j] - V[i][j]) / del;
            double temp2 = (V[i][j+1] - V[i][j]) / del;
            sum += (del * del) * (0.5 * temp1 * temp1 + 0.5 * temp2 * temp2 - den_values[i][j] * V[i][j]);
        }
    }
    return sum;
}

void global_calc(paramet params, FILE* file_sum, FILE* file_result, FILE* file_sigma)
{
    double** Vs = matrix_alloc(params.xn+1, params.yn+1);
    double** Vn = matrix_alloc(params.xn+1, params.yn+1);
    double** den_values = matrix_alloc(params.xn+1, params.yn+1);

    for (int i = 0; i < params.xn + 1; i++) {
        Vs[i][0] = Vn[i][0] = params.V1;
        Vs[i][params.yn] = Vn[i][params.yn] = params.V2;
        for (int j = 0; j < params.yn + 1; j++) {
            den_values[i][j] = den(i * params.del, j * params.del, params.max_x, params.max_y);
        }
    }

    int iter_count = 0;
    double prev_S, curr_S;
    curr_S = sum_calc(Vs, params.xn, params.yn, params.del, den_values);
    do {
        prev_S = curr_S;
        for (int i = 1; i < params.xn; i++)
        {
            for (int j = 1; j < params.yn; j++)
            {
                Vn[i][j] = 0.25 * (Vs[i+1][j] + Vs[i-1][j] + Vs[i][j+1] + Vs[i][j-1] + params.del * params.del
                        / params.eps * den_values[i][j]);
            }
        }

        for (int j = 0; j < params.yn + 1; j++)
        {
            Vn[0][j] = Vn[1][j];
            Vn[params.xn][j] = Vn[params.xn-1][j];
        }

        for (int i = 0; i < params.xn + 1; i++)
        {
            for (int j = 0; j < params.yn + 1; j++)
            {
                Vs[i][j] = (1.0 - params.ome) * Vs[i][j] + params.ome * Vn[i][j];
            }
        }
        curr_S = sum_calc(Vs, params.xn, params.yn, params.del, den_values);
        fprintf(file_sum, "%d %f\n", iter_count++, curr_S);
    }
    while (fabs((curr_S - prev_S) / prev_S) > params.TOL);

    for (int i = 0; i < params.xn + 1; ++i) {
        for (int j = 0; j < params.yn + 1; ++j) {
            fprintf(file_result, "%d %d %f\n", i, j, Vs[i][j]);
        }
        fprintf(file_result, "\n");
    }

    for (int i = 1; i < params.xn; ++i) {
        for (int j = 1; j < params.yn; ++j) {
            double temp1 = (Vs[i+1][j] - 2.0*Vs[i][j] + Vs[i-1][j]) / (params.del * params.del);
            double temp2 = (Vs[i][j+1] - 2.0*Vs[i][j] + Vs[i][j-1]) / (params.del * params.del);
            fprintf(file_sigma, "%d %d %f\n", i, j, temp1 + temp2 + den_values[i][j] / params.eps);
        }
        fprintf(file_sigma, "\n");
    }

    matrix_free(Vs, params.xn + 1);
    matrix_free(Vn, params.xn + 1);
    matrix_free(den_values, params.xn + 1);
}

void local_calc(paramet params, FILE* file)
{
    double **V = matrix_alloc(params.xn + 1, params.yn + 1);
    double **den_values = matrix_alloc(params.xn + 1, params.yn + 1);
    for (int i = 0; i < params.xn + 1; i++) {
        V[i][0] = params.V1;
        V[i][params.yn] = params.V2;
        for (int j = 0; j < params.yn + 1; j++) {
            den_values[i][j] = den(i * params.del, j * params.del, params.max_x, params.max_y);
        }
    }

    int iter_count = 0;
    double prev_S, curr_S;
    curr_S = sum_calc(V, params.xn, params.yn, params.del, den_values);
    do {
        for (int i = 1; i < params.xn; i++) {
            for (int j = 1; j < params.yn; j++) {
                V[i][j] = (1.0 - params.ome) * V[i][j] + params.ome / 4.0 * (V[i+1][j] + V[i-1][j] + V[i][j+1] + V[i][j-1]
                                                                   + params.del * params.del / params.eps * den_values[i][j]);
            }
        }
        for (int j = 0; j < params.yn + 1; j++) {
            V[0][j] = V[1][j];
            V[params.xn][j] = V[params.xn - 1][j];
        }
        prev_S = curr_S;
        curr_S = sum_calc(V, params.xn, params.yn, params.del, den_values);
        fprintf(file, "%d %f\n", ++iter_count, curr_S);
    } while (fabs((curr_S - prev_S) / prev_S) > params.TOL);

    matrix_free(V, params.xn + 1);
    matrix_free(den_values, params.xn + 1);
}