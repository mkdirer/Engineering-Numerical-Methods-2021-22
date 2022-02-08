#include <cmath>
#include <iostream>
#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

constexpr double delta = 1.0;
constexpr int nx = 40;
constexpr int ny = 40;
constexpr int n = (nx + 1) * (ny + 1);
constexpr double dt = 1.0;
constexpr int TA = 40;
constexpr int TB = 0;
constexpr int TC = 30;
constexpr int TD = 0;
constexpr double KB {0.1};
constexpr double KD {0.6};
constexpr int IT_MAX = 2000;


int l(const int i, const int j) {
    return i + j * (nx + 1);
}

double x(const int i) {
    return i * delta;
}

double y(const int j) {
    return j * delta;
}

void matrix_a(gsl_matrix* a) {

    // Wnętrze obszaru (szary obszar na rysunku)
    for (int i = 1; i <= nx - 1; i++) {
        for (int j = 1; j <= ny - 1; j++) {
            gsl_matrix_set(a, l(i, j), l(i, j) - nx - 1, dt / (2 * std::pow(delta, 2)));
            gsl_matrix_set(a, l(i, j), l(i, j) - 1, dt / (2 * std::pow(delta, 2)));
            gsl_matrix_set(a, l(i, j), l(i, j) + 1, dt / (2 * std::pow(delta, 2)));
            gsl_matrix_set(a, l(i, j), l(i, j) + nx + 1, dt / (2 * std::pow(delta, 2)));
            gsl_matrix_set(a, l(i, j), l(i, j), (-(2 * dt) / std::pow(delta, 2)) - 1);
        }
    }

    // WB Dirichleta (lewy i prawy brzeg)
    for (int j = 0; j <= ny; j++) {
        gsl_matrix_set(a, l(0, j), l(0, j), 1);
    }

    for (int j = 0; j <= ny; j++) {
        gsl_matrix_set(a, l(nx, j), l(nx, j), 1);
    }

    // WB von Neumanna na górnym brzegu dla chwili n + 1

    for (int i = 1; i <= nx - 1; i++) {
        gsl_matrix_set(a, l(i, ny), l(i, ny) - nx - 1, -1 / (KB * delta));
        gsl_matrix_set(a, l(i, ny), l(i, ny), 1 + 1 / (KB * delta));
    }

    // WB von Neumanna na dolnym brzegu dla chwili n + 1
    for (int i = 1; i <= nx - 1; i++) {
        gsl_matrix_set(a, l(i, 0), l(i, 0), 1 + 1 / (KD * delta));
        gsl_matrix_set(a, l(i, 0), l(i, 0) + nx + 1, -1 / (KD * delta));
    }
    // for(int i=0; i<nx;i++){
    //     for(int j=0; j<ny;j++){
    //         printf("%f ", gsl_matrix_get(a, i, j));
    //     }
    //     printf("\n");
    // }

}

void matrix_b(gsl_matrix* b) {
    // Wnętrze obszaru (szary obszar na rysunku)
    for (int i = 1; i <= nx - 1; i++) {
        for (int j = 1; j <= ny - 1; j++) {
            gsl_matrix_set(b, l(i, j), l(i, j) - nx - 1, -dt / (2 * std::pow(delta, 2)));
            gsl_matrix_set(b, l(i, j), l(i, j) - 1, -dt / (2 * std::pow(delta, 2)));
            gsl_matrix_set(b, l(i, j), l(i, j) + 1, -dt / (2 * std::pow(delta, 2)));
            gsl_matrix_set(b, l(i, j), l(i, j) + nx + 1, -dt / (2 * std::pow(delta, 2)));
            gsl_matrix_set(b, l(i, j), l(i, j), (2 * dt) / std::pow(delta, 2) - 1);
        }
    }
    // WB Dirichleta (lewy i prawy brzeg)
    for (int j = 0; j <= ny; j++) {
        gsl_matrix_set(b, l(0, j), l(0, j), 1);
    }

    for (int j = 0; j <= ny; j++) {
        gsl_matrix_set(b, l(nx, j), l(nx, j), 1);
    }
    // WB von Neumanna na górnym brzegu dla chwili n + 1
    for (int k = 0; k < n; k++) {
        for (int i = 1; i <= nx - 1; i++) {
            gsl_matrix_set(b, l(i, ny), k, 0);
        }
    }
    // WB von Neumanna na dolnym brzegu dla chwili n + 1
    for (int k = 0; k < n; k++) {
        for (int i = 1; i <= nx - 1; i++) {
            gsl_matrix_set(b, l(i, 0), k, 0);
        }
    }
}

void vector_c(gsl_vector* c) {
// WB Dirichleta (lewy i prawy brzeg)

    for (int j = 0; j <= ny; j++) {
        gsl_vector_set(c, l(0, j), 0);
    }

    for (int j = 0; j <= ny; j++) {
        gsl_vector_set(c, l(nx, j), 0);
    }

// WB von Neumanna na górnym brzegu dla chwili n + 1
    for (int i = 1; i <= nx - 1; i++) {
        gsl_vector_set(c, l(i, ny), TB);
    }
// WB von Neumanna na dolnym brzegu dla chwili n + 1
    for (int i = 1; i <= nx - 1; i++) {
        gsl_vector_set(c, l(i, 0), TD);
    }
}

void vector_T(gsl_vector* T) {
    // lewy brzeg
    for (int j = 0; j <= ny; j++) {
        gsl_vector_set(T, l(0, j), TA);
    }
// prawy brzeg
    for (int j = 0; j <= ny; j++) {
        gsl_vector_set(T, l(nx, j), TC);
    }
// reszta
    for (int i = 1; i <= nx - 1; i++) {
        for (int j = 0; j <= ny; j++) {
            gsl_vector_set(T, l(i, j), 0);
        }
    }
}

double d2T(gsl_vector* T, int k) {
    return ((gsl_vector_get(T, k + 1) - 2 * gsl_vector_get(T, k) + gsl_vector_get(T, k - 1)) / std::pow(delta, 2))
        + ((gsl_vector_get(T, k + nx + 1) - 2 * gsl_vector_get(T, k) + gsl_vector_get(T, k - nx - 1)) / std::pow(delta, 2));
}

int main() {

    gsl_matrix *a = gsl_matrix_calloc(n, n);
    gsl_matrix *b = gsl_matrix_calloc(n, n);
    gsl_vector *c = gsl_vector_calloc(n);

    gsl_vector *d = gsl_vector_calloc(n);
    gsl_vector *T = gsl_vector_calloc(n);

    gsl_permutation* p = gsl_permutation_calloc(n);
    

    matrix_a(a);
    matrix_b(b);
    vector_c(c);
    vector_T(T);

    int signum = 0;
    gsl_linalg_LU_decomp(a, p, &signum);

    for( int it=0;  it <= IT_MAX; it++ ){

        gsl_blas_dgemv( CblasNoTrans, 1, b, T, 0, d );
        gsl_blas_daxpy( 1, c, d );

        gsl_linalg_LU_solve( a, p, d, T );
        
        if( it == 100 || it == 200 || it == 500 || it == 1000 || it == 2000 ){
            FILE* f1 = fopen(("data/out_T_" + std::to_string(it) + ".dat").c_str(), "w");
            FILE* f2 = fopen(("data/out_d2T_" + std::to_string(it) + ".dat").c_str(), "w");
            for( int i=1; i<=nx-1; i++ ){
                for( int j=1; j<=ny-1; j++ ){
                    fprintf(f1, "%f %f %f\n", x(i), y(j), gsl_vector_get(T, l(i, j)));
                    fprintf(f2, "%f %f %f\n", x(i), y(j), d2T(T, l(i, j)));
                }
                fprintf(f1, "\n");
                fprintf(f2, "\n");
            }
            fclose(f1);
            fclose(f2);
        }
    }

    gsl_matrix_free(a);
    gsl_matrix_free(b);
    gsl_vector_free(c);
    gsl_vector_free(T);
    gsl_permutation_free(p);

    return 0;
}