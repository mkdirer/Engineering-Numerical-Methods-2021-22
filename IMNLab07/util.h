#pragma once
#include "stdlib.h"

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

void zero_matrix(double** matrix, int x, int y) {
    for (int i = 0; i < x; ++i) {
        for (int j = 0; j < y; ++j) {
            matrix[i][j] = 0.0;
        }
    }
}

void matrix_to_file(double** matrix, int i_max, int j_max, FILE* file) {
    for (int i = 1; i < i_max; ++i) {
        for (int j = 1; j < j_max; ++j) {
            fprintf(file, "%d %d %f\n", i, j, matrix[i][j]);
        }
        fprintf(file, "\n");
    }
}

