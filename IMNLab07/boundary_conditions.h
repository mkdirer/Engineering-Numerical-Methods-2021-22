#pragma once
#include <math.h>
#include "constants.h"

double Q_out(double q_in) {
    return q_in * (pow(YNY, 3) - pow(YJ1, 3) - 3.0 * YJ1 * YNY*YNY + 3.0 * YJ1*YJ1 * YNY) / pow(YNY, 3);
}

double bcA_psi(int j, double q_in) {
    double y = DELTA * j;
    double res = q_in / (2.0 * MU) * (pow(y, 3) / 3.0 - y*y / 2.0 * (YJ1 + YNY) + y * YJ1 * YNY);
    return res;
}

double bcC_psi(int j, double q_in) {
    double y = DELTA * j;
    double q_out = Q_out(q_in);
    return q_out / (2.0 * MU) * (pow(y, 3) / 3.0 - y*y / 2.0 * YNY) + (q_in * YJ1*YJ1 * (3.0 * YNY - YJ1)) / (12.0 * MU);
}

double bcA_zeta(int j, double q_in) {
    double y = DELTA * j;
    return q_in / (2.0 * MU) * (2.0 * y - YJ1 - YNY);
}

double bcC_zeta(int j, double q_in) {
    double y = DELTA * j;
    double q_out = Q_out(q_in);
    return q_out / (2.0 * MU) * (2.0 * y - YNY);
}

double bcB_zeta(double** psi, int i) {
    return 2.0 / (DELTA*DELTA) * (psi[i][NY-1] - psi[i][NY]);
}

double bcD_zeta(double** psi, int i) {
    return 2.0 / (DELTA*DELTA) * (psi[i][1] - psi[i][0]);
}

double bcE_zeta(double** psi, int j) {
    return 2.0 / (DELTA*DELTA) * (psi[I1+1][j] - psi[I1][j]);
}

double bcF_zeta(double** psi, int i) {
    return 2.0 / (DELTA*DELTA) * (psi[i][J1+1] - psi[i][J1]);
}

double bcEF_corner_zeta(double** zeta) {
    return 0.5 * (zeta[I1-1][J1] + zeta[I1][J1-1]);
}

void update_psi_bounds(double** psi, double q_in) {
    // A 
    for (int j = J1; j <= NY; ++j) {
        psi[0][j] = bcA_psi(j, q_in);
    }
    // B 
    for (int i = 1; i < NX; ++i) {
        psi[i][NY] = psi[0][NY];
    }
    // C 
    for (int j = 0; j <= NY; ++j) {
        psi[NX][j] = bcC_psi(j, q_in);
    }
    // D
    for (int i = I1; i < NX; ++i) {
        psi[i][0] = psi[0][J1];
    }
    // E
    for (int j = 1; j <= J1; ++j) {
        psi[I1][j] = psi[0][J1];
    }
    // F
    for (int i = 1; i <= I1; ++i) {
        psi[i][J1] = psi[0][J1];
    }
}

void update_zeta_bounds(double** psi, double** zeta, double q_in) {
    // A 
    for (int j = J1; j <= NY; ++j) {
        zeta[0][j] = bcA_zeta(j, q_in);
    }
    // C
    for (int j = 0; j <= NY; ++j) {
        zeta[NX][j] = bcC_zeta(j, q_in);
    }
    // B
    for (int i = 1; i < NX; ++i) {
        zeta[i][NY] = bcB_zeta(psi, i);
    }
    // D
    for (int i = I1+1; i < NX; ++i) {
        zeta[i][0] = bcD_zeta(psi, i);
    }
    // E
    for (int j = 1; j < J1; ++j) {
        zeta[I1][j] = bcE_zeta(psi, j);
    }
    // F
    for (int i = 1; i <= I1; ++i) {
        zeta[i][J1] = bcF_zeta(psi, i);
    }

    zeta[I1][J1] = bcEF_corner_zeta(zeta);
}
