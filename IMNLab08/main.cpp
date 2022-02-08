#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <fstream>


const double PI = 3.14159265358979323846;

const int nx = 400;
const int ny = 90;
const int i_1 = 200;
const int i_2 = 210;
const int j_1 = 50;
const double delta = 0.01;
const double sigma = 10 * delta;
const double xA = 0.45;
const double yA = 0.45;
const int ITMAX = 10000;

void pole_predkosci(const double psi[nx+1][ny+1], double (&vx)[nx+1][ny+1], double (&vy)[nx+1][ny+1], double &vmax, double &deltaT, FILE *file1, FILE *file2, const bool write) {
    for (int i = 1; i <= nx - 1; i++) {
        for (int j = 1; j <= ny-1; j++) {
            vx[i][j] = (psi[i][j+1] - psi[i][j-1]) / (2.0 * delta);
            vy[i][j] = -(psi[i+1][j] - psi[i-1][j]) / (2.0 * delta);
        }
    }

    //na zastawce
    for (int i = i_1; i <= i_2; i++) {
        for (int j = 0; j <= j_1; j++) {
            vx[i][j] = 0.0;
            vy[i][j] = 0.0;
        }
    }
    // na dolnym i gornym brzegu
    for (int i = 1; i <= nx-1; i++) {
        vx[i][0] = 0.0;
        vy[i][ny] = 0.0;
    }
    //na lewym i prawym
    for (int j = 0; j <= ny; j++) {
        vx[0][j] = vx[1][j];
        vx[nx][j] = vx[nx-1][j];
    } 

    //zapis do pliku
    if (write) {
        for (int i = 0; i <= nx; i++) {
            for (int j = 0; j <= ny; j++) {
                fprintf(file1, "%g ", vx[i][j]);
                fprintf(file2, "%g ", vy[i][j]);
            }
            fprintf(file1, "\n");
            fprintf(file2, "\n");
        }
    }


    for (int i = 0; i <= nx; i++) {
        for (int j = 0; j <= ny; j++) {
            if(sqrt (pow(vx[i][j], 2) + pow(vy[i][j], 2)) > vmax){
                vmax = sqrt (pow(vx[i][j], 2) + pow(vy[i][j], 2));
            }
        }
    }

   deltaT = delta / (4.0 * vmax);
}

void AD(const double D, const double psi[nx+1][ny+1], FILE *file1, FILE *file2, FILE *file3, FILE *file4, const bool write) {

    double u0[nx+1][ny+1];
    double u1[nx+1][ny+1];
    double vx[nx+1][ny+1] = {};
    double vy[nx+1][ny+1] = {};

    double vmax = 0.0;
    double deltaT = 0.0;

    pole_predkosci(psi, vx, vy, vmax, deltaT,  file1, file2, write);

    printf("vmax = %g deltaT = %g \n", vmax, deltaT);

    for (int i = 0; i <= nx; i++) {
        for (int j = 0; j <= ny; j++) {
            u0[i][j] = (1.0 / (2.0 * PI * pow(sigma, 2))) * exp(-(pow((delta*i) - xA, 2) + pow((delta*j) - yA, 2)) / (2.0 * pow(sigma, 2)));
        }
    }

    int licznik = 0;
    for (int it = 1; it <= ITMAX; it++) {

        for (int i = 0; i <= nx; i++){
			for (int j = 0; j <= ny; j++) {
				u1[i][j] = u0[i][j];
			}
        }

        for (int k = 1; k <= 20; k++) {
            for (int i = 0; i <= nx; i++) {
                for (int j = 1; j <= ny - 1; j++) {
                    if (i < i_1 || i > i_2 || j > j_1) {
                        if (i==0) {
                           u1[i][j] =  ( 1.0/( 1.0+( (2.0*D*deltaT) / pow(delta, 2)) ) ) * ( u0[i][j] - (deltaT/2.0) * vx[i][j] *
                            ( ( (u0[i+1][j] - u0[nx][j])/(2.0*delta) ) + (u1[i+1][j] - u1[nx][j])/(2.0*delta) ) - (deltaT / 2.0) * vy[i][j] * 
                            ( ( u0[i][j+1] - u0[i][j-1] )/(2.0*delta) + (u1[i][j+1] - u1[i][j-1])/(2.0*delta) ) + (deltaT/2.0) * D * 
                            ( ( u0[i+1][j] + u0[nx][j] + u0[i][j+1] + u0[i][j-1] - 4*u0[i][j] )/pow(delta,2) + ( u1[i+1][j] + u1[nx][j] + u1[i][j+1] + u1[i][j-1] )/pow(delta,2) ));
                        } else if (i==nx) {
                           u1[i][j] =  ( 1.0/( 1.0+( (2.0*D*deltaT) / pow(delta, 2)) ) ) * ( u0[i][j] - (deltaT/2.0) * vx[i][j] *
                            ( ( (u0[0][j] - u0[i-1][j])/(2.0*delta) ) + (u1[0][j] - u1[i-1][j])/(2.0*delta) ) - (deltaT / 2.0) * vy[i][j] * 
                            ( ( u0[i][j+1] - u0[i][j-1] )/(2.0*delta) + (u1[i][j+1] - u1[i][j-1])/(2.0*delta) ) + (deltaT/2.0) * D * 
                            ( ( u0[0][j] + u0[i-1][j] + u0[i][j+1] + u0[i][j-1] - 4*u0[i][j] )/pow(delta,2) + ( u1[0][j] + u1[i-1][j] + u1[i][j+1] + u1[i][j-1] )/pow(delta,2) ));                            
                        } else {
                           u1[i][j] =  ( 1.0/( 1.0+( (2.0*D*deltaT) / pow(delta, 2)) ) ) * ( u0[i][j] - (deltaT/2.0) * vx[i][j] *
                            ( ( (u0[i+1][j] - u0[i-1][j])/(2.0*delta) ) + (u1[i+1][j] - u1[i-1][j])/(2.0*delta) ) - (deltaT / 2.0) * vy[i][j] * 
                            ( ( u0[i][j+1] - u0[i][j-1] )/(2.0*delta) + (u1[i][j+1] - u1[i][j-1])/(2.0*delta) ) + (deltaT/2.0) * D * 
                            ( ( u0[i+1][j] + u0[i-1][j] + u0[i][j+1] + u0[i][j-1] - 4*u0[i][j] )/pow(delta,2) + ( u1[i+1][j] + u1[i-1][j] + u1[i][j+1] + u1[i][j-1] )/pow(delta,2) ));                             
                        }
                    }
                }
            }
        }

        for (int i = 0; i <= nx; i++){
			for (int j = 0; j <= ny; j++) {
				u0[i][j] = u1[i][j];
			}
        }

        double c = 0.0;
        double xsr = 0.0;

        for (int i = 0; i <= nx; i++){
			for (int j = 0; j <= ny; j++) {
                c+= u0[i][j];
                xsr += (delta * i) * u0[i][j];
            }
        }

        c*= pow(delta, 2);
        xsr *= pow(delta, 2);


        fprintf(file3, "%g %g %g \n",it*deltaT, c, xsr);

        if (it%700== 0 && licznik < 5) {
            for (int i = 0; i <= nx;i++) {
		        for (int j = 0; j <= ny; j++) {
			        fprintf(file4, "%g ",u1[i][j]);	
	    	    }
		        fprintf(file4, "\n");
            }
             fprintf(file4, "\n\n");
             licznik++;
             printf("%d \n", licznik);
        }
       
    }
    fprintf(file3, "\n\n");	
}




int main() {

    std::ifstream inFile;
    FILE *file1;
	    file1 = fopen("vx.dat", "w");
    FILE *file2;
	    file2 = fopen("vy.dat", "w");
    FILE *file3;
	    file3 = fopen("c_xsr.dat", "w");
    FILE *file4;
	    file4 = fopen("mapa_0.dat", "w");
    FILE *file5;
	    file5 = fopen("mapa_0.1.dat", "w");

    inFile.open("psi.dat");

    double psi[nx+1][ny+1] = {};
    const double D1 = 0.0;
    const double D2 = 0.1;

    int i_psi;
    int j_psi;

    for (int i = 0; i <= nx; i++) {
        for(int j = 0; j <= ny; j++) {
            inFile>>i_psi>>j_psi>>psi[i][j];
        }
    }

    AD(D1, psi, file1, file2, file3, file4, true);
    AD(D2, psi, file1, file2, file3, file5, false);


    inFile.close();
    fclose(file1);
    fclose(file2);
    fclose(file3);
    fclose(file4);
    fclose(file5);

    return 0;
}