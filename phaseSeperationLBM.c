#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#define nx 100
#define ny 100

int main()
{
    int i, j, a, ia, ja, q = 9;
    int ex[q], ey[q], kb[q];
    double ***f, ***ft;
    double **ux, **uy, **rho;
    double tau = 1.0, **psi, **fx, **fy; // tau is the relaxation time; psi is the pseudopotential
    double num, den;
    double source[q];
    double w[q];
    double rho0 = 1.0; 
    double temp1, temp2, feq[9];
    int ts = 0, time = 10000;
    double wm[q];
    double G = -0.44 / rho0;
    double err = 1, eps = 1.0e-03;
    int count = 0, flag = 0;
    double cont = 1.0 - 0.5 / tau;
    double visc = (tau - 0.5) / 3.0;
    FILE *soln;
    double umax, rho_max, rho_min;
    double log2 = log(2.0);
    char name[10] = {"solu"}, name1[20];
    char name2[20], name3[5] = {".dat"}, str_u[3] = {"_"};

    ex[0] = 0; ey[0] = 0;
    ex[1] = 1; ey[1] = 0;
    ex[2] = 0; ey[2] = 1;
    ex[3] = -1; ey[3] = 0;
    ex[4] = 0; ey[4] = -1;
    ex[5] = 1; ey[5] = 1;
    ex[6] = -1; ey[6] = 1;
    ex[7] = -1; ey[7] = -1;
    ex[8] = 1; ey[8] = -1;

    kb[0] = 0; kb[1] = 3; kb[2] = 4; kb[3] = 1; kb[4] = 2; kb[5] = 7; kb[6] = 8; kb[7] = 5; kb[8] = 6;

    for (a = 0; a < 9; a++) {
        if (a == 0) { w[a] = 4.0 / 9.0; wm[a] = 1.0;
        }
        else if (a >= 1 && a <= 4) {
        w[a] = 1.0/9.0; wm[a] = 1.0;
        }
        else {
            w[a] = 1.0 / 36.0; wm[a] = 0.25;
            }
    }

    f = (double ***)malloc(q * sizeof(double **));
    ft = (double ***)malloc(q * sizeof(double **));

    for (a = 0; a < q; a++) {
        f[a] = (double **)malloc(nx * sizeof(double *));
        ft[a] = (double **)malloc(nx * sizeof(double *));
        for (i = 0; i < nx; i++) {
            f[a][i] = (double *)malloc(ny * sizeof(double));
            ft[a][i] = (double *)malloc(ny * sizeof(double));
        }
    }

    ux = (double **)malloc(nx * sizeof(double *));
    uy = (double **)malloc(nx * sizeof(double *));
    fx = (double **)malloc(nx * sizeof(double *));
    fy = (double **)malloc(nx * sizeof(double *));
    rho = (double **)malloc(nx * sizeof(double *));
    psi = (double **)malloc(nx * sizeof(double *));

    for (i = 0; i < nx; i++) {
        ux[i] = (double *)malloc(ny * sizeof(double));
        uy[i] = (double *)malloc(ny * sizeof(double));
        fx[i] = (double *)malloc(ny * sizeof(double));
        fy[i] = (double *)malloc(ny * sizeof(double));
        rho[i] = (double *)malloc(ny * sizeof(double));
        psi[i] = (double *)malloc(ny * sizeof(double));
    }

    for (a = 0; a < q; a++) {
        for (i = 0; i < nx; i++) {
            for (j = 0; j < ny; j++) {
                if (i == (nx/2) && j == (ny/2)) {
                    f[a][i][j] = w[a] * rho0 * log2 * (1.0 + eps);
                }
                else {
                    f[a][i][j] = w[a] * rho0 * log2;
                }

                ux[i][j] = 0.0;
                uy[i][j] = 0.0;
            }
        }
    }

    for (ts = 1; ts <= time; ts++) {
        count++;
        num = 0.0; den = 0.0;
        for (i = 0; i < nx; i++) {
            for (j = 0; j < ny; j++) {
                rho[i][j] = 0.0; ux[i][j] = 0.0; uy[i][j] = 0.0;
                for (a = 0; a < q; a++) {
                    rho[i][j] += f[a][i][j];
                    ux[i][j] += ex[a] * f[a][i][j];
                    uy[i][j] += ey[a] * f[a][i][j];
                }
            }
        }

        // Calculate interparticle force
        for (i = 0; i < nx; i++) {
            for (j = 0; j < ny; j++) {
                psi[i][j] = rho0 * (1.0 - exp(-rho[i][j] / rho0));
            }
        }

        for (i = 0; i < nx; i++) {
            for (j = 0; j < ny; j++) {
                temp1 = 0.0; temp2 = 0.0;
                for (a = 0; a < q; a++) {
                    ia = i + ex[a];
                    ja = j + ey[a];

                    if (ia > nx - 1) { ia = 0;}
                    if (ia < 0) { ia = nx - 1;}
                    if (ja > ny - 1) { ja = 0;}
                    if (ja < 0) { ja = ny - 1;}

                    temp1 += wm[a] * ex[a] * psi[ia][ja];
                    temp2 += wm[a] * ey[a] * psi[ia][ja];

                }
                fx[i][j] = - G * psi[i][j] * temp1;
                fy[i][j] = - G * psi[i][j] * temp2;
            }
        }

        for (i = 0; i < nx; i++){
            for (j = 0; j < ny; j++) {
                ux[i][j] += 0.5 * fx[i][j]; ux[i][j] /= rho[i][j];
                uy[i][j] += 0.5 * fy[i][j]; uy[i][j] /= rho[i][j];
                for (a = 0; a < q; a++) {
                    temp1 = ex[a] * ux[i][j] + ey[a] * uy[i][j];
                    temp2 = ux[i][j] * ux[i][j] + uy[i][j] * uy[i][j];
                    feq[a] = w[a] * rho[i][j] * (1.0 + 3.0 * temp1 + 4.5 * temp1 * temp1 - 1.5 * temp2);
                    temp1 = 3.0 * (ex[a] - ux[i][j]) * fx[i][j] + 3.0 * (ey[a] - uy[i][j]) * fy[i][j];
                    temp2 = 9.0 * (ex[a] * ux[i][j] + ey[a] * uy[i][j]) * (ex[a] * fx[i][j] + ey[a] * fx[i][j]);
                    source[a] = cont * w[a] * (temp1 + temp2);
                    ft[a][i][j] = f[a][i][j] - (f[a][i][j] - feq[a]) / tau + source[a];
                }
            }
        }

        for (i = 0; i < nx; i++) {
            for (j = 0; j < ny; j++) {
                for (a = 0; a < q; a++) {
                    ia = i + ex[a];
                    ja = j + ey[a];
                    if (ia < 0) {
                        ia = nx - 1;
                    }
                    if (ia > nx - 1) {
                        ia = 0;
                    }
                    if (ja < 0) {
                        ja = ny - 1;
                    }
                    if (ja > ny - 1) {
                        ja = 0;
                    }
                    f[a][ia][ja] = ft[a][i][j];
                }
            }
        }

        rho_max = rho[0][0]; rho_min = rho[0][0];
        for (i = 0; i < nx; i++) {
            for (j = 0; j < ny; j++) {
                if (rho_max < rho[i][j]) {
                    rho_max = rho[i][j];
                }
                if (rho_min > rho[i][j]) {
                    rho_min = rho[i][j];
                }
            }
        }
        if (ts % 100 == 0) {
            printf("ts = %d\t, Maximum density = %8.4f\t, Minimum density = %8.4f\t, Ratio = %8.4f\n", ts, rho_max, rho_min, rho_max / rho_min);

            strcpy(name1, name2);
            int ai = sprintf(name2, "%d", ts);
            strcat(name1, name2);
            strcat(name1, name3);
            soln = fopen(name1, "w");
            fprintf(soln, "Variable = x, y, u, v, rho\n");
            fprintf(soln, "Zone I = %d, J = %d\n\n", nx, ny);
            for (j = 0; j < ny; j++) {
                for (i = 0; i < nx; i++) {
                    fprintf(soln, "%d %d %8.6f %8.6f %8.6f", i, j, ux[i][j], uy[i][j], rho[i][j]);
                }
            }
            fclose(soln);
        }

    }

    return 0;
}