#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>

#define ARGS_NUM 9

int Nx, Ny, MaxI;
double Hx, Hy, W;
double UDivisor;
double Pipi = M_PI * M_PI; // Pipi = Pi * Pi

void print_errno(void);
void getParams(int argc, char* argv[], FILE *fp);
double f(double x, double y);
void u(int i, int j, double **u, double **f);
double residueNorm(double *r);
