#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <time.h>

#define ARGS_NUM 9

int Nx, Ny, MaxI;
double Hx, Hy, W;
double UDivisor;

double timestamp(void);
void getParams(int argc, char* argv[]);
double f(int n);
double calcU(int n, double *u);
double subsRow(int n, double *u);
void sor(double *x, double *r, double *timeSor, double *timeResNorm);