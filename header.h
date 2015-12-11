#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <time.h>
#include <likwid.h>

#define ARGS_NUM 9
#define BLOCK_SIZE 3

double timestamp(void);
FILE* getParams(int argc, char* argv[], double *hx, double *hy, int *maxI);
double f(int x, int y);
double calcU(int n, double *u, double *fMem, double uDivisor, double hx, double hy, int nx, double coef1, double coef2, double coef3, double coef4);
double subsRow(int n, double *u, double uDivisor, double hx, double hy, int nx, double coef1, double coef2, double coef3, double coef4);
void sor(double *x, double *r, double *fMem, double *timeSor, double *timeResNorm, double w, double divided, double hx, double hy, int nx, int ny, int maxI);
