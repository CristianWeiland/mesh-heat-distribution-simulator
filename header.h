#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <time.h>

#define ARGS_NUM 9

double timestamp(void);
void getParams(int argc, char* argv[], double *hx, double *hy, int *maxI);
double f(int n, double hx, double hy, int nx);
//double calcU(int n, double *u, double *f, double uDivisor, double hx, double hy, int nx) {
double calcU(int n, double *u, double uDivisor, double hx, double hy, int nx);
double subsRow(int n, double *u, double uDivisor, double hx, double hy, int nx);
//void sor(double *x, double *r, double *f, double *timeSor, double *timeResNorm, double w, double uDivisor, double hx, double hy, int nx, int ny, int maxI) {
void sor(double *x, double *r, double *timeSor, double *timeResNorm, double w, double uDivisor, double hx, double hy, int nx, int ny, int maxI);
