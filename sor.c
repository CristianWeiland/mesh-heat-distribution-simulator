#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#define SIZE 9

int N = SIZE;
int MaxI = 10;

inline int in(int i, int j) {
	return i*N + j;
}

void generate_matrix(double *A, double *b) {
/*	A[in(0,0)] = 10;
	A[in(0,1)] = 4;
	A[in(0,2)] = 2;
	A[in(1,0)] = -2;
	A[in(1,1)] = 15;
	A[in(1,2)] = 3;
	A[in(2,0)] = 2;
	A[in(2,1)] = -2;
	A[in(2,2)] = 13;
	b[0] = 5;
	b[1] = 3;
	b[2] = 4;*/
	int i,j;
	for(i=0; i<N*N; ++i) {
		A[i] = 0;
	}
	for(i=0; i<N; ++i) {
		for(j=0; j<N; ++j) {
			if(i == j)
				A[in(i,j)] = -4;
		}
	}
	A[in(0,1)] = 1;
	A[in(0,3)] = 1;
	A[in(1,0)] = 1;
	A[in(1,2)] = 1;
	A[in(1,4)] = 1;
	A[in(2,1)] = 1;
	A[in(2,5)] = 1;
	A[in(3,0)] = 1;
	A[in(3,4)] = 1;
	A[in(3,6)] = 1;
	A[in(4,1)] = 1;
	A[in(4,3)] = 1;
	A[in(4,5)] = 1;
	A[in(4,7)] = 1;
	A[in(5,2)] = 1;
	A[in(5,4)] = 1;
	A[in(5,8)] = 1;
	A[in(6,3)] = 1;
	A[in(6,7)] = 1;
	A[in(7,4)] = 1;
	A[in(7,6)] = 1;
	A[in(7,8)] = 1;
	A[in(8,5)] = 1;
	A[in(8,7)] = 1;
	b[0] = -100;
	b[1] = -20;
	b[2] = -20;
	b[3] = -80;
	b[4] = 0;
	b[5] = 0;
	b[6] = -260;
	b[7] = -180;
	b[8] = -180;
}

void print_vector(double *x) {
	int i;
	for(i=0; i<N; ++i) {
		printf("%f ", x[i]);
	}
	printf("\n");
}

void print_matrix(double *A) {
	int i,j;
	for(i=0; i<N; ++i) {
		for(j=0; j<N; ++j) {
			printf("%.0f ", A[in(i,j)]);
		}
		printf("\n");
	}
	printf("\n\n\n");
}

void sor(double *A, double *b, double *x, double w) {
	int i,j,k;
	double sigma;
	for(k=0; k<MaxI; k++) {
		for(i=0; i<N; ++i) {
			sigma = 0;
			for(j=0; j<N; ++j) {
				sigma = sigma + A[in(i,j)] * x[j];
			}
			sigma = sigma - A[in(i,i)] * x[i];
			x[i] = x[i] + w * ((b[i] - sigma)/A[in(i,i)] - x[i]);
		}
	}
}

int main(int argc, char* argv[]) {
	double *A,*b,*x,w;
	int i,j;
	A = malloc(N * N * sizeof(double));
	b = malloc(N * sizeof(double));
	x = malloc(N * sizeof(double));
	generate_matrix(A,b);
	for(i=0; i<N; ++i) {
		x[i] = 1;
	}
	print_vector(x);
	print_matrix(A);
	print_vector(b);
	w = 1.25;
	sor(A,b,x,w);

	print_vector(x);
}