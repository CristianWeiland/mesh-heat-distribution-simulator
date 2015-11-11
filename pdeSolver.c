#include "header.h"

/* Final version of simplified equation:
u(i,j) = f(x,y) + (u(i+1,j) + u(i-1,j))/Δx² + (u(i,j+1) + u(i,j-1))/Δy² + (-u(i+1,j)+u(i-1,j))/2Δx + (-u(i,j+1)+u(i,j-1))/2Δy
         --------------------------------------------------------------------------------------------------------------------
                                                          2/Δx² + 2/Δy² + 4π²

*/

void print_errno(void) {
    printf("%s",strerror(errno));
}

void getParams(int argc, char* argv[], FILE *fp) {
    if(argc != ARGS_NUM) {
        fprintf(stderr,"Wrong number of arguments.\n");
        exit(-1);
    }
    int i;
    for(i=1; i<ARGS_NUM; i+=2) {
        if(strcmp(argv[i],"-hx") == 0) {
            Hx = atof(argv[i+1]);
        } else if(strcmp(argv[i],"-hy") == 0) {
            Hy = atof(argv[i+1]);
        } else if(strcmp(argv[i],"-i") == 0) {
            MaxI = atoi(argv[i+1]);
        } else if(strcmp(argv[i],"-o") == 0) {
            fp = fopen(argv[i+1],"w");
        } else {
            fprintf(stderr,"Incorrect parameter.\n");
            exit(-1);
        }
    }
}

double f(double x, double y) { // If stored in memory, it will need Nx * Ny * 8 bytes of memory. Should we do it?
	return 4 * M_PI * ( ( sin(2 * M_PI * x) ) * ( sinh(M_PI * y) ) + ( sin(2 * Pipi - M_PI * x) ) * ( sinh(Pipi - M_PI * y) ));
}

inline int in(int i, int j) {
// Calcula o indice do vetor, como se fosse uma matriz.
	return i*Ny + j;
}

void calcU(int i, int j, double *u) {
	double res = 0;
	res += f(i,j) + (u[ in(i+1,j) ] + u[ in(i-1,j) ] ) / Hx * Hx + (u[ in(i,j+1) ] + u[ in(i,j-1) ]) / Hy * Hy;
	res += (u[ in(i-1,j) ] - u[ in(i+1,j) ]) / 2 * Hx + (u[ in(i,j-1) ] - u[ in(i,j+1) ]) / 2 * Hy;
	res = res / UDivisor;
	u[ in(i,j) ] = res;
}

double residueNorm(double *r) {
// u is only a row and r is the last residue
	double res = 0;
	int i;
	for(i=0; i<Ny; i++) {
		res += r[i];
	}
	return sqrt(res);
}

int main(int argc, char *argv[]) {
	int i, j, k;
	FILE *fpExit;

	getParams(argc,argv,fpExit);

	Nx = (1/Hx);
    Nx++;
	Ny = (1/Hy);
	Ny++;
	W = (2 - (Hx + Hy)) / 2;
	UDivisor = (2 / Hx * Hx) + (2 / Hy * Hy) + 4 * Pipi;

	double sigma;
	double *u = calloc(sizeof(double), Nx * Ny);

	sigma = sinh(M_PI * M_PI);
	for(i=0; i<Nx; i++) { // Calculate side extremities (top and bottom are always 0).
		u[ in(i,0) ] = sin(2 * M_PI * (M_PI - (i * Hx))) * sigma;
		u[ in(i,Nx) ] = sin(2 * M_PI * (i * Hx));
	}

	// Nx columns and Ny rows.

	/* // SOR
    for(k=0; k<MaxI; k++) {
		for(i=1; i<Ny; i++) {
			sigma = 0;
			for(j=1; j<Nx; j++) {
				sigma += a[i][j] * fi[i];
			}
			sigma -= a[i][i] * fi[i];
			//fi[i] = (1 - W) * fi[i] + W/a[i][i] * (b[i] - sigma);
			fi[i] = fi[i] + W * ((b[i] - sigma) / a[i][i] - fi[i])
		}
	}*/

	return 0;
}

/*
u(i,j) = f(x,y) + (u(i+1,j) + u(i-1,j))/Δx² + (u(i,j+1) + u(i,j-1))/Δy² + (-u(i+1,j)+u(i-1,j))/2Δx + (-u(i,j+1)+u(i,j-1))/2Δy
         --------------------------------------------------------------------------------------------------------------------
                                                          2/Δx² + 2/Δy² + 4π²


Formato:

Ax = B

A [a11 a12 a13] * X [x1] = B [b1] -> a11*x1 + a12*x2 + a13*x3 = b1 -> x1 = (b1 - a12*x2 - a13*x3)/a11
  [a21 a22 a23]     [x2]     [b2] -> a11*x1 + a12*x2 + a13*x3 = b2 -> x2 = (b2 - a11*x1 - a13*x3)/a12
  [a31 a32 a33]     [x3]     [b3] -> a11*x1 + a12*x2 + a13*x3 = b3 -> x3 = (b3 - a11*x1 - a12*x2)/a13

se x1 = u(1,1), sei calcular x1 da forma double u(int i, int j) {}
Mesmo vale pra x2 e x3. E valeria pra x4, x5, x6, etc.


Se u for uma matriz de 6 elementos. Ou seja, Hx = 0.25 e Hy = 0.3

u  u   u   u
u  x5  x6  u
u  x3  x4  u
u  x1  x2  u
u  u   u   u

5 linhas, 4 colunas.
Ny = 3, Nx = 2.

x1 = u11
x2 = u12
x3 = u13
x4 = u21
x5 = u22
x6 = u23

Matriz A, olhando pela equação u(i,j) obtida, seria algo assim:
A = [1        ]
    [  1      ]
    [    1    ]
    [      1  ]
    [        1]
*/