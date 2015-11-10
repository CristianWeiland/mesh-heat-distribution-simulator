#include "header.h"

/* Final version of simplified equation:
u(i,j) = f(x,y) + (u(i+1,j) + u(i-1,j))/Δx² + (u(i,j+1) + u(i,j-1))/Δy² + (-u(i+1,j)+u(i-1,j))/2Δx + (-u(i,j+1)+u(i,j-1))/2Δy
         --------------------------------------------------------------------------------------------------------------------
                                                          2/Δx² + 2/Δy² + 4π²

*/

//#define MAX_SIZE 100000000 // 100 MB.
//int inMemory = 0;

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
	Ny = (1/Hy);
	Nx++;
	Ny++;
	if(Nx != Ny) {
		puts("Not a square matrix. Result will be undetermined or impossible.");
		exit(-1);
	}
	W = (2 - (Hx + Hy)) / 2;
	if(W < 1) {
		puts("Relaxation factor less than 1 - program may not work as you expect.");
	}
	UDivisor = (2 / Hx * Hx) + (2 / Hy * Hy) + 4 * Pipi;
	/*if((Nx * Ny * 8) <= MAX_SIZE) { // Should we save a matrix with f(x,y) values?
		inMemory = 1;
	}*/

	double sigma;
	//double **u,**f;

	/*for(i=0; i<Nx; i++) {
		for(j=0; i<Ny; i++) {
			u[i][j] = calloc(sizeof(double));
			f[i][j] = malloc(sizeof(double));
			f[i][j] = f(i*Hx,j*Hy);
		}
	}*/

	double *u = calloc(sizeof(double), Nx * Ny);

	sigma = sinh(M_PI * M_PI);
	for(i=0; i<Nx; i++) { // Calculate side extremities (top and bottom are always 0).
		u[ in(i,0) ] = sin(2 * M_PI * (M_PI - (i * Hx))) * sigma;
		u[ in(i,Nx) ] = sin(2 * M_PI * (i * Hx));
	}

	// Nx columns and Ny rows.

	/*for(k=0; k<MaxI; k++) {
		for(i=1; i<Ny; i++) {
			sigma = 0;
			for(j=1; j<Nx; j++) {
				/*if(j != i) {
					sigma += a[i][j] * fi[i];
				}*/
				/*sigma += a[i][j] * fi[i];
			}
			sigma -= a[i][i] * fi[i];
			//fi[i] = (1 - W) * fi[i] + W/a[i][i] * (b[i] - sigma);
			fi[i] = fi[i] + W * ((b[i] - sigma) / a[i][i] - fi[i])
		}
	}*/

	return 0;
}