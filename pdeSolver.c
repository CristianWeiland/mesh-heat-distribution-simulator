#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* Final version of simplified equation:
u(i,j) = f(x,y) + (u(i+1,j) + u(i-1,j))/Δx² + (u(i,j+1) + u(i,j-1))/Δy² + (-u(i+1,j)+u(i-1,j))/2Δx + (-u(i,j+1)+u(i,j-1))/2Δy
         --------------------------------------------------------------------------------------------------------------------
                                                          2/Δx² + 2/Δy² + 4π²

*/

//#define MAX_SIZE 100000000 // 100 MB.
//int inMemory = 0;
int Nx, Ny, MaxI;
double Hx, Hy, W;
double UDivisor;
double Pipi = M_PI * M_PI; // Pipi = Pi * Pi

double f(double x, double y) { // If stored in memory, it will need Nx * Ny * 8 bytes of memory. Should we do it?
	return 4 * M_PI * ( ( sin(2 * M_PI * x) ) * ( sinh(M_PI * y) ) + ( sin(2 * Pipi - M_PI * x) ) * ( sinh(Pipi - M_PI * y) ));
}

void u(int i, int j, double **u, double **f) {
	double res = 0;
	res += f[i][j] + (u[i+1][j] + u[i-1][j]) / Hx * Hx + (u[i][j+1] + u[i][j-1]) / Hy * Hy;
	res += (u[i-1][j] - u[i+1][j]) / 2 * Hx + (u[i][j-1] - u[i][j+1]) / 2 * Hy;
	res = res / UDivisor;
	u[i][j] = res;
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

	fpExit = fopen(argv[5],"w");

	if(argc == 9) {
		if((strcmp(argv[1],"-hx") != 0) || (strcmp(argv[3],"-hy") != 0) || (strcmp(argv[5],"-i") != 0) || (strcmp(argv[7],"-o") != 0)) {
			puts("Please use the following format: './pdeSolver -hx <Hx> -hy <Hy> -i <maxIter> -o arquivo_saida.'");
			exit(-1);
		}
	} else {
		puts("Please use the following format: './pdeSolver -hx <Hx> -hy <Hy> -i <maxIter> -o arquivo_saida.'");
	}

	Hx = atof(argv[2]);
	Hy = atof(argv[4]);
	MaxI = atoi(argv[6]);

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
	double **u,**f;

	for(i=0; i<Nx; i++) {
		for(j=0; i<Ny; i++) {
			u[i][j] = calloc(sizeof(double));
			f[i][j] = malloc(sizeof(double));
			f[i][j] = f(i*Hx,j*Hy);
		}
	}
	sigma = sinh(M_PI * M_PI);
	for(i=0; i<Nx; i++) { // Calculate side extremities (top and bottom are always 0).
		u[i][0] = sin(2 * M_PI * (M_PI - (i * Hx))) * sigma;
		u[i][Nx] = sin(2 * M_PI * (i * Hx));
	}

	// Nx columns and Ny rows.

	for(k=0; k<MaxI; k++) {
		for(i=1; i<Ny; i++) {
			sigma = 0;
			for(j=1; j<Nx; j++) {
				/*if(j != i) {
					sigma += a[i][j] * fi[i];
				}*/
				sigma += a[i][j] * fi[i];
			}
			sigma -= a[i][i] * fi[i];
			//fi[i] = (1 - W) * fi[i] + W/a[i][i] * (b[i] - sigma);
			fi[i] = fi[i] + W * ((b[i] - sigma) / a[i][i] - fi[i])
		}
	}

	return 0;
}