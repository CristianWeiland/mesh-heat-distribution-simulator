#include "header.h"
/*
Banda de Memória: utilizar o grupo MEM do likwid, e apresentar o resultado de "Memory bandwidth [MBytes/s]";
likwid-perfctr -f -C 1 -g MEM -m ./pdeSolver -hx 0.05 -hy 0.05 -i 20 -o graph | grep "Memory bandwidth" | cut -c 39-50

Cache miss: utilizar o grupo CACHE do likwid, e apresentar o resultado de "data cache miss ratio";
likwid-perfctr -f -C 1 -g CACHE -m ./pdeSolver -hx 0.05 -hy 0.05 -i 20 -o graph

Operações aritméticas: utilizar o grupo FLOPS_DP do likwid, e apresentar o resultado de "DP MFLOP/s";
likwid-perfctr -f -C 1 -g FLOPS_DP -m ./pdeSolver -hx 0.05 -hy 0.05 -i 20 -o graph

---------

10 Pontos
likwid-perfctr -f -C 1 -g MEM -m ./pdeSolver -hx 0.34 -hy 0.34 -i 20 -o graph
likwid-perfctr -f -C 1 -g CACHE -m ./pdeSolver -hx 0.34 -hy 0.34 -i 20 -o graph
likwid-perfctr -f -C 1 -g FLOPS_DP -m ./pdeSolver -hx 0.34 -hy 0.34 -i 20 -o graph

100 Pontos
likwid-perfctr -f -C 1 -g MEM -m ./pdeSolver -hx 0.0317 -hy 0.0317 -i 20 -o graph
likwid-perfctr -f -C 1 -g CACHE -m ./pdeSolver -hx 0.0317 -hy 0.0317 -i 20 -o graph
likwid-perfctr -f -C 1 -g FLOPS_DP -m ./pdeSolver -hx 0.0317 -hy 0.0317 -i 20 -o graph

127 Pontos
likwid-perfctr -f -C 1 -g MEM -m ./pdeSolver -hx 0.025 -hy 0.025 -i 20 -o graph
likwid-perfctr -f -C 1 -g CACHE -m ./pdeSolver -hx 0.025 -hy 0.025 -i 20 -o graph
likwid-perfctr -f -C 1 -g FLOPS_DP -m ./pdeSolver -hx 0.025 -hy 0.025 -i 20 -o graph

128 Pontos
likwid-perfctr -f -C 1 -g MEM -m ./pdeSolver -hx 0.0248 -hy 0.0248 -i 20 -o graph
likwid-perfctr -f -C 1 -g CACHE -m ./pdeSolver -hx 0.0248 -hy 0.0248 -i 20 -o graph
likwid-perfctr -f -C 1 -g FLOPS_DP -m ./pdeSolver -hx 0.0248 -hy 0.0248 -i 20 -o graph

200 Pontos
likwid-perfctr -f -C 1 -g MEM -m ./pdeSolver -hx 0.0158 -hy 0.0158 -i 20 -o graph
likwid-perfctr -f -C 1 -g CACHE -m ./pdeSolver -hx 0.0158 -hy 0.0158 -i 20 -o graph
likwid-perfctr -f -C 1 -g FLOPS_DP -m ./pdeSolver -hx 0.0158 -hy 0.0158 -i 20 -o graph

255 Pontos
likwid-perfctr -f -C 1 -g MEM -m ./pdeSolver -hx 0.01235 -hy 0.01235 -i 20 -o graph
likwid-perfctr -f -C 1 -g CACHE -m ./pdeSolver -hx 0.01235 -hy 0.01235 -i 20 -o graph
likwid-perfctr -f -C 1 -g FLOPS_DP -m ./pdeSolver -hx 0.01235 -hy 0.01235 -i 20 -o graph

256 Pontos
likwid-perfctr -f -C 1 -g MEM -m ./pdeSolver -hx 0.0123 -hy 0.0123 -i 20 -o graph
likwid-perfctr -f -C 1 -g CACHE -m ./pdeSolver -hx 0.0123 -hy 0.0123 -i 20 -o graph
likwid-perfctr -f -C 1 -g FLOPS_DP -m ./pdeSolver -hx 0.0123 -hy 0.0123 -i 20 -o graph

500 Pontos
likwid-perfctr -f -C 1 -g MEM -m ./pdeSolver -hx 0.0063 -hy 0.0063 -i 20 -o graph
likwid-perfctr -f -C 1 -g CACHE -m ./pdeSolver -hx 0.0063 -hy 0.0063 -i 20 -o graph
likwid-perfctr -f -C 1 -g FLOPS_DP -m ./pdeSolver -hx 0.0063 -hy 0.0063 -i 20 -o graph

511 Pontos
likwid-perfctr -f -C 1 -g MEM -m ./pdeSolver -hx 0.00616 -hy 0.00616 -i 20 -o graph
likwid-perfctr -f -C 1 -g CACHE -m ./pdeSolver -hx 0.00616 -hy 0.00616 -i 20 -o graph
likwid-perfctr -f -C 1 -g FLOPS_DP -m ./pdeSolver -hx 0.00616 -hy 0.00616 -i 20 -o graph

512 Pontos
likwid-perfctr -f -C 1 -g MEM -m ./pdeSolver -hx 0.00615 -hy 0.00615 -i 20 -o graph
likwid-perfctr -f -C 1 -g CACHE -m ./pdeSolver -hx 0.00615 -hy 0.00615 -i 20 -o graph
likwid-perfctr -f -C 1 -g FLOPS_DP -m ./pdeSolver -hx 0.00615 -hy 0.00615 -i 20 -o graph

1000 Pontos
likwid-perfctr -f -C 1 -g MEM -m ./pdeSolver -hx 0.003145 -hy 0.003145 -i 20 -o graph
likwid-perfctr -f -C 1 -g CACHE -m ./pdeSolver -hx 0.003145 -hy 0.003145 -i 20 -o graph
likwid-perfctr -f -C 1 -g FLOPS_DP -m ./pdeSolver -hx 0.003145 -hy 0.003145 -i 20 -o graph

1023 Pontos
likwid-perfctr -f -C 1 -g MEM -m ./pdeSolver -hx 0.003073 -hy 0.003073 -i 20 -o graph
likwid-perfctr -f -C 1 -g CACHE -m ./pdeSolver -hx 0.003073 -hy 0.003073 -i 20 -o graph
likwid-perfctr -f -C 1 -g FLOPS_DP -m ./pdeSolver -hx 0.003073 -hy 0.003073 -i 20 -o graph

1024 Pontos
likwid-perfctr -f -C 1 -g MEM -m ./pdeSolver -hx 0.00307 -hy 0.00307 -i 20 -o graph
likwid-perfctr -f -C 1 -g CACHE -m ./pdeSolver -hx 0.00307 -hy 0.00307 -i 20 -o graph
likwid-perfctr -f -C 1 -g FLOPS_DP -m ./pdeSolver -hx 0.00307 -hy 0.00307 -i 20 -o graph

2000 Pontos
likwid-perfctr -f -C 1 -g MEM -m ./pdeSolver -hx 0.0015715 -hy 0.0015715 -i 20 -o graph
likwid-perfctr -f -C 1 -g CACHE -m ./pdeSolver -hx 0.0015715 -hy 0.0015715 -i 20 -o graph
likwid-perfctr -f -C 1 -g FLOPS_DP -m ./pdeSolver -hx 0.0015715 -hy 0.0015715 -i 20 -o graph

2047 Pontos
likwid-perfctr -f -C 1 -g MEM -m ./pdeSolver -hx 0.0015355 -hy 0.0015355 -i 20 -o graph
likwid-perfctr -f -C 1 -g CACHE -m ./pdeSolver -hx 0.0015355 -hy 0.0015355 -i 20 -o graph
likwid-perfctr -f -C 1 -g FLOPS_DP -m ./pdeSolver -hx 0.0015355 -hy 0.0015355 -i 20 -o graph

2048 Pontos
likwid-perfctr -f -C 1 -g MEM -m ./pdeSolver -hx 0.001535 -hy 0.001535 -i 20 -o graph
likwid-perfctr -f -C 1 -g CACHE -m ./pdeSolver -hx 0.001535 -hy 0.001535 -i 20 -o graph
likwid-perfctr -f -C 1 -g FLOPS_DP -m ./pdeSolver -hx 0.001535 -hy 0.001535 -i 20 -o graph
*/

double timestamp(void) {
    struct timeval tp;
    gettimeofday(&tp, NULL);
    return((double)(tp.tv_sec + tp.tv_usec/1000000.0));
}

FILE* getParams(int argc, char* argv[], double *hx, double *hy, int *maxI) {
    if(argc != ARGS_NUM) {
        fprintf(stderr,"Wrong number of arguments.\n");
        exit(-1);
    }

    int i;
    FILE *fp;

    for(i=1; i<ARGS_NUM; i+=2) {
        if(strcmp(argv[i],"-hx") == 0) {
            *hx = atof(argv[i+1]);
            if(*hx <= 0.0f) {
            	fprintf(stderr,"Distance should not be less than 0.\n");
            	exit(-4);
            }
        } else if(strcmp(argv[i],"-hy") == 0) {
            *hy = atof(argv[i+1]);
            if(*hy <= 0.0f) {
            	fprintf(stderr,"Distance should not be less than 0.\n");
            	exit(-4);
            }
        } else if(strcmp(argv[i],"-i") == 0) {
            *maxI = atoi(argv[i+1]);
            if(*maxI < 1) {
            	fprintf(stderr,"You need to do at least 1 iteration.\n");
            	exit(-3);
            }
        } else if(strcmp(argv[i],"-o") == 0) {
			if((fp = fopen(argv[i+1],"w")) == NULL) {
				fprintf(stderr,"Could not open file.");
				exit(-6);
			}
            fprintf(fp,"splot \"solution.dat\"\npause -1");
        } else {
            fprintf(stderr,"Incorrect parameter.\n");
            exit(-2);
        }
    }

    return fp;
}

inline double f(int x, int y) {
/* f(x,y) = 4π²[ sin(2πx)sinh(πy) + sin(2π(π−x))sinh(π(π−y)) ] */
	//double x = j * hx, y = i * hy;
	return (4*M_PI*M_PI * ( (sin(2*M_PI*x)) * (sinh(M_PI*y)) + (sin(2*M_PI*(M_PI-x))) * (sinh(M_PI*(M_PI-y))) ));
}

inline double calcU(int n, double *u, double *fMem, double divided, double hx, double hy, int nx, double coef1, double coef2, double coef3, double coef4) {
/*
Final version of simplified equation:
                            coef1                    coef2                      coef3                    coef4
u(i,j) = f(x,y) + u(i+1,j)*(1/Δx²-1/2Δx) + u(i-1,j)*(1/Δx²+1/2Δx) + u(i,j+1)*1/(1/Δy²-1/2Δy) + u(i,j-1)*(1/Δy²+1/2Δy)
         --------------------------------------------------------------------------------------------------------------
                                                        2/Δx² + 2/Δy² + 4π²
*/
	return ((fMem[n] + u[n+nx] * coef1 + u[n-nx] * coef2 + u[n+1] * coef3 + u[n-1] * coef4) * divided);
}

inline double subsRow(int n, double *u, double uDivisor, double hx, double hy, int nx, double coef1, double coef2, double coef3, double coef4) {
/*
f(x,y) =
(2/Δx²+2/Δy²+4π²)*u(i,j) - ( (u(i+1,j)+u(i-1,j))/Δx² + (u(i,j+1)+u(i,j-1))/Δy² + (-u(i+1,j)+u(i-1,j))/2Δx + (-u(i,j+1)+u(i,j-1))/2Δy) )
*/
/*
2/Δx²+2/Δy²+4π² * u(i,j) = f(x,y) + u(i+1,j) * 1/(Δx(Δx-2)) + u(i-1,j) * 1/(Δx(Δx+2)) + u(i,j+1) * 1/(Δy(Δy-2)) + u(i,j-1) * 1/(Δy(Δy+2))
f(x,y) = 2/Δx²+2/Δy²+4π² * u(i,j) - (u(i+1,j) * 1/(Δx(Δx-2)) + u(i-1,j) * 1/(Δx(Δx+2)) + u(i,j+1) * 1/(Δy(Δy-2)) + u(i,j-1) * 1/(Δy(Δy+2)))
*/
    return (uDivisor * u[n] + (-1 * (u[n+nx]*coef1 + u[n-nx]*coef2 + u[n+1]*coef3 + u[n-1]*coef4)));
}

void sor(double *x, double *r, double *fMem, double *timeSor, double *timeResNorm, double w, double uDivisor, double hx, double hy, int nx, int ny, int maxI) {
	int i, j, k, l, m, row, col, index, inx;
	double now, res, tRes, maxRes = 0, divided; // tRes is total residue in this iteration, maxRes is the biggest residue.
    double coef1, coef2, coef3, coef4;

    coef1 = (1/(hx*hx)) - (1/(2*hx)); // u(i+1,j)
    coef2 = (1/(hx*hx)) + (1/(2*hx)); // u(i-1,j)
    coef3 = (1/(hy*hy)) - (1/(2*hy)); // u(i,j+1)
    coef4 = (1/(hy*hy)) + (1/(2*hy)); // u(i,j-1)
    divided = 1 / uDivisor;

	for(k=0; k<maxI; ++k) {
		now = timestamp(); // Starting iteration time counter.
		LIKWID_MARKER_START("sor");

        for(i=1; i<ny-1; i+=BLOCK_SIZE) {
            inx = i * nx;
            for(j=1; j<nx-1; j+=BLOCK_SIZE) {
                for(l=0; l<BLOCK_SIZE && (l+i)<ny-1; ++l) {
                    row = inx + l*nx;
                    for(m=0; m<BLOCK_SIZE && (m+j)<nx-1; ++m) {
                        index = row+j+m;
                        x[index] = x[index] + w * (calcU(index,x,fMem,divided,hx,hy,nx,coef1,coef2,coef3,coef4) - x[index]);
                    }
                }
            }
        }

        LIKWID_MARKER_STOP("sor");
		*timeSor += timestamp() - now; // Get iteration time.
		now = timestamp(); // Start residue norm time counter.
		LIKWID_MARKER_START("residue");

		tRes = 0.0f;

	    for(i=1; i<ny-1; ++i) { // Ignoring borders.
            inx = i * nx;
	        for(j=1; j<nx-1; ++j) { // Ignoring borders as well.
	            res = fMem[inx+j] - subsRow(inx+j,x,uDivisor,hx,hy,nx,coef1,coef2,coef3,coef4);
				if(res > maxRes)
					maxRes = res;
				tRes = tRes + res * res; // Adds res² to the total residue of this iteration.
	        }
	    }

		r[k] = sqrt(tRes); // Store the norm of the residue in a vector (r).

        LIKWID_MARKER_STOP("residue");
		*timeResNorm += timestamp() - now; // Get residue norm time.
	}

	*timeSor = *timeSor / maxI; // Get average values.
	*timeResNorm = *timeResNorm / maxI;
}

int main(int argc, char *argv[]) {
	int i, j, nx, ny, maxI, alpha, inx;
	double hx, hy, w, beta, sigma, uDivisor, *x, *r, *fMem, timeSor, timeResNorm, y;
	FILE *fpExit, *fpData;

	fpExit = getParams(argc,argv,&hx,&hy,&maxI);

    LIKWID_MARKER_INIT;

	nx = (round(M_PI/hx)) + 1;
	ny = (round(M_PI/hy)) + 1;
	printf("Nx = %d, Ny = %d\n",nx,ny);
	w = 2 - ((hx + hy) / 2);
	uDivisor = (2 / (hx * hx)) + (2 / (hy * hy)) + 4 * M_PI * M_PI;

	if((x = malloc(nx * ny * sizeof(double))) == NULL) {
		fprintf(stderr,"Could not allocate memory.");
		exit(-5);
	}
	if((r = malloc(maxI * sizeof(double))) == NULL) {
		fprintf(stderr,"Could not allocate memory.");
		exit(-5);
	}

    if((fMem = malloc(nx * ny * sizeof(double))) == NULL) {
        fprintf(stderr,"Could not allocate memory.");
        exit(-5);
    }

	if((fpData = fopen("solution.dat","w")) == NULL) {
		fprintf(stderr,"Could not open file.");
		exit(-6);
	}

    timeSor = 0.0f;
    timeResNorm = 0.0f;

	for(i = 1; i < ny - 1; ++i) { // This 'for' has to ignore borders.
        inx = i*nx;
		for(j = 0; j < nx; ++j) { // This 'for' cant ignore borders.
			x[inx+j] = 0.0f;
		}
	}

    sigma = sinh(M_PI * M_PI);
    alpha = nx * ny - nx;
    beta = 2 * M_PI * hx;

	for(i=0; i<nx; ++i) { // Creating borders
        x[i] = sin(2 * M_PI * (M_PI - (i * hx))) * sigma;
        x[alpha+i] = sin(beta * i) * sigma;
	}

	// Initializing f(x,y)
    for(i=1; i<ny-1; ++i) { // Ignoring borders.
        inx = i * nx;
        y = i * hy;
        for(j=1; j<nx-1; ++j) { // Ignoring borders as well.
            fMem[inx+j] = f(j*hx,y);
        }
    }

    sor(x,r,fMem,&timeSor,&timeResNorm,w,uDivisor,hx,hy,nx,ny,maxI);

	fprintf(fpExit,"\n\n\n###########\n# Tempo Método SOR: %lf\n# Tempo Resíduo: %lf\n\n# Norma do Resíduo\n",timeSor,timeResNorm);

	for(i=0; i<maxI; ++i) {
		fprintf(fpExit,"# i=%d: %.15lf\n",i,r[i]);
	}
	fprintf(fpExit,"###########\n");

	for(i = 1; i < ny - 1; ++i) { // This 'for' has to ignore borders.
        beta = i * hy;
        inx = i * nx;
		for(j = 1; j < nx - 1; ++j) {
			fprintf(fpData,"%.15lf %.15lf %.15lf\n",j*hx,beta,x[inx+j]);
		}
	}

	fclose(fpExit);
	fclose(fpData);

    LIKWID_MARKER_CLOSE;

	return 0;
}

/*
6 A30 A31 A32 A33 A34
5 A25 A26 A27 A28 A29
4 A20 A21 A22 A23 A24
3 A15 A16 A17 A18 A19
2 A10 A11 A12 A13 A14
1  A5  A6  A7  A8  A9
0  A0  A1  A2  A3  A4
    0   1   2   3   4

Nx = 5
Ny = 7

*/
//f(x,y) + (u(i+1,j) + u(i-1,j))/Δx² + (u(i,j+1) + u(i,j-1))/Δy² + (-u(i+1,j)+u(i-1,j))/2Δx + (-u(i,j+1)+u(i,j-1))/2Δy
/*
Final version of simplified equation:
u(i,j) = f(x,y) + (u(i+1,j) + u(i-1,j))/Δx² + (u(i,j+1) + u(i,j-1))/Δy² + (-u(i+1,j)+u(i-1,j))/2Δx + (-u(i,j+1)+u(i,j-1))/2Δy
         --------------------------------------------------------------------------------------------------------------------
                                                          2/Δx² + 2/Δy² + 4π²
Juntando u(i+1,j) com u(i+1,j), etc.

(2/Δx² + 2/Δy² + 4π²) * u(i,j) - f(x,y) =

u(i+1,j)/Δx² + u(i-1,j)/Δx² + u(i,j+1)/Δy² + u(i,j-1)/Δy² -u(i+1,j)/2Δx + u(i-1,j)/2Δx -u(i,j+1)2Δy + u(i,j-1)/2Δy =

u(i+1,j)/Δx² -u(i+1,j)/2Δx + u(i-1,j)/Δx² + u(i-1,j)/2Δx + u(i,j+1)/Δy² -u(i,j+1)2Δy + u(i,j-1)/Δy² + u(i,j-1)/2Δy =

u(i+1,j) * (1/Δx² - 1/2Δx)

u(i+1,j) * 1/(Δx²-2Δx) + u(i-1,j) * 1/(Δx²+2Δx) + u(i,j+1) * 1/(Δy²-2Δy) + u(i,j-1) * 1/(Δy²+2Δy)
           A                        B                        C                        D

u(i+1,j) * 1/(Δx(Δx-2)) + u(i-1,j) * 1/(Δx(Δx+2)) + u(i,j+1) * 1/(Δy(Δy-2)) + u(i,j-1) * 1/(Δy(Δy+2))
           A                        B                        C                        D

A = 1/(hx * (hx - 2))
B = 1/(hx * (hx + 2))
C = 1/(hy * (hy - 2))
D = 1/(hy * (hy + 2))
*/
