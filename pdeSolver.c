#include "header.h"

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

double f(int i, int j, double hx, double hy) {
/* f(x,y) = 4π²[ sin(2πx)sinh(πy) + sin(2π(π−x))sinh(π(π−y)) ] */
	double x = j * hx, y = i * hy;
	return (4*M_PI*M_PI * ( (sin(2*M_PI*x)) * (sinh(M_PI*y)) + (sin(2*M_PI*(M_PI-x))) * (sinh(M_PI*(M_PI-y))) ));
}

inline double calcU(int n, double *u, double *fMem, double uDivisor, double hx, double hy, int nx, double coef1, double coef2, double coef3, double coef4) {
/*
Final version of simplified equation:
                            coef1                    coef2                      coef3                    coef4
u(i,j) = f(x,y) + u(i+1,j)*(1/Δx²-1/2Δx) + u(i-1,j)*(1/Δx²+1/2Δx) + u(i,j+1)*1/(1/Δy²-1/2Δy) + u(i,j-1)*(1/Δy²+1/2Δy)
         --------------------------------------------------------------------------------------------------------------
                                                        2/Δx² + 2/Δy² + 4π²
*/
	return ((fMem[n] + u[n+nx] * coef1 + u[n-nx] * coef2 + u[n+1] * coef3 + u[n-1] * coef4) / uDivisor);
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

void sor(double *x, double *r, double *fMem, double *timeSor, double *timeResNorm, double w, double uDivisor, double hx, double hy, int nx, int ny, int maxI, int e) {
	int i, j, k, l, m, row, col;
	double now, res, tRes, maxRes = 0; // tRes is total residue in this iteration, maxRes is the biggest residue.
    double coef1, coef2, coef3, coef4;

    coef1 = (1/(hx*hx)) - (1/(2*hx)); // u(i+1,j)
    coef2 = (1/(hx*hx)) + (1/(2*hx)); // u(i-1,j)
    coef3 = (1/(hy*hy)) - (1/(2*hy)); // u(i,j+1)
    coef4 = (1/(hy*hy)) + (1/(2*hy)); // u(i,j-1)

	for(k=0; k<maxI; ++k) {
		now = timestamp(); // Starting iteration time counter.

        for(i=1; i<ny-1; i+=BLOCK_SIZE) {
            for(j=1; j<nx-1; j+=BLOCK_SIZE) {
                for(l=0; l<BLOCK_SIZE && l+i<ny-1; ++l) {
                    for(m=0; m<BLOCK_SIZE && m+j<nx-1; ++m) {
                        row = i+l;
                        col = j+m;
                        x[row*(nx+e)+col] = x[row*(nx+e)+col] + w * (calcU(row*(nx+e)+col,x,fMem,uDivisor,hx,hy,nx+e,coef1,coef2,coef3,coef4) - x[row*(nx+e)+col]);
                    }
                }
            }
        }

		*timeSor += timestamp() - now; // Get iteration time.
		now = timestamp(); // Start residue norm time counter.

		tRes = 0.0f;

	    for(i=1; i<ny-1; ++i) { // Ignoring borders.
	        for(j=1; j<nx-1; ++j) { // Ignoring borders as well.
	            res = fMem[i*(nx+e)+j] - subsRow(i*(nx+e)+j,x,uDivisor,hx,hy,nx+e,coef1,coef2,coef3,coef4);
				if(res > maxRes)
					maxRes = res;
				tRes += res * res; // Adds res² to the total residue of this iteration.
	        }
	    }

		r[k] = sqrt(tRes); // Store the norm of the residue in a vector (r).

		*timeResNorm += timestamp() - now; // Get residue norm time.
	}

	*timeSor = *timeSor / maxI; // Get average values.
	*timeResNorm = *timeResNorm / maxI;
}

int main(int argc, char *argv[]) {
	int e, i, j, nx, ny, maxI, alpha, n, a;
	double hx, hy, w, beta, sigma, uDivisor, *x, *r, *fMem, timeSor, timeResNorm;
	FILE *fpExit, *fpData;

	fpExit = getParams(argc,argv,&hx,&hy,&maxI);

	nx = (round(M_PI/hx)) + 1;
	ny = (round(M_PI/hy)) + 1;
	w = 2 - ((hx + hy) / 2);
	uDivisor = (2 / (hx * hx)) + (2 / (hy * hy)) + 4 * M_PI * M_PI;

    for(e = 1, n = nx, a = 1, i = 0; i < 32; ++i){
        if((n & a != n) && (n & a != 0))
            e = 0;
        a = a << 1;
    }
/*
    while(pot > 1 && ((pot % 2) == 0))
        pot = pot / 2;
    e = (pot == 1) ? 1 : 0; // Trying to avoid cache trashing.
*/
    printf("e = %d\n");

	if((x = malloc((nx + e) * ny * sizeof(double))) == NULL) {
		fprintf(stderr,"Could not allocate memory.");
		exit(-5);
	}
	if((r = malloc(maxI * sizeof(double))) == NULL) {
		fprintf(stderr,"Could not allocate memory.");
		exit(-5);
	}

    if((fMem = malloc((nx + e) * ny * sizeof(double))) == NULL) {
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
		for(j = 0; j < nx; ++j) { // This 'for' cant ignore borders.
			x[i*(nx+e)+j] = 0.0f;
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
        for(j=1; j<nx-1; ++j) { // Ignoring borders as well.
            fMem[i*(nx+e)+j] = f(i,j,hx,hy);
        }
    }

    sor(x,r,fMem,&timeSor,&timeResNorm,w,uDivisor,hx,hy,nx,ny,maxI,e);

	fprintf(fpExit,"\n\n\n###########\n# Tempo Método SOR: %lf\n# Tempo Resíduo: %lf\n\n# Norma do Resíduo\n",timeSor,timeResNorm);

	for(i=0; i<maxI; ++i) {
		fprintf(fpExit,"# i=%d: %.15lf\n",i,r[i]);
	}
	fprintf(fpExit,"###########\n");

	for(i = 1; i < ny - 1; ++i) { // This 'for' has to ignore borders.
		for(j = 1; j < nx - 1; ++j) {
			fprintf(fpData,"%.15lf %.15lf %.15lf\n",j*hx,i*hy,x[i*(nx+e)+j]);
		}
	}

	fclose(fpExit);
	fclose(fpData);

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
