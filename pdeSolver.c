#include "header.h"

double timestamp(void) {
    struct timeval tp;
    gettimeofday(&tp, NULL);
    return((double)(tp.tv_sec + tp.tv_usec/1000000.0));
}

void getParams(int argc, char* argv[], double *hx, double *hy, int *maxI) {
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
			fclose(fp);
        } else {
            fprintf(stderr,"Incorrect parameter.\n");
            exit(-2);
        }
    }
}

double f(int n, double hx, double hy, int nx) {
/*
f(x,y) = 4π²[ sin(2πx)sinh(πy) + sin(2π(π−x))sinh(π(π−y)) ]
*/
	int i = n % nx, j = n / nx;
	double x = i * hx, y = j * hy;
	return (4*M_PI*M_PI * ( (sin(2*M_PI*x)) * (sinh(M_PI*y)) + (sin(2*M_PI*(M_PI-x))) * (sinh(M_PI*(M_PI-y))) ));
}

//double calcU(int n, double *u, double *f, double uDivisor, double hx, double hy, int nx) {
double calcU(int n, double *u, double uDivisor, double hx, double hy, int nx) {
/*
Final version of simplified equation:
u(i,j) = f(x,y) + (u(i+1,j) + u(i-1,j))/Δx² + (u(i,j+1) + u(i,j-1))/Δy² + (-u(i+1,j)+u(i-1,j))/2Δx + (-u(i,j+1)+u(i,j-1))/2Δy
         --------------------------------------------------------------------------------------------------------------------
                                                          2/Δx² + 2/Δy² + 4π²
*/
	double res = 0;
    //res += f[n] + (u[n+nx] + u[n-nx] ) / (hx * hx) + (u[n+1] + u[n-1]) / (hy * hy);
	res += f(n,hx,hy,nx) + (u[n+nx] + u[n-nx] ) / (hx * hx) + (u[n+1] + u[n-1]) / (hy * hy);
	res += (u[n-nx] - u[n+nx]) / (2 * hx) + (u[n-1] - u[n+1]) / (2 * hy);
	res = res / uDivisor;
	return res;
}

double subsRow(int n, double *u, double uDivisor, double hx, double hy, int nx) {
/*
f(x,y) =
(2/Δx²+2/Δy²+4π²)*u(i,j) - ( (u(i+1,j)+u(i-1,j))/Δx² + (u(i,j+1)+u(i,j-1))/Δy² + (-u(i+1,j)+u(i-1,j))/2Δx + (-u(i,j+1)+u(i,j-1))/2Δy) )
*/
	double res = 0;
	res = uDivisor * u[n];
	res -= ((u[n+nx] + u[n-nx]) / (hx * hx) + (u[n+1] + u[n-1]) / (hy * hy) + (u[n-nx] - u[n+nx]) / (2 * hx) + (u[n-1] - u[n+1]) / (2 * hy));
	return res;
}

//void sor(double *x, double *r, double *f, double *timeSor, double *timeResNorm, double w, double uDivisor, double hx, double hy, int nx, int ny, int maxI) {
void sor(double *x, double *r, double *timeSor, double *timeResNorm, double w, double uDivisor, double hx, double hy, int nx, int ny, int maxI) {
	int i, j, k;
	double sigma, now, fxy, res, maxRes = 0, tRes = 0; // maxRes is the biggest residue, tRes is total residue in this iteration.

	for(k=0; k<maxI; ++k) {
		now = timestamp(); // Starting iteration time counter.
		for(i = 1 + nx; i < nx * ny - nx - 1; ++i) { // Start at u[1][1], which means u[ny+1] and do not calculate borders
			//x[i] = x[i] + w * (calcU(i,x,f,uDivisor,hx,hy,nx) - x[i]);
            x[i] = x[i] + w * (calcU(i,x,uDivisor,hx,hy,nx) - x[i]);
		}
		*timeSor += timestamp() - now; // Get iteration time.
		now = timestamp(); // Start residue norm time counter.
		for(i = 1 + nx; i < nx * ny - nx - 1; ++i) {
            //res = f(i);
			res = f(i,hx,hy,nx); // res = f(x,y)
			res -= subsRow(i,x,uDivisor,hx,hy,nx); // res = f(x,y) - (a0 * x0 + a1 * x1 + ... + an * xn)
			if(res > maxRes)
				maxRes = res;
			tRes += res * res; // Adds res² to the total residue of this iteration.
		}
		r[k] = sqrt(tRes); // Store the norm of the residue in a vector (r).
		tRes = 0;
		*timeResNorm += timestamp() - now; // Get residue norm time.
	}

	*timeSor = *timeSor / maxI; // Get average values.
	*timeResNorm = *timeResNorm / maxI;
}

int main(int argc, char *argv[]) {
	int i, j, k, nx, ny, maxI, alpha;
	double hx, hy, w, beta, gama, sigma, uDivisor, *x, *r, *f, *timeSor, *timeResNorm;
	FILE *fpExit;

	getParams(argc,argv,&hx,&hy,&maxI);

	nx = (round(M_PI/hx)) + 1;
	ny = (round(M_PI/hy)) + 1;
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
	if((timeSor = calloc(1,sizeof(double))) == NULL) {
		fprintf(stderr,"Could not allocate memory.");
		exit(-5);
	}
	if((timeResNorm = calloc(1,sizeof(double))) == NULL) {
		fprintf(stderr,"Could not allocate memory.");
		exit(-5);
	}
	if((fpExit = fopen("solution.dat","w")) == NULL) {
		fprintf(stderr,"Could not open file.");
		exit(-6);
	}
    /*
    if((f = malloc(nx * ny * sizeof(double))) == NULL) {
        fprintf(stderr,"Could not allocate memory.");
        exit(-5);
    }
    for(i=1; i<ny-1; i++) { // Ignoring borders.
        for(j=1; j<nx-1; j++) { // Ignoring borders as well.
            f[i*nx+j] = f(i*nx+j); // Please check those indexes. (looking to lines 181-185, it looks OK).
        }
    }
    */



    sigma = sinh(M_PI * M_PI);
    alpha = nx * ny - nx;
    beta = 2 * M_PI * hx;
    gama = 2 * M_PI * M_PI;

	for(i = nx; i < alpha; ++i) { // Initialize central points (with left and right borders) as 0.
		x[i] = 0.0f;
	}

	for(i=0; i<nx; ++i) { // Creating borders
        x[i] = sin(gama - (i * hx)) * sigma;
        x[alpha+i] = sin(beta * i) * sigma;
	}

    //sor(x,r,f,timeSor,timeResNorm,w,uDivisor,hx,hy,nx,ny,maxI);
	sor(x,r,timeSor,timeResNorm,w,uDivisor,hx,hy,nx,ny,maxI);

	fprintf(fpExit,"###########\n# Tempo Método SOR: %lf\n# Tempo Resíduo: %lf\n\n# Norma do Resíduo\n",*timeSor,*timeResNorm);
	for(i=0; i<maxI; ++i) {
		fprintf(fpExit,"# i=%d: %.15lf\n",i,r[i]);
	}
	fprintf(fpExit,"###########\n");

	for(i = 1; i < ny - 1; ++i) { // This 'for' has to ignore borders.
		for(j = 1; j < nx - 1; ++j) {
			fprintf(fpExit,"%.15lf %.15lf %.15lf\n",j*hx,i*hy,x[i*nx+j]);
		}
	}

	fclose(fpExit);

	return 0;
}
