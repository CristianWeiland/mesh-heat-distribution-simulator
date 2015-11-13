#include "header.h"

/* Final version of simplified equation:
u(i,j) = f(x,y) + (u(i+1,j) + u(i-1,j))/Δx² + (u(i,j+1) + u(i,j-1))/Δy² + (-u(i+1,j)+u(i-1,j))/2Δx + (-u(i,j+1)+u(i,j-1))/2Δy
         --------------------------------------------------------------------------------------------------------------------
                                                          2/Δx² + 2/Δy² + 4π²
Residue:
f(x,y) = (2/Δx² + 2/Δy² + 4π²) * u(i,j) - ((u(i+1,j) + u(i-1,j))/Δx² + (u(i,j+1) + u(i,j-1))/Δy² + (-u(i+1,j)+u(i-1,j))/2Δx + (-u(i,j+1)+u(i,j-1))/2Δy)
*/

void print_vector(double *x) {
	int i;
	for(i=0; i<Nx * Ny; ++i) {
		printf("%f ", x[i]);
	}
	printf("\n");
}

void print_matrix(double *A) {
	int i,j;
	for(i=0; i<Nx * Ny; ++i) {
		for(j=0; j<Nx * Ny; ++j) {
			printf("%.0f ", A[in(i,j)]);
		}
		printf("\n");
	}
	printf("\n\n\n");
}

double timestamp(void) {
    struct timeval tp;
    gettimeofday(&tp, NULL);
    return((double)(tp.tv_sec + tp.tv_usec/1000000.0));
}

void getParams(int argc, char* argv[], FILE **fp) {
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
            *fp = fopen(argv[i+1],"w");
        } else {
            fprintf(stderr,"Incorrect parameter.\n");
            exit(-1);
        }
    }
}

double f(int n) {
	int i = n / Nx, j = n % Nx;
	double x = i * Hx, y = j * Hy;
	return (4 * M_PI * ( ( sin(2 * M_PI * x) ) * ( sinh(M_PI * y) ) + ( sin(2 * Pipi - M_PI * x) ) * ( sinh(Pipi - M_PI * y) )));
}

inline int in(int i, int j) { // Calculate vector index, like its a matrix.
	return i*Ny + j;
}

double calcU(int n, double *u) {
/*
u(i,j) = f(x,y) + (u(i+1,j) + u(i-1,j))/Δx² + (u(i,j+1) + u(i,j-1))/Δy² + (-u(i+1,j)+u(i-1,j))/2Δx + (-u(i,j+1)+u(i,j-1))/2Δy
         --------------------------------------------------------------------------------------------------------------------
                                                          2/Δx² + 2/Δy² + 4π²
*/
	double res = 0;
	res += f(n) + (u[n+Ny] + u[n-Ny] ) / Hx * Hx + (u[n+1] + u[n-1]) / Hy * Hy;
	res += (u[n-Ny] - u[n+Ny]) / 2 * Hx + (u[n-1] - u[n+1]) / 2 * Hy;
	res = res / UDivisor;
	return res;
}

double subsRow(int n, double *u) {
//f(x,y) = (2/Δx² + 2/Δy² + 4π²) * u(i,j) - ((u(i+1,j) + u(i-1,j))/Δx² + (u(i,j+1) + u(i,j-1))/Δy² + (-u(i+1,j)+u(i-1,j))/2Δx + (-u(i,j+1)+u(i,j-1))/2Δy)
	double res = 0;
	res = UDivisor * u[n];
	res -= ((u[n+Ny] + u[n-Ny] ) / Hx * Hx + (u[n+1] + u[n-1]) / Hy * Hy + (u[n-Ny] - u[n+Ny]) / 2 * Hx + (u[n-1] - u[n+1]) / 2 * Hy);
	return res;
}

void sor(double *b, double *x, double *r, double *timeSor, double *timeResNorm) {
	int i,j,k;
	double sigma, now, fxy, res, maxRes = 0, tRes = 0; // maxRes is the biggest residue, tRes is total residue in this iteration.
	for(k=0; k<MaxI; k++) { // Iterate MaxI times.
		now = timestamp(); // Starting iteration time counter.
		for(i = 1 + Ny; i < Nx * Ny - Ny - 1; ++i) { // Start at u[1][1], which means u[Ny+1] and do not calculate last row/column (they are borders)
			x[i] = x[i] + W * (calcU(i,x) - x[i]);
		}
		*timeSor += timestamp() - now; // Get iteration time.
		now = timestamp(); // Start residue norm time counter.
		for(i = 1 + Ny; i < Nx * Ny - Ny - 1; ++i) {
			res = f(i); // res = f(x,y)
			res -= subsRow(i,x); // res = f(x,y) - (a0 * x0 + a1 * x1 + ... + an * xn)
			if(res > maxRes)
				maxRes = res;
			tRes += res * res; // Adds res² to the total residue of this iteration.
		}
		r[k] = sqrt(tRes); // Store the norm of the residue in a vector (r).
		tRes = 0;
		*timeResNorm += timestamp() - now; // Get residue norm time.
	}
	*timeSor = *timeSor / MaxI; // Get average values.
	*timeResNorm = *timeResNorm / MaxI;
}

int main(int argc, char *argv[]) {
	int i, j, k;
	FILE *fpExit;

	getParams(argc,argv,&fpExit);

	//Nx = (M_PI/Hx);
	//Ny = (M_PI/Hy);
	//Nx++;
	//Ny++;
	Nx = (round(M_PI/Hx)) + 1;
	Ny = (round(M_PI/Hy)) + 1;
	W = 2 - ((Hx + Hy) / 2);
	UDivisor = (2 / Hx * Hx) + (2 / Hy * Hy) + 4 * Pipi;

	double sigma;
	double *b, *x, *r, *timeSor, *timeResNorm;

	b = malloc(Nx * Ny * sizeof(double));
	x = malloc(Nx * Ny * sizeof(double));
	r = malloc(MaxI * sizeof(double));
	timeSor = calloc(1,sizeof(double));
	timeResNorm = calloc(1,sizeof(double));

	sigma = sinh(Pipi);

	for(i=0; i<Nx; i++) { // Creating borders
		x[ in(i,0) ] = sin(2 * M_PI * (M_PI - (i * Hx))) * sigma;
		x[ in(i,Nx) ] = sin(2 * M_PI * (i * Hx));
	}

	sor(b,x,r,timeSor,timeResNorm);

	printf("TimeSor: %lf\nTimeResNorm: %lf\n\nNorma do Resíduo\n",*timeSor,*timeResNorm);
	for(i=0; i<MaxI; ++i) {
		printf("# i=%d: %lf\n",i,r[i]);
	}

	fprintf(fpExit,"###########\n# Tempo Método SOR: %lf\n# Tempo Resíduo: %lf\n\n# Norma do Resíduo\n",*timeSor,*timeResNorm);
	for(i=0; i<MaxI; ++i) {
		fprintf(fpExit,"# i=%d: %lf\n",i,r[i]);
	}
	fprintf(fpExit,"###########\n");
	//FILE *expected = fopen("expected.txt","w");
	for(i=0; i<Ny; i++) { // Print bottom border.
		fprintf(fpExit, "%lf %lf %lf\n", 0.0f, i*Hy, subsRow(i, x));
		//fprintf(expected, "%lf %lf %lf\n", 0.0f, i*Hy, f(i));
	}
	for(i=1; i<Nx; i++) { // Since im using subsRow, I cant start at position 0.
		fprintf(fpExit,"%lf %lf %lf\n", i*Hx, 0.0f, subsRow(i*Ny,x)); // Border.
		//fprintf(expected,"%lf %lf %lf\n", i*Hx, 0.0f, f(i)); // Border.
		for(j=1; j<Ny; j++) {
			fprintf(fpExit,"%lf %lf %lf\n",i*Hx,j*Hy,subsRow(i*Ny+j, x)); // Should we change i*Hx with j*Hy? (Columns by rows?)
			//fprintf(expected,"%lf %lf %lf\n",i*Hx,j*Hy,f(i*Ny+j));
		}
	}
	for(i=0; i<Ny; i++) { // Print top border. We should check it.
		fprintf(fpExit, "%lf %lf %lf\n", Nx*Hx, i*Hy, subsRow(Nx*Ny - Ny + i,x));
		//fprintf(expected, "%lf %lf %lf\n", Nx*Hx, i*Hy, f(Nx*Ny - Ny + i));
	}
	// Nx columns and Ny rows.
	fclose(fpExit);
	//fclose(expected);

	return 0;
}
