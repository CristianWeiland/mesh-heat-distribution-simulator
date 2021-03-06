#include "header.h"

/*
Final version of simplified equation:
u(i,j) = f(x,y) + (u(i+1,j) + u(i-1,j))/Δx² + (u(i,j+1) + u(i,j-1))/Δy² + (-u(i+1,j)+u(i-1,j))/2Δx + (-u(i,j+1)+u(i,j-1))/2Δy
         --------------------------------------------------------------------------------------------------------------------
                                                          2/Δx² + 2/Δy² + 4π²
Residue:
f(x,y) =
(2/Δx² + 2/Δy² + 4π²)*u(i,j) - ( (u(i+1,j)+u(i-1,j))/Δx² + (u(i,j+1)+u(i,j-1))/Δy² + (-u(i+1,j)+u(i-1,j))/2Δx + (-u(i,j+1)+u(i,j-1))/2Δy )
*/

double timestamp(void) {
    struct timeval tp;
    gettimeofday(&tp, NULL);
    return((double)(tp.tv_sec + tp.tv_usec/1000000.0));
}

void getParams(int argc, char* argv[]) {
    if(argc != ARGS_NUM) {
        fprintf(stderr,"Wrong number of arguments.\n");
        exit(-1);
    }
    int i;
    FILE *fp;
    for(i=1; i<ARGS_NUM; i+=2) {
        if(strcmp(argv[i],"-hx") == 0) {
            Hx = atof(argv[i+1]);
        } else if(strcmp(argv[i],"-hy") == 0) {
            Hy = atof(argv[i+1]);
        } else if(strcmp(argv[i],"-i") == 0) {
            MaxI = atoi(argv[i+1]);
        } else if(strcmp(argv[i],"-o") == 0) {
            fp = fopen(argv[i+1],"w");
            fprintf(fp,"splot \"solution.dat\"\npause -1");
        } else {
            fprintf(stderr,"Incorrect parameter.\n");
            exit(-1);
        }
    }
}

double f(int n) {
	int i = n / Nx, j = n % Nx;
	double x = i * Hx, y = j * Hy;
	return (4 * M_PI * M_PI * ( (sin(2 * M_PI * x)) * (sinh(M_PI * y)) + (sin(2 * M_PI * (M_PI - x))) * (sinh(M_PI * (M_PI - y))) ));
}

double calcU(int n, double *u) {
/*
u(i,j) = f(x,y) + (u(i+1,j) + u(i-1,j))/Δx² + (u(i,j+1) + u(i,j-1))/Δy² + (-u(i+1,j)+u(i-1,j))/2Δx + (-u(i,j+1)+u(i,j-1))/2Δy
         --------------------------------------------------------------------------------------------------------------------
                                                          2/Δx² + 2/Δy² + 4π²
*/
	double res = 0;
	res += f(n) + (u[n+Nx] + u[n-Nx] ) / (Hx * Hx) + (u[n+1] + u[n-1]) / (Hy * Hy);
	res += (u[n-Nx] - u[n+Nx]) / (2 * Hx) + (u[n-1] - u[n+1]) / (2 * Hy);
	res = res / UDivisor;
	return res;
}

double subsRow(int n, double *u) {
//f(x,y) = (2/Δx² + 2/Δy² + 4π²) * u(i,j) - ((u(i+1,j) + u(i-1,j))/Δx² + (u(i,j+1) + u(i,j-1))/Δy² + (-u(i+1,j)+u(i-1,j))/2Δx + (-u(i,j+1)+u(i,j-1))/2Δy)
	double res = 0;
	res = UDivisor * u[n];
	res -= ((u[n+Nx] + u[n-Nx]) / (Hx * Hx) + (u[n+1] + u[n-1]) / (Hy * Hy) + (u[n-Nx] - u[n+Nx]) / (2 * Hx) + (u[n-1] - u[n+1]) / (2 * Hy));
	return res;
}

void sor(double *x, double *r, double *timeSor, double *timeResNorm) {
	int i,j,k;
	double sigma, now, fxy, res, maxRes = 0, tRes = 0; // maxRes is the biggest residue, tRes is total residue in this iteration.
	for(k=0; k<MaxI; k++) { // Iterate MaxI times.
		now = timestamp(); // Starting iteration time counter.
		for(i = 1 + Nx; i < Nx * Ny - Nx - 1; ++i) { // Start at u[1][1], which means u[Ny+1] and do not calculate last row/column (they are borders)
			x[i] = x[i] + W * (calcU(i,x) - x[i]);
		}
		*timeSor += timestamp() - now; // Get iteration time.
		now = timestamp(); // Start residue norm time counter.
		for(i = 1 + Nx; i < Nx * Ny - Nx - 1; ++i) {
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
	double sigma, *x, *r, *timeSor, *timeResNorm;
	FILE *fpExit;

	getParams(argc,argv);

	Nx = (round(M_PI/Hx)) + 1;
	Ny = (round(M_PI/Hy)) + 1;
	W = 2 - ((Hx + Hy) / 2);
	UDivisor = (2 / (Hx * Hx)) + (2 / (Hy * Hy)) + 4 * M_PI * M_PI;

	x = malloc(Nx * Ny * sizeof(double));
	r = malloc(MaxI * sizeof(double));
	timeSor = calloc(1,sizeof(double));
	timeResNorm = calloc(1,sizeof(double));
	fpExit = fopen("solution.dat","w");

	sigma = sinh(M_PI * M_PI);

	for(i = Nx; i < Nx*Ny - Nx; ++i) { // Initialize central points (with left and right borders) as 0.
		x[i] = 0.0f;
	}

	for(i=0; i<Nx; ++i) { // Creating borders
		x[i] = sin(2 * M_PI * (M_PI - (i * Hx))) * sigma;
		x[Nx*Ny-Nx+i] = sin(2 * M_PI * (i * Hx)) * sigma;
	}

	sor(x,r,timeSor,timeResNorm);

	fprintf(fpExit,"###########\n# Tempo Método SOR: %lf\n# Tempo Resíduo: %lf\n\n# Norma do Resíduo\n",*timeSor,*timeResNorm);
	for(i=0; i<MaxI; ++i) {
		fprintf(fpExit,"# i=%d: %lf\n",i,r[i]);
	}
	fprintf(fpExit,"###########\n");
	FILE *expected = fopen("expected.txt","w");
/*
	for(i = 0; i < Nx; ++i) { // Print bottom border.
		printf("1- Salvando no indice %d %d, ponto %f %f\n",i,0,i*Hx,0.0f);
		fprintf(fpExit, "%lf %lf %lf\n", i*Hx, 0.0f, x[i]);
		fprintf(expected, "%lf %lf %lf\n", i*Hx, 0.0f, f(i));
	}*/
	for(i = 1; i < Ny - 1; ++i) { // Since im using subsRow, I cant start at position 0.
		/*printf("2- Salvando no indice %d %d, ponto %f %f\n",0,i,0.0f,i*Hy);
		fprintf(fpExit,"%lf %lf %lf\n", 0.0f, i*Hy, x[i*Nx]); // Left border.
		fprintf(expected,"%lf %lf %lf\n", 0.0f, i*Hy, f(i*Nx)); // Left border.*/
		for(j = 1; j < Nx - 1; ++j) {
			//printf("3- Salvando no indice %d %d, ponto %f %f\n",j,i,j*Hx,i*Hy);
			fprintf(fpExit,"%.15lf %.15lf %.15lf\n",j*Hx,i*Hy,subsRow(i*Nx+j, x));
			fprintf(fpExit,"%.15lf %.15lf %.15lf\n",j*Hx,i*Hy,f(i*Nx+j));
		}
		/*printf("4- Salvando no indice %d %d, ponto %f %f\n",Ny,i,M_PI,i*Hy);
		fprintf(fpExit,"%lf %lf %lf\n", M_PI, i*Hy, x[i*Nx+Nx-1]); // Right border.
		fprintf(expected,"%lf %lf %lf\n", M_PI, i*Hy, f(i*Nx+Nx-1)); // Right border.*/
	}/*
	for(i = 0; i < Nx; ++i) { // Print top border. We should check it.
		printf("5- Salvando no indice %d %d, ponto %f %f\n",i,Ny,i*Hx,M_PI);
		fprintf(fpExit, "%lf %lf %lf\n", i*Hx, M_PI, x[Nx*Ny - Nx + i]); // Top border.
		fprintf(expected, "%lf %lf %lf\n", i*Hx, M_PI, f(Nx*Ny - Nx + i)); // Top border.
	}*/

	fclose(fpExit);
	fclose(expected);

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