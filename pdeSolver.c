#include "header.h"

/* Final version of simplified equation:
u(i,j) = f(x,y) + (u(i+1,j) + u(i-1,j))/Δx² + (u(i,j+1) + u(i,j-1))/Δy² + (-u(i+1,j)+u(i-1,j))/2Δx + (-u(i,j+1)+u(i,j-1))/2Δy
         --------------------------------------------------------------------------------------------------------------------
                                                          2/Δx² + 2/Δy² + 4π²
Resíduo:
f(x,y) = (2/Δx² + 2/Δy² + 4π²) * u(i,j) - ((u(i+1,j) + u(i-1,j))/Δx² + (u(i,j+1) + u(i,j-1))/Δy² + (-u(i+1,j)+u(i-1,j))/2Δx + (-u(i,j+1)+u(i,j-1))/2Δy)
*/

//#define MAX_SIZE 100000000 // 100 MB.
//int inMemory = 0;
int N = 9;

void print_vector(double *x) {
	int i;
	for(i=0; i<N; ++i) {
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
/*
double f(double x, double y) { // If stored in memory, it will need Nx * Ny * 8 bytes of memory. Should we do it?
	return 4 * M_PI * ( ( sin(2 * M_PI * x) ) * ( sinh(M_PI * y) ) + ( sin(2 * Pipi - M_PI * x) ) * ( sinh(Pipi - M_PI * y) ));
}
*/
double f(int n) { // If stored in memory, it will need Nx * Ny * 8 bytes of memory. Should we do it?
	int i = n / Nx, j = n % Nx;
	double x = i * Hx, y = j * Hy;
	return 4 * M_PI * ( ( sin(2 * M_PI * x) ) * ( sinh(M_PI * y) ) + ( sin(2 * Pipi - M_PI * x) ) * ( sinh(Pipi - M_PI * y) ));
}

inline int in(int i, int j) {
// Calcula o indice do vetor, como se fosse uma matriz.
	return i*Ny + j;
}
/*
void calcU(int i, int j, double *u) {
	double res = 0;
	res += f(i,j) + (u[ in(i+1,j) ] + u[ in(i-1,j) ] ) / Hx * Hx + (u[ in(i,j+1) ] + u[ in(i,j-1) ]) / Hy * Hy;
	res += (u[ in(i-1,j) ] - u[ in(i+1,j) ]) / 2 * Hx + (u[ in(i,j-1) ] - u[ in(i,j+1) ]) / 2 * Hy;
	res = res / UDivisor;
	u[ in(i,j) ] = res;
}
*/
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
	//u[n] = res;
	return res;
}

double calcRes(int n, double *u) {
//f(x,y) = (2/Δx² + 2/Δy² + 4π²) * u(i,j) - ((u(i+1,j) + u(i-1,j))/Δx² + (u(i,j+1) + u(i,j-1))/Δy² + (-u(i+1,j)+u(i-1,j))/2Δx + (-u(i,j+1)+u(i,j-1))/2Δy)
	double res = 0;
	//res = (2/Hx*Hx + 2/Hy*Hy + 4*M_PI*M_PI) * u[n];
	res = UDivisor * u[n];
	res -= ((u[n+Ny] + u[n-Ny] ) / Hx * Hx + (u[n+1] + u[n-1]) / Hy * Hy + (u[n-Ny] - u[n+Ny]) / 2 * Hx + (u[n-1] - u[n+1]) / 2 * Hy);
	return res;
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

/*
a11 * x1 + a12 * x2 + a13 * x3 = b1 -> x1 = (b1 - (a12 * x2 + a13 * x3)) / a11
x1 = u11;
x2 = u12;
x3 = u13;
x4 = u21;
x5 = u22;
x6 = u23;

Res[1] = b1 - (a11 * x1 + a12 * x2 + a13 * x3)
*/

void sor(double *b, double *x, double *r, double *timeSor, double *timeResNorm) {
	int i,j,k; // Ax = b -> ?u = f
	double sigma, now, fxy, res, maxRes = 0, tRes = 0; // maxRes is the biggest residue, tRes is total residue in this iteration.
	for(k=0; k<MaxI; k++) {
		now = timestamp();
		for(i=0; i<Nx * Ny; ++i) {
			//b[i] = f();
			/*sigma = 0;
			for(j=0; j<N; ++j) {
				sigma = sigma + A[in(i,j)] * x[j];
			}
			sigma = sigma - A[in(i,i)] * x[i];*/
			//                  b[1] - somatorio / a11
			//x[i] = x[i] + W * ((b[i] - sigma)/A[in(i,i)] - x[i]);
			x[i] = x[i] + W * (calcU(i,x) - x[i]);
			//res = b[i] - x[i];
		}
		*timeSor += timestamp() - now;
		now = timestamp();
		res = f(k);
		printf("Res = %lf - ",res);
		/*for(i = 0; i < Nx * Ny; ++i)
			res -= A[in(i,i)] * x[i];*/
		res -= calcRes(k,x);
		printf("Res = %lf\n",res);
		if(res > maxRes)
			maxRes = res;
		tRes += res * res;
		//residueNorm(r);
		r[k] = sqrt(tRes);
		tRes = 0;
		*timeResNorm += timestamp() - now;
	}
	*timeSor = *timeSor / MaxI;
	*timeResNorm = *timeResNorm / MaxI;
}

int main(int argc, char *argv[]) {
	int i, j, k;
	FILE *fpExit;

	getParams(argc,argv,fpExit);

	Nx = (M_PI/Hx);
	Ny = (M_PI/Hy);
	Nx++;
	Ny++;
	W = 2 - ((Hx + Hy) / 2);
	UDivisor = (2 / Hx * Hx) + (2 / Hy * Hy) + 4 * Pipi;

	double sigma;
	//double *A = calloc(sizeof(double), Nx * Ny * Nx * Ny);
	double *b, *x, *r, *timeSor, *timeResNorm;

	b = malloc(Nx * Ny * sizeof(double));
	x = malloc((Nx + 2) * (Ny + 2) * sizeof(double));
	r = malloc(MaxI * sizeof(double));
	timeSor = calloc(1,sizeof(double));
	timeResNorm = calloc(1,sizeof(double));

	sigma = sinh(M_PI * M_PI);
	/*for(i=0; i<Nx; ++i) { // Calculate side extremities (top and bottom are always 0, so calloc already took care of it).
		A[ in(i,0) ] = sin(2 * M_PI * (M_PI - (i * Hx))) * sigma;
		A[ in(i,Nx) ] = sin(2 * M_PI * (i * Hx));
	}*/

	// Creating borders
	for(i=0; i<Nx; i++) {
		x[ in(i,0) ] = sin(2 * M_PI * (M_PI - (i * Hx))) * sigma;
		x[ in(i,Nx+1) ] = sin(2 * M_PI * (i * Hx));
	}

	//for(i=0; i<Nx+2; i++)

	//print_matrix(A);

	sor(b,x,r,timeSor,timeResNorm);

	print_vector(x);

	printf("TimeSor: %lf\nTimeResNorm: %lf\n\nNorma do Resíduo\n",*timeSor,*timeResNorm);
	for(i=0; i<MaxI; ++i) {
		printf("i=%d: %lf\n",i,r[i]);
	}
	// Nx columns and Ny rows.

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
