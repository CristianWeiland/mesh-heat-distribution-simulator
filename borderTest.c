#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 100000

double timestamp(void) {
    struct timeval tp;
    gettimeofday(&tp, NULL);
    return((double)(tp.tv_sec + tp.tv_usec/1000000.0));
}

int main() {
    int i, nx, ny, count, alpha;
    double beta, gama, sigma, hx, *x, used, begin;

    nx = 100;
    ny = 100;
    x = malloc(nx * ny * sizeof(double));
    hx = 0.01;

    used = 0.0f;

    // Versão anterior, sem armazenar cálculos repetidos em variáveis (exceto sigma)

    for(count=0; count<N; ++count) {
        begin = timestamp();

        for(i = nx; i < nx * ny - nx; ++i) {
            x[i] = 0.0f;
        }

        for(i=0; i<nx; ++i) {
            x[i] = sin(2 * M_PI * (M_PI - (i * hx))) * sigma;
            x[nx*ny-nx+i] = sin(2 * M_PI * (i * hx)) * sigma;
        }

        used += timestamp() - begin;
    }

    // Versão atual, armazenando.
/*
    for(count=0; count<N; ++count) {
        begin = timestamp();

        sigma = sinh(M_PI * M_PI);
        alpha = nx * ny - nx;
        beta = 2 * M_PI * hx;
        gama = 2 * M_PI * M_PI;

        for(i = nx; i < alpha; ++i) {
            x[i] = 0.0f;
        }

        for(i=0; i<nx; ++i) {
            x[i] = sin(gama - (i * hx)) * sigma;
            x[alpha+i] = sin(beta * i) * sigma;
        }

        used += timestamp() - begin;
    }
*/
    // Deixe uma das duas versões comentadas e teste a outra.
    printf("Time used: %.15lf\n", used/N);

    return 1;
}
