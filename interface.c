#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#define ARGS_NUM 9

double Hx,Hy;
int MaxIter;
FILE *fp;

void print_errno(void) {
    printf("%s",strerror(errno));
}

void getParams(int argc, char* argv[]) {
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
            MaxIter = atoi(argv[i+1]);
        } else if(strcmp(argv[i],"-o") == 0) {
            fp = fopen(argv[i+1],"w");
        } else {
            fprintf(stderr,"Incorrect parameter.\n");
            exit(-1);
        }
    }
}

int main(int argc, char* argv[]) {
    getParams(argc,argv);
}