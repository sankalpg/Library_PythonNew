
#include <stdio.h>
#include "../basicDSPFuncs/basicDSPCFuncs.h"
#include <math.h>
#include <stdlib.h>
#include <time.h>

#define M_PI 3.14159265358979323846

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

FILE *fp;

int main (int argc , char *argv[])
{
    
    double *sine, *sine2, *sineIntC,*sineIntQ, *sineIntL, factor=1.1;
    float  *indexes;
    int N = 2000, M=1000, f = 100;
    sine = (double *)malloc(sizeof(double)*N);
    sine2 = (double *)malloc(sizeof(double)*N);
    sineIntC = (double *)malloc(sizeof(double)*N);
    sineIntL = (double *)malloc(sizeof(double)*N);
    sineIntQ = (double *)malloc(sizeof(double)*N);
    indexes = (float *)malloc(sizeof(float)*N);
    
    for (int ii=0;ii<N;ii++)
    {
        sine[ii] = sin(2*M_PI*(float)ii*f/(float)N);
        sine2[ii] = sin(2*M_PI*(float)ii*f*factor/(float)N);
        indexes[ii] = factor*ii;
    }
    
    
    cubicInterpolate(sine, sineIntC, indexes, M);
    linearInterpolateV1(sine, sineIntL, indexes, M);
    quadraticInterpolate(sine, sineIntQ, indexes, M);
    
    
    fp = fopen("sine.txt","w");
    for(int ii=0;ii<M;ii++)
    {
        fprintf(fp, "%f\t%f\t%f\t%f\t%f\n",sine[ii], sine2[ii],sineIntL[ii], sineIntQ[ii], sineIntC[ii]);
    }
    fclose(fp);
    
    
}