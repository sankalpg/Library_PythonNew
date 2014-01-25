
#include <stdio.h>
#include "../basicDSPFuncs/basicDSPCFuncs.h"
#include <math.h>
#include <stdlib.h>
#include <time.h>

#define M_PI 3.14159265358979323846

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

FILE *fp1, *fp2, *fp3;

int runningMinMaxBruteforce(double *data, double *U, double* L, int len, int winLen)
{
    int ii,jj, str,end;
    double minV, maxV;
    
    
    for (ii=0;ii<len;ii++)
    {
        str = MAX(0,ii-winLen);
        end = MIN(len-1,ii+winLen);
        minV = data[str];
        maxV = data[str];
        for(jj=str;jj<=end;jj++)
        {
            if(data[jj]>maxV)
            {
                maxV=data[jj];
            }
            if(data[jj]<minV)
            {
                minV=data[jj];
            }
        }
        U[ii]=maxV;
        L[ii]=minV;
    }
}


int main (int argc , char *argv[])
{
    double *data, *U1, *U2, *L1, *L2;
    int len = 500;
    int ii,jj;
    double t1,t2;
    
    data = (double*)malloc(sizeof(double)*len);
    for(ii=0;ii<len;ii++)
    {
        data[ii]= sin(2*M_PI*(double)ii/(double)len);
    }
    
    U1 = (double*)malloc(sizeof(double)*len);
    U2 = (double*)malloc(sizeof(double)*len);
    L1 = (double*)malloc(sizeof(double)*len);
    L2 = (double*)malloc(sizeof(double)*len);
    
    t1 = clock();
    for (ii=0;ii<100000;ii++)
    {
        computeRunningMinMax(data, U1, L1, len, 50);
    }
    t2=clock();
    printf("Time taken to do smart way :%f\n",(t2-t1)/CLOCKS_PER_SEC);
    
    t1 = clock();
    for (ii=0;ii<100000;ii++)
    {
        runningMinMaxBruteforce(data, U2, L2, len, 50);
    }
    t2=clock();
    printf("Time taken to do brute force way :%f\n",(t2-t1)/CLOCKS_PER_SEC);
    
    
    
    fp1 = fopen("sine.txt","w");
     for(ii=0;ii<len;ii++)
    {
        fprintf(fp1, "%f\t%f\t%f\t%f\t%f\n",data[ii], U1[ii],L1[ii],U2[ii],L2[ii]);
    }
    fclose(fp1);
    
    
}