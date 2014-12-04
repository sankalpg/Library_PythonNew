/*
 * This implementation of running min and max is proposed by "Daniel Lemire" in 
 * STREAMING MAXIMUM-MINIMUM FILTER USING NO
MORE THAN THREE COMPARISONS PER ELEMENT
 * This is supposed to be quite a fast implementation of running min and max
 *
 * Author: Sankalp Gulati
 * No part of this code is taken from anywhere and is solely based on the pseudocode in above paper
 * NOTE: this code can be made memory efficient by using push pop kind of things rather than allocating big memory chunk.
 * 
 * I made it computationally efficient without worrying much about memory inefficiency
 * Note that min and max is returned for a range of 2*winLen+1 indices
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "basicDSPCFuncs.h"

int computeRunningMinMax(double *data, double *U, double *L, int lenData, int winLen)
{
    
    int UstartInd, UnumPoints, LstartInd, LnumPoints, ii;
    int *utemp, *ltemp;
    
    
    //without worrying much about memory just allocating memory for temporary monotonic wedges used for this implementation
    utemp = (int*)malloc(sizeof(int)*lenData);
    ltemp = (int*)malloc(sizeof(int)*lenData);
    
    //Initializing
    UstartInd=0;
    utemp[UstartInd] = 0;
    UnumPoints=1;
    
    LstartInd=0;
    ltemp[LstartInd] = 0;
    LnumPoints=1;
    
    for(ii=1;ii<lenData;ii++)
    {
        if (ii>winLen)
        {
            U[ii-winLen-1] = data[utemp[UstartInd]];
            L[ii-winLen-1] = data[ltemp[LstartInd]];
        }
        
        if(data[ii]>data[ii-1])
        {
            UnumPoints-=1;
            
            while ((UnumPoints>0)&&(data[ii]>data[utemp[UstartInd+UnumPoints-1]]))
            {
                UnumPoints-=1;
            }
        }
        else
        {
            LnumPoints-=1;
            
            while ((LnumPoints>0)&&(data[ii]< data[ltemp[LstartInd+LnumPoints-1]]))
            {
                LnumPoints-=1;
            }
            
        }
        UnumPoints+=1;
        LnumPoints+=1;
        utemp[UstartInd+UnumPoints-1]=ii;
        ltemp[LstartInd+LnumPoints-1]=ii;
        
        if (ii == 2*winLen+1 + utemp[UstartInd])
        {
            UstartInd+=1;
            UnumPoints-=1;
        }        
        else if (ii == 2*winLen+1 + ltemp[LstartInd])
        {
            LstartInd+=1;
            LnumPoints-=1;
        }
            
        
    }
    
    for (ii = lenData; ii < lenData+winLen+1; ii++) 
    {

        U[ii-winLen-1] = data[utemp[UstartInd]];
        L[ii-winLen-1] = data[ltemp[LstartInd]];

        if (ii-utemp[UstartInd] >= 2*winLen+1)           
        {
            UstartInd+=1;
            UnumPoints-=1;
        }

        if (ii-ltemp[LstartInd] >= 2*winLen+1)           
        {
            LstartInd+=1;
            LnumPoints-=1;
        }
    }
    
    free(utemp);
    free(ltemp);
    
    return 1;
    
}

void linearInterpolateV1(double *dataInp, double *dataOut, float *indInt, int N)
{
        int ii =0, ind;
        double x0, x1, x2,x3,x4,t, a0, a1,a2,a3;
        
        for(ii=0;ii<N;ii++)
        {
            ind = floor(indInt[ii]);
            t = indInt[ii] - ind;
            
            dataOut[ii] = dataInp[ind]*(1-t)  + dataInp[ind+1]*t ; 
        }
        
}


/*
 * READ this reference for a nice comparison
 * reference : http://radio.feld.cvut.cz/AES/atp2002/proc/paper07.pdf
 */

/*Not that this is not cubic spline interpolation*/
void cubicInterpolate(double *dataInp, double *dataOut, float *indInt, int N)
{
    //NOTE you should pass N such that dataInp doesn't get out of index
    
        int ii =0, ind, max_ind;
        double x0, x1, x2,x3,x4,t, a0, a1,a2,a3;
        max_ind = ceil(indInt[N-1]);
        
        for(ii=0;ii<N;ii++)
        {
            ind = floor(indInt[ii]);
            t = indInt[ii] - ind;
            x0 = dataInp[max(0, ind-1)];
            x1 = dataInp[ind];
            x2 = dataInp[ind+1];
            x3 = dataInp[min(max_ind, ind+2)];
            dataOut[ii] = x1 + (t/6.0)*(((-1.0*(t*t) + (3.0*t) -2 )*x0)  + 3*(((t*t) - (2.0*t) - 1)*x1) + 3.0*((-1.0*(t*t) + (t) + 2)*x2) + (((t*t) - 1)*x3) );
            
        }

}
/*This is faster than cubic and not much worse in performance*/
void quadraticInterpolate(double *dataInp, double *dataOut, float *indInt, int N)
{
        int ii =0, ind, max_ind;
        double x0, x1, x2,x3,x4,t, a0, a1,a2,a3;
        
        max_ind = ceil(indInt[N-1]);
        
        for(ii=0;ii<N;ii++)
        {
            ind = floor(indInt[ii]);
            t = indInt[ii] - ind;
            x0 = dataInp[max(0, ind-1)];
            x1 = dataInp[ind];
            x2 = dataInp[ind+1];
            x3 = dataInp[min(max_ind, ind+2)];
            dataOut[ii] = x1 + (t/2.0)*(t*(x1-2.0*x2+x3) -3.0*x1+4.0*x2-x3);
        }

}

double quantizePitch(double pitchCents, int binsPDiv)
{
    
    double pitchOut;
    
    pitchOut = floor((pitchCents/binsPDiv)+0.5)*binsPDiv;
    
    return pitchOut;
}

double computeMean(double* data, int len)
{
    int ii;
    double m = 0;
    for(ii=0;ii<len;ii++)
    {
        m+=data[ii];
    }
    return m/len;
}

double computeSTD(double *data, int len, double mean)
{
    double diffSq = 0;
    int ii;
    double diff;
    
    for(ii=0;ii<len;ii++)
    {
        diff = (data[ii]-mean);
        diffSq+=diff*diff;
    }
    return sqrt(diffSq/len);
}


int compare (const void * a, const void * b)
{
  return ( *(double*)a - *(double*)b );
}

double computeMedian(double* data, int len)
{
    int ii;
    int median;

    qsort (data, len, sizeof(double), compare);
    
    if (len%2==0)
    {
        median = (data[(int)floor(len/2)-1] + data[(int)ceil(len/2)-1])/2;
    }
    else
    {
        median = data[(int)ceil(len/2)-1];
    }
    
    return median;
}

double computeMAD(double* data, int len, double median)
{
    int ii;
    int mad;
    
    for(ii=0; ii<len; ii++)
    {
        data[ii] = abs(data[ii]-median);
    }

    qsort (data, len, sizeof(double), compare);
    
    if (len%2==0)
    {
        mad = (data[(int)floor(len/2)-1] + data[(int)ceil(len/2)-1])/2;
    }
    else
    {
        mad = data[(int)ceil(len/2)-1];
    }
    
    return mad;
}

