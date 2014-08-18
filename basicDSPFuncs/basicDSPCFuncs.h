
#ifndef BASICDSPFUNC_H


#define BASICDSPFUNC_H
#ifndef max
    #define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

#ifndef min
    #define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif


int computeRunningMinMax(double *data, double *U, double *L, int lenData, int winLen);

void cubicInterpolate(double *dataInp, double *dataOut, float *indInt, int N);

void linearInterpolateV1(double *dataInp, double *dataOut, float *indInt, int N);

void quadraticInterpolate(double *dataInp, double *dataOut, float *indInt, int N);

double quantizePitch(double pitchCents, int binsPDiv);

double computeMean(double* data, int len);

double computeSTD(double *data, int len, double mean);

double computeMAD(double* data, int len, double median);

int compare (const void * a, const void * b);

double computeMedian(double* data, int len);

#endif //BASICDSPFUNC_H

