
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

int zNorm(double *data, int len);

int meanNorm(double *data, int len);

int medianNorm(double *data, int len);

int MADNorm(double *data, int len);

int noNorm(double *data, int len);

int tonicNorm(double *data, int len);

int normalizePASAPA(double *data1, int len1, double *data2, int len2);

double computeInflectionPoints2(double *data, int len);

double computeInflectionPoints1(double *data, int len);

double measureGlobalComplexity2(double *data, int len);

double measureGlobalComplexity1(double *data, int len);


#endif //BASICDSPFUNC_H

