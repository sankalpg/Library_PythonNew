
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


#endif //BASICDSPFUNC_H