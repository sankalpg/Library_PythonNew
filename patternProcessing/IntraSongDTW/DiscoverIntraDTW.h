
#ifndef DiscoverIntraDTW_H

#define DiscoverIntraDTW_H


//TYPE DEFS
#define DATATYPE       double
#define DISTTYPE        double
#define INDTYPE        long long

//preprocessor defines
//#define DEBUG
//#define ENABLE_EA


// other defines
#define INF 999999999999.0
#define LOG2  0.693147180559945
#define EPS 0.0000000000000001
#define computeLBkimFL(a,b,c,d) ((a-b)*(a-b)) + ((c-d)*(c-d))


#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <fcntl.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <signal.h>
#include <stdint.h>
#include <string.h>
#include "../../similarityMeasures/dtw/dtw.h"
#include "../../basicDSPFuncs/basicDSPCFuncs.h"
#include <float.h>





// structures
typedef struct motifInfo
{
  
    DISTTYPE dist;
    INDTYPE ind1;
    INDTYPE ind2;
    
}motifInfo;

//functions
DISTTYPE manageTopKMotifs(motifInfo *topKmotifs, float *tStamps, int K, INDTYPE ind1 , INDTYPE ind2, DISTTYPE dist, float blackDur);

long long getNumLines(const char *file);

void linearlyInterpolate(float *array, int size, float val1, float val2);

#endif //#ifndef DiscoverIntraDTW_H