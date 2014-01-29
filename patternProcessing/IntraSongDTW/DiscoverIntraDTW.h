
#ifndef DiscoverIntraDTW_H

#define DiscoverIntraDTW_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <signal.h>
#include <stdint.h>
#include <string.h>
#include "../../similarityMeasures/dtw/dtw.h"
#include "../../basicDSPFuncs/basicDSPCFuncs.h"
#include <float.h>

//preprocessor defines!!
//#define FIXPDATA
#define FLOATDATA


#ifdef FIXPDATA
#define DATATYPE       uint32_t
#else
#define DATATYPE       double
#endif

#define DISTTYPE        double
#define INDTYPE        long long

#define INF 999999999999.0
# define computeLBkimFL(a,b,c,d) ((a-b)*(a-b)) + ((c-d)*(c-d))



typedef struct motifInfo
{
  
    DISTTYPE dist;
    INDTYPE ind1;
    INDTYPE ind2;
    
}motifInfo;

#endif //#ifndef DiscoverIntraDTW_H