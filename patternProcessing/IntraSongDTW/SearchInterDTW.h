
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
//#define computeLBkimFL(a,b,c,d) ((a-b)*(a-b)) + ((c-d)*(c-d))


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
    
}motifInfo_t;

typedef struct segInfo
{
    float str;
    float end;
}segInfo_t;

typedef struct segInfoInterp
{
    float str;
    float end;
    float endInterpH;
    float endInterpL;
    
}segInfoInterp_t;

typedef struct procLogs
{
    char *commitID;
    
    float timeDataLoad;
    float timeGenSubs;
    float timeRemBlacklist; 
    float timeGenEnvelops;
    float timeDiscovery;
    
    long long totalPitchSamples;
    long long totalPitchNonSilSamples;
    long long totalSubsGenerated;
    long long totalSubsBlacklisted;
    long long totalSubsInterpolated;
    
    long long totalFLDone;
    long long totalLBKeoghEQ;
    long long totalLBKeoghEC;
    long long totalDTWComputations;
    long long totalPriorityUpdates;
}procLogs_t;

typedef struct procParams
{
    int binsPOct;
    int dsFactor;
    int removeTaniSegs;
    float minPossiblePitch;
    float allowedSilDur;
    float varDur;
    float threshold;
    float flatThreshold;
    float maxPauseDur;
    float factorLow;
    float factorHigh;
    float durMotif;
    float blackDur;
}procParams_t;

typedef struct fileExts
{
    char *pitchExt;
    char *tonicExt;
    char *segExt;
    char *seedMotifExt;
    char *logExt;
    char *searchExt;
    char *motifExt;
    char *mappExt;
}fileExts_t;



//functions
DISTTYPE manageTopKMotifs(motifInfo *topKmotifs, segInfo_t *tStamps1, segInfo_t *tStamps2, int K, INDTYPE ind1 , INDTYPE ind2, DISTTYPE dist, float blackDur);

#endif //#ifndef DiscoverIntraDTW_H


