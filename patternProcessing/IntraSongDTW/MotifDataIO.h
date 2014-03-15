
#ifndef MOTIFDATAIO_H

#define MOTIFDATAIO_H

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
#include <float.h>
#include <string.h>

#include "../../similarityMeasures/dtw/dtw.h"
#include "../../basicDSPFuncs/basicDSPCFuncs.h"
#include "MotifDataIO.h"
#include "patternProcStructs.h"



// structures
typedef struct procLogs
{
    char *commitID;
    
    float timeDataLoad;
    float timeGenSubs;
    float timeRemBlacklist; 
    float timeGenEnvelops;
    float timeDiscovery;
    float timeWriteData;
    float timeTotal;
    
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
    char *paramOutExt;
}fileExts_t;


typedef struct mappInfo
{
    INDTYPE last_line;
    
}mappInfo_t;


INDTYPE getNumLines(const char *file);

INDTYPE readPreProcessGenDB(DATATYPE ***d, segInfo_t **t, int *motifLen, char *baseName, fileExts_t *myFileExts, procParams_t *myProcParams, procLogs_t *myProcLogs, int verbos);
INDTYPE loadSeedMotifSequence(DATATYPE ***d, segInfo_t **t, int *motifLen, char *baseName, fileExts_t *myFileExts, procParams_t *myProcParams, procLogs_t *myProcLogs, int maxNMotifsPairs, int verbos);

void dumpSearchMotifInfo(char *motifFile, motifInfo** topKmotifs, segInfo_t *tStampsInterpSeed, segInfo_t *tStampsInterp, int NSeeds, INDTYPE K, int nInterFact, int verbos);
void dumpDiscoveryLogs(char *logFile, procLogs_t myProcLogs, int verbos);
void dumpDiscoveredMotifInfo(char *motifFile, motifInfo *topKmotifs, segInfo_t *tStampsInterp, int K, int verbos);

void dumpParameterValuesUsed(char *paramOutFile, procParams_t *myProcParams);

#endif //MOTIFDATAIO_H

