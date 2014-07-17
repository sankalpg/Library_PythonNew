

#ifndef TSA_DATASTRUCTS_H

#define TSA_DATASTRUCTS_H

#include "TSAhashDefs.h"

typedef struct TSAsam
{
    TSADATA value;
    float tStamp;
}TSAsam_t;


typedef struct TSAsubSeq
{
    TSADATA *pData;
    float *pTStamps;
    TSAIND id;
    float sTime;
    float eTime;
    int len;
    int fileId;
}TSAsubSeq_t;

typedef struct TSAmPair
{
    TSAIND id1;
    TSAIND id2;
    TSADIST dist;
}TSAmPair_t;

typedef struct TSAPattern
{
    float sTime;
    float eTime;
    int fileId;
    int id;
}TSAPattern_t;

typedef struct procParams
{
    //params used for preprocessing
    float durMotif;
    float blackDur;
    int dsFactor;
    int binsPOct;    
    float minPossiblePitch;
    float varDur;
    float threshold;
    float flatThreshold;
    float maxPauseDur;
    float DTWBand;
    float interpFac[MAXNTEMPOFACTORS];
    int nInterpFac; 
    int motifLengths[MAXNTEMPOFACTORS];
    int motifLengthsM1[MAXNTEMPOFACTORS];
    int indexMotifLenReal;
    int indexMotifLenLongest;
    int **combMTX;
    int *simMeasureRankRefinement;
    int nSimMeasuresUsed;
    int SimMeasuresUsed;
    int verbos;
    int removeTaniSegs;
    int dumpLogs;
    
    
}procParams_t;

typedef struct fileExts
{
    char tsFileExt[MAX_FEXT_CHARS];
    char tonicExt[MAX_FEXT_CHARS];
    char blackTimeExt[MAX_FEXT_CHARS];
    char logFileExt[MAX_FEXT_CHARS];
    char paramsDumpExt[MAX_FEXT_CHARS];
    char outFileExt[MAX_FEXT_CHARS];
}fileExts_t;

typedef struct procLogs
{
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


#endif  // TSA_DATASTRUCTS_H
