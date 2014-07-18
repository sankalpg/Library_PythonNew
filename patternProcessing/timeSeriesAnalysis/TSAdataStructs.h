

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

typedef struct TSAmotifInfo
{
    TSAIND ind1;
    TSAIND ind2;
    TSADIST dist;
}TSAmotifInfo_t;

typedef struct TSAPattern
{
    float sTime;
    float eTime;
    int fileId;
    int id;
}TSAPattern_t;

typedef struct TSAmotifDataStorage
{
  
    TSAIND patternID;
    TSADIST *data;
    float strTime;
    float endTime;
    
    
    
}TSAmotifDataStorage_t;

typedef struct TSAmotifInfoExt
{
  
    TSADIST dist;
    TSAIND id1;
    TSAIND id2;
    int searchFileID;
    TSAIND patternID;
    TSAmotifDataStorage_t *storagePtr;

}TSAmotifInfoExt_t;


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
    fileExts
}fileExts_t;

typedef struct procLogs
{
    float tLoadTS;
    float tProcTS;
    float tGenEnv; 
    float tGenUniScal;
    float tTotal;
    float tDump;
    
    long long lenTS;
    long long lenProcTS;
    long long nSubSeqs;
    long long nSubSeqsBL;
    long long nProcSubSeq;
    
    long long nLB_KIM_FL;
    long long nLB_Keogh_EQ;
    long long nLB_Keogh_EC;
    long long nDTW_EA;
    long long nPriorityUpdates;
    
}procLogs_t;


#endif  // TSA_DATASTRUCTS_H