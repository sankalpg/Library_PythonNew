

#ifndef TSA_DATASTRUCTS_H

#define TSA_DATASTRUCTS_H

#include "TSAhashDefs.h"

typedef struct TSAseg
{
    float sTime;
    float eTime;
    TSAIND id;
}TSAseg_t;


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
    float mean;
    float std;
    
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
    TSADATA *data;
    int len;
    float sTime;
    float eTime;
    float mean;
    float std;    
}TSAmotifDataStorage_t;

typedef struct TSAmotifInfoExt
{
  
    TSADIST dist;
    TSAIND ind1;
    TSAIND ind2;
    int searchFileID;
    TSAIND patternID;
    TSAIND patternID1;
    TSAIND patternID2;
    float sTime;
    float eTime;
    TSAmotifDataStorage_t *storagePtr;

}TSAmotifInfoExt_t;


typedef struct sortElem
{
    TSAIND index;
    TSADIST value;
}sortElem_t;

typedef struct TSADistParams
{
    int distType;
    float DTWBand;
    int DTWType;
    int rankRefDistType;
    int distNormType;   

}TSADistParams_t;

typedef struct TSARepParams
{
    int TSRepType;
    int quantSize;
    int sampleRate;
    int normType;
    int binsPOct;
    float minPossiblePitch;
    int removeTaniSegs;
    float varDur;
    float threshold;
    float flatThreshold;
    float maxPauseDur;
    int dsFactor;
    float pitchHop;
    int useFlatNoteFile;//flat to indicate how the qubseq flatness filtering is performed

}TSARepParams_t;

typedef struct TSAPattParams
{
    float durMotif;
    int maxNMotifsPairs;
    int nInterpFac;
    int subSeqLen;

}TSAPattParams_t;


typedef struct procParams
{
    //params used for preprocessing
    TSADistParams_t distParams;
    TSARepParams_t  repParams;
    TSAPattParams_t pattParams;
    
    float interpFac[MAXNTEMPOFACTORS];
    int motifLengths[MAXNTEMPOFACTORS];
    int motifLengthsM1[MAXNTEMPOFACTORS];
    int indexMotifLenReal;
    int indexMotifLenLongest;
    int **combMTX;
    int *simMeasureRankRefinement;
    int nSimMeasuresUsed;
    int verbos;
    int dumpLogs;
    int methodVariant;
    int complexityMeasure;
    
}procParams_t;

typedef struct fileExts
{
    char tsFileExt[MAX_FEXT_CHARS];
    char tonicExt[MAX_FEXT_CHARS];
    char blackTimeExt[MAX_FEXT_CHARS];
    char subSeqFileExt[MAX_FEXT_CHARS];
    char subSeqTNFileExt[MAX_FEXT_CHARS]; //this is subSeqs which are already tonic normalized
    char subSeqInfoFileExt[MAX_FEXT_CHARS];
    char srchListExt[MAX_FEXT_CHARS];
    char flatNoteFileExt[MAX_FEXT_CHARS];
    char srchMappFileExt[MAX_FEXT_CHARS];
    char srchQueryFileExt[MAX_FEXT_CHARS];
    char disOutFileExt[MAX_FEXT_CHARS];
    char srchOutFileExt[MAX_FEXT_CHARS];
    char knnOutFileExt[MAX_FEXT_CHARS];
    char distOutFileExt[MAX_FEXT_CHARS];
    char disLogFileExt[MAX_FEXT_CHARS];
    char srchLogFileExt[MAX_FEXT_CHARS];
    char dumpLogFileExt[MAX_FEXT_CHARS];
    char outFileExt[MAX_FEXT_CHARS];

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

