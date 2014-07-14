

#ifndef TSA_DATASTRUCTS_H

#define TSA_DATASTRUCTS_H

#include "hashDefs.h"

typedef struct TSAsam
{
    TSADATA value;
    float tStamp;
}TSAsam_t;


typedef struct TSAsubSeq
{
    TSADATA *pData;
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



typedef struct procParams
{
    //params used for preprocessing
    
    int binsPOct;
    int dsFactor;
    int removeTaniSegs;
    float minPossiblePitch;
    float varDur;
    float threshold;
    float flatThreshold;
    float maxPauseDur;
    float durMotif;
    float blackDur;
    float DTWBand;
    INDTYPE nPitchSamples;
    float interpFac[MAXNTEMPOFACTORS];
    int nInterpFac;
    int motifLengths[MAXNTEMPOFACTORS];
    int motifLengthsM1[MAXNTEMPOFACTORS];
    int indexMotifLenReal;
    int indexMotifLenLongest;
    int **combMTX;
    int *simMeasureRankRefinement;
    int nSimMeasuresUsed;
    
    
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
    char *paramOutExt;
}fileExts_t;



#endif  // TSA_DATASTRUCTS_H
