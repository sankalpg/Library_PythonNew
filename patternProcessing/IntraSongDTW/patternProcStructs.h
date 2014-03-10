

#ifndef PATTERNPROCSTRUCTURES_H

#define PATTERNPROCSTRUCTURES_H

#include "hashDefs.h"

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


typedef struct motifInfo
{
  
    DISTTYPE dist;
    INDTYPE ind1;
    INDTYPE ind2;
    
}motifInfo_t;



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
    float DTWBand;
}procParams_t;



#endif //PATTERNPROCSTRUCTURES_H