
#ifndef DiscoverIntraDTW_H

#define DiscoverIntraDTW_H


#include "MotifDataIO.h"

//functions
DISTTYPE manageTopKMotifs(motifInfo *topKmotifs, segInfo_t *tStamps1, segInfo_t *tStamps2, int K, INDTYPE ind1 , INDTYPE ind2, DISTTYPE dist, float blackDur, int searchFileID);
int manageTopKMotifsData(motifInfo *topKmotifs, longTermDataStorage_t *longTermDataStorage, DATATYPE** dataInterp, segInfo_t *tStampsInterp, INDTYPE *patternID, int *emptySpaceInd, int lenMotifReal, int K, int searchFileID);
#endif //#ifndef DiscoverIntraDTW_H




