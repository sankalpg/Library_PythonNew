
#ifndef DiscoverIntraDTW_H

#define DiscoverIntraDTW_H


#include "MotifDataIO.h"

//functions
DISTTYPE manageTopKMotifs(motifInfo *topKmotifs, segInfo_t *tStamps1, segInfo_t *tStamps2, int K, INDTYPE ind1 , INDTYPE ind2, DISTTYPE dist, float blackDur);

#endif //#ifndef DiscoverIntraDTW_H



