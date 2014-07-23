#ifndef TSAPOOL_H

#define TSAPOOL_H

#include "TSAhashDefs.h"
#include "TSAdataStructs.h"


class TSApool
{
public:
    
    int K;
    float blackDur;
    TSAmotifInfo_t *priorityQDisc;
    TSAmotifInfoExt_t **priorityQSear;
    
    TSApool();
    TSApool(int K, float blackDur);
    int initPriorityQDisc();
    int initPriorityQSear(TSAIND nQueries);
    
    TSADIST managePriorityQDisc(TSAsubSeq_t *subSeqPtr, TSAIND ind1, TSAIND ind2, TSADIST dist);
    TSADIST managePriorityQSear(TSAIND queryInd, TSAsubSeq_t *subSeqPtr, TSAIND ind1, TSAIND ind2, TSADIST dist, int searchFileID);
    
    
};


#endif  //TSAPOOL_H