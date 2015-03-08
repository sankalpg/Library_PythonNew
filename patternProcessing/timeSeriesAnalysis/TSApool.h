#ifndef TSAPOOL_H

#define TSAPOOL_H

#include "TSAhashDefs.h"
#include "TSAdataStructs.h"


class TSApool
{
public:
    
    int K;
    TSAIND numQueries;
    int discOrSear;
    int useLTStorage;
    TSAmotifInfo_t *priorityQDisc;
    TSAmotifInfoExt_t **priorityQSear;
    TSAIND *pattsPerQ;  //for mainting queues for < distance threshold type of situation, we need to keep track of how many elements are added to each queue since it will be different for each query or queue.
    TSAIND *spaceAllocQ;
    int isQKNN;
    int isQDist;
    TSAmotifDataStorage_t **longTermDataStorage;
    TSAIND patternID;
    int *emptySpaceInd;
    
    TSApool();
    ~TSApool();
    TSApool(int K);
    int initPriorityQDisc();
    int initPriorityQSear(TSAIND nQueries);
    int initPriorityQSearDist(TSAIND nQueries);
    int extendPriorityQSearDist(TSAIND queryInd);
    
    TSADIST managePriorityQDisc(TSAsubSeq_t *subSeqPtr, TSAIND ind1, TSAIND ind2, TSADIST dist, float blackDur);
    TSADIST managePriorityQSear(TSAIND queryInd, TSAsubSeq_t *subSeqPtr, TSAIND ind1, TSAIND ind2, TSADIST dist, int searchFileID, float blackDur);
    int     managePriorityQSearDIST(TSAIND queryInd, TSAsubSeq_t *subSeqPtr, TSAIND ind1, TSAIND ind2, TSADIST dist, int searchFileID);
    
    int     updatePattStorageData(TSAIND queryInd, TSAsubSeq_t *subSeqPtr, int lenMotifReal, int searchFileID);
    int     initPattStorage(TSAIND nQueries, int lenMotifReal);
    
    int     sortQSearch(TSAIND queryInd);
    
    
};


#endif  //TSAPOOL_H