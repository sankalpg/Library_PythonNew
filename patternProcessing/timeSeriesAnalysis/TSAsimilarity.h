
#ifndef TSASIMILARITY_H

#define TSASIMILARITY_H

#include "TSAhashDefs.h"
#include "TSAdataStructs.h"

class TSAdtwSimilarity
{
public:
    
    int bandDTW;
    
    TSAsubSeq_t *subSeqQueryPtr;
    TSAsubSeq_t *subSeqCandPtr;
    
    TSADATA **envLQueryPtr;
    TSADATA **envUQueryPtr;
    TSAIND nQuery;
    int lenQuery;
    
    TSADATA **envLCandPtr;
    TSADATA **envUCandPtr;
    TSAIND nCand;
    int lenCand;
    
    TSADIST *accLB_Keogh_EQ;
    TSADIST *accLB_Keogh_EC;
    
    TSADIST **costMTX;
    
    TSADIST *bsfArray;
    
    TSAdtwSimilarity();
    
    int     configureTSASimilarity(int lenQ, int lenC, float globalConst);
    int     setQueryPtr(TSAsubSeq_t *qPtr, TSAIND nQ);
    int     setCandPtr(TSAsubSeq_t *cPtr, TSAIND nC);
    int     computeQueryEnvelops();
    int     computeCandEnvelops();
    int     copyQueryEnv2Cand();
    int     initArrayBSF(TSAIND len);
    
    
};

#endif //TSASIMILARITY_H