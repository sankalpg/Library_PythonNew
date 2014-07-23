
#include "TSAsimilarity.h"

TSAdtwSimilarity::TSAdtwSimilarity()
{
    nQuery=-1;
    lenQuery=-1;
    
    nCand=-1;
    lenCand=-1;
    
}

int TSAdtwSimilarity::configureTSASimilarity(int lenQ, int lenC, float globalConst)
{
    
    if ((lenQuery!=-1)||(lenCand!=-1))
    {
        //free memory before reconfigure
        free(accLB_Keogh_EQ);
        free(accLB_Keogh_EC);
        for(int ii=0; ii<lenQuery;ii++)
        {
            free(costMTX[ii]);
        }
        free(costMTX);
        
      }
    
    lenQuery = lenQ;
    lenCand = lenC;
    accLB_Keogh_EQ = (TSADATA*)malloc(sizeof(TSADATA)*lenQuery);
    accLB_Keogh_EC = (TSADATA*)malloc(sizeof(TSADATA)*lenCand);
    
    costMTX = (TSADIST **)malloc(sizeof(TSADIST *)*lenQuery);
    for(int ii=0; ii<lenQuery;ii++)
    {
        costMTX[ii] = (TSADATA *)malloc(sizeof(TSADATA)*lenCand);
        for(int jj=0; jj<lenCand;jj++)
        {
            costMTX[ii][jj]=FLT_MAX;
        }
    }
    
    bandDTW = (int)floor(lenQuery*globalConst);
    
}

int TSAdtwSimilarity::setQueryPtr(TSAsubSeq_t *qPtr, TSAIND nQ)
{
    subSeqQueryPtr = qPtr;
    nQuery = nQ;
    
}
int TSAdtwSimilarity::setCandPtr(TSAsubSeq_t *cPtr, TSAIND nC)
{
    subSeqCandPtr = cPtr;
    nCand = nC;
}

int TSAdtwSimilarity::computeQueryEnvelops()
{
    //assign memory to store envelops
    envLQueryPtr = (TSADATA **)malloc(sizeof(TSADATA *)*nQuery);
    envUQueryPtr = (TSADATA **)malloc(sizeof(TSADATA *)*nQuery);
    
    for (TSAIND ii=0; ii < nQuery; ii++)
    {
        envLQueryPtr[ii] = (TSADATA *)malloc(sizeof(TSADATA)*subSeqQueryPtr[ii].len);
        envUQueryPtr[ii] = (TSADATA *)malloc(sizeof(TSADATA)*subSeqQueryPtr[ii].len);
        computeRunningMinMax(subSeqQueryPtr[ii].pData, envUQueryPtr[ii], envLQueryPtr[ii], subSeqQueryPtr[ii].len, bandDTW);
    }
    return 1;
}

int TSAdtwSimilarity::computeCandEnvelops()
{
    //assign memory to store envelops
    envLCandPtr = (TSADATA **)malloc(sizeof(TSADATA *)*nCand);
    envUCandPtr = (TSADATA **)malloc(sizeof(TSADATA *)*nCand);
    
    for (TSAIND ii=0; ii < nCand; ii++)
    {
        envLCandPtr[ii] = (TSADATA *)malloc(sizeof(TSADATA)*subSeqCandPtr[ii].len);
        envUCandPtr[ii] = (TSADATA *)malloc(sizeof(TSADATA)*subSeqCandPtr[ii].len);
        computeRunningMinMax(subSeqCandPtr[ii].pData, envUCandPtr[ii], envLCandPtr[ii], subSeqCandPtr[ii].len, bandDTW);
    }
    return 1;
}

int TSAdtwSimilarity::copyQueryEnv2Cand()
{
    envLCandPtr = envLQueryPtr;
    envUCandPtr = envUQueryPtr;
    nCand = nQuery;
    
}




