
#include "TSAsimilarity.h"


int compareSortElems(const void *a, const void *b)
{
    if (((sortElem_t*)a)->value > ((sortElem_t*)b)->value)
    {
        return 1;
    }
    else if (((sortElem_t*)a)->value < ((sortElem_t*)b)->value)
    {
        return -1;
    }
    return 0;
}

TSAdtwSimilarity::TSAdtwSimilarity(procLogs_t *LogsPtr)
{
    nQuery=-1;
    lenQuery=-1;
    
    nCand=-1;
    lenCand=-1;
    
    isBSFArrayInit=-1;
    isCopyByRef=-1;
    
    procLogsPtr = LogsPtr;
    
}

TSAdtwSimilarity::~TSAdtwSimilarity()
{
    if ((lenQuery!=-1))
    {
        //free memory before reconfigure
        free(accLB_Keogh_EQ);
        free(accLB_Keogh_EC);
        for(int ii=0; ii<lenQuery;ii++)
        {
            free(costMTX[ii]);
        }
        free(costMTX);
        lenQuery=-1;
    }
    if (nQuery!=-1)    
    {
        for(TSAIND ii=0; ii< nQuery; ii++)
        {
            free(envLQueryPtr[ii]);
            free(envUQueryPtr[ii]);
        }
        free(envLQueryPtr);
        free(envUQueryPtr);
        nQuery=-1;
    }
    if ((nCand!=-1)&&(isCopyByRef!=1))
    {
        deleteCandEnvMem();
    }
    if(isBSFArrayInit==1)
    {
        free(bsfArray);
    }
    
}

int TSAdtwSimilarity::deleteCandEnvMem()
{
    for(TSAIND ii=0; ii< nCand; ii++)
    {
        free(envLCandPtr[ii]);
        free(envUCandPtr[ii]);
    }
    free(envLCandPtr);
    free(envUCandPtr);
    nCand=-1;
}


int TSAdtwSimilarity::configureTSASimilarity(int lenQ, int lenC, float globalConst)
{
    
    lenQuery = lenQ;
    lenCand = lenC;
    accLB_Keogh_EQ = (TSADATA*)malloc(sizeof(TSADATA)*lenQuery);
    accLB_Keogh_EC = (TSADATA*)malloc(sizeof(TSADATA)*lenCand);
    
    for(int ii=0; ii< lenQuery; ii++)
    {
        accLB_Keogh_EQ[ii]=0;
        accLB_Keogh_EC[ii]=0;
    }
    
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

int TSAdtwSimilarity::initArrayBSF(TSAIND len)
{
    bsfArray = (TSADIST *)malloc(sizeof(TSADIST)*len);
    isBSFArrayInit=1;
    for (TSAIND ii=0;ii<len;ii++)
    {
        bsfArray[ii] = INF;
    }
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
    float t1,t2;
    t1 = clock();
    //assign memory to store envelops
    envLQueryPtr = (TSADATA **)malloc(sizeof(TSADATA *)*nQuery);
    envUQueryPtr = (TSADATA **)malloc(sizeof(TSADATA *)*nQuery);
    
    for (TSAIND ii=0; ii < nQuery; ii++)
    {
        envLQueryPtr[ii] = (TSADATA *)malloc(sizeof(TSADATA)*subSeqQueryPtr[ii].len);
        envUQueryPtr[ii] = (TSADATA *)malloc(sizeof(TSADATA)*subSeqQueryPtr[ii].len);
        computeRunningMinMax(subSeqQueryPtr[ii].pData, envUQueryPtr[ii], envLQueryPtr[ii], subSeqQueryPtr[ii].len, bandDTW);
    }
    t2=clock();
    procLogsPtr->tGenEnv += (t2-t1)/CLOCKS_PER_SEC;
    
    return 1;
}

int TSAdtwSimilarity::computeCandEnvelops()
{
    float t1,t2;
    t1 = clock();
    //assign memory to store envelops
    envLCandPtr = (TSADATA **)malloc(sizeof(TSADATA *)*nCand);
    envUCandPtr = (TSADATA **)malloc(sizeof(TSADATA *)*nCand);
    
    for (TSAIND ii=0; ii < nCand; ii++)
    {
        envLCandPtr[ii] = (TSADATA *)malloc(sizeof(TSADATA)*subSeqCandPtr[ii].len);
        envUCandPtr[ii] = (TSADATA *)malloc(sizeof(TSADATA)*subSeqCandPtr[ii].len);
        computeRunningMinMax(subSeqCandPtr[ii].pData, envUCandPtr[ii], envLCandPtr[ii], subSeqCandPtr[ii].len, bandDTW);
    }
    
    t2=clock();
    procLogsPtr->tGenEnv += (t2-t1)/CLOCKS_PER_SEC;
    
    return 1;
}

int TSAdtwSimilarity::copyQueryEnv2Cand()
{
    envLCandPtr = envLQueryPtr;
    envUCandPtr = envUQueryPtr;
    nCand = nQuery;
    isCopyByRef=1;
    
}




