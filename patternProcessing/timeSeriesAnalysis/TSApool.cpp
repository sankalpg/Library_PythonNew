
#include "TSApool.h"


TSApool::TSApool()
{
    K=-1;
    patternID = 0;
}
TSApool::TSApool(int n, float bDur)
{
    K=n;
    blackDur = bDur;
}

int TSApool::initPriorityQDisc()
{
    if(K==-1)
    {
        printf("Initialize the number of elements in priority Queue first");
        exit(0);
    }
    priorityQDisc = (TSAmotifInfo_t *)malloc(sizeof(TSAmotifInfo_t)*K);
    for(int ii=0;ii<K;ii++)
    {
        priorityQDisc[ii].dist=INF;
        priorityQDisc[ii].ind1=0;
        priorityQDisc[ii].ind2=0;
        
    }
    
    return 1;    
}
int TSApool::initPriorityQSear(TSAIND nQueries)
{
    if(K==-1)
    {
        printf("Initialize the number of elements in priority Queue first");
        exit(0);
    }
    
    priorityQSear = (TSAmotifInfoExt_t **)malloc(sizeof(TSAmotifInfoExt_t*)*nQueries);
    for (int ii=0; ii < nQueries; ii++)
    {
        priorityQSear[ii] = (TSAmotifInfoExt_t *)malloc(sizeof(TSAmotifInfoExt_t)*K);
        for(int jj=0;jj<K;jj++)
        {
            priorityQSear[ii][jj].dist=INF;
            priorityQSear[ii][jj].ind1=0;
            priorityQSear[ii][jj].ind2=0;
            priorityQSear[ii][jj].sTime = -1;
            priorityQSear[ii][jj].eTime = -1;
            priorityQSear[ii][jj].patternID=PID_DEFAULT2;
            priorityQSear[ii][jj].searchFileID = FID_DEFAULT1;
        }
    }
    
    return 1;
}

int TSApool::initPattStorage(TSAIND nQueries, int lenMotifReal)
{
    longTermDataStorage = (TSAmotifDataStorage_t **)malloc(sizeof(TSAmotifDataStorage_t*)*nQueries);
    for(TSAIND jj=0;jj<nQueries;jj++)
    {
        longTermDataStorage[jj] = (TSAmotifDataStorage_t *)malloc(sizeof(TSAmotifDataStorage_t)*K);
        
        for(TSAIND ii=0;ii<K;ii++)
        {
            longTermDataStorage[jj][ii].data = (DATATYPE *)malloc(sizeof(DATATYPE)*lenMotifReal);
            longTermDataStorage[jj][ii].patternID = PID_DEFAULT1;
            longTermDataStorage[jj][ii].len = lenMotifReal;
        }
        
    }
    emptySpaceInd = (int *)malloc(sizeof(int)*K);
    
}

int TSApool::updatePattStorageData(TSAmotifInfoExt_t *priorityQSear, TSAsubSeq_t *subSeqPtr, int lenMotifReal, int searchFileID)
{
    int match_found, *emptySpacePtr;
    int emptySpaceCnt=0;
    int pp,ss;
    
    //Lets find out what patterns are removed from the priority list after processing one seed over entire search file. if patterns are removed, store there location is emptySapceInd buffer which we later use to copy new data (data which is newly added to priority list)
    match_found=0;
    for(pp=0;pp<K;pp++)
    {
        match_found=0;
        for (ss=0;ss<K;ss++)
        {
            if (longTermDataStorage[pp].patternID == priorityQSear[ss].patternID)
            {
                match_found=1;
                break;
            }
        }
        if(match_found==0)
        {
            emptySpaceInd[emptySpaceCnt] = pp;
            emptySpaceCnt++;
        }
    }
    
    emptySpacePtr = emptySpaceInd;
    //After processing every seed motif over a file, check what are the new patterns added in the priority list and store data corresponding to them in the longTermStorage structure and assign them pattern ID so that we can track these patterns in the future
    for(pp=0;pp<K;pp++)
    {
        //newly added patterns will be the ones added from the current search file so just search only in that domain
        if ((priorityQSear[pp].searchFileID == searchFileID)&&(priorityQSear[pp].patternID==PID_DEFAULT3))
        {
            priorityQSear[pp].patternID = patternID;
            patternID++;
            
            //Also in such case add data to the long term storage;
            if(emptySpaceCnt>0)
            {
                memcpy(longTermDataStorage[emptySpacePtr[0]].data, subSeqPtr[priorityQSear[pp].ind2].pData, sizeof(DATATYPE)*lenMotifReal);
                longTermDataStorage[emptySpacePtr[0]].patternID = priorityQSear[pp].patternID;
                longTermDataStorage[emptySpacePtr[0]].sTime = subSeqPtr[priorityQSear[pp].ind2].sTime;
                longTermDataStorage[emptySpacePtr[0]].eTime = subSeqPtr[priorityQSear[pp].ind2].eTime;
                priorityQSear[pp].storagePtr = &longTermDataStorage[emptySpacePtr[0]];
                emptySpacePtr++;
                emptySpaceCnt--;
            }
            else
            {
                printf("SOMETHING TERRIBLE HAPPENED");
                return -1;
            }
        }
    }
    
    return 1;
}



TSADIST TSApool::managePriorityQDisc(TSAsubSeq_t *subSeqPtr, TSAIND ind1, TSAIND ind2, TSADIST dist)
{
    int sortInd=-1, matchInd=-1;
    
    for(int ii=0;ii<K;ii++)
    {
        if ((priorityQDisc[ii].dist > dist)&&(sortInd==-1))
        {
            sortInd=ii;
        }
        // searching if we already have a motif in out top K list which is near to the currently good match
        if ((fabs(subSeqPtr[priorityQDisc[ii].ind1].sTime - subSeqPtr[ind1].sTime) < blackDur) || (fabs(subSeqPtr[priorityQDisc[ii].ind2].sTime -subSeqPtr[ind1].sTime) < blackDur) || (fabs(subSeqPtr[priorityQDisc[ii].ind1].sTime-subSeqPtr[ind2].sTime) < blackDur) || (fabs(subSeqPtr[priorityQDisc[ii].ind2].sTime-subSeqPtr[ind2].sTime) < blackDur))
        {
            matchInd=ii;
            break;
        }
    }
    if (sortInd==-1)//we couldn't get any satisfactory replacement before we get a close neighbour
    {
        return priorityQDisc[K-1].dist;
    }
    //There are three possibilities
    //1) There is no match found in the existing top motifs, simplest
    if (matchInd==-1)
    {
        memmove(&priorityQDisc[sortInd+1], &priorityQDisc[sortInd], sizeof(TSAmotifInfo_t)*(K-(sortInd+1)));
        priorityQDisc[sortInd].dist = dist;
        priorityQDisc[sortInd].ind1 = ind1;
        priorityQDisc[sortInd].ind2 = ind2;
    }
    else if (sortInd == matchInd)
    {
        priorityQDisc[sortInd].dist = dist;
        priorityQDisc[sortInd].ind1 = ind1;
        priorityQDisc[sortInd].ind2 = ind2;
    }
    else if (sortInd < matchInd)
    {
        memmove(&priorityQDisc[sortInd+1], &priorityQDisc[sortInd], sizeof(TSAmotifInfo_t)*(matchInd-sortInd));
        priorityQDisc[sortInd].dist = dist;
        priorityQDisc[sortInd].ind1 = ind1;
        priorityQDisc[sortInd].ind2 = ind2;
    }
    
    return priorityQDisc[K-1].dist;
    
}



TSADIST TSApool::managePriorityQSear(TSAIND queryInd, TSAsubSeq_t *subSeqPtr, TSAIND ind1, TSAIND ind2, TSADIST dist, int searchFileID)
{
    int sortInd=-1, matchInd=-1;
    
    for(int ii=0;ii<K;ii++)
    {
        if ((priorityQSear[queryInd][ii].dist > dist)&&(sortInd==-1))
        {
            sortInd=ii;
        }
    }
    
    
    for(int ii=0;ii<K;ii++)
    {
        if ((priorityQSear[queryInd][ii].dist > dist)&&(sortInd==-1))
        {
            sortInd=ii;
        }

        // searching if we already have a motif in out top K list which is near to the currently good match
        if ((priorityQSear[queryInd][ii].searchFileID==searchFileID) && (fabs(priorityQSear[queryInd][ii].sTime -subSeqPtr[ind2].sTime) < blackDur))
        {
            matchInd=ii;
            break;
        }
    }
    
    if (sortInd==-1)//we couldn't get any satisfactory replacement before we get a close neighbour
    {
        return priorityQSear[queryInd][K-1].dist;
    }
    //There are threbsfe possibilities
    //1) There is no match found in the existing top motifs, simplest
    if (matchInd==-1)
    {
        memmove(&priorityQSear[queryInd][sortInd+1], &priorityQSear[queryInd][sortInd], sizeof(TSAmotifInfoExt_t)*(K-(sortInd+1)));
        priorityQSear[queryInd][sortInd].dist = dist;
        priorityQSear[queryInd][sortInd].ind1 = ind1;
        priorityQSear[queryInd][sortInd].ind2 = ind2;
        priorityQSear[queryInd][sortInd].sTime = subSeqPtr[ind2].sTime;
        priorityQSear[queryInd][sortInd].eTime = subSeqPtr[ind2].eTime;
        priorityQSear[queryInd][sortInd].searchFileID = searchFileID;
        priorityQSear[queryInd][sortInd].patternID = PID_DEFAULT3;
    }
    else if (sortInd == matchInd)
    {
        priorityQSear[queryInd][sortInd].dist = dist;
        priorityQSear[queryInd][sortInd].ind1 = ind1;
        priorityQSear[queryInd][sortInd].ind2 = ind2;
        priorityQSear[queryInd][sortInd].sTime = subSeqPtr[ind2].sTime;
        priorityQSear[queryInd][sortInd].eTime = subSeqPtr[ind2].eTime;
        priorityQSear[queryInd][sortInd].searchFileID = searchFileID;
        priorityQSear[queryInd][sortInd].patternID = PID_DEFAULT3;
    }
    else if (sortInd < matchInd)
    {
        memmove(&priorityQSear[queryInd][sortInd+1], &priorityQSear[queryInd][sortInd], sizeof(TSAmotifInfoExt_t)*(matchInd-sortInd));
        priorityQSear[queryInd][sortInd].dist = dist;
        priorityQSear[queryInd][sortInd].ind1 = ind1;
        priorityQSear[queryInd][sortInd].ind2 = ind2;
        priorityQSear[queryInd][sortInd].sTime = subSeqPtr[ind2].sTime;
        priorityQSear[queryInd][sortInd].eTime = subSeqPtr[ind2].eTime;
        priorityQSear[queryInd][sortInd].searchFileID = searchFileID;
        priorityQSear[queryInd][sortInd].patternID = PID_DEFAULT3;
    }
    
    return priorityQSear[queryInd][K-1].dist;
    
}

