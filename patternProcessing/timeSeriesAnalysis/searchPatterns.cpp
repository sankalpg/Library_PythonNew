
#include "TSAdataIO.h"
#include "TSAsimilarity.h"
#include "TSApool.h"
#include "TSAlogs.h"


using namespace std;


int main( int argc , char *argv[])
{
    procParams_t *myProcParamsPtr;
    fileExts_t *myFileExtsPtr;
    int Err=0;
    int verbos;
    
    TSAparamHandle paramHand;
    TSAlogs logs;
    fileNameHandler fHandle;
    
    
    //checking if the number of input arguments are correct 
    if(argc < 6 || argc > 7)
    {
        printf("\nInvalid number of arguments!!!\n");
        exit(1);
    }
    //reading commandline parameters
    char *baseName = argv[1];
    char *paramFile = argv[2];
    char *fileExtFile = argv[3];
    int kNN = atoi(argv[4]);
    TSADIST distTsld = atof(argv[5]); 
    if( argc == 7 ){verbos = atoi(argv[6]);}
    
    //read params from the paramFile
    paramHand.readParamsFromFile(paramFile);
    myProcParamsPtr = paramHand.getParamPtr();
    
    //read file extensions from the file
    paramHand.readFileExtsInfoFile(fileExtFile);
    myFileExtsPtr = paramHand.getExtPtr();
    
    
    fHandle.initialize(baseName, myFileExtsPtr);
    
    //create a data handler object
    TSAdataHandler *TSData1 = new TSAdataHandler(baseName, &logs.procLogs, myFileExtsPtr, myProcParamsPtr);
    TSData1->loadMotifDataTemplate1();
    
    TSAdtwSimilarity dtwUCR;
    
    int lenMotifReal = TSData1->procParams.motifLengths[TSData1->procParams.indexMotifLenReal];
    int nInterFact = TSData1->procParams.nInterpFac;
    dtwUCR.configureTSASimilarity(lenMotifReal, lenMotifReal, TSData1->procParams.DTWBand);
    
    dtwUCR.setQueryPtr(TSData1->subSeqPtr, TSData1->nSubSeqs);
    dtwUCR.computeQueryEnvelops();
    dtwUCR.initArrayBSF(ceil(TSData1->nSubSeqs/nInterFact));
    //dtwUCR.copyQueryEnv2Cand();
    
    TSADIST LB_kim_FL, LB_Keogh_EQ, LB_Keogh_EC, realDist, bsf=FLT_MAX;
    TSApool pool(kNN, myProcParamsPtr->blackDur);
    pool.initPriorityQSear(ceil(TSData1->nSubSeqs/nInterFact));
    pool.initPattStorage(ceil(TSData1->nSubSeqs/nInterFact), lenMotifReal);
    
    // iterating over all files 
    FILE *fp2 = fopen(fHandle.getMappFileName(),"w");
    int searchFileID=0;
    
    fHandle.loadSearchFileList();
    TSAIND queryInd=0;
    
    for(TSAIND ss=0; ss < fHandle.nSearchFiles; ss++)
    {
        TSAdataHandler *TSData2 = new TSAdataHandler(fHandle.searchFileNames[ss], &logs.procLogs, myFileExtsPtr, myProcParamsPtr);
        TSData2->genTemplate1SubSeqs();
        dtwUCR.setCandPtr(TSData2->subSeqPtr, TSData2->nSubSeqs);
        dtwUCR.computeCandEnvelops();
        
        searchFileID = ss;
        
        for(TSAIND jj=0;jj< TSData2->nSubSeqs;jj++)
        {
            for(TSAIND ii=0;ii< TSData1->nSubSeqs;ii++)
            {
                queryInd = (TSAIND)floor(ii/nInterFact);
                
                if (paramHand.procParams.combMTX[ii%nInterFact][jj%nInterFact]==0)
                    continue;

                if ((strcmp(baseName, fHandle.searchFileNames[ss])==0)&& (fabs(TSData1->subSeqPtr[ii].sTime-TSData2->subSeqPtr[jj].sTime)< TSData1->procParams.blackDur))
                    //beware that basename and searchFile name should both have either full path or relative path.
                {
                    continue;
                }
                LB_kim_FL = computeLBkimFL(TSData1->subSeqPtr[ii].pData[0], TSData2->subSeqPtr[jj].pData[0], TSData1->subSeqPtr[ii].pData[lenMotifReal-1], TSData2->subSeqPtr[jj].pData[lenMotifReal-1], SqEuclidean);
                if (LB_kim_FL< dtwUCR.bsfArray[queryInd]) 
                {
                    LB_Keogh_EQ = computeKeoghsLB(dtwUCR.envUQueryPtr[ii],dtwUCR.envLQueryPtr[ii],dtwUCR.accLB_Keogh_EQ, TSData2->subSeqPtr[jj].pData,lenMotifReal, dtwUCR.bsfArray[queryInd], SqEuclidean);
                    if(LB_Keogh_EQ < dtwUCR.bsfArray[queryInd])
                    {
                        LB_Keogh_EC = computeKeoghsLB(dtwUCR.envUCandPtr[jj],dtwUCR.envLCandPtr[jj],dtwUCR.accLB_Keogh_EC, TSData1->subSeqPtr[ii].pData,lenMotifReal, dtwUCR.bsfArray[queryInd], SqEuclidean);
                        if(LB_Keogh_EC < dtwUCR.bsfArray[queryInd])
                        {
                            realDist = dtw1dBandConst(TSData1->subSeqPtr[ii].pData, TSData2->subSeqPtr[jj].pData, lenMotifReal, lenMotifReal, dtwUCR.costMTX, SqEuclidean, dtwUCR.bandDTW, dtwUCR.bsfArray[queryInd], dtwUCR.accLB_Keogh_EQ);
                            if (realDist <= dtwUCR.bsfArray[queryInd])
                            {
                                dtwUCR.bsfArray[queryInd] = pool.managePriorityQSear(queryInd, TSData2->subSeqPtr, ii, jj, realDist, searchFileID);
                            }
                        }
                    }
                }
            }
        }
        for(TSAIND jj=0;jj< ceil(TSData1->nSubSeqs/nInterFact);jj++)
        {
            pool.updatePattStorageData(jj, TSData2->subSeqPtr, lenMotifReal, searchFileID);
        }
        delete TSData2;
        dtwUCR.deleteCandEnvMem();
        
    }
    fclose(fp2);
    
    //lets do rank refinement
    //############## Rank Refinement using sophisticated DTW ########################
    for (int mm=0;mm<paramHand.procParams.nSimMeasuresUsed;mm++)
    {
            // Since there can be multiple similarity measure used for rank refinement (mainly during experiment phase) this rank refinement step should be in loop, no need to loop rest of the steps
            
            //recomputing the distance
            for(TSAIND ii=0;ii<ceil(TSData1->nSubSeqs/nInterFact);ii++)
            {
                for(TSAIND jj=0;jj<pool.K;jj++)
                {
                    if (pool.priorityQSear[ii][jj].dist < INF)   //do refinement only for a valid top entry, leave the infinites!!
                    {
                        pool.priorityQSear[ii][jj].dist = dtw1dBandConst_localConst(TSData1->subSeqPtr[pool.priorityQSear[ii][jj].ind1].pData, pool.priorityQSear[ii][jj].storagePtr->data, lenMotifReal, lenMotifReal, dtwUCR.costMTX, paramHand.procParams.simMeasureRankRefinement[mm], dtwUCR.bandDTW, INF, dtwUCR.accLB_Keogh_EQ);
                    }
                    else
                    {
                        if (verbos)
                        {
                            printf("There is some serious problem in rank refinement step %lld,%lld",jj,ii);
                        }
                    }
                    
                }
                //sorting the priority list
                pool.sortQSearch(ii);
            }   
                
            
            TSData1->dumpSearMotifInfo(fHandle.getOutFileNamePostRR(paramHand.procParams.simMeasureRankRefinement[mm]), pool.priorityQSear, TSData1->nSubSeqs/nInterFact, pool.K, verbos);
    
    }    
    //TSData1->dumpDiscMotifInfo(TSData1->fHandle.getOutFileName(), pool.priorityQDisc, pool.K, verbos);
    //TSData1->dumpSearMotifInfo(fHandle.getOutFileName(), pool.priorityQSear, TSData1->nSubSeqs/nInterFact, pool.K, verbos);
    //generate pattern sub sequences
    
    
    //prepare lower bounding data
    
    //enter into nested loop
    
    //free memories
    
    //dump the generate data
    delete TSData1;
    
    
    if (verbos){printf("Processing done!\n");}
    return 1;
}

