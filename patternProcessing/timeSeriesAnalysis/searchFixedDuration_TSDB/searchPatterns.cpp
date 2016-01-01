
#include "../TSAdataIO.h"
#include "../TSAsimilarity.h"
#include "../TSApool.h"
#include "../TSAlogs.h"


using namespace std;

float t1,t2, ti,tf;


int main( int argc , char *argv[])
{
    TSAparamHandle paramHand;
    TSAlogs logs;
    fileNameHandler fHandle;
    procParams_t *myProcParamsPtr;
    fileExts_t *myFileExtsPtr;
    int Err=0, verbos;
    float offset=0, mean2;
    double *offseted_data;
    
    //checking if the number of input arguments are correct 
    if(argc < 6 || argc > 7)
    {
        printf("\nInvalid number of arguments!!!\n");
        exit(1);
    }
    
    ti = clock();
    
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
    
    TSAdtwSimilarity dtwUCR( &logs.procLogs);
    
    int lenMotifReal = TSData1->procParams.motifLengths[TSData1->procParams.indexMotifLenReal];
    int nInterFact = TSData1->procParams.pattParams.nInterpFac;
    dtwUCR.configureTSASimilarity(lenMotifReal, lenMotifReal, TSData1->procParams.distParams.DTWBand);
    
    dtwUCR.setQueryPtr(TSData1->subSeqPtr, TSData1->nSubSeqs);
    dtwUCR.computeQueryEnvelops();
    dtwUCR.initArrayBSF(ceil(TSData1->nSubSeqs/nInterFact));
    //dtwUCR.copyQueryEnv2Cand();
    
    TSADIST LB_kim_FL, LB_Keogh_EQ, LB_Keogh_EC, realDist, bsf=FLT_MAX;
    TSApool pool(kNN);
    pool.initPriorityQSear(ceil(TSData1->nSubSeqs/nInterFact));
    pool.initPattStorage(ceil(TSData1->nSubSeqs/nInterFact), lenMotifReal);
    
    // iterating over all files 
    FILE *fp2 = fopen(fHandle.getMappFileName(),"w");
    int searchFileID=0;
    
    fHandle.loadSearchFileList();
    TSAIND queryInd=0;

    offseted_data = (double*)malloc(sizeof(double)*lenMotifReal);
    
    for(TSAIND ss=0; ss < fHandle.nSearchFiles; ss++)
    {
        searchFileID = ss;
        fprintf(fp2, "%d\t%s\n",searchFileID, fHandle.searchFileNames[ss]);
        TSAdataHandler *TSData2 = new TSAdataHandler(fHandle.searchFileNames[ss], &logs.procLogs, myFileExtsPtr, myProcParamsPtr);
        TSData2->genTemplate1SubSeqs();
        dtwUCR.setCandPtr(TSData2->subSeqPtr, TSData2->nSubSeqs);
        dtwUCR.computeCandEnvelops();
        
        
        
        for(TSAIND jj=0;jj< TSData2->nSubSeqs;jj++)
        {
            for(TSAIND ii=0;ii< TSData1->nSubSeqs;ii++)
            {
                queryInd = (TSAIND)floor(ii/nInterFact);
                
                if (paramHand.procParams.combMTX[ii%nInterFact][jj%nInterFact]==0)
                    continue;
                //Criterion for determining overlapping phrases is changed after flat note compression integration
                if ((strcmp(baseName, fHandle.searchFileNames[ss])==0)&& ((TSData1->subSeqPtr[ii].eTime- TSData2->subSeqPtr[jj].sTime)*(TSData2->subSeqPtr[jj].eTime-TSData1->subSeqPtr[ii].sTime) > 0))
                    //beware that basename and searchFile name should both have either full path or relative path.
                {
                    continue;
                }
                offset = TSData1->estimateOffset(TSData1->subSeqPtr[ii].mean, TSData2->subSeqPtr[jj].mean);
                LB_kim_FL = computeLBkimFL(TSData1->subSeqPtr[ii].pData[0], TSData2->subSeqPtr[jj].pData[0] + offset, TSData1->subSeqPtr[ii].pData[lenMotifReal-1], TSData2->subSeqPtr[jj].pData[lenMotifReal-1] + offset, SqEuclidean);
                logs.procLogs.nLB_KIM_FL++;
                if (LB_kim_FL< dtwUCR.bsfArray[queryInd]) 
                {
                    LB_Keogh_EQ = computeKeoghsLB(dtwUCR.envUQueryPtr[ii],dtwUCR.envLQueryPtr[ii],dtwUCR.accLB_Keogh_EQ, TSData2->subSeqPtr[jj].pData,lenMotifReal, dtwUCR.bsfArray[queryInd], SqEuclidean, offset);
                    logs.procLogs.nLB_Keogh_EQ++;
                    if(LB_Keogh_EQ < dtwUCR.bsfArray[queryInd])
                    {
                        LB_Keogh_EC = computeKeoghsLB(dtwUCR.envUCandPtr[jj],dtwUCR.envLCandPtr[jj],dtwUCR.accLB_Keogh_EC, TSData1->subSeqPtr[ii].pData,lenMotifReal, dtwUCR.bsfArray[queryInd], SqEuclidean, -1*offset);
                        logs.procLogs.nLB_Keogh_EC++;
                        if(LB_Keogh_EC < dtwUCR.bsfArray[queryInd])
                        {
                            if ((offset !=0)&&(TSData1->procParams.repParams.normType == OCTAVE_NORM)){
                                for (int zz=0;zz<lenMotifReal;zz++){offseted_data[zz] = TSData2->subSeqPtr[jj].pData[zz]+offset;}    //saving the offsetted copy
                                realDist = dtw1dBandConst(TSData1->subSeqPtr[ii].pData, offseted_data, lenMotifReal, lenMotifReal, dtwUCR.costMTX, SqEuclidean, dtwUCR.bandDTW, dtwUCR.bsfArray[queryInd], dtwUCR.accLB_Keogh_EQ);
                            }
                            else{
                                realDist = dtw1dBandConst(TSData1->subSeqPtr[ii].pData, TSData2->subSeqPtr[jj].pData, lenMotifReal, lenMotifReal, dtwUCR.costMTX, SqEuclidean, dtwUCR.bandDTW, dtwUCR.bsfArray[queryInd], dtwUCR.accLB_Keogh_EQ);
                            }
                            logs.procLogs.nDTW_EA++;
                            if (realDist <= dtwUCR.bsfArray[queryInd])
                            {
                                dtwUCR.bsfArray[queryInd] = pool.managePriorityQSear(queryInd, TSData2->subSeqPtr, ii, jj, realDist, searchFileID);
                                logs.procLogs.nPriorityUpdates++;
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
                        mean2 = computeMean(pool.priorityQSear[ii][jj].storagePtr->data, lenMotifReal);
                        offset = TSData1->estimateOffset(TSData1->subSeqPtr[pool.priorityQSear[ii][jj].ind1].mean, mean2);
                        if ((offset !=0)&&(TSData1->procParams.repParams.normType == OCTAVE_NORM)){
                                for (int zz=0;zz<lenMotifReal;zz++){offseted_data[zz] = pool.priorityQSear[ii][jj].storagePtr->data[zz]+offset;}    //saving the offsetted copy
                                pool.priorityQSear[ii][jj].dist = dtw1dBandConst_localConst(TSData1->subSeqPtr[pool.priorityQSear[ii][jj].ind1].pData, offseted_data, lenMotifReal, lenMotifReal, dtwUCR.costMTX, paramHand.procParams.simMeasureRankRefinement[mm], dtwUCR.bandDTW, INF, dtwUCR.accLB_Keogh_EQ);
                        }
                        else{
                            pool.priorityQSear[ii][jj].dist = dtw1dBandConst_localConst(TSData1->subSeqPtr[pool.priorityQSear[ii][jj].ind1].pData, pool.priorityQSear[ii][jj].storagePtr->data, lenMotifReal, lenMotifReal, dtwUCR.costMTX, paramHand.procParams.simMeasureRankRefinement[mm], dtwUCR.bandDTW, INF, dtwUCR.accLB_Keogh_EQ);
                        }
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
            t1 = clock();

            TSData1->dumpSearMotifInfo(fHandle.getOutFileNamePostRR(paramHand.procParams.simMeasureRankRefinement[mm]), pool.priorityQSear, TSData1->nSubSeqs/nInterFact, pool.K, verbos);
            
            t2 = clock();
            logs.procLogs.tDump += (t2-t1)/CLOCKS_PER_SEC;
    
    }
    
    tf = clock();
    
    logs.procLogs.tTotal += (tf-ti)/CLOCKS_PER_SEC;
    
    
    logs.dumpProcLogs(TSData1->fHandle.getFileName(TSData1->fHandle.fileExtPtr->srchLogFileExt), verbos);
    
    delete TSData1;
    
    
    if (verbos){printf("Processing done!\n");}
    return 1;
}

