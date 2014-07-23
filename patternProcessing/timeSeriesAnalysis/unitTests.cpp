
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
    TSAdataHandler TSData1(baseName, &logs.procLogs, myFileExtsPtr, myProcParamsPtr);
    
    TSData1.loadMotifDataTemplate1();
    printf("Hello10\n");
    
    TSAsubSeq_t *subSeqPtr = TSData1.subSeqPtr;
    
    
    TSAdtwSimilarity dtwUCR;
    
    int lenMotifReal = TSData1.procParams.motifLengths[TSData1.procParams.indexMotifLenReal];
    dtwUCR.configureTSASimilarity(lenMotifReal, lenMotifReal, TSData1.procParams.DTWBand);
    
    printf("Hello11\n");
    
    dtwUCR.setQueryPtr(TSData1.subSeqPtr, TSData1.nSubSeqs);
    dtwUCR.computeQueryEnvelops();
    //dtwUCR.copyQueryEnv2Cand();
    
    printf("Hello12\n");
    
    int nInterFact = TSData1.procParams.nInterpFac;
    TSADIST LB_kim_FL, LB_Keogh_EQ, LB_Keogh_EC, realDist, bsf=FLT_MAX;
    
    TSApool pool(kNN, myProcParamsPtr->blackDur);
    pool.initPriorityQSear(ceil(TSData1.nSubSeqs/nInterFact));
    printf("you have chosen this KNN %d\n", pool.K);
    
    printf("Hello13\n");
    
    // iterating over all files 
    FILE *fp2 = fopen(fHandle.getMappFileName(),"w");
    int searchFileID=0;
    
    fHandle.loadSearchFileList();
    TSAIND queryInd=0;
    
    for(TSAIND ss=0; ss < fHandle.nSearchFiles; ss++)
    {
        TSAdataHandler TSData2(fHandle.searchFileNames[ss], &logs.procLogs, myFileExtsPtr, myProcParamsPtr);
        TSData2.genTemplate1SubSeqs();
        dtwUCR.setCandPtr(TSData2.subSeqPtr, TSData2.nSubSeqs);
        dtwUCR.computeCandEnvelops();
        
        searchFileID = ss;
        
        for(TSAIND jj=0;jj< TSData2.nSubSeqs;jj++)
        {
            for(TSAIND ii=0;ii< TSData1.nSubSeqs;ii++)
            {
                queryInd = (TSAIND)floor(queryInd/nInterFact);
                if (TSData1.procParams.combMTX[ii%nInterFact][jj%nInterFact]==0)
                    continue;
                if (fabs(subSeqPtr[ii].sTime-subSeqPtr[jj].sTime)< TSData1.procParams.blackDur)
                {
                    continue;
                }
                if ((strcmp(baseName, fHandle.searchFileNames[ss])==0)&& (fabs(TSData1.subSeqPtr[ii].sTime-TSData2.subSeqPtr[jj].sTime)< TSData1.procParams.blackDur))
                    //beware that basename and searchFile name should both have either full path or relative path.
                {
                    continue;
                }
                LB_kim_FL = computeLBkimFL(TSData1.subSeqPtr[ii].pData[0], TSData2.subSeqPtr[jj].pData[0], TSData1.subSeqPtr[ii].pData[lenMotifReal-1], TSData2.subSeqPtr[jj].pData[lenMotifReal-1], SqEuclidean);
                if (LB_kim_FL< bsf) 
                {
                    LB_Keogh_EQ = computeKeoghsLB(dtwUCR.envUQueryPtr[ii],dtwUCR.envLQueryPtr[ii],dtwUCR.accLB_Keogh_EQ, TSData2.subSeqPtr[jj].pData,lenMotifReal, bsf, SqEuclidean);
                    if(LB_Keogh_EQ < bsf)
                    {
                        LB_Keogh_EC = computeKeoghsLB(dtwUCR.envUCandPtr[jj],dtwUCR.envLCandPtr[jj],dtwUCR.accLB_Keogh_EC, TSData1.subSeqPtr[ii].pData,lenMotifReal, bsf, SqEuclidean);
                        if(LB_Keogh_EC < bsf)
                        {
                            realDist = dtw1dBandConst(TSData1.subSeqPtr[ii].pData, TSData2.subSeqPtr[jj].pData, lenMotifReal, lenMotifReal, dtwUCR.costMTX, SqEuclidean, dtwUCR.bandDTW, bsf, dtwUCR.accLB_Keogh_EQ);
                            if (realDist <= bsf)
                            {
                                bsf = pool.managePriorityQSear(queryInd, subSeqPtr, ii, jj, realDist, searchFileID);
                            }
                        }
                    }
                }
            }
        }
    }
    fclose(fp2);
        
    
    //TSData1.dumpDiscMotifInfo(TSData1.fHandle.getOutFileName(), pool.priorityQDisc, pool.K, verbos);
    TSData1.dumpSearMotifInfo(fHandle.getOutFileName(), pool.priorityQSear, TSData1.nSubSeqs/nInterFact, pool.K, verbos);
    //generate pattern sub sequences
    
    
    //prepare lower bounding data
    
    //enter into nested loop
    
    //free memories
    
    //dump the generate data
    
    
    if (verbos){printf("Processing done!\n");}
    return 1;
}

