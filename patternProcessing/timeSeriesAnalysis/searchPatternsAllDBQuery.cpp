
#include "TSAdataIO.h"
#include "TSAsimilarity.h"
#include "TSApool.h"
#include "TSAlogs.h"


using namespace std;


int main( int argc , char *argv[])
{
    int Err=0;
    int verbos;
    TSADIST realDist, LB_kim_FL, LB_Keogh_EQ, LB_Keogh_EC;
    
    
    
    procParams_t *myProcParamsPtr;
    fileExts_t *myFileExtsPtr;
    TSAparamHandle paramHand;
    TSAlogs logs;
    fileNameHandler fHandleTemp, fHandle;
    
    
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
    
    //lets use this to read the search file names
    fHandleTemp.initialize(baseName, myFileExtsPtr);
    fHandleTemp.loadSearchFileList();
    
    for(int ff=0; ff< fHandleTemp.nSearchFiles; ff++)
    {
        for(int ff=0; ff< fHandleTemp.nSearchFiles; ff++)
        {
            
            
        }
        
    }
    
    fHandle.searchFileNames[ii]
    
    //create a data handler object
    TSAdataHandler *TSData1 = new TSAdataHandler(baseName, &logs.procLogs, myFileExtsPtr, myProcParamsPtr);
    
    //preparing data1
    TSData1->readTSData(fHandle.getTSFileName());
    TSData1->readHopSizeTS(fHandle.getTSFileName());
    TSData1->downSampleTS();
    //TSData1->filterSamplesTS();
    TSData1->convertHz2Cents(fHandle.getTonicFileName());
    TSData1->readQueryTimeStamps(fHandle.getQueryFileName(), VIGNESH_MOTIF_ANNOT_FORMAT);
    
    fHandle.loadSearchFileList();
    TSAIND queryInd=0;
    
    TSApool pool(kNN, myProcParamsPtr->blackDur);
    pool.initPriorityQSear(TSData1->nQueries);
    
    // iterating over all files 
    FILE *fp2 = fopen(fHandle.getMappFileName(),"w");
    for(int ii=0;ii<fHandle.nSearchFiles;ii++)
    {
        fprintf(fp2, "%d\t%s\n", ii, fHandle.searchFileNames[ii]);
    }
    fclose(fp2);
    int searchFileID=0;
    
    FILE *fp = fopen(TSData1->fHandle.getOutFileName(), "w");
    fclose(fp);
    
    
    printf("Hello1\n");
    for (int qq=0; qq < TSData1->nQueries; qq++ )
    {
        printf("%lld\t%d\n", TSData1->nQueries, qq);
        printf("Hello1.1\n");
        TSData1->genSubSeqsWithTStamps(&TSData1->queryTStamps [qq], 1);
        printf("Hello1.2\n");
        TSData1->genUniScaledSubSeqs();
        printf("Hello1.3\n");
        printf("Hello1.4\n");
        TSData1->procParams.durMotif = TSData1->subSeqPtr[TSData1->procParams.indexMotifLenReal].eTime-TSData1->subSeqPtr[TSData1->procParams.indexMotifLenReal].sTime;
        printf("%f\t%f\n", TSData1->pHop, TSData1->procParams.durMotif);
        printf("Hello1.5\n");
        TSData1->calculateDiffMotifLengths();
        int lenMotifReal = TSData1->procParams.motifLengths[TSData1->procParams.indexMotifLenReal];
        int nInterFact = TSData1->procParams.nInterpFac;
        
        printf("Hello1.6\n");
        TSAdtwSimilarity *dtwUCR = new TSAdtwSimilarity();
        printf("Hello1.7\n");
        printf("motif len %d\n", lenMotifReal);
        dtwUCR->configureTSASimilarity(lenMotifReal, lenMotifReal, TSData1->procParams.DTWBand);
        printf("Hello1.8\n");
        dtwUCR->setQueryPtr(TSData1->subSeqPtr, TSData1->nSubSeqs);
        printf("Hello1.9\n");
        dtwUCR->computeQueryEnvelops();
        dtwUCR->initArrayBSF(ceil(TSData1->nSubSeqs/nInterFact));
        
        printf("Hello2\n");
        
        
        
        for(TSAIND ss=0; ss < fHandle.nSearchFiles; ss++)
        {
        
            TSAdataHandler *TSData2 = new TSAdataHandler(fHandle.searchFileNames[ss], &logs.procLogs, myFileExtsPtr, &(TSData1->procParams));
            //read the time series data    
            TSData2->readTSData(TSData2->fHandle.getTSFileName());
            TSData2->readHopSizeTS(TSData2->fHandle.getTSFileName());
            TSData2->downSampleTS();
            //remove silence pitch regions
            //TSData2->filterSamplesTS();
            TSData2->convertHz2Cents(TSData2->fHandle.getTonicFileName());
            TSData2->calculateDiffMotifLengths();
            printf("motif len %d\n", TSData2->procParams.motifLengths[TSData2->procParams.indexMotifLenReal]);
            TSData2->genSlidingWindowSubSeqs();
            TSData2->genUniScaledSubSeqs();
            
            dtwUCR->setCandPtr(TSData2->subSeqPtr, TSData2->nSubSeqs);
            dtwUCR->computeCandEnvelops();
            
            searchFileID = ss;
            printf("Hello3\n");
       
            for(TSAIND jj=0;jj< TSData2->nSubSeqs;jj++)
            {
                for(TSAIND ii=0;ii< TSData1->nSubSeqs;ii++)
                {
                    queryInd = (TSAIND)floor(ii/nInterFact);
                    
                    if (paramHand.procParams.combMTX[ii%nInterFact][jj%nInterFact]==0)
                        continue;

                    if (TSData2->subSeqPtr[jj].sTime>77.0)
                    {
                        TSData2->subSeqPtr[jj].sTime=TSData2->subSeqPtr[jj].sTime;
                    }
                    
                    if ((strcmp(baseName, fHandle.searchFileNames[ss])==0)&& (fabs(TSData1->subSeqPtr[ii].sTime-TSData2->subSeqPtr[jj].sTime)< TSData1->procParams.durMotif))
                        //beware that basename and searchFile name should both have either full path or relative path.
                    {
                        continue;
                    }

                    
                    LB_kim_FL = computeLBkimFL(TSData1->subSeqPtr[ii].pData[0], TSData2->subSeqPtr[jj].pData[0], TSData1->subSeqPtr[ii].pData[lenMotifReal-1], TSData2->subSeqPtr[jj].pData[lenMotifReal-1], SqEuclidean);
                    if (LB_kim_FL< dtwUCR->bsfArray[queryInd]) 
                    {
                        LB_Keogh_EQ = computeKeoghsLB(dtwUCR->envUQueryPtr[ii],dtwUCR->envLQueryPtr[ii],dtwUCR->accLB_Keogh_EQ, TSData2->subSeqPtr[jj].pData,lenMotifReal, dtwUCR->bsfArray[queryInd], SqEuclidean);
                        if(LB_Keogh_EQ < dtwUCR->bsfArray[queryInd])
                        {
                            LB_Keogh_EC = computeKeoghsLB(dtwUCR->envUCandPtr[jj],dtwUCR->envLCandPtr[jj],dtwUCR->accLB_Keogh_EC, TSData1->subSeqPtr[ii].pData,lenMotifReal, dtwUCR->bsfArray[queryInd], SqEuclidean);
                            if(LB_Keogh_EC < dtwUCR->bsfArray[queryInd])
                            {
                                realDist = dtw1dBandConst(TSData1->subSeqPtr[ii].pData, TSData2->subSeqPtr[jj].pData, lenMotifReal, lenMotifReal, dtwUCR->costMTX, SqEuclidean, dtwUCR->bandDTW, dtwUCR->bsfArray[queryInd], dtwUCR->accLB_Keogh_EQ);
                                if (realDist <= dtwUCR->bsfArray[queryInd])
                                {
                                    dtwUCR->bsfArray[queryInd] = pool.managePriorityQSear(qq, TSData2->subSeqPtr, ii, jj, realDist, searchFileID, TSData1->procParams.durMotif);
                                }
                            }
                        }
                    }
                }
            }
            dtwUCR->deleteCandEnvMem();
            /*if(qq==0)
            {FILE *fp;
            fp = fopen("tempData.bin", "wb");
            
            fwrite(TSData1->subSeqPtr[2].pData,sizeof(TSADATA), lenMotifReal, fp);
            fwrite(TSData2->subSeqPtr[2920].pData,sizeof(TSADATA), lenMotifReal, fp);
            fwrite(TSData2->subSeqPtr[24664].pData,sizeof(TSADATA), lenMotifReal, fp);
            fclose(fp);}*/
            delete TSData2;
            printf("Hello4\n");
        }
        {
            FILE *fp;
            fp = fopen(TSData1->fHandle.getOutFileName(), "ab");
            
            for(TSAIND ii=0;ii<kNN;ii++)
            {
                fprintf(fp, "%d\t%d\t%f\t%f\t%lld\t%f\n", qq, pool.priorityQSear[qq][ii].searchFileID, pool.priorityQSear[qq][ii].sTime, pool.priorityQSear[qq][ii].eTime,pool.priorityQSear[qq][ii].ind2, pool.priorityQSear[qq][ii].dist);
            }
            fclose(fp);
            
            
            
        }
        
        
        
        delete dtwUCR;
        if (qq <TSData1->nQueries-1)
        {TSData1->freeSubSeqsMem();}
        printf("Hello5\n");
    }
    printf("Hello6\n");
    
    delete TSData1;
    
    if (verbos){printf("Processing done!\n");}
    return 1;
}

