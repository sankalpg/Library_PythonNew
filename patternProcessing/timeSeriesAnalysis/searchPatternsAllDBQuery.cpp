
#include "TSAdataIO.h"
#include "TSAsimilarity.h"
#include "TSApool.h"
#include "TSAlogs.h"


using namespace std;

typedef double (*distFunc)(double*, double*, int, int, double**, int, int, double, double*);

int main( int argc , char *argv[])
{
    int Err=0;
    int verbos;
    TSADIST realDist, LB_kim_FL, LB_Keogh_EQ, LB_Keogh_EC;
    int distType=-1;
    
    distFunc myDist[3]={NULL};
    myDist[0] = &euclideanSeq;
    myDist[1] = &dtw1dBandConst;
    myDist[2] = &dtw1dBandConst_localConst;
    
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
    
    FILE *fp10;
    fp10 = fopen(fHandleTemp.getOutFileName(), "w");
    fclose(fp10);
    
    int nInterFact = myProcParamsPtr->pattParams.nInterpFac;
    float blackDur=0;
    
    
    if (myProcParamsPtr->distParams.distType==0)
    {
        distType=0;
    }
    else if (myProcParamsPtr->distParams.distType==1)
    {
        if (myProcParamsPtr->distParams.DTWType==0)
        {
            distType=1;
        }
        else if (myProcParamsPtr->distParams.DTWType==1)
        {
            distType=2;
        }

    }
    
    for(int ff1=0; ff1< fHandleTemp.nSearchFiles; ff1++)
    {
        //create a data handler object
        TSAdataHandler *TSData1 = new TSAdataHandler(fHandleTemp.searchFileNames[ff1], &logs.procLogs, myFileExtsPtr, myProcParamsPtr);
		TSData1->fHandle.loadSearchFileList();
        TSData1->readTSData(TSData1->fHandle.getTSFileName());
        TSData1->readHopSizeTS(TSData1->fHandle.getTSFileName());
        TSData1->downSampleTS();
        TSData1->convertHz2Cents(TSData1->fHandle.getTonicFileName());
        TSData1->quantizeSampleTS(TSData1->procParams.repParams.TSRepType);
        TSData1->readQueryTimeStamps(TSData1->fHandle.getQueryFileName(), MY_MOTIF_ANNOT_FORMAT);
        
        TSAIND qq=0;
        int searchFileID=0;
        TSApool *pool = new TSApool(kNN);
        pool->initPriorityQSear(TSData1->nQueries);
        
        TSAdtwSimilarity *dtwUCRTemp = new TSAdtwSimilarity(&logs.procLogs);    //only for managing bsf array
        dtwUCRTemp->configureTSASimilarity(1, 1, myProcParamsPtr->distParams.DTWBand);
        dtwUCRTemp->initArrayBSF(TSData1->nQueries);
        
        for(int ff2=0; ff2< TSData1->fHandle.nSearchFiles; ff2++)
        {
            searchFileID = ff2;
            
            TSAdataHandler *TSData2 = new TSAdataHandler(TSData1->fHandle.searchFileNames[ff2], &logs.procLogs, myFileExtsPtr, myProcParamsPtr);
            TSData2->readTSData(TSData2->fHandle.getTSFileName());
            TSData2->readHopSizeTS(TSData2->fHandle.getTSFileName());
            TSData2->downSampleTS();
            TSData2->convertHz2Cents(TSData2->fHandle.getTonicFileName());
            TSData2->quantizeSampleTS(TSData2->procParams.repParams.TSRepType);
            
            for(int qq=0; qq < TSData1->nQueries; qq++)
            {
                TSData1->genSubSeqsWithTStamps(&TSData1->queryTStamps[qq], 1);
                TSData1->genUniScaledSubSeqs();
                TSData1->normalizeSubSeqs(TSData1->procParams.repParams.normType);
                
                TSData2->procParams.pattParams.durMotif = TSData1->procParams.pattParams.durMotif;
                blackDur = TSData1->procParams.pattParams.durMotif*TSData1->procParams.pattParams.blackDurFact; //blackDur is black dur factor w.r.t. motif length
                
                TSData2->calculateDiffMotifLengths();
                TSData2->genSlidingWindowSubSeqs();
                TSData2->genUniScaledSubSeqs();
                TSData2->normalizeSubSeqs(TSData2->procParams.repParams.normType);
                
                
                int lenMotifReal = TSData1->procParams.motifLengths[TSData1->procParams.indexMotifLenReal];
                
                TSAdtwSimilarity *dtwUCR = new TSAdtwSimilarity(&logs.procLogs);
                dtwUCR->configureTSASimilarity(lenMotifReal, lenMotifReal, myProcParamsPtr->distParams.DTWBand);
                
                dtwUCR->setQueryPtr(TSData1->subSeqPtr, TSData1->nSubSeqs);
                dtwUCR->setCandPtr(TSData2->subSeqPtr, TSData2->nSubSeqs);
                dtwUCR->computeQueryEnvelops();
                dtwUCR->computeCandEnvelops();
                
                for(TSAIND jj=0;jj< TSData2->nSubSeqs;jj++)
                {
                    for(TSAIND ii=0;ii< TSData1->nSubSeqs;ii++)
                    {

                        /*if (paramHand.procParams.combMTX[ii%nInterFact][jj%nInterFact]==0)
                            continue;
                        */
                        if ((strcmp(fHandleTemp.searchFileNames[ff1], TSData1->fHandle.searchFileNames[ff2])==0)&& (fabs(TSData1->subSeqPtr[ii].sTime-TSData2->subSeqPtr[jj].sTime)< blackDur))
                            //beware that basename and searchFile name should both have either full path or relative path.
                        {
                            continue;
                        }
                        
                        realDist = myDist[distType](TSData1->subSeqPtr[ii].pData, TSData2->subSeqPtr[jj].pData, lenMotifReal, lenMotifReal, dtwUCR->costMTX, SqEuclidean, dtwUCR->bandDTW, dtwUCRTemp->bsfArray[qq], dtwUCR->accLB_Keogh_EQ);
						dtwUCRTemp->bsfArray[qq] = pool->managePriorityQSear(qq, TSData2->subSeqPtr, ii, jj, realDist, searchFileID, blackDur);
                        
                    }
                }
                delete dtwUCR;
                TSData1->freeSubSeqsMem();
                TSData2->freeSubSeqsMem();
            }
            delete TSData2;
        }
        for (int qq=0; qq < TSData1->nQueries; qq++ )
        {
            {
                FILE *fp;
                fp = fopen(fHandleTemp.getOutFileName(), "ab");
                
                for(TSAIND ii=0;ii<kNN;ii++)
                {
                    fprintf(fp, "%d\t%d\t%d\t%f\t%f\t%f\n", ff1, qq, pool->priorityQSear[qq][ii].searchFileID, pool->priorityQSear[qq][ii].sTime, pool->priorityQSear[qq][ii].eTime, pool->priorityQSear[qq][ii].dist);
                }
                fclose(fp);
            }
        } 
        delete TSData1;
        delete pool;
        delete dtwUCRTemp;
    }
    
    if (verbos){printf("Processing done!\n");}
    return 1;
}

