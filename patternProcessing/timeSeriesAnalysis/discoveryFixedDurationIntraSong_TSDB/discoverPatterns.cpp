
#include "../TSAdataIO.h"
#include "../TSAsimilarity.h"
#include "../TSApool.h"
#include "../TSAlogs.h"


using namespace std;
float t1,t2, ti,tf;


int main( int argc , char *argv[])
{
    procParams_t *myProcParamsPtr;
    fileExts_t *myFileExtsPtr;
    int Err=0;
    int verbos;
    int rVal=0;
    float offset=0;
    double *offseted_data;
    
    TSAparamHandle paramHand;
    TSAlogs logs;
    
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
    
    //create a data handler object
    TSAdataHandler *TSData1 = new TSAdataHandler(baseName, &logs.procLogs, myFileExtsPtr, myProcParamsPtr);
    
    rVal = TSData1->genTemplate1SubSeqs();
    if (rVal==0)
    {
        return 0;
    }
    
    TSAsubSeq_t *subSeqPtr = TSData1->subSeqPtr;
    
    
    TSAdtwSimilarity dtwUCR( &logs.procLogs);
    
    int lenMotifReal = TSData1->procParams.motifLengths[TSData1->procParams.indexMotifLenReal];
    dtwUCR.configureTSASimilarity(lenMotifReal, lenMotifReal, TSData1->procParams.distParams.DTWBand);
    
    dtwUCR.setQueryPtr(TSData1->subSeqPtr, TSData1->nSubSeqs);
    dtwUCR.computeQueryEnvelops();
    dtwUCR.copyQueryEnv2Cand();
    
    int nInterFact = TSData1->procParams.pattParams.nInterpFac;
    TSADIST LB_kim_FL, LB_Keogh_EQ, LB_Keogh_EC, realDist, bsf=FLT_MAX;
    TSADATA **U, **L, *accLB1, *accLB2;
    
    U = dtwUCR.envUQueryPtr;
    L = dtwUCR.envLQueryPtr;
    accLB1 = dtwUCR.accLB_Keogh_EQ;
    accLB2 = dtwUCR.accLB_Keogh_EC;
    
    TSApool *pool = new TSApool(kNN);
    pool->initPriorityQDisc();

    offseted_data = (double*)malloc(sizeof(double)*lenMotifReal);
    
    for(TSAIND ii=0;ii< TSData1->nSubSeqs;ii++)
    {
        for(TSAIND jj=ii+1;jj< TSData1->nSubSeqs;jj++)
        {
            if (paramHand.procParams.combMTX[ii%nInterFact][jj%nInterFact]==0)
                continue;
            if (fabs(subSeqPtr[ii].sTime-subSeqPtr[jj].sTime)< TSData1->procParams.pattParams.blackDur)
            {
                continue;
            }
            offset = TSData1->estimateOffset(subSeqPtr[ii].mean, subSeqPtr[jj].mean);
            //printf("Mean1, mean2 and offset are: %f\t%f\t%f\n", subSeqPtr[ii].mean, subSeqPtr[jj].mean, offset);
            LB_kim_FL = computeLBkimFL(subSeqPtr[ii].pData[0], subSeqPtr[jj].pData[0] + offset, subSeqPtr[ii].pData[lenMotifReal-1], subSeqPtr[jj].pData[lenMotifReal-1]+ offset, SqEuclidean);
            logs.procLogs.nLB_KIM_FL++;
            if (LB_kim_FL< bsf) 
            {
                LB_Keogh_EQ = computeKeoghsLB(U[ii],L[ii],accLB1, subSeqPtr[jj].pData,lenMotifReal, bsf, SqEuclidean, offset);
                logs.procLogs.nLB_Keogh_EQ++;
                if(LB_Keogh_EQ < bsf)
                {
                    LB_Keogh_EC = computeKeoghsLB(U[jj],L[jj],accLB2, subSeqPtr[ii].pData,lenMotifReal, bsf, SqEuclidean, -1*offset);
                    logs.procLogs.nLB_Keogh_EC++;
                    if(LB_Keogh_EC < bsf)
                    {
                        if ((offset !=0)&&(TSData1->procParams.repParams.normType == OCTAVE_NORM)){
                            for (int zz=0;zz<lenMotifReal;zz++){offseted_data[zz] = subSeqPtr[jj].pData[zz]+offset;}    //saving the offsetted copy
                            realDist = dtw1dBandConst(subSeqPtr[ii].pData, offseted_data, lenMotifReal, lenMotifReal, dtwUCR.costMTX, SqEuclidean, dtwUCR.bandDTW, bsf, accLB1);
                        }
                        else{
                            realDist = dtw1dBandConst(subSeqPtr[ii].pData, subSeqPtr[jj].pData, lenMotifReal, lenMotifReal, dtwUCR.costMTX, SqEuclidean, dtwUCR.bandDTW, bsf, accLB1);
                        }


                        logs.procLogs.nDTW_EA++;
                        if (realDist <= bsf)
                        {
                            bsf = pool->managePriorityQDisc(subSeqPtr, ii, jj, realDist, TSData1->procParams.pattParams.blackDur);
                            logs.procLogs.nPriorityUpdates++;
                        }
                    }
                }
            }
        }
    }
    tf = clock();
    
    logs.procLogs.tTotal += (tf-ti)/CLOCKS_PER_SEC;
    
    t1 = clock();
    
    TSData1->dumpDiscMotifInfo(TSData1->fHandle.getOutFileName(), pool->priorityQDisc, pool->K, verbos);
    
    t2 = clock();
    logs.procLogs.tDump += (t2-t1)/CLOCKS_PER_SEC;
    
    
    
    //generate pattern sub sequences
    
    
    //prepare lower bounding data
    
    //enter into nested loop
    
    //free memories
    
    //dump the generate data
    
    logs.dumpProcLogs(TSData1->fHandle.getLogFileName(), verbos);
    
    delete TSData1;
    delete pool;
    if (verbos){printf("Processing done!\n");}
    
    
    
    
    return 1;
}


