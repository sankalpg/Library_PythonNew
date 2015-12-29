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
    TSAparamHandle paramHand;
    TSAlogs logs;
    TSADIST LB_kim_FL, LB_Keogh_EQ, LB_Keogh_EC, realDist, bsf=FLT_MAX;
    TSADATA **U, **L, *accLB1, *accLB2;
    int Err=0, verbos, rVal=0;
    float offset=0;
    double *offseted_data;
    

    //checking if the number of input arguments are correct 
    if(argc < 5 || argc > 6)
    {
        printf("\nInvalid number of arguments!!!\n");
        exit(1);
    }
    //registering starting time
    ti = clock();
    
    //reading commandline parameters
    char *baseName = argv[1];
    char *paramFile = argv[2];
    char *fileExtFile = argv[3];
    int nPairs = atoi(argv[4]);
    if( argc == 6 ){verbos = atoi(argv[5]);}
    
    //read params from the paramFile
    paramHand.readParamsFromFile(paramFile);
    myProcParamsPtr = paramHand.getParamPtr();
    
    //read file extensions from the file
    paramHand.readFileExtsInfoFile(fileExtFile);
    myFileExtsPtr = paramHand.getExtPtr();
    
    //create a data handler object
    TSAdataHandler *TSData1 = new TSAdataHandler(baseName, &logs.procLogs, myFileExtsPtr, myProcParamsPtr);

    //generating subsequences from the time series data
    rVal = TSData1->genTemplate1SubSeqs();
    if (rVal==0)
    {
        return 0;
    }
    
    TSAsubSeq_t *subSeqPtr = TSData1->subSeqPtr;
    //initializing object for handling lower bounding processing
    TSAdtwSimilarity dtwUCR(&logs.procLogs);
    //initializing object for handling storage of top nPairs closest motif pairs
    //int nPairs = myProcParamsPtr->maxNMotifsPairs;
    TSApool *pool = new TSApool(nPairs);
    
    int lenMotifReal = TSData1->procParams.motifLengths[TSData1->procParams.indexMotifLenReal];
    int nInterFact = TSData1->procParams.pattParams.nInterpFac;

    // configuring object to initialize lower bounding pointers/arrays
    dtwUCR.configureTSASimilarity(lenMotifReal, lenMotifReal, TSData1->procParams.distParams.DTWBand);
    dtwUCR.setQueryPtr(TSData1->subSeqPtr, TSData1->nSubSeqs);
    dtwUCR.computeQueryEnvelops();
    dtwUCR.copyQueryEnv2Cand();

    U = dtwUCR.envUQueryPtr;
    L = dtwUCR.envLQueryPtr;
    accLB1 = dtwUCR.accLB_Keogh_EQ;
    accLB2 = dtwUCR.accLB_Keogh_EC;

    pool->initPriorityQDisc();

    //allocating memory for storing offseted data (this is used for implementing octave normalization)
    offseted_data = (double*)malloc(sizeof(double)*lenMotifReal);

    //Starting nested loop for motif discovery
    for(TSAIND ii=0;ii< TSData1->nSubSeqs;ii++)
    {
        for(TSAIND jj=ii+1;jj< TSData1->nSubSeqs;jj++)
        {
            //For interpolated subSeqs we don't compute distance for every combination of interpolation values. But only the ones that allow for all unique combinations.
            if (paramHand.procParams.combMTX[ii%nInterFact][jj%nInterFact]==0)
                continue;
            //NOTE: this condition is to avoid overlapping candidates.
            //OLD condition when flat note compression wasn't there: if (fabs(subSeqPtr[ii].sTime-subSeqPtr[jj].sTime)< TSData1->procParams.pattParams.blackDur)
            if((subSeqPtr[ii].eTime-subSeqPtr[jj].sTime)*(subSeqPtr[jj].eTime-subSeqPtr[ii].sTime) > 0)
            {
                continue;
            }
            //Computing offset between two subSeqs, to know if they are octave transposed or not
            offset = TSData1->estimateOffset(subSeqPtr[ii].mean, subSeqPtr[jj].mean);

            //Computing FL lower bound
            LB_kim_FL = computeLBkimFL(subSeqPtr[ii].pData[0], subSeqPtr[jj].pData[0] + offset, subSeqPtr[ii].pData[lenMotifReal-1], subSeqPtr[jj].pData[lenMotifReal-1]+ offset, SqEuclidean);
            logs.procLogs.nLB_KIM_FL++;

            if (LB_kim_FL< bsf) 
            {
                //Computing LB_Keogh EQ lower bound
                LB_Keogh_EQ = computeKeoghsLB(U[ii],L[ii],accLB1, subSeqPtr[jj].pData,lenMotifReal, bsf, SqEuclidean, offset);
                logs.procLogs.nLB_Keogh_EQ++;

                if(LB_Keogh_EQ < bsf)
                {
                    //Computing EB_Keogh EC lower bound
                    LB_Keogh_EC = computeKeoghsLB(U[jj],L[jj],accLB2, subSeqPtr[ii].pData,lenMotifReal, bsf, SqEuclidean, -1*offset);
                    logs.procLogs.nLB_Keogh_EC++;

                    if(LB_Keogh_EC < bsf)
                    {
                        //If the normalizatino type selected is octave normalization and the computed offset is not zero (i.e. we have detected a potential octave transposition)
                        if ((offset !=0)&&(TSData1->procParams.repParams.normType == OCTAVE_NORM)){
                            //creating a copy of the data with offseted samples
                            for (int zz=0;zz<lenMotifReal;zz++){offseted_data[zz] = subSeqPtr[jj].pData[zz]+offset;}    //saving the offsetted copy

                            //Computing the actual DTW distance now!!
                            realDist = dtw1dBandConst(subSeqPtr[ii].pData, offseted_data, lenMotifReal, lenMotifReal, dtwUCR.costMTX, SqEuclidean, dtwUCR.bandDTW, bsf, accLB1);
                        }
                        else{

                            //Computing the actual DTW distance now!!
                            realDist = dtw1dBandConst(subSeqPtr[ii].pData, subSeqPtr[jj].pData, lenMotifReal, lenMotifReal, dtwUCR.costMTX, SqEuclidean, dtwUCR.bandDTW, bsf, accLB1);
                        }
                        logs.procLogs.nDTW_EA++;
                        if (realDist <= bsf)
                        {
                            // if the distance found is less than the best so far lets update the queue
                            bsf = pool->managePriorityQDisc(subSeqPtr, ii, jj, realDist);
                            logs.procLogs.nPriorityUpdates++;
                        }
                    }
                }
            }
        }
    }
    //computing the final time
    tf = clock();
    logs.procLogs.tTotal += (tf-ti)/CLOCKS_PER_SEC;
    t1 = clock();
    TSData1->dumpDiscMotifInfo(TSData1->fHandle.getFileName(TSData1->fHandle.fileExtPtr->disOutFileExt), pool->priorityQDisc, pool->K, verbos);
    t2 = clock();
    logs.procLogs.tDump += (t2-t1)/CLOCKS_PER_SEC;
    logs.dumpProcLogs(TSData1->fHandle.getFileName(TSData1->fHandle.fileExtPtr->disLogFileExt), verbos);

    //clearing memory
    delete TSData1;
    delete pool;
    delete offseted_data;
    if (verbos){printf("Processing done!\n");}

    //All set, lets exit!
    return 1;
}


