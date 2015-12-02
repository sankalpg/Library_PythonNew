#include "../TSAdataIO.h"
#include "../TSAsimilarity.h"
#include "../TSApool.h"
#include "../TSAlogs.h"

using namespace std;

int main( int argc , char *argv[])
{
    int Err=0;
    int verbos;
    TSADIST realDist, LB_kim_FL, LB_Keogh_EQ, LB_Keogh_EC, bsf_local;
    TSAIND searchFileID;

    
    procParams_t *myProcParamsPtr;
    fileExts_t *myFileExtsPtr;
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
    if( argc == 6 ){verbos = atoi(argv[5]);}    
    
    //read params from the paramFile
    paramHand.readParamsFromFile(paramFile);
    myProcParamsPtr = paramHand.getParamPtr();
    
    //read file extensions from the file
    paramHand.readFileExtsInfoFile(fileExtFile);
    myFileExtsPtr = paramHand.getExtPtr();
    
    fHandle.initialize(baseName, myFileExtsPtr);
    fHandle.loadSearchFileList();
    
    //create a data handler object
    TSAdataHandler *TSData1 = new TSAdataHandler(baseName, &logs.procLogs, myFileExtsPtr, myProcParamsPtr);
    TSAIND NPatternsFile1 = TSData1->getNumLines(TSData1->fHandle.getSubSeqInfoFileName());
    printf("Number of lines in file %lld\n",NPatternsFile1);
    //reading data from stored subsequences
    TSData1->calculateDiffMotifLengths();   //note that this give lengths of patterns before any downsampling. Need it for estimating subseqlen
    int subSeqLen = TSData1->procParams.motifLengths[TSData1->procParams.indexMotifLenLongest];
    TSData1->procParams.pattParams.subSeqLen = subSeqLen;
    TSData1->procParams.pattParams.subSeqLen = subSeqLen;
    TSData1->readSubSeqData(TSData1->fHandle.getSubSeqFileName(), NPatternsFile1);
    TSData1->readSubSeqInfo(TSData1->fHandle.getSubSeqInfoFileName());
    TSData1->setSubSeqLengthsTStamps();
    TSData1->downSampleSubSeqs();
    TSData1->calculateDiffMotifLengths();   //this will give us the motif lengths that we should be using (i.e. lengths after downsampling)
    int lenMotifReal = TSData1->procParams.motifLengths[TSData1->procParams.indexMotifLenReal];
    int nInterFact = TSData1->procParams.pattParams.nInterpFac;
    TSData1->setSubSeqLengthsFIX(lenMotifReal);
    TSData1->genUniScaledSubSeqsVarLen();
    TSData1->normalizeSubSeqs(TSData1->procParams.repParams.normType);
    
    
    
    TSAdtwSimilarity dtwUCR( &logs.procLogs);
    dtwUCR.configureTSASimilarity(lenMotifReal, lenMotifReal, TSData1->procParams.distParams.DTWBand);
    dtwUCR.setQueryPtr(TSData1->subSeqPtr, TSData1->nSubSeqs);
    dtwUCR.computeQueryEnvelops();
    dtwUCR.initArrayBSF(NPatternsFile1);
    
    TSApool pool(kNN);
    pool.initPriorityQSear(NPatternsFile1);
    
    fHandle.loadSearchFileList();
    TSAIND queryInd=0;
    TSAIND ind1, ind2;
    //printf("Hello1");
    //iterating over all the files
    for(TSAIND ss=0; ss < fHandle.nSearchFiles; ss++)
    {
        searchFileID = ss;
        //printf("Processing searchfile %s\n",fHandle.searchFileNames[ss]);
        TSAdataHandler *TSData2 = new TSAdataHandler(fHandle.searchFileNames[ss], &logs.procLogs, myFileExtsPtr, myProcParamsPtr);
        TSAIND NPatternsFile2 = TSData2->getNumLines(TSData2->fHandle.getSubSeqInfoFileName());
        TSData2->procParams.pattParams.subSeqLen = subSeqLen;
        //reading data from stored subsequences
        TSData2->procParams.pattParams.subSeqLen = subSeqLen;
        TSData2->readSubSeqData(TSData2->fHandle.getSubSeqFileName(), NPatternsFile2);
        TSData2->readSubSeqInfo(TSData2->fHandle.getSubSeqInfoFileName());
        TSData2->setSubSeqLengthsTStamps();
        TSData2->downSampleSubSeqs();
        TSData2->setSubSeqLengthsFIX(lenMotifReal);
        TSData2->genUniScaledSubSeqsVarLen();
        TSData2->normalizeSubSeqs(TSData1->procParams.repParams.normType);
       
        dtwUCR.setCandPtr(TSData2->subSeqPtr, TSData2->nSubSeqs);
        dtwUCR.computeCandEnvelops();
        
        for(TSAIND ii=0;ii< NPatternsFile1;ii++)
        {
            //printf("%lld",ii);
            for(TSAIND jj=0;jj< NPatternsFile2;jj++)
            {
                //printf("%lld",jj);
                queryInd = ii;
                ind1 = ii*nInterFact;
                ind2 = jj*nInterFact;
                
                //we have to avoid overlapping patterns being detected as neighbors
                if ((strcmp(baseName, fHandle.searchFileNames[ss])==0)&& (fabs(TSData1->subSeqPtr[ind1].sTime-TSData2->subSeqPtr[ind2].sTime)< TSData1->procParams.pattParams.blackDur))
                    //beware that basename and searchFile name should both have either full path or relative path but should be in the same format always.
                {
                    continue;
                }
                bsf_local = dtwUCR.bsfArray[queryInd];
                for (int pp = 0; pp < nInterFact; pp++)
                {
                    for (int mm = 0; mm < nInterFact; mm++)
                    {                
                        if (paramHand.procParams.combMTX[pp][mm]==0)
                            continue;
                        
                        LB_kim_FL = computeLBkimFL(TSData1->subSeqPtr[ind1+pp].pData[0], TSData2->subSeqPtr[ind2+mm].pData[0], TSData1->subSeqPtr[ind1+pp].pData[lenMotifReal-1], TSData2->subSeqPtr[ind2+mm].pData[lenMotifReal-1], SqEuclidean);
                        logs.procLogs.nLB_KIM_FL++;
                        if (LB_kim_FL< bsf_local) 
                        {
                            LB_Keogh_EQ = computeKeoghsLB(dtwUCR.envUQueryPtr[ind1+pp],dtwUCR.envLQueryPtr[ind1+pp],dtwUCR.accLB_Keogh_EQ, TSData2->subSeqPtr[ind2+mm].pData,lenMotifReal, dtwUCR.bsfArray[queryInd], SqEuclidean, 0.0);
                            logs.procLogs.nLB_Keogh_EQ++;
                            if(LB_Keogh_EQ < bsf_local)
                            {
                                LB_Keogh_EC = computeKeoghsLB(dtwUCR.envUCandPtr[ind2+mm],dtwUCR.envLCandPtr[ind2+mm],dtwUCR.accLB_Keogh_EC, TSData1->subSeqPtr[ind1+pp].pData,lenMotifReal, dtwUCR.bsfArray[queryInd], SqEuclidean, 0.0);
                                logs.procLogs.nLB_Keogh_EC++;
                                if(LB_Keogh_EC < bsf_local)
                                {
                                    realDist = dtw1dBandConst(TSData1->subSeqPtr[ind1+pp].pData, TSData2->subSeqPtr[ind2+mm].pData, lenMotifReal, lenMotifReal, dtwUCR.costMTX, SqEuclidean, dtwUCR.bandDTW, dtwUCR.bsfArray[queryInd], dtwUCR.accLB_Keogh_EQ);
                                    logs.procLogs.nDTW_EA++;
                                    if(realDist < bsf_local)
                                    {
                                        bsf_local = realDist;
                                    }
                                }
                            }
                        }
                    }
                }
                if (bsf_local <= dtwUCR.bsfArray[queryInd])
                {
                    if (bsf_local==0)
                    {
                        continue;
                    }
                    dtwUCR.bsfArray[queryInd] = pool.managePriorityQSear(queryInd, TSData2->subSeqPtr,  ind1, ind2, bsf_local, searchFileID, TSData1->procParams.pattParams.blackDur);
                    logs.procLogs.nPriorityUpdates++;
                }
            }
            //printf("%f\n",dtwUCR.bsfArray[queryInd]);
        }
        
        delete TSData2;
        dtwUCR.deleteCandEnvMem(); 
        
        
    }
    TSData1->dumpPatternKNNInfo(TSData1->fHandle.getOutFileName(), pool.priorityQSear, NPatternsFile1, kNN, verbos);
    delete TSData1;
    
    
}
    
    
   
    
    
    