
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
    
    TSAparamHandle paramHand;
    TSAlogs logs;
    fileNameHandler fHandle;
    
    
    //checking if the number of input arguments are correct 
    if(argc < 4 || argc > 5)
    {
        printf("\nInvalid number of arguments!!!\n");
        exit(1);
    }
    
    ti = clock();
    
    //reading commandline parameters
    char *baseName = argv[1];
    char *paramFile = argv[2];
    char *fileExtFile = argv[3];
    if( argc == 5 ){verbos = atoi(argv[4]);}
    
    //read params from the paramFile
    paramHand.readParamsFromFile(paramFile);
    myProcParamsPtr = paramHand.getParamPtr();
    
    //read file extensions from the file
    paramHand.readFileExtsInfoFile(fileExtFile);
    myFileExtsPtr = paramHand.getExtPtr();
    
    
    fHandle.initialize(baseName, myFileExtsPtr);
    
    //create a data handler object
    TSAdataHandler *TSData1 = new TSAdataHandler(baseName, &logs.procLogs, myFileExtsPtr, myProcParamsPtr);
    
    //TSData1->loadMotifDataTemplate1();
    TSData1->readTSData(fHandle.getTSFileName());
    
    TSData1->readHopSizeTS(fHandle.getTSFileName());
    
    TSData1->downSampleTS();
    
    TSData1->filterSamplesTS();
    
    TSData1->convertHz2Cents(fHandle.getTonicFileName());
    
    //calculate different motif lengths before doing sliding window candidate generation
    TSData1->calculateDiffMotifLengths();
    
    TSData1->readQueryTimeStamps(fHandle.getSubSeqInfoFileName(), PATTERNS_PER_FILE_DUMP);
    
    TSData1->genSubSeqsWithTStarts(TSData1->queryTStamps, TSData1->nQueries);
    
    int lenRawMotifData = TSData1->procParams.motifLengths[TSData1->procParams.indexMotifLenLongest];
    
    FILE *fp = fopen(fHandle.getSubSeqFileName(), "wb");
    
    for (int ii=0;ii<TSData1->nQueries; ii++)
    {
        fwrite(TSData1->subSeqPtr[ii].pData, sizeof(TSADATA), lenRawMotifData, fp);
    }
    delete TSData1;
    
    if (verbos){printf("Processing done!\n");}
    return 1;
}

