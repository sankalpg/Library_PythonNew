
#include "TSAdataIO.h"


using namespace std;


int main( int argc , char *argv[])
{
    procParams_t *myProcParamsPtr;
    fileExts_t *myFileExtsPtr;
    int Err=0;
    int verbos;
    procLogs_t myProcLogs;
    
    
    TSAparamHandle paramHand;
    
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
    
    //initialize the log counts
    initializeLogCounts(&myProcLogs);
    
    //create a data handler object
    TSAdataHandler TSData1(baseName, &myProcLogs, myFileExtsPtr, myProcParamsPtr);
    
    TSData1.genTemplate1SubSeqs();
    
    
    
    
    //generate pattern sub sequences
    
    
    //prepare lower bounding data
    
    //enter into nested loop
    
    //free memories
    
    //dump the generate data
    
    
    if (verbos){printf("Processing done!\n");}
    return 1;
}

