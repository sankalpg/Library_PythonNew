
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
    float resolution = 0.01;
    
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
    
    TSData1->readQueryTimeStamps(fHandle.getSubSeqInfoFileName(), PATTERNS_PER_FILE_DUMP);
    
    TSData1->readKNNPatternDump(fHandle.getPatternKNNFileName(), MOTIFID1_MOTIFID2_DIST);
    
    if ( TSData1->nQueries != TSData1->nPatternPairs)
    {
        if (verbos ==1)
        {
            printf("There is some problem with the info and knn file for the file %s\n",baseName);
            return 0;
        }
    }
    TSData1->initializeBlackList(TSData1->nQueries);
    
    sortElem_t *sortArray = (sortElem_t*)malloc(sizeof(sortElem_t)*TSData1->nQueries);
    
    for (int ii=0;ii<TSData1->nQueries;ii++)
    {
        sortArray[ii].value = TSData1->patternPairs[ii].dist;
        sortArray[ii].index = ii;
    }
 
    //qsort(sortArray, TSData1->nQueries, sizeof(sortElem_t), compareSortElems);
    qsort(sortArray, TSData1->nQueries, sizeof(sortElem_t), compareSortElems);
    
    // finding out the biggest value of the time stamp of a pattern in the current file
    float timeMax = -INF;
    for(int jj=0; jj< TSData1->nPatternPairs; jj++)
    {
        if (TSData1->queryTStamps[jj].eTime > timeMax)
        {
            timeMax = TSData1->queryTStamps[jj].eTime;
        }
    }
    TSAIND nTimeSamples = (TSAIND)ceil(timeMax/resolution);
    int *overlapArray = (int*)malloc(sizeof(int)*nTimeSamples);
    for(int jj=0; jj< nTimeSamples; jj++)
    {
        overlapArray[jj]=-1;
    }
    int hasOverlap=0, OverlapIndex=-1;
    float str, end;
    int jj_sort;
    for (int jj=0;jj<TSData1->nPatternPairs;jj++)
    {
        jj_sort = sortArray[jj].index;
        hasOverlap=0;
        OverlapIndex=-1;
        if (TSData1->queryTStamps[jj_sort].id!=TSData1->patternPairs[jj_sort].ind1)
        {
            if (verbos==1)
            {
                printf("There is a mismatch in Info and Dist file for filename %s\n",baseName);
                return 0;
            }
            break;
        }
            
            
        str = (int)floor(TSData1->queryTStamps[jj_sort].sTime/resolution);
        end = (int)floor(TSData1->queryTStamps[jj_sort].eTime/resolution);
        
        // first check if there is an overlap
        for(TSAIND mm=str;mm<=end; mm++)
        {
            if (overlapArray[mm]!=-1)
            {
                hasOverlap =1;
                OverlapIndex = overlapArray[mm];
                break;
            }
        }
        if (hasOverlap==0)
        {
            for(TSAIND mm=str;mm<=end; mm++)
            {
                overlapArray[mm] = jj_sort;
            }
            TSData1->blacklist[jj_sort]=0;
            
        }
        else
        {
             TSData1->blacklist[jj_sort]=1;
        }
        
    }

    FILE *fp1 = fopen(fHandle.getBlackListSegFileName(), "w");
    for(int jj=0;jj<TSData1->nPatternPairs;jj++)
    {
        fprintf(fp1, "%d\n",TSData1->blacklist[jj]);
    }
    fclose(fp1);
    
    free(sortArray);
    free(overlapArray);
    
    delete TSData1;
    
    if (verbos){printf("Processing done!\n");}
    return 1;
}

