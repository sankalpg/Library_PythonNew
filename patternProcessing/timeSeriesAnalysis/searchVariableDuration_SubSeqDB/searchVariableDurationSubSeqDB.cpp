#include "../TSAdataIO.h"
#include "../TSAsimilarity.h"
#include "../TSApool.h"
#include "../TSAlogs.h"


using namespace std;

typedef double (*distFunc)(double*, double*, int, int, double**, int, int, double, double*);
typedef int (*normFunc)(double*, int);
typedef int (*pathLenFunc)(double **, int, int);
typedef double (*complexityFunc)(double *, int);

TSADIST mu, median, stdev, mad;


/*These are the functions which are needed specifically for this file and are not generic enough to be moved to common files
 */

typedef struct couplet
{
    TSADIST dist;
    TSAIND ind;
}couplet_t;


TSAIND countQueries(char *fileName)
{
    FILE *fp;  
    fp =fopen(fileName,"r");
    int tempi[3];
    float tempf[2];
    
    if (fp==NULL)
    {
        printf("Error opening file %s\n", fileName);
        return 0;
    }
    TSAIND cnt = 0;ICASSP2015_Experiment_O3_CUR
    while (fscanf(fp, "%f\t%f\t%d\t%d\t%d\n",&tempf[0], &tempf[1], &tempi[0], &tempi[1], &tempi[2])!=EOF)    //read till the end of the file
    {
        if (tempi[1]!=-1)
        {
            cnt++;
        }
    }
    fclose(fp);
    
    return cnt;
    
}
int compareDistPatts(const void *a, const void *b)
{
    if (((couplet_t*)a)->dist > ((couplet_t*)b)->dist)
    {
        return 1;
    }
    else if (((couplet_t*)a)->dist < ((couplet_t*)b)->dist)
    {
        return -1;
    }
    return 0;
}

int copyLinStrechedBuffer(TSADATA *dst, TSADATA *src, int lenSrc, int lenDst)
{
    TSAIND ii,jj;
    if (lenSrc==lenDst)
    {
        memcpy(dst, src, sizeof(TSAIND)*lenDst);
    }
    else
    {
        for(ii=0;ii<lenDst;ii++)
        {
            dst[ii] = src[(int)floor(ii*lenSrc/lenDst)];
        }
    }
}

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


int main( int argc , char *argv[])
{
    int Err=0;
    int verbos;
    TSADIST realDist, LB_kim_FL, LB_Keogh_EQ, LB_Keogh_EC;
    int distType=-1;
    float complexity1, complexity2;
    
    //declaring function pointers (I prefer to use function pointers to respond to different parma values rather than if else condition because its much much faster specially if its in a loop. Its doesn't even involve any comparison.)    
    distFunc myDist[3]={NULL};
    normFunc myNorm[10]={NULL};
    pathLenFunc myPathLen[3]={NULL};
    complexityFunc myComplexity[4]={};
    
    //Filling in all the method variants for computing distances
    myDist[0] = &euclideanSeq;
    myDist[1] = &dtw1dBandConst;
    myDist[2] = &dtw1dBandConst_localConst;
    
    //variants for different kind of path lengths. Note that even though path doesn't make sense in the case of euclidean distance, for an efficient implementation I make a dummy wrapper.    
    myPathLen[0] = &path_Euclidean;
    myPathLen[1] = &path_11;
    myPathLen[2] = &path_12;
    
    //variants for normalization of melody sequence (Here also several functions are dummp and they retun the input signal as it is. But they are maded so that they fit the template of function pointers to avoid comparison in loop)
    myNorm[0] = &noNorm;
    myNorm[1] = &tonicNorm;
    myNorm[2] = &zNorm;
    myNorm[3] = &meanNorm;
    myNorm[4] = &medianNorm;
    myNorm[5] = &MADNorm;
    
    //variants for complexity measures (for carnatic music some of these measures can be benefitial)
    myComplexity[0] = &measureGlobalComplexity1;
    myComplexity[1] = &measureGlobalComplexity2;
    myComplexity[2] = &computeInflectionPoints1;
    myComplexity[3] = &computeInflectionPoints2;
    
    //declaring some other useful variables
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
    
    //since we dump output in appending mode. We open that file once in 'w' mode at the beginning to clear out all the contents
    FILE *fp10;
    fp10 = fopen(fHandleTemp.getOutFileName(), "w");
    fclose(fp10);
    
    
    int nInterFact = myProcParamsPtr->pattParams.nInterpFac;
    float blackDur=0;
    
    //reading distance type based on parameters read from param file
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
    
    //loading the data from the subsequence database
    TSAdataHandler *TSData1 = new TSAdataHandler(fHandleTemp.searchFileNames[0], &logs.procLogs, myFileExtsPtr, myProcParamsPtr);
    TSAIND nSubs = TSData1->getNumLines(TSData1->fHandle.getSubSeqInfoFileName());
    //if user has chosen tonic norm or pa sa norm, in that case we load tonic normalized subsequences
    if ((myProcParamsPtr->repParams.normType == TONIC_NORM) || (myProcParamsPtr->repParams.normType == TONIC_NORM_PASAPA))
    {
        TSData1->readSubSeqData(TSData1->fHandle.getSubSeqTNFileName(), nSubs);
    }
    else
    {
        TSData1->readSubSeqData(TSData1->fHandle.getSubSeqFileName(), nSubs);
    }
    
    
    char *tempChar = TSData1->fHandle.getSubSeqInfoFileName();
    //count number of queries in the subseq database
    TSAIND nQueries = countQueries(tempChar);    
    //reading pattern lengths from info file, downsampling, generating multiple time-stretched copies and quantizing pitch values.
    TSData1->readSubSeqLengths(tempChar);
    TSData1->downSampleSubSeqs();
    TSData1->genUniScaledSubSeqsVarLen();
    TSData1->quantizeSampleSubSeqs(TSData1->procParams.repParams.TSRepType);
    
    // computing total number of candidates
    TSAIND nCands = (int)(TSData1->nSubSeqs/nInterFact);
    
    couplet_t *pArray = (couplet_t *)malloc(sizeof(couplet_t)*nSubs);
    
    TSAIND ii, jj, ss, kk;
    TSAIND ind1, ind2;
    int pattLenFinal=0, bandDTW, nSamplesCandidate;
    TSADIST *temp101, min_dist;
    int pathLen =0;
    
    //length of the subsequences stores. Note that this is not the length fo the pattern. pattern length should be < this length. We need some extra sample to do time -compacting. 
    int subSeqLen = TSData1->procParams.pattParams.subSeqLen;
    
    //temp buffers to store the linearly streched version of a subseq to make lengths equal
    TSADATA *buff1= (TSADATA *)malloc(sizeof(TSADATA)*subSeqLen);
    TSADATA *buff2= (TSADATA *)malloc(sizeof(TSADATA)*subSeqLen);    
    
    //Initializing costMTX for the DTW computation
    TSADIST **costMTX = (TSADIST **)malloc(sizeof(TSADIST*)*subSeqLen);
    for(ii=0;ii<subSeqLen;ii++)
    {
        costMTX[ii] = (TSADIST *)malloc(sizeof(TSADIST)*subSeqLen);
    }
    
    //in case we chose a distance normalization where all patterns are brought to the same length, we have to compute max length of the patterns over entire dataset.
    if (TSData1->procParams.distParams.distNormType >= MAXLEN_NO_NORM)
    {   
        pattLenFinal =0;
        for(int tt=0; tt<nCands;tt++)
            if (pattLenFinal < TSData1->subSeqPtr[tt].len)
            {
                pattLenFinal = TSData1->subSeqPtr[tt].len;
            }
    }
    
    //Looping for every query pattern
    for(ii=0; ii<nQueries; ii++)
    {
        //since there are nInterFact number of time stretched version, computing the starting index of the new query
		ind1 = ii*nInterFact;
        
        //for every query nullifying the array where we store the distance with the candidates
        for(ss=0; ss<nCands;ss++)
        {
            pArray[ss].dist = INF;
            pArray[ss].ind = -1;
        }
        
        // Looping for every candidate pattern
        for(jj=0; jj<nCands; jj++)
        {
            //if the candidate happens to be the quey itself mark the distance as INF
            if(ii==jj)
            {
                pArray[jj].ind = jj;
                pArray[jj].dist = INF;
                continue;
            }
            
            //starting index of the candidate pattern because every candidate also has nInterFact number of time-stretched versions
            ind2 = jj*nInterFact;
            
            //Var1: both query len and candidate lengths are stretched to same length. Hence number of samples taken from candidate patterns is its length
            //Var2: Candidates are taken to be the same length as th query. So number of samples to be considered for candidates is the length of the query
            if (TSData1->procParams.methodVariant == Var1)
                {nSamplesCandidate =    TSData1->subSeqPtr[ind2].len;}
                else if (TSData1->procParams.methodVariant == Var2)
                {nSamplesCandidate =    TSData1->subSeqPtr[ind1].len;}
                else
                {
                    printf("Please provide a valid method variant type\n");
                    return -1;
                }

            //in Var1 pattLenFinal is the maximum of the lengths of the query and candidate. Whereas in Var2 it is the length of the query always.
            if (TSData1->procParams.distParams.distNormType < MAXLEN_NO_NORM)
            {
                if (TSData1->procParams.methodVariant == Var1)
                {
                    pattLenFinal = max(TSData1->subSeqPtr[ind1].len, TSData1->subSeqPtr[ind2].len);
                }
                else if (TSData1->procParams.methodVariant == Var2)
                {
                    pattLenFinal = TSData1->subSeqPtr[ind1].len;    
                }
                else
                {
                    printf("Please provide a valid method variant type\n");
                    return -1;
                }
				
            }
            bandDTW = (int)floor(pattLenFinal*myProcParamsPtr->distParams.DTWBand);
            
            //Since we are dealing with variable length queries and candidates, we need to reinitialize dtw CostMTX in every iteration. Note that typically this operation is done in DTW function. But in cases where 
            //query length remains the same, we dont have to do this in every iteration (aka also not in DTW function) because number of pixels filled by DTW function in each iteration are exactly the same. 
            for(ss=0;ss<pattLenFinal;ss++)
            {
                for(kk=0;kk<pattLenFinal;kk++)
                {
                    costMTX[ss][kk]=INF;
                }
            }
            
            
            min_dist = INF;
            //for one query and one candidate we iterate over all possible combinations of time-stretched verisons and take the minimum distance.
            for(ss=0; ss< nInterFact; ss++)
            {
                for(kk=0; kk< nInterFact; kk++)
                {
                    //initializing complexity of both query and candidate to 1. This way if we dont compute the complexity (if specified in the param file) we get the same distance as computed from the DTW/Euclidean distance.
                    complexity1 =1;
                    complexity2 =1;
                    
                    //if we didn't have to compute the minimum overall possible combinations we could have chosen set of combinations for which we do this computation, like done for the unsupervised anlaysis. Thats why next line is commented because here we dont do it.!
                    //if (paramHand.procParams.combMTX[ss][kk]==0)
                    //continue;
                    
                    //Performing linear time stretch by step sampling and making two patterns of different lengths the same length
                    copyLinStrechedBuffer(buff1, TSData1->subSeqPtr[ind1+ss].pData, TSData1->subSeqPtr[ind1].len, pattLenFinal);
                    copyLinStrechedBuffer(buff2, TSData1->subSeqPtr[ind2+kk].pData, nSamplesCandidate, pattLenFinal);
                    
                    //If user chooses to weight the distance measure with the 
                    if (TSData1->procParams.complexityMeasure>=0)
                    {
                        complexity1 = myComplexity[TSData1->procParams.complexityMeasure](buff1, pattLenFinal-1);
                        complexity2 = myComplexity[TSData1->procParams.complexityMeasure](buff2, pattLenFinal-1);
                    }
                    
                    if (myProcParamsPtr->repParams.normType == TONIC_NORM_PASAPA)
                    {
                      normalizePASAPA(buff1, pattLenFinal, buff2, pattLenFinal);
                    }
                    else
                    {
                      myNorm[myProcParamsPtr->repParams.normType](buff1, pattLenFinal);
                      myNorm[myProcParamsPtr->repParams.normType](buff2, pattLenFinal);
                      
                    }
                    
                    realDist = myDist[distType](buff1, buff2, pattLenFinal, pattLenFinal, costMTX, SqEuclidean, bandDTW, -1, temp101);
                    realDist = realDist*(max(complexity1, complexity2)/min(complexity1, complexity2));                    
                    if ((TSData1->procParams.distParams.distNormType == NORM_1) || (TSData1->procParams.distParams.distNormType==MAXLEN_NO_NORM))
                    {
                        pathLen = 1;    
                    }
                    else if ((TSData1->procParams.distParams.distNormType == PATH_LEN) || (TSData1->procParams.distParams.distNormType==MAXLEN_PATH_LEN))
                    {
                        pathLen = myPathLen[distType](costMTX, pattLenFinal, pattLenFinal);
                    }
                    else
                    {
                        pathLen =  pattLenFinal;
                    }
                    realDist= realDist/pathLen;
                    
                    if(realDist < min_dist)
                    {
                        min_dist=realDist;
                    }
                }
            }
            
            pArray[jj].ind = jj;
            pArray[jj].dist = min_dist;
        }
        qsort(pArray, nCands, sizeof(couplet_t), compareDistPatts);
        FILE *fp;
        fp = fopen(fHandleTemp.getOutFileName(), "ab");
        
        for(TSAIND pp=0;pp<kNN;pp++)
        {
            fprintf(fp, "%lld\t%lld\t%f\n", ii, pArray[pp].ind, pArray[pp].dist);
        }
        fclose(fp);
        
    }
    
    
    for(int ii=0;ii<TSData1->nSubSeqs;ii++)
    {
        free(TSData1->subSeqPtr[ii].pData);
    }
    free(TSData1->subSeqPtr);
    for(ii=0;ii<subSeqLen;ii++)
    {
        free(costMTX[ii]);
    }
    free(costMTX);
    free(buff1);
    free(buff2);
    free(pArray);
    if (verbos){printf("Processing done!\n");}
    return 1;
}

