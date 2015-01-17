
#include "TSAdataIO.h"
#include "TSAsimilarity.h"
#include "TSApool.h"
#include "TSAlogs.h"


using namespace std;

typedef double (*distFunc)(double*, double*, int, int, double**, int, int, double, double*);
typedef int (*normFunc)(double*, int);
typedef int (*pathLenFunc)(double **, int, int);

TSADIST mu, median, stdev, mad;
TSADATA *tempArr = (TSADATA *)malloc(sizeof(TSADATA)*5000);

int zNormorm(TSADATA *data, int len)
{
    mu = computeMean(data, len);
    stdev = computeSTD(data, len, mu);
    for (int jj=0;jj<len;jj++)
    {
        data[jj] = (data[jj]-mu)/stdev;
    }
    return 1;
}

int meanNorm(TSADATA *data, int len)
{
    mu = computeMean(data, len);
    for (int jj=0;jj<len;jj++)
    {
        data[jj] = data[jj]-mu;
    }
    
}

int medianNorm(TSADATA *data, int len)
{
    memcpy(tempArr, data, sizeof(TSADATA)*len);
    median = computeMedian(tempArr, len);
    for (int jj=0;jj<len;jj++)
    {
        data[jj] = data[jj]-median;
    } 
}

int MADNorm(TSADATA *data, int len)
{
    memcpy(tempArr, data, sizeof(TSADATA)*len);
    median = computeMedian(tempArr, len);
    mad = computeMAD(tempArr, len, median);
    for (int jj=0;jj<len;jj++)
    {
        data[jj] = (data[jj]-median)/mad;
    } 
    
}

int noNorm(TSADATA *data, int len)
{
    return 1;
}
int tonicNorm(TSADATA *data, int len)
{
    return 1;
}

typedef struct couplet
{
    TSADIST dist;
    TSAIND ind;
}couplet_t;


TSAIND countQueries(char *fileName)
{
    FILE *fp;  
    fp =fopen(fileName,"r");
    int tempi[2];
    float tempf[2];
    
    if (fp==NULL)
    {
        printf("Error opening file %s\n", fileName);
        return 0;
    }
    TSAIND cnt = 0;
    while (fscanf(fp, "%f\t%f\t%d\t%d\n",&tempf[0], &tempf[1], &tempi[0], &tempi[1])!=EOF)    //read till the end of the file
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

double measureGlobalComplexity(TSADATA *data, int len)
{
  int elems =0;
  double accum=0;
  TSADATA tempDiff=0;
  for (int ii=0; ii<len-1;ii++)
  {
    tempDiff = fabs(data[ii]-data[ii+1]);
    if (tempDiff<=600)
    {
      accum+=tempDiff*tempDiff;
      elems+=1;
    }
    
  }
  if (elems>0)
  {
    return sqrt(accum/float(elems));
  }
  return EPS;
}


int main( int argc , char *argv[])
{
    int Err=0;
    int verbos;
    TSADIST realDist, LB_kim_FL, LB_Keogh_EQ, LB_Keogh_EC;
    int distType=-1;
    
    distFunc myDist[3]={NULL};
    normFunc myNorm[10]={NULL};
    pathLenFunc myPathLen[3]={NULL};
    
    typedef int (*pathLenFunc)(double **, int, int);
    
    myDist[0] = &euclideanSeq;
    myDist[1] = &dtw1dBandConst;
    myDist[2] = &dtw1dBandConst_localConst;
    
    myPathLen[0] = &path_Euclidean;
    myPathLen[1] = &path_11;
    myPathLen[2] = &path_12;
    
    myNorm[0] = &noNorm;
    myNorm[1] = &tonicNorm;
    myNorm[2] = &zNormorm;
    myNorm[3] = &meanNorm;
    myNorm[4] = &medianNorm;
    myNorm[5] = &MADNorm;


    
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
    
    //loading the data
    TSAdataHandler *TSData1 = new TSAdataHandler(fHandleTemp.searchFileNames[0], &logs.procLogs, myFileExtsPtr, myProcParamsPtr);
    TSAIND nSubs = TSData1->getNumLines(TSData1->fHandle.getSubSeqInfoFileName());
    if (myProcParamsPtr->repParams.normType == TONIC_NORM)
    {
        TSData1->readSubSeqData(TSData1->fHandle.getSubSeqTNFileName(), nSubs);
    }
    else
    {
        TSData1->readSubSeqData(TSData1->fHandle.getSubSeqFileName(), nSubs);
    }
    char *tempChar = TSData1->fHandle.getSubSeqInfoFileName();
    TSAIND nQueries = countQueries(tempChar);
    TSData1->readSubSeqLengths(tempChar);
    TSData1->downSampleSubSeqs();
    TSData1->genUniScaledSubSeqsVarLen();
    TSData1->quantizeSampleSubSeqs(TSData1->procParams.repParams.TSRepType);
    
    TSAIND nCands = (int)(TSData1->nSubSeqs/nInterFact);
    
    couplet_t *pArray = (couplet_t *)malloc(sizeof(couplet_t)*nSubs);
    
    TSAIND ii, jj, ss, kk;
    TSAIND ind1, ind2;
    int pattLenFinal=0, bandDTW;
    TSADIST *temp101, min_dist;
    int pathLen =0;
    
    int subSeqLen = TSData1->procParams.pattParams.subSeqLen;
    
    //temp buffers to store the linearly streched version of a subseq to make lengths equal
    TSADATA *buff1= (TSADATA *)malloc(sizeof(TSADATA)*subSeqLen);
    TSADATA *buff2= (TSADATA *)malloc(sizeof(TSADATA)*subSeqLen);    
    
    
    TSADIST **costMTX = (TSADIST **)malloc(sizeof(TSADIST*)*subSeqLen);
    for(ii=0;ii<subSeqLen;ii++)
    {
        costMTX[ii] = (TSADIST *)malloc(sizeof(TSADIST)*subSeqLen);
    }
    
    float complexity1, complexity2;
    for(ii=0; ii<nQueries; ii++)
    {
        for(ss=0; ss<nCands;ss++)
        {
            pArray[ss].dist = INF;
            pArray[ss].ind = -1;
        }
        for(jj=0; jj<nCands; jj++)
        {
            if(ii==jj)
            {
                pArray[jj].ind = jj;
                pArray[jj].dist = INF;
                continue;
            }
            ind1 = ii*nInterFact;
            ind2 = jj*nInterFact;

           
            if (TSData1->procParams.distParams.distNormType >=MAXLEN_NO_NORM)
            {   
                pattLenFinal =0;
                for(int tt=0; tt<nCands;tt++)
                    if (pattLenFinal < TSData1->subSeqPtr[tt].len)
                    {
                        pattLenFinal = TSData1->subSeqPtr[tt].len;
                    }
            }
            else
            {
                pattLenFinal = max(TSData1->subSeqPtr[ind1].len, TSData1->subSeqPtr[ind2].len);
            }
            

            bandDTW = (int)floor(pattLenFinal*myProcParamsPtr->distParams.DTWBand);
            for(ss=0;ss<pattLenFinal;ss++)
            {
                for(kk=0;kk<pattLenFinal;kk++)
                {
                    costMTX[ss][kk]=INF;
                }
            }
            min_dist = INF;
            for(ss=0; ss< nInterFact; ss++)
            {
                for(kk=0; kk< nInterFact; kk++)
                {
                    //if (paramHand.procParams.combMTX[ss][kk]==0)
                    //continue;
                    copyLinStrechedBuffer(buff1, TSData1->subSeqPtr[ind1+ss].pData, TSData1->subSeqPtr[ind1].len, pattLenFinal);
                    copyLinStrechedBuffer(buff2, TSData1->subSeqPtr[ind2+kk].pData, TSData1->subSeqPtr[ind2].len, pattLenFinal);
                    complexity1 = measureGlobalComplexity(buff1, pattLenFinal-1);
                    complexity2 = measureGlobalComplexity(buff2, pattLenFinal-1);
                    myNorm[myProcParamsPtr->repParams.normType](buff1, pattLenFinal);
                    myNorm[myProcParamsPtr->repParams.normType](buff2, pattLenFinal);
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

