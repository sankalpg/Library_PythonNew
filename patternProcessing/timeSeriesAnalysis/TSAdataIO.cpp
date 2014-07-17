
#include "TSAdataIO.h"

using namespace std;

TSAparamHandle::TSAparamHandle()
{
    
    memset(fileExts.tsFileExt, '\0', sizeof(char)*MAX_FEXT_CHARS);
    memset(fileExts.tonicExt, '\0', sizeof(char)*MAX_FEXT_CHARS);
    memset(fileExts.blackTimeExt, '\0', sizeof(char)*MAX_FEXT_CHARS);
    memset(fileExts.logFileExt, '\0', sizeof(char)*MAX_FEXT_CHARS);
    memset(fileExts.paramsDumpExt, '\0', sizeof(char)*MAX_FEXT_CHARS);
    memset(fileExts.outFileExt, '\0', sizeof(char)*MAX_FEXT_CHARS);
}
 
int TSAparamHandle::readParamsFromFile(char *paramFile)
{
    FILE *fp;
    char tempFilename[400]={'\0'};
    char field[100]={'\0'};
    char value[100]={'\0'};
    
    fp = fopen(paramFile, "r");
    if (fp == NULL){printf("Error opening file %s\n", paramFile);exit(0);}
    while(fgets(tempFilename, 400, fp))
    {
        sscanf(tempFilename, "%s %s\n", field, value);
        if (strcmp(field, "durMotif:")==0){procParams.durMotif=atof(value);}
        if (strcmp(field, "blackDurFactor:")==0){procParams.blackDur=atof(value)*procParams.durMotif;}
        if (strcmp(field, "dsFactor:")==0){procParams.dsFactor=atoi(value);}
        if (strcmp(field, "binsPOct:")==0){procParams.binsPOct=atoi(value);}
        if (strcmp(field, "minPossiblePitch:")==0){procParams.minPossiblePitch=atof(value);}
        if (strcmp(field, "varDur:")==0){procParams.varDur=atof(value);}
        if (strcmp(field, "threshold:")==0){procParams.threshold=atof(value);}
        if (strcmp(field, "flatThreshold:")==0){procParams.flatThreshold=atof(value);}
        if (strcmp(field, "maxPauseDur:")==0){procParams.maxPauseDur=atof(value);}
        if (strcmp(field, "DTWBand:")==0){procParams.DTWBand=atof(value);}
        if (strcmp(field, "nInterpFac:")==0){procParams.nInterpFac=atoi(value);}
        if (strcmp(field, "SimMeasuresUsed:")==0){procParams.SimMeasuresUsed=atoi(value);}
        if (strcmp(field, "removeTaniSegs:")==0){procParams.removeTaniSegs=atoi(value);}
        if (strcmp(field, "dumpLogs:")==0){procParams.dumpLogs=atoi(value);}
    }
    fclose(fp);
    
    return 1;
}

int TSAparamHandle::readFileExtsInfoFile(char *fileExtsFile)
{
    FILE *fp;
    char tempFilename[400]={'\0'};
    char field[100]={'\0'};
    char value[100]={'\0'};
    
    fp = fopen(fileExtsFile, "r");
    if (fp == NULL){printf("Error opening file %s\n", fileExtsFile);exit(0);}
    while(fgets(tempFilename, 400, fp))
    {
        sscanf(tempFilename, "%s %s\n", field, value);
        if (strcmp(field, "tsFileExt:")==0){strcat(fileExts.tsFileExt, value);}
        if (strcmp(field, "tonicExt:")==0){strcat(fileExts.tonicExt, value);}
        if (strcmp(field, "blackTimeExt:")==0){strcat(fileExts.blackTimeExt, value);}
        if (strcmp(field, "logFileExt:")==0){strcat(fileExts.logFileExt, value);}
        if (strcmp(field, "paramsDumpExt:")==0){strcat(fileExts.paramsDumpExt, value);}
        if (strcmp(field, "outFileExt:")==0){strcat(fileExts.outFileExt, value);}
    }
    fclose(fp);
    
    return 1;
}
procParams_t* TSAparamHandle::getParamPtr()
{
    return &procParams;
};
fileExts_t* TSAparamHandle::getExtPtr()
{
    return &fileExts;
};

int initializeLogCounts(procLogs_t *myProcLogs)
{
    myProcLogs->timeDataLoad=0;
    myProcLogs->timeGenSubs=0;
    myProcLogs->timeRemBlacklist=0; 
    myProcLogs->timeGenEnvelops=0;
    myProcLogs->timeDiscovery=0;
    myProcLogs->timeWriteData=0;
    myProcLogs->timeTotal=0;
    myProcLogs->totalPitchSamples=0;
    myProcLogs->totalPitchNonSilSamples=0;
    myProcLogs->totalSubsGenerated=0;
    myProcLogs->totalSubsBlacklisted=0;
    myProcLogs->totalSubsInterpolated=0;
    myProcLogs->totalFLDone=0;
    myProcLogs->totalLBKeoghEQ=0;
    myProcLogs->totalLBKeoghEC=0;
    myProcLogs->totalDTWComputations=0;
    myProcLogs->totalPriorityUpdates=0;
}



TSAdataHandler::TSAdataHandler(char *bName, procLogs_t *procLogs, fileExts_t *fileExts, procParams_t *pParams)
{
    
    procLogPtr = procLogs;
    fileExtPtr = fileExts;
    procParams = *pParams;
    baseName = bName;
    
    fHandle.initialize(baseName, fileExtPtr);
	
}

int TSAdataHandler::genTemplate1SubSeqs()
{
    //read the time series data    
    readTSData(fHandle.getTSFileName());
    readHopSizeTS(fHandle.getTSFileName());
    
    printf("%lld\n", lenTS);
    
    //downsample
    downSampleTS();
    
    printf("%lld\n", lenTS);
    
    //remove silence pitch regions
    filterSamplesTS();
    
    printf("%lld\n", lenTS);
    
    convertHz2Cents(fHandle.getTonicFileName());
	
	//calculate different motif lengths before doing sliding window candidate generation
	calculateDiffMotifLengths();
    
    genSlidingWindowSubSeqs();
    
    initializeBlackList();
    
    updateBLDurThsld();
    
    updateBLStdThsld();
    
    updateBLInvalidSegment(fHandle.getBlackListSegFileName());
    
    
    
    
}

int updateBLInvalidSegment(char *segmentFile)
{
    TSAIND ii,jj,nLines;
    FILE *fp;
    float temp[4] = {0};
    TSAPattern_t *invalidSegs;
    
    nLines = getNumLines(segmentFile);
    invalidSegs = (segInfo_t *)malloc(sizeof(segInfo_t)*nLines);
    fp =fopen(segmentFile,"r");
    if (fp==NULL)
    {
        printf("Error opening file %s\n", segmentFile);
        return 0;
    }
    //reading segments of the tani sections
    ii=0;
    while (fscanf(fp, "%f %f\n",&temp[0],&temp[1])!=EOF)
    {
        invalidSegs[ii].sTime = temp[0];
        invalidSegs[ii].eTime = temp[1];
        ii++;
    }
    fclose(fp);
    //blacklisting the indexes whose time stamps lie between these segments
    for(ii=0;ii<nSubSeqs;ii++)
    {
        for(jj=0;jj<nLines;jj++)
        {
            if ((subSeqPtr[ii].eTime >= invalidSegs[jj].sTime) && (subSeqPtr[ii].sTime <= invalidSegs[jj].eTime))
                blacklist[ii]=1;
        }
        
    }
    
    free(invalidSegs);
    
    return 1;
}

/*
 * This function computes a running variance of the input data
 */
int TSAdataHandler::computeStdTSLocal(float **std, int varSam)
{
    float *mean, *stdVec, temp[2]={0};
    TSAIND ii;
    int N;
    
    
    mean = (float *)malloc(sizeof(float)*lenTS);
    memset(mean,0,sizeof(float)*lenTS);
    stdVec = (float *)malloc(sizeof(float)*lenTS);
    memset(stdVec,0,sizeof(float)*lenTS);
    
    N = 2*varSam+1 ;    
    for(ii=0;ii<N;ii++)
    {
        mean[varSam] +=  samPtr[ii].value ;
    }    
    for(ii=varSam+1;ii<lenTS-varSam;ii++)
    {
        mean[ii] =  mean[ii-1] - samPtr[ii-(varSam+1)].value + samPtr[ii+varSam].value ;  //running mean
    }
    for(ii=0;ii<lenTS;ii++)
    {
        mean[ii] =  mean[ii]/N;  //running mean
    }
    
    //computing running variance
    for(ii=0;ii<N;ii++)
    {
        temp[0] = (samPtr[ii].value- mean[ii]);
        stdVec[varSam] += temp[0]*temp[0];
    } 
    for(ii=varSam+1;ii<lenTS-varSam;ii++)
    {
        temp[0] = samPtr[ii-(varSam+1)].value -mean[ii-(varSam+1)];
        stdVec[ii] =  stdVec[ii-1] - temp[0]*temp[0] ; 
        temp[1] = (samPtr[ii+varSam].value -mean[ii+varSam]);
        stdVec[ii] += temp[1]*temp[1] ;
    }
    for(ii=0;ii<lenTS;ii++)
    {
        stdVec[ii] =  sqrt(stdVec[ii]/N);
    }    
    
    *std = stdVec;
    free(mean);
    
}


int updateBLStdThsld()
{
    //computing standard deviation
    float *stdVec
    int varSam = (int)round(procParams.varDur/pHop);
    computeStdTSLocal(&stdVec, varSam);
    
    // Assigning weather a point belongs to a flat region or non flar region based on a threhsold
    for(int ii=varSam;ii<nPitchSamples-varSam;ii++)
    {
        if (stdVec[ii]>myProcParams->threshold)
        {
            stdVec[ii] = 1;
        }
        else
        {
            stdVec[ii] = 0;
        }            
            
    }
    
     //blackCnt=lenTS;
    TSAIND ind=0;
    float ex=0;
    int lenRawMotifData = procParams.motifLengths[procParams.indexMotifLenLongest];
    int lenRawMotifDataM1 = lenRawMotifData-1;
    
    //############################## generating subsequences ##########################################
    while(ind<lenTS)
    {
        ex+=stdVec[ind];
        if (ind >= lenRawMotifDataM1)
        {
            if (blacklist[ind-lenRawMotifDataM1]==0)
            {
                if(ex < procParams.flatThreshold*(float)lenRawMotifData)
                {
                    
                    blacklist[ind-lenRawMotifDataM1]=1;
                }
                
            }
            ex -= stdVec[ind-lenRawMotifDataM1];
        }
       ind++;        
    }
    
    return 1;
}

int TSAdataHandler::updateBLDurThsld()
{
    float start, end, maxMotifDur;
    
    maxMotifDur = (procParams.durMotif*procParams.interpFac[procParams.indexMotifLenLongest]) + procParams.maxPauseDur ; 
    
    for(int ind=0; ind<nSubSeqs; ind++)
    {
        start = samPtr[ind].tStamp;
        end = samPtr[ind+procParams.motifLengths[procParams.indexMotifLenLongest]].tStamp;
        
        if (fabs(start-end) > maxMotifDur)       //allow 200 ms pauses in total not more than that
        {
            blacklist[ind]=1;
        }
    }
    
    return 1;
}



int TSAdataHandler::initializeBlackList()
{
    blacklist = (int *)malloc(sizeof(int)*nSubSeqs);
    
    return 1;
}

int TSAdataHandler::genSlidingWindowSubSeqs()
{
    
    
    int lenRawMotifData = procParams.motifLengths[procParams.indexMotifLenLongest];
    int lenRawMotifDataM1 = lenRawMotifData-1;
    
    
    subSeqPtr = (TSAsubSeq_t *)malloc(sizeof(TSAsubSeq_t)*lenTS);
    
    TSAIND ind=0;
    float ex=0;
    
    //############################## generating subsequences ##########################################
    while(ind<lenTS)
    {
        subSeqPtr[ind].sTime  = samPtr[ind].tStamp;
        subSeqPtr[ind].eTime  = samPtr[ind+procParams.motifLengths[procParams.indexMotifLenReal]].tStamp;
        
        if (ind < lenTs - lenRawMotifDataM1)
        {
            subSeqPtr[ind].pData = (TSADATA *)malloc(sizeof(TSADATA)*lenRawMotifData);
        }
        
        for(TSAIND ll = min(ind, lenTS - lenRawMotifDataM1-1) ; ll >= max(0,ind-lenRawMotifDataM1) ; ll--)
        {
            subSeqPtr[ll].pData[ind-ll] = samPtr[ind].value; 
        }
        
       ind++;        
    }
    
    nSubSeqs = ind;
    
    return 1;

}

int TSAdataHandler::calculateDiffMotifLengths()
{
    
    if (procParams.nInterpFac==1)
    {
        procParams.interpFac[0]=1.0;
        
        procParams.combMTX = (int **)malloc(sizeof(int*)*procParams.nInterpFac);
        for(ii=0;ii<procParams.nInterpFac;ii++)
        {
            procParams.combMTX[ii] =  (int *)malloc(sizeof(int)*procParams.nInterpFac);
            for(jj=0;jj<procParams.nInterpFac;jj++)
            {
                procParams.combMTX[ii][jj] = 1;
            }
        }
    }
    else if (procParams.nInterpFac==3)
    {
        procParams.interpFac[0]=0.9;
        procParams.interpFac[1]=1.0;
        procParams.interpFac[2]=1.1;
        
        procParams.combMTX = (int **)malloc(sizeof(int*)*procParams.nInterpFac);
        for(ii=0;ii<procParams.nInterpFac;ii++)
        {
            procParams.combMTX[ii] =  (int *)malloc(sizeof(int)*procParams.nInterpFac);
            for(jj=0;jj<procParams.nInterpFac;jj++)
            {
                procParams.combMTX[ii][jj] = combAllwd_3[ii][jj];
            }
        }
        
        
    }
    else if (procParams.nInterpFac==5)
    {
        procParams.interpFac[0]=0.9;
        procParams.interpFac[1]=0.95;
        procParams.interpFac[2]=1.0;
        procParams.interpFac[3]=1.05;
        procParams.interpFac[4]=1.1;
        
        procParams.combMTX = (int **)malloc(sizeof(int*)*procParams.nInterpFac);
        for(ii=0;ii<procParams.nInterpFac;ii++)
        {
            procParams.combMTX[ii] =  (int *)malloc(sizeof(int)*procParams.nInterpFac);
            for(jj=0;jj<procParams.nInterpFac;jj++)
            {
                procParams.combMTX[ii][jj] = combAllwd_5[ii][jj];
            }
        }
    }
    
    
    //######### computing all the motif lengths for several interpolation factors ##############
    int max_factor=-1;
    for (ii=0;ii<procParams.nInterpFac; ii++)
    {
        procParams.motifLengths[ii] = (int)ceil((procParams.durMotif*procParams.interpFac[ii])/pHop);
        if (procParams.interpFac[ii]==1)
        {
            procParams.indexMotifLenReal = ii;
        }
        if (procParams.interpFac[ii]>max_factor)
        {
            max_factor = procParams.interpFac[ii];
            procParams.indexMotifLenLongest = ii;
        }
    }
    
    //CRUCIAL POINT !!! since cubic interpolation needs 4 points (2 ahead) just store
    if (procParams.nInterpFac >1)
    {
        procParams.motifLengths[procParams.indexMotifLenLongest]+=1;
    }
    
}

int TSAdataHandler::convertHz2Cents(char *tonicFileName)
{
    float tonic;
    
    //Opening tonic file
    FILE *fp =fopen(tonicFileName,"r");
    if (fp==NULL)
    {
        printf("Error opening file %s\n", tonicFileName);
        return 0;
    }
    int nRead = fscanf(fp, "%f\n",&tonic);
    fclose(fp);
    
    float temp1 = ((float)procParams.binsPOct)/LOG2;
    
    for(int ii=0;ii<lenTS;ii++)
    {
        samPtr[ii].value = (TSADATA)(temp1*log((samPtr[ii].value+EPS)/tonic));
    }
    
    return 1;
}


int TSAdataHandler::filterSamplesTS()
{
    
    // There can be many filtering criterion. Currently what is implemented is to filter out silence regions of pitch sequence. 
    
    TSAIND nValidSamples = 0;
    TSAIND *validSampleInd = (TSAIND *)malloc(sizeof(TSAIND)*lenTS);
    //store all locations which have valid samples
    for(int ii=0;ii<lenTS; ii++)
    {
        if(samPtr[ii].value > procParams.minPossiblePitch)
        {
            validSampleInd[nValidSamples]=ii;
            nValidSamples++;
        }
    }
    // just copy valid samples to another location
    TSAIND lenTS_new = nValidSamples;
    TSAsam_t * samPtr_new = (TSAsam_t *)malloc(sizeof(TSAsam_t)*lenTS_new);
    for (int ii=0;ii<lenTS_new; ii++)
    {
        samPtr_new[ii] = samPtr[validSampleInd[ii]];
    }
    
    //free old memory
    free(samPtr);
    free(validSampleInd);
    samPtr = samPtr_new;
    lenTS = lenTS_new;
    
    return 1;
}

int TSAdataHandler::readHopSizeTS(char *fileName)
{
    FILE *fp;
    int nRead;
    float tsData;
    float tsTime1, tsTime2;
    
    // Opening time series file (JUST TO OBTAIN HOP SIZE)
    fp =fopen(fileName,"r");
    if (fp==NULL)
    {
        printf("Error opening file %s\n", fileName);
        return 0;
    }
    //reading just first two lines, in order to obtain hopsize//
    nRead = fscanf(fp, "%f\t%f\n",&tsData,&tsTime1);
    nRead = fscanf(fp, "%f\t%f\n",&tsData,&tsTime2);
    pHop = (tsTime2-tsTime1);
    fclose(fp);
    return 1;
}
int TSAdataHandler::downSampleTS()
{
    if (procParams.dsFactor <=0)
    {
        return 1;
    }
    
    TSAIND lenTS_new = ceil(lenTS/procParams.dsFactor);
    TSAsam_t * samPtr_new = (TSAsam_t *)malloc(sizeof(TSAsam_t)*lenTS_new);
    for (int ii=0;ii<lenTS_new; ii++)
    {
        samPtr_new[ii] = samPtr[ii*procParams.dsFactor];
    }
    
    //correcting hop size as well
    pHop = pHop*procParams.dsFactor;
    
    //free old memory
    free(samPtr);
    samPtr = samPtr_new;
    lenTS = lenTS_new;
    
    return 1;
}

int TSAdataHandler::readTSData(char *fileName)
{
    FILE *fp;
    float tsData;
    float tsTime;
    TSAIND cnt;
    
    
    //count number of lines in the file specified
    nLinesFile = getNumLines(fHandle.getTSFileName());
    procLogPtr->totalPitchSamples += nLinesFile;
    lenTS = nLinesFile;
    
    //allocate memory to store the time series data 
    samPtr = (TSAsam_t *)malloc(sizeof(TSAsam_t)*lenTS);
    fp =fopen(fileName,"r");
    if (fp==NULL)
    {
        printf("Error opening file %s\n", fileName);
        return 0;
    }
    cnt = 0;
    while (fscanf(fp, "%f\t%f\n",&tsData,&tsTime)!=EOF)    //read till the end of the file
    {
        samPtr[cnt].value  = tsData;
        samPtr[cnt].tStamp  = tsTime;
        cnt++;
    }
    fclose(fp);
    
    return 1;
}
    
    

/*
 * This function quickly returns number of lines in a text file
 * It uses memory mapp methods to perform this task. Its super fast as compared to reading and obtaining number of lines. 
 * By reading number of lines beforehand we can know how many lines are there in pitch file and we can then pre-allocate the memory to store pitch samples. Otherwise adaptively its really slow to read data (because that requires rellocation of buffers).
 */
TSAIND TSAdataHandler::getNumLines(const char *file)
{
    int fp;
    struct stat fs;
    char *buf;
    TSAIND line=0, ii;
    
    fp = open(file,O_RDONLY);
    if (fp == -1) 
    {
                printf("Error opening file1 %s\n",file);
                return 0;
    }
    if (fstat(fp, &fs) == -1)
    {
                printf("Error opening file2 %s\n",file);
                return 0;
    }
    buf = (char*)mmap(0, fs.st_size, PROT_READ, MAP_SHARED, fp, 0);
    
    for (ii=0;ii<fs.st_size;ii++)
    {
        if (buf[ii]=='\n')
            line++;
    }
    
    munmap(buf, fs.st_size);
    close(fp);
    
    return line;
    
}

fileNameHandler::fileNameHandler()
{

}
int fileNameHandler::initialize(char *bName, fileExts_t *fExtPtr)
{
    baseName = bName;
    fileExtPtr = fExtPtr;
    
    return 1;
}

char* fileNameHandler::getTSFileName()
{
    memset(fileName, '\0', sizeof(char)*MAX_FNAME_CHARS);
    strcat(fileName, baseName);
    strcat(fileName, fileExtPtr->tsFileExt);
    
    return fileName;
    
}
char* fileNameHandler::getTonicFileName()
{
    memset(fileName, '\0', sizeof(char)*MAX_FNAME_CHARS);
    strcat(fileName, baseName);
    strcat(fileName, fileExtPtr->tonicExt);
    
    return fileName;
    
}
char* fileNameHandler::getBlackListSegFileName()
{
    memset(fileName, '\0', sizeof(char)*MAX_FNAME_CHARS);
    strcat(fileName, baseName);
    strcat(fileName, fileExtPtr->blackTimeExt);
    
    return fileName;
    
}
char* fileNameHandler::getLogFileName()
{
    memset(fileName, '\0', sizeof(char)*MAX_FNAME_CHARS);
    strcat(fileName, baseName);
    strcat(fileName, fileExtPtr->logFileExt);
    
    return fileName;
    
}
char* fileNameHandler::getParamDumpFileName()
{
    memset(fileName, '\0', sizeof(char)*MAX_FNAME_CHARS);
    strcat(fileName, baseName);
    strcat(fileName, fileExtPtr->paramsDumpExt);
    
    return fileName;
    
}
char* fileNameHandler::getOutFileName()
{
    memset(fileName, '\0', sizeof(char)*MAX_FNAME_CHARS);
    strcat(fileName, baseName);
    strcat(fileName, fileExtPtr->outFileExt);
    
    return fileName;
    
}
  

