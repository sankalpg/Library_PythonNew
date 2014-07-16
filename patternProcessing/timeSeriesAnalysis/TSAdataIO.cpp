
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
  

