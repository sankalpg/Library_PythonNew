
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
    memset(fileExts.mappFileExt, '\0', sizeof(char)*MAX_FEXT_CHARS);
    memset(fileExts.searchListExt, '\0', sizeof(char)*MAX_FEXT_CHARS);
    memset(fileExts.queryFileExt, '\0', sizeof(char)*MAX_FEXT_CHARS);
    memset(fileExts.subSeqFileExt, '\0', sizeof(char)*MAX_FEXT_CHARS);
    memset(fileExts.subSeqTNFileExt, '\0', sizeof(char)*MAX_FEXT_CHARS);
    memset(fileExts.subSeqInfoFileExt, '\0', sizeof(char)*MAX_FEXT_CHARS);
    memset(fileExts.patternKNNExt, '\0', sizeof(char)*MAX_FEXT_CHARS);
    memset(fileExts.FNFileExt, '\0', sizeof(char)*MAX_FEXT_CHARS);

    //I need to initialize other structure members as well, #TODO
    procParams.distParams.distNormType = PATH_LEN;
    procParams.methodVariant=Var1;
    procParams.complexityMeasure=-1;

}

TSAparamHandle::~TSAparamHandle()
{
    free(procParams.simMeasureRankRefinement);
    for(int ii=0; ii< procParams.pattParams.nInterpFac; ii++)
    {
        free(procParams.combMTX[ii]);
    }
    free(procParams.combMTX);
}

int TSAparamHandle::readParamsFromFile(char *paramFile)
{
    FILE *fp;
    char tempFilename[400]={'\0'};
    char field[100]={'\0'};
    char value[100]={'\0'};

    //Initilizations
    procParams.repParams.normType=0;



    fp = fopen(paramFile, "r");
    if (fp == NULL){printf("Error opening file %s\n", paramFile);exit(0);}
    while(fgets(tempFilename, 400, fp))
    {
        sscanf(tempFilename, "%s %s\n", field, value);
        /*Distane related parameters*/
        if (strcmp(field, "distType:")==0){procParams.distParams.distType=atoi(value);}
        if (strcmp(field, "DTWBand:")==0){procParams.distParams.DTWBand=atof(value);}
        if (strcmp(field, "DTWType:")==0){procParams.distParams.DTWType=atoi(value);}
        if (strcmp(field, "rankRefDistType:")==0){procParams.distParams.rankRefDistType=atoi(value);}
        if (strcmp(field, "distNormType:")==0){procParams.distParams.distNormType=atoi(value);}
        /*representation related parameters*/
        if (strcmp(field, "TSRepType:")==0){procParams.repParams.TSRepType=atoi(value);}
        if (strcmp(field, "quantSize:")==0){procParams.repParams.quantSize=atoi(value);}
        if (strcmp(field, "normType:")==0){procParams.repParams.normType=atoi(value);}
        if (strcmp(field, "sampleRate:")==0){procParams.repParams.sampleRate=atof(value);}
        if (strcmp(field, "binsPOct:")==0){procParams.repParams.binsPOct=atoi(value);}
        if (strcmp(field, "minPossiblePitch:")==0){procParams.repParams.minPossiblePitch=atof(value);}
        if (strcmp(field, "removeTaniSegs:")==0){procParams.repParams.removeTaniSegs=atoi(value);}
        if (strcmp(field, "varDur:")==0){procParams.repParams.varDur=atof(value);}
        if (strcmp(field, "threshold:")==0){procParams.repParams.threshold=atof(value);}
        if (strcmp(field, "flatThreshold:")==0){procParams.repParams.flatThreshold=atof(value);}
        if (strcmp(field, "maxPauseDur:")==0){procParams.repParams.maxPauseDur=atof(value);}
        /*pattern processing parameters*/
        if (strcmp(field, "durMotif:")==0){procParams.pattParams.durMotif=atof(value);}
        if (strcmp(field, "blackDurFact:")==0){procParams.pattParams.blackDurFact=atof(value);}
        if (strcmp(field, "maxNMotifsPairs:")==0){procParams.pattParams.maxNMotifsPairs=atoi(value);}
        if (strcmp(field, "nInterpFac:")==0){procParams.pattParams.nInterpFac=atoi(value);}
        if (strcmp(field, "pitchHop:")==0){procParams.repParams.pitchHop=atof(value);}
        if (strcmp(field, "subSeqLen:")==0){procParams.pattParams.subSeqLen=atoi(value);}
        if (strcmp(field, "methodVariant:")==0){procParams.methodVariant=atoi(value);}
        if (strcmp(field, "complexityMeasure:")==0){procParams.complexityMeasure=atoi(value);}

        if (strcmp(field, "dumpLogs:")==0){procParams.dumpLogs=atoi(value);}


        if (procParams.pattParams.durMotif>0)
        {
            procParams.pattParams.blackDur=procParams.pattParams.blackDurFact*procParams.pattParams.durMotif;
        }

    }
    fclose(fp);

    if (procParams.distParams.rankRefDistType>0)
    {
        procParams.simMeasureRankRefinement = (int*)malloc(sizeof(int)*1);
        procParams.simMeasureRankRefinement[0] = procParams.distParams.rankRefDistType;
        procParams.nSimMeasuresUsed =1;
    }
    else    //In such case use 4 different similarity measures needed for experimentation
    {
        procParams.simMeasureRankRefinement = (int*)malloc(sizeof(int)*4);
        procParams.simMeasureRankRefinement[0] = SqEuclidean;
        procParams.simMeasureRankRefinement[1] = CityBlock;
        procParams.simMeasureRankRefinement[2] = ShiftCityBlock;
        procParams.simMeasureRankRefinement[3] = ShiftLinExp;
        procParams.nSimMeasuresUsed =4;
    }

    int ii, jj;

    if (procParams.pattParams.nInterpFac==1)
    {
        procParams.interpFac[0]=1.0;

        procParams.combMTX = (int **)malloc(sizeof(int*)*procParams.pattParams.nInterpFac);
        for(ii=0;ii<procParams.pattParams.nInterpFac;ii++)
        {
            procParams.combMTX[ii] =  (int *)malloc(sizeof(int)*procParams.pattParams.nInterpFac);
            for(jj=0;jj<procParams.pattParams.nInterpFac;jj++)
            {
                procParams.combMTX[ii][jj] = 1;
            }
        }
    }
    else if (procParams.pattParams.nInterpFac==3)
    {
        procParams.interpFac[0]=0.9;
        procParams.interpFac[1]=1.0;
        procParams.interpFac[2]=1.1;

        procParams.combMTX = (int **)malloc(sizeof(int*)*procParams.pattParams.nInterpFac);
        for(ii=0;ii<procParams.pattParams.nInterpFac;ii++)
        {
            procParams.combMTX[ii] =  (int *)malloc(sizeof(int)*procParams.pattParams.nInterpFac);
            for(jj=0;jj<procParams.pattParams.nInterpFac;jj++)
            {
                procParams.combMTX[ii][jj] = combAllwd_3[ii][jj];
            }
        }


    }
    else if (procParams.pattParams.nInterpFac==5)
    {
        procParams.interpFac[0]=0.9;
        procParams.interpFac[1]=0.95;
        procParams.interpFac[2]=1.0;
        procParams.interpFac[3]=1.05;
        procParams.interpFac[4]=1.1;

        procParams.combMTX = (int **)malloc(sizeof(int*)*procParams.pattParams.nInterpFac);
        for(ii=0;ii<procParams.pattParams.nInterpFac;ii++)
        {
            procParams.combMTX[ii] =  (int *)malloc(sizeof(int)*procParams.pattParams.nInterpFac);
            for(jj=0;jj<procParams.pattParams.nInterpFac;jj++)
            {
                procParams.combMTX[ii][jj] = combAllwd_5[ii][jj];
            }
        }
    }

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
        if (strcmp(field, "mappFileExt:")==0){strcat(fileExts.mappFileExt, value);}
        if (strcmp(field, "searchListExt:")==0){strcat(fileExts.searchListExt, value);}
        if (strcmp(field, "queryFileExt:")==0){strcat(fileExts.queryFileExt, value);}
        if (strcmp(field, "subSeqFileExt:")==0){strcat(fileExts.subSeqFileExt, value);}
        if (strcmp(field, "subSeqTNFileExt:")==0){strcat(fileExts.subSeqTNFileExt, value);}
        if (strcmp(field, "subSeqInfoFileExt:")==0){strcat(fileExts.subSeqInfoFileExt, value);}
        if (strcmp(field, "patternKNNExt:")==0){strcat(fileExts.patternKNNExt, value);}
        if (strcmp(field, "FNFileExt:")==0){strcat(fileExts.FNFileExt, value);}

        memset(tempFilename, '\0', sizeof(char)*400);
        memset(field, '\0', sizeof(char)*100);
        memset(value, '\0', sizeof(char)*100);
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

TSAdataHandler::TSAdataHandler(char *bName, procLogs_t *procLogs, fileExts_t *fileExts, procParams_t *pParams)
{

    procLogPtr = procLogs;
    fileExtPtr = fileExts;
    procParams = *pParams;
    baseName = bName;

    fHandle.initialize(baseName, fileExtPtr);

    pHop = procParams.repParams.pitchHop; //TODO: this pHop should be replaced by procParams.repParams.pitchHop everywhere

    isBlackListAlloc=-1;
    nQueries=-1;
    nSubSeqs=-1;
    lenTS =-1;
    nPatternPairs=-1;

}
TSAdataHandler::~TSAdataHandler()
{
    if(lenTS!=-1)
    {
        free(samPtr);
    }
    if(nSubSeqs!=-1)
    {
        freeSubSeqsMem();
    }

    if (isBlackListAlloc==1)
    {free(blacklist);}

    if(nQueries>0)
    {free(queryTStamps);}

    if(nPatternPairs>0)
    {free(patternPairs);}

}

int TSAdataHandler::freeSubSeqsMem()
{
    for(TSAIND ii=0; ii< nSubSeqs; ii++)
    {
        free(subSeqPtr[ii].pData);
        free(subSeqPtr[ii].pTStamps);
    }
    free(subSeqPtr);
    nSubSeqs=-1;
}


int TSAdataHandler::loadMotifDataTemplate1()
{
     //read the time series data
    readTSData(fHandle.getTSFileName());
    readHopSizeTS(fHandle.getTSFileName());

    //downsample
    downSampleTS();

    //remove silence pitch regions
    filterSamplesTS();

    convertHz2Cents(fHandle.getTonicFileName());

    //calculate different motif lengths before doing sliding window candidate generation
    calculateDiffMotifLengths();

    readQueryTimeStamps(fHandle.getQueryFileName(), MOTIFPAIR_DUMP_FORMAT);

    genSubSeqsWithTStarts(queryTStamps, nQueries);

    genUniScaledSubSeqs();

    computeMeanSTDSubSeqs(procParams.motifLengths[procParams.indexMotifLenReal]);

    normalizeSubSeqs(procParams.repParams.normType);


}

int TSAdataHandler::genSubSeqsWithTStamps(TSAseg_t *qTStamps, TSAIND nQueries)
{
    TSAIND *seedMotifStrInd = (TSAIND*)malloc(sizeof(TSAIND)*nQueries);
    TSAIND *seedMotifEndInd = (TSAIND*)malloc(sizeof(TSAIND)*nQueries);

    for(TSAIND ii=0;ii<nQueries;ii++)
    {
        TSADIST min_valS = INF;
        TSADIST min_valE = INF;
        for(TSAIND jj=0;jj<lenTS;jj++)
        {
            if (fabs(samPtr[jj].tStamp-qTStamps[ii].sTime) < min_valS)
            {
                min_valS = fabs(samPtr[jj].tStamp-qTStamps[ii].sTime);
                seedMotifStrInd[ii] = jj;
            }
            if (fabs(samPtr[jj].tStamp-qTStamps[ii].eTime) < min_valE)
            {
                min_valE = fabs(samPtr[jj].tStamp-qTStamps[ii].eTime);
                seedMotifEndInd[ii] = jj;
            }
        }
    }
    subSeqPtr = (TSAsubSeq_t *)malloc(sizeof(TSAsubSeq_t)*nQueries);

    for(TSAIND ii=0;ii<nQueries;ii++)
    {

        //whats the length of this motif
        procParams.pattParams.durMotif = samPtr[seedMotifEndInd[ii]].tStamp - samPtr[seedMotifStrInd[ii]].tStamp;
        calculateDiffMotifLengths();
        int lenRawMotifData = procParams.motifLengths[procParams.indexMotifLenLongest];

        subSeqPtr[ii].pData = (TSADATA *)malloc(sizeof(TSADATA)*lenRawMotifData);
        subSeqPtr[ii].pTStamps = (float *)malloc(sizeof(float)*lenRawMotifData);

        //at this point in time just store the original length start end without uniform scaling
        subSeqPtr[ii].sTime = samPtr[seedMotifStrInd[ii]].tStamp;
        subSeqPtr[ii].eTime = samPtr[seedMotifEndInd[ii]].tStamp;

        for(TSAIND jj=0;jj<lenRawMotifData;jj++)
        {
            subSeqPtr[ii].pData[jj]=samPtr[seedMotifStrInd[ii]+jj].value;
            subSeqPtr[ii].pTStamps[jj]=samPtr[seedMotifStrInd[ii]+jj].tStamp;
        }
    }
    nSubSeqs = nQueries;
    free(seedMotifStrInd);
    free(seedMotifEndInd);

    return 1;


}

int TSAdataHandler::genSubSeqsWithTStarts(TSAseg_t *qTStamps, TSAIND nQueries)
{
    int lenRawMotifData = procParams.motifLengths[procParams.indexMotifLenLongest];

    TSAIND *seedMotifInd = (TSAIND*)malloc(sizeof(TSAIND)*nQueries);

    for(TSAIND ii=0;ii<nQueries;ii++)
    {
        TSADIST min_val = INF;
        for(TSAIND jj=0;jj<lenTS;jj++)
        {
            if (fabs(samPtr[jj].tStamp-qTStamps[ii].sTime)<min_val)
            {
                min_val = fabs(samPtr[jj].tStamp-qTStamps[ii].sTime);
                seedMotifInd[ii] = jj;
            }
        }
    }

    subSeqPtr = (TSAsubSeq_t *)malloc(sizeof(TSAsubSeq_t)*nQueries);

    for(TSAIND ii=0;ii<nQueries;ii++)
    {

        subSeqPtr[ii].pData = (TSADATA *)malloc(sizeof(TSADATA)*lenRawMotifData);
        subSeqPtr[ii].pTStamps = (float *)malloc(sizeof(float)*lenRawMotifData);

        subSeqPtr[ii].sTime = samPtr[seedMotifInd[ii]].tStamp;
        subSeqPtr[ii].eTime = samPtr[seedMotifInd[ii]+lenRawMotifData-1].tStamp;

        for(TSAIND jj=0;jj<lenRawMotifData;jj++)
        {
            subSeqPtr[ii].pData[jj]=samPtr[seedMotifInd[ii]+jj].value;
            subSeqPtr[ii].pTStamps[jj]=samPtr[seedMotifInd[ii]+jj].tStamp;
        }
    }
    nSubSeqs = nQueries;
    free(seedMotifInd);

    return 1;


}

int TSAdataHandler::readKNNPatternDump(char *patternKNNFile, int format)
{
    FILE *fp;
    TSAmotifInfo_t *pPairs;
    TSAIND nPairs;

    if (format == MOTIFID1_MOTIFID2_DIST)
    {
        float temp1;
        TSAIND temp[2]={0};
        TSAIND ii=0;
        fp = fopen(patternKNNFile,"r");
        if (fp==NULL)
        {
            printf("Error opening file %s\n", patternKNNFile);
            return 0;
        }
        int nLines = getNumLines(patternKNNFile);
        pPairs = (TSAmotifInfo_t*)malloc(sizeof(TSAmotifInfo_t)*nLines);
        while(fscanf(fp,"%lld\t%lld\t%f\n", &temp[0], &temp[1],&temp1)!=EOF)
        {
            if(ii>=nLines)
                break;
            pPairs[ii].ind1 = temp[0];
            pPairs[ii].ind2 = temp[1];
            pPairs[ii].dist = temp1;
            ii++;
        }
        fclose(fp);
        nPatternPairs=ii;
        patternPairs = pPairs;
    }

}

int TSAdataHandler::readQueryTimeStamps(char *queryFileName, int format)
{
    FILE *fp;
    TSAseg_t *qTStamps;

    if (format == MOTIFPAIR_DUMP_FORMAT)
    {
        float temp[10]={0};
        TSAIND ii=0, jj=0;
        TSADIST temp4;
        fp = fopen(queryFileName,"r");
        if (fp==NULL)
        {
            printf("Error opening file %s\n", queryFileName);
            return 0;
        }
        //reading number of lines in the file
        int nLines = getNumLines(queryFileName);
        qTStamps = (TSAseg_t*)malloc(sizeof(TSAseg_t)*nLines*2);

        while(fscanf(fp,"%f\t%f\t%f\t%f\t%lf\t%f\t%f\n", &temp[0], &temp[1], &temp[2], &temp[3], &temp4, &temp[5], &temp[6])!=EOF)
        {
            if(ii>=nLines)
                break;

            if(temp4<INF)
            {
                qTStamps[jj].sTime=temp[0];
                qTStamps[jj].eTime=temp[1];
                jj++;

                qTStamps[jj].sTime=temp[2];
                qTStamps[jj].eTime=temp[3];
                jj++;
            }
            ii++;
        }
        fclose(fp);
        nQueries=jj;
        queryTStamps = qTStamps;
    }
    else if (format == VIGNESH_MOTIF_ANNOT_FORMAT)
    {
        float temp[10]={0};
        TSAIND ii=0, jj=0;
        char temp4[10]={'\0'};
        fp = fopen(queryFileName,"r");
        if (fp==NULL)
        {
            printf("Error opening file %s\n", queryFileName);
            return 0;
        }
        //reading number of lines in the file
        int nLines = getNumLines(queryFileName);
        qTStamps = (TSAseg_t*)malloc(sizeof(TSAseg_t)*nLines*1);

        while(fscanf(fp,"%f %f %s\n", &temp[0], &temp[1], temp4)!=EOF)
        {
            qTStamps[jj].sTime=temp[0];
            qTStamps[jj].eTime=temp[1]+qTStamps[jj].sTime;
            jj++;
        }
        fclose(fp);
        nQueries=jj;
        queryTStamps = qTStamps;
    }
    else if (format == MY_MOTIF_ANNOT_FORMAT)
    {
        float temp[10]={0};
        TSAIND ii=0, jj=0;
        int temp4;
        fp = fopen(queryFileName,"r");
        if (fp==NULL)
        {
            printf("Error opening file %s\n", queryFileName);
            return 0;
        }
        //reading number of lines in the file
        int nLines = getNumLines(queryFileName);
        qTStamps = (TSAseg_t*)malloc(sizeof(TSAseg_t)*nLines*1);

        while(fscanf(fp,"%f\t%f\t%d\n", &temp[0], &temp[1], &temp4)!=EOF)
        {
            qTStamps[jj].sTime=temp[0];
            qTStamps[jj].eTime=temp[1];
            jj++;
        }
        fclose(fp);
        nQueries=jj;
        queryTStamps = qTStamps;
    }
    else if (format == PATTERNS_PER_FILE_DUMP)
    {
        float temp[10]={0};
        TSAIND ii=0, jj=0;
        long long int temp4;
        fp = fopen(queryFileName,"r");
        if (fp==NULL)
        {
            printf("Error opening file %s\n", queryFileName);
            return 0;
        }
        //reading number of lines in the file
        int nLines = getNumLines(queryFileName);
        qTStamps = (TSAseg_t*)malloc(sizeof(TSAseg_t)*nLines*1);

        while(fscanf(fp,"%lld\t%f\t%f\n", &temp4, &temp[0], &temp[1])!=EOF)
        {
            qTStamps[jj].sTime=temp[0];
            qTStamps[jj].eTime=temp[1];
            qTStamps[jj].id = temp4;
            jj++;
        }
        fclose(fp);
        nQueries=jj;
        queryTStamps = qTStamps;
    }


    return 1;

}

int TSAdataHandler::genTemplate1SubSeqs()
{
    int rVal = 0;
    //read the time series data
    rVal = readTSData(fHandle.getTSFileName());
    if (rVal==0)
    {return 1;}

    readHopSizeTS(fHandle.getTSFileName());

    //downsample
    downSampleTS();
    
    //remove silence pitch regions
    filterSamplesTS();

    procLogPtr->lenProcTS+=lenTS;

    convertHz2Cents(fHandle.getTonicFileName());

    //calculate different motif lengths before doing sliding window candidate generation
    calculateDiffMotifLengths();

    genSlidingWindowSubSeqs();

    procLogPtr->nSubSeqs+=nSubSeqs;

    initializeBlackList(nSubSeqs);

    updateBLDurThsld();
    
    //updateBLStdThsld();
    // Now we have a new function to find blacklisted subs which have a lot of flat segment in it.
    updateBlackListFlatSeqs(fHandle.getFlatNoteFileName());

    updateBLInvalidSegment(fHandle.getBlackListSegFileName());

    filterBlackListedSubSeqs();

    procLogPtr->nSubSeqsBL+=procLogPtr->nSubSeqs - nSubSeqs;

    genUniScaledSubSeqs();

    computeMeanSTDSubSeqs(procParams.motifLengths[procParams.indexMotifLenReal]);

    normalizeSubSeqs(procParams.repParams.normType);

    procLogPtr->nProcSubSeq+=nSubSeqs;

}
/*
This function computes an offset between a pair of subsequences based on the normalization type.
For example: if we choose to perform octave normalization then if the two subseqs are detection to be octave transposed
the offset would be N*1200. Similarly another normalization strategy can be PASA norm.
*/
float TSAdataHandler::estimateOffset(float mean1, float mean2)
{
    float diff, dev, rounded_diff;
    if (procParams.repParams.normType == OCTAVE_NORM)
    {
        diff = mean1 - mean2;
        rounded_diff = round(diff/procParams.repParams.binsPOct)*procParams.repParams.binsPOct;
        dev = diff -rounded_diff;
        if (fabs(dev) <= procParams.repParams.binsPOct/8.0)
        {  //note that 1/8 bins per octave is 150 cents equivalent
            return rounded_diff;
        }
        else
        {
            return 0.0;
        }
    }

    return 0.0;
}


/*
This function as of now is used by patternDistancesFixedDuration.cpp for a small amount of data (as it is applied in the context
where there is no lower bounding). The idea is to replace this optional step which is controlled by #define to a parameter
for octave and PASA normalization. That would be by integrating this procedurein the lower bounding also.
*/
int TSAdataHandler::copyAndNormalizeSubSeqsPASA(TSADATA *out1, TSADATA *out2, TSADATA *inp1, TSADATA *inp2, float mean1, float mean2, int len)
{
    memcpy(out1, inp1, sizeof(TSADATA)*len);
    memcpy(out2, inp2, sizeof(TSADATA)*len);

    if (procParams.repParams.normType != TONIC_NORM_PASAPA)
    {
        return 1;
    }

    float diff = mean1-mean2;
    float sign=1;
    if (diff!=0)
    { sign= diff/fabs(diff);}

    float offset=0;
    if ((fabs(diff)>550) && (fabs(diff)<850))
    {
        offset=sign*700;
    }
    else if  ((fabs(diff)>1050) && (fabs(diff)<1350))
    {
        offset=sign*1200;
    }
    else if((fabs(diff)>1750) && (fabs(diff)<2050))
    {
        offset=sign*1900;
    }

    for(int ii=0;ii<len;ii++)
    {
        out1[ii]-=offset;
    }

    return 1;
}

int TSAdataHandler::genUniScaledSubSeqsVarLen()
{
    float t1,t2;
    t1 = clock();

    int nInterpFac = procParams.pattParams.nInterpFac;
    TSAIND nSubSeqs_new = nSubSeqs*nInterpFac;

    TSAsubSeq_t *subSeqPtr_new = (TSAsubSeq_t *)malloc(sizeof(TSAsubSeq_t)*nSubSeqs_new);

    int max_pos_len = (int)floor(procParams.pattParams.subSeqLen/procParams.interpFac[nInterpFac-1]);

    //we do interpolation as well
    float **indInterp = (float**)malloc(sizeof(float*)*nInterpFac);
    for (int jj=0;jj<nInterpFac;jj++)
    {
        indInterp[jj] = (float *)malloc(sizeof(float)*max_pos_len);
        for (int ii=0;ii<max_pos_len; ii++)
        {
            indInterp[jj][ii] = procParams.interpFac[jj]*ii ;
        }
    }
    TSAIND ind=0;
    for (TSAIND ii=0;ii<nSubSeqs;ii++)
    {
        for (int jj=0;jj<nInterpFac;jj++)
        {
            if (procParams.interpFac[jj]==1)
            {
                subSeqPtr_new[ind] = subSeqPtr[ii];
                subSeqPtr_new[ind].len = subSeqPtr[ii].len;
                subSeqPtr_new[ind].id = subSeqPtr[ii].id;
                subSeqPtr_new[ind].sTime = subSeqPtr[ii].sTime;
                subSeqPtr_new[ind].eTime = subSeqPtr[ii].eTime;
                ind++;
            }
            else
            {
                subSeqPtr_new[ind].pData = (TSADATA *)malloc(sizeof(TSADATA)*max_pos_len);
                subSeqPtr_new[ind].pTStamps = (float *)malloc(sizeof(float)*1); //just a dummy allocation
                 //just a dummy allocation
                subSeqPtr_new[ind].len = subSeqPtr[ii].len;
                subSeqPtr_new[ind].id = subSeqPtr[ii].id;
                subSeqPtr_new[ind].sTime = subSeqPtr[ii].sTime;
                subSeqPtr_new[ind].eTime = subSeqPtr[ii].eTime;
                //subSeqPtr_new[ind].sTime  = subSeqPtr[ii].pTStamps[0];
                //subSeqPtr_new[ind].eTime  = subSeqPtr[ii].pTStamps[procParams.motifLengths[jj]-1];
                cubicInterpolate(subSeqPtr[ii].pData, subSeqPtr_new[ind].pData, indInterp[jj], max_pos_len);
                ind++;
            }
        }
    }

    for (int jj=0;jj<nInterpFac;jj++)
    {
        free(indInterp[jj]);
    }
    free(indInterp);

    free(subSeqPtr);
    subSeqPtr = subSeqPtr_new;
    nSubSeqs = nSubSeqs_new;
    procParams.pattParams.subSeqLen = max_pos_len;

    t2 = clock();
    procLogPtr->tGenUniScal += (t2-t1)/CLOCKS_PER_SEC;

    return 1;

}


int TSAdataHandler::genUniScaledSubSeqs()
{
    float t1,t2;
    t1 = clock();

    int nInterpFac = procParams.pattParams.nInterpFac;
    int lenMotifReal = procParams.motifLengths[procParams.indexMotifLenReal];

    TSAIND nSubSeqs_new = nSubSeqs*nInterpFac;

    TSAsubSeq_t *subSeqPtr_new = (TSAsubSeq_t *)malloc(sizeof(TSAsubSeq_t)*nSubSeqs_new);


    //we do interpolation as well
    float **indInterp = (float**)malloc(sizeof(float*)*nInterpFac);
    for (int jj=0;jj<nInterpFac;jj++)
    {
        indInterp[jj] = (float *)malloc(sizeof(float)*lenMotifReal);
        for (int ii=0;ii<lenMotifReal; ii++)
        {
            indInterp[jj][ii] = procParams.interpFac[jj]*ii ;
        }
    }
    TSAIND ind=0;
    for (TSAIND ii=0;ii<nSubSeqs;ii++)
    {
        for (int jj=0;jj<nInterpFac;jj++)
        {
            if (procParams.interpFac[jj]==1)
            {
                subSeqPtr_new[ind] = subSeqPtr[ii];
                subSeqPtr_new[ind].len = lenMotifReal;
                subSeqPtr_new[ind].sTime  = subSeqPtr[ii].pTStamps[0];
                subSeqPtr_new[ind].eTime  = subSeqPtr[ii].pTStamps[lenMotifReal-1];
                ind++;
            }
            else
            {
                subSeqPtr_new[ind].pData = (TSADATA *)malloc(sizeof(TSADATA)*lenMotifReal);
                subSeqPtr_new[ind].pTStamps = (float *)malloc(sizeof(float)*1); //just a dummy allocation
                subSeqPtr_new[ind].len = lenMotifReal;
                cubicInterpolate(subSeqPtr[ii].pData, subSeqPtr_new[ind].pData, indInterp[jj], lenMotifReal);
                subSeqPtr_new[ind].sTime  = subSeqPtr[ii].pTStamps[0];
                subSeqPtr_new[ind].eTime  = subSeqPtr[ii].pTStamps[procParams.motifLengths[jj]-1];
                ind++;
            }
        }
    }

    for (int jj=0;jj<nInterpFac;jj++)
    {
        free(indInterp[jj]);
    }
    free(indInterp);

    free(subSeqPtr);
    subSeqPtr = subSeqPtr_new;
    nSubSeqs = nSubSeqs_new;

    t2 = clock();
    procLogPtr->tGenUniScal += (t2-t1)/CLOCKS_PER_SEC;

    return 1;

}


int TSAdataHandler::filterBlackListedSubSeqs()
{

    //counting number of valid segments
    TSAIND cnt=0;
    for (TSAIND ii=0; ii<nSubSeqs; ii++)
    {
        if(blacklist[ii]==0)
        {
            cnt++;
        }
    }

    TSAIND nSubSeqs_new = cnt;

    TSAsubSeq_t *subSeqPtr_new = (TSAsubSeq_t *)malloc(sizeof(TSAsubSeq_t)*nSubSeqs_new);
    cnt=0;
    for(TSAIND ii=0;ii<nSubSeqs;ii++)
    {
        if(blacklist[ii]==0)
        {
            subSeqPtr_new[cnt] = subSeqPtr[ii];
            cnt++;
        }
        else
        {
            free(subSeqPtr[ii].pData);
            free(subSeqPtr[ii].pTStamps);
        }
    }
    free(subSeqPtr);
    subSeqPtr = subSeqPtr_new;
    nSubSeqs = nSubSeqs_new;

}

int TSAdataHandler::updateBLInvalidSegment(char *segmentFile)
{
    TSAIND ii,jj,nLines;
    FILE *fp;
    float temp[4] = {0};
    TSAPattern_t *invalidSegs;

    nLines = getNumLines(segmentFile);
    invalidSegs = (TSAPattern_t *)malloc(sizeof(TSAPattern_t)*nLines);
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

/*
 * This function computes mean of all the subsequences
 */
int TSAdataHandler::computeMeanSTDSubSeqs(int len)
{
    for(TSAIND ii=0; ii < nSubSeqs; ii++)
    {
        subSeqPtr[ii].mean  = computeMean(subSeqPtr[ii].pData, len);
        subSeqPtr[ii].std  = computeSTD(subSeqPtr[ii].pData, len, subSeqPtr[ii].mean);
    }

}
/*
 * This function updates blacklist array based on flat note information. If subseqs has too much of flat region the blacklist flag is set.
 */
int TSAdataHandler::updateBlackListFlatSeqs(char *FNFile){
    
    //reading the flatnotesegment file
    FILE *fp;
    float temp[4] = {0};
    TSAPattern_t *flatSegs;
    int nLines = getNumLines(FNFile);
    flatSegs = (TSAPattern_t *)malloc(sizeof(TSAPattern_t)*nLines);
    fp =fopen(FNFile,"r");
    if (fp==NULL)
    {
        printf("Error opening file %s\n", FNFile);
        return 0;
    }
    //reading segments of the tani sections
    int ii=0;
    while (fscanf(fp, "%f %f\n",&temp[0],&temp[1])!=EOF)
    {
        flatSegs[ii].sTime = temp[0];
        flatSegs[ii].eTime = temp[1];
        ii++;
    }
    fclose(fp);
    
    int *flatSamples;
    flatSamples = (int *)malloc(sizeof(int)*lenTS);
    memset(flatSamples, 0, sizeof(int)*lenTS);
    TSAIND flat_note_index = 0;
    int cnt=0;
    for(int ii=1; ii< lenTS; ii++){
        if((samPtr[ii].tStamp >= flatSegs[flat_note_index].sTime)&&(samPtr[ii].tStamp <= flatSegs[flat_note_index].eTime)){
            flatSamples[ii] = 1;
        }
        else if(flatSamples[ii-1]==1){
            flat_note_index+=1;
            if((samPtr[ii].tStamp >= flatSegs[flat_note_index].sTime)&&(samPtr[ii].tStamp <= flatSegs[flat_note_index].eTime)){
                flatSamples[ii] = 1;
            }
        }
    }
    TSAIND ind=0;
    float ex=0;
    int lenRawMotifData = procParams.motifLengths[procParams.indexMotifLenLongest];
    int lenRawMotifDataM1 = lenRawMotifData-1;
    while(ind<lenTS)
    {
        ex+=flatSamples[ind];
        if (ind >= lenRawMotifDataM1)
        {
            if (blacklist[ind-lenRawMotifDataM1]==0)
            {
                if(ex < procParams.repParams.flatThreshold*(float)lenRawMotifData)
                {
                    blacklist[ind-lenRawMotifDataM1]=1;
                    cnt+=1;
                }
            }
            ex -= flatSamples[ind-lenRawMotifDataM1];
        }
       ind++;
    }
    free(flatSamples);
    return 1;
}


int TSAdataHandler::updateBLStdThsld()
{
    //computing standard deviation
    float *stdVec;
#ifdef DEBUG    
    float *stdVec_temp;
#endif    
    int varSam = (int)round(procParams.repParams.varDur/pHop);
    computeStdTSLocal(&stdVec, varSam);
#ifdef DEBUG        
    stdVec_temp = (float*)malloc(sizeof(float)*lenTS);
#endif    
    
    // Assigning weather a point belongs to a flat region or non flar region based on a threhsold
    for(int ii=varSam;ii<lenTS-varSam;ii++)
    {
#ifdef DEBUG            
        stdVec_temp[ii] = stdVec[ii];
#endif        
        if (stdVec[ii]>procParams.repParams.threshold)
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
#ifdef DEBUG    
    FILE *fp1;
    fp1 =fopen("dumpDebug.txt","w");
#endif    

    //############################## generating subsequences ##########################################
    while(ind<lenTS)
    {
        ex+=stdVec[ind];
        
        if (ind >= lenRawMotifDataM1)
        {
            
            if (blacklist[ind-lenRawMotifDataM1]==0)
            {
                if(ex < procParams.repParams.flatThreshold*(float)lenRawMotifData)
                {

                    blacklist[ind-lenRawMotifDataM1]=1;
                }

            }
#ifdef DEBUG                
            printf("indexMotifLenReal %d\n",lenRawMotifData);
            fprintf(fp1, "%f\t%f\t%f\t%f\t%d\n", samPtr[ind-lenRawMotifDataM1].tStamp, stdVec_temp[ind - lenRawMotifDataM1], stdVec[ind - lenRawMotifDataM1], ex, blacklist[ind-lenRawMotifDataM1]);
#endif            
            ex -= stdVec[ind-lenRawMotifDataM1];
        }
       ind++;
    }
#ifdef DEBUG        
    fclose(fp1);
#endif    
    free(stdVec);
    return 1;
}

int TSAdataHandler::updateBLDurThsld()
{
    int nSilSamples=0;
    for(int ind=0; ind<nSubSeqs; ind++)
    {
        nSilSamples=0;
        for (int jj=ind+1; jj<ind+procParams.motifLengths[procParams.indexMotifLenLongest]; jj++){
            //NOTE: ind+1 because the first sample shouldn't carry from previous samples the number of silence frames. Imagine this to be the first sample after silence
            nSilSamples+=nSamSil[jj];
        }
        if (nSilSamples*pHop > procParams.repParams.maxPauseDur){
            //printf("Encuontered silence region starttime: %f, endtime: %f, duration of silence: %f\n", start, end, nSilSamples*pHop );
            blacklist[ind]=1;
        }
    }
    return 1;
}

int TSAdataHandler::loadBlackList(char *blackListFile)
{
    FILE *fp;
    int temp;
    TSAIND ii=0;

    fp = fopen(blackListFile,"r");
    if (fp==NULL)
    {
        printf("Error opening file %s\n", blackListFile);
        return 0;
    }
     while(fscanf(fp,"%d\n", &temp)!=EOF)
    {
        blacklist[ii] = temp;
        ii++;
    }
    fclose(fp);

}

int TSAdataHandler::initializeBlackList(TSAIND N)
{
    blacklist = (int *)malloc(sizeof(int)*N);
    for(TSAIND ii=0; ii<N; ii++)
    {
        blacklist[ii]=0;
    }
    isBlackListAlloc=1;
    return 1;
}

int TSAdataHandler::genSlidingWindowSubSeqs()
{
    float t1,t2;

    t1 = clock();

    int lenRawMotifData = procParams.motifLengths[procParams.indexMotifLenLongest];
    int lenRawMotifDataM1 = lenRawMotifData-1;

    subSeqPtr = (TSAsubSeq_t *)malloc(sizeof(TSAsubSeq_t)*lenTS);

    TSAIND ind=0;
    float ex=0;

    //############################## generating subsequences ##########################################
    while(ind<lenTS)
    {

        if (ind < lenTS - lenRawMotifDataM1)
        {
            subSeqPtr[ind].sTime  = samPtr[ind].tStamp;
            subSeqPtr[ind].eTime  = samPtr[ind+procParams.motifLengths[procParams.indexMotifLenReal]-1].tStamp;

            subSeqPtr[ind].pData = (TSADATA *)malloc(sizeof(TSADATA)*lenRawMotifData);
            subSeqPtr[ind].len = lenRawMotifData;
            subSeqPtr[ind].pTStamps = (float *)malloc(sizeof(float)*lenRawMotifData);
        }

        for(TSAIND ll = min(ind, lenTS - lenRawMotifDataM1-1) ; ll >= max(0,ind-lenRawMotifDataM1) ; ll--)
        {
            subSeqPtr[ll].pData[ind-ll] = samPtr[ind].value;
            subSeqPtr[ll].pTStamps[ind-ll] = samPtr[ind].tStamp;
        }

       ind++;
    }
    nSubSeqs = lenTS - lenRawMotifDataM1;

    t2 = clock();
    procLogPtr->tGenUniScal += (t2-t1)/CLOCKS_PER_SEC;

    return 1;

}

int TSAdataHandler::calculateDiffMotifLengths()
{
    int ii,jj;
    float t1,t2;

    t1 = clock();
    //######### computing all the motif lengths for several interpolation factors ##############
    int max_factor=-1;
    for (ii=0;ii<procParams.pattParams.nInterpFac; ii++)
    {
        procParams.motifLengths[ii] = (int)ceil((procParams.pattParams.durMotif*procParams.interpFac[ii])/pHop);
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
    if (procParams.pattParams.nInterpFac >1)
    {
        procParams.motifLengths[procParams.indexMotifLenLongest]+=1;
    }
    t2 = clock();
    procLogPtr->tProcTS += (t2-t1)/CLOCKS_PER_SEC;

    return 1;

}

int TSAdataHandler::convertHz2Cents(char *tonicFileName)
{
    float tonic, t1,t2;

    t1 = clock();

    //Opening tonic file
    FILE *fp =fopen(tonicFileName,"r");
    if (fp==NULL)
    {
        printf("Error opening file %s\n", tonicFileName);
        return 0;
    }
    int nRead = fscanf(fp, "%f\n",&tonic);
    fclose(fp);

    float temp1 = ((float)procParams.repParams.binsPOct)/LOG2;

    for(int ii=0;ii<lenTS;ii++)
    {
        samPtr[ii].value = (TSADATA)(temp1*log((samPtr[ii].value+EPS)/tonic));
    }

    t2 = clock();
    procLogPtr->tProcTS += (t2-t1)/CLOCKS_PER_SEC;

    return 1;
}


int TSAdataHandler::filterSamplesTS()
{
    float t1,t2;

    t1 = clock();

    // There can be many filtering criterion. Currently what is implemented is to filter out silence regions of pitch sequence.
    nSamSil = (int *)malloc(sizeof(int)*lenTS);
    memset(nSamSil, 0, sizeof(int)*lenTS);
    
    TSAIND nValidSamples = 0;
    TSAIND nSilSamples = 0;
    TSAIND *validSampleInd = (TSAIND *)malloc(sizeof(TSAIND)*lenTS);
    //store all locations which have valid samples
    for(TSAIND ii=0;ii<lenTS; ii++)
    {
        if(samPtr[ii].value > procParams.repParams.minPossiblePitch)
        {
            nSamSil[nValidSamples] = nSilSamples;
            validSampleInd[nValidSamples]=ii;
            nValidSamples++;
            nSilSamples=0;
        }
        else
        {
         nSilSamples++;   
        }
    }
    // just copy valid samples to another location
    TSAIND lenTS_new = nValidSamples;
    TSAsam_t * samPtr_new = (TSAsam_t *)malloc(sizeof(TSAsam_t)*lenTS_new);
    for (TSAIND ii=0;ii<lenTS_new; ii++)
    {
        samPtr_new[ii] = samPtr[validSampleInd[ii]];
    }

    //free old memory
    free(samPtr);
    free(validSampleInd);
    samPtr = samPtr_new;
    lenTS = lenTS_new;

    t2 = clock();
    procLogPtr->tProcTS += (t2-t1)/CLOCKS_PER_SEC;

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
    nRead = fscanf(fp, "%f\t%f\n",&tsTime1, &tsData);
    nRead = fscanf(fp, "%f\t%f\n",&tsTime2, &tsData);

    pHop = (tsTime2-tsTime1);

    procParams.repParams.dsFactor = (int)floor((procParams.repParams.sampleRate/(pHop*1000))+0.5);

    fclose(fp);
    return 1;
}

int TSAdataHandler::quantizeSampleSubSeqs(int quantizationType)
{
    float t1,t2;
    t1 = clock();
    float temp = floor(procParams.repParams.binsPOct/procParams.repParams.quantSize);

    if (quantizationType == NO_QUANT)
    {
        return 1;
    }

    if (quantizationType == NEAREST_NOTE_QUANT)
    {
        for (int ii=0;ii<nSubSeqs; ii++)
        {
            for(int jj=0;jj<procParams.pattParams.subSeqLen;jj++)
            {
                subSeqPtr[ii].pData[jj] = quantizePitch(subSeqPtr[ii].pData[jj], temp);

            }
        }
    }

    t2 = clock();
    procLogPtr->tProcTS += (t2-t1)/CLOCKS_PER_SEC;

    return 1;
}


int TSAdataHandler::quantizeSampleTS(int quantizationType)
{
    float t1,t2;
    t1 = clock();
    float temp = floor(procParams.repParams.binsPOct/procParams.repParams.quantSize);

    if (quantizationType == NO_QUANT)
    {
        return 1;
    }

    if (quantizationType == NEAREST_NOTE_QUANT)
    {
        for (int ii=0;ii<lenTS; ii++)
        {
            samPtr[ii].value = quantizePitch(samPtr[ii].value, temp);
        }
    }

    t2 = clock();
    procLogPtr->tProcTS += (t2-t1)/CLOCKS_PER_SEC;

    return 1;
}

int TSAdataHandler::normalizeSubSeqs(int normType)
{
    float t1,t2;
    t1 = clock();

    if (normType==NO_NORM)
    {
        return 1;
    }
    else if(normType==Z_NORM)
    {
        float mu, std;
        for(int ii=0; ii<nSubSeqs; ii++)
        {
            mu = subSeqPtr[ii].mean;
            std = subSeqPtr[ii].std;

            for (int jj=0;jj<subSeqPtr[ii].len;jj++)
            {
                subSeqPtr[ii].pData[jj] = (subSeqPtr[ii].pData[jj]-mu)/std;
            }
             subSeqPtr[ii].mean = 0;
             subSeqPtr[ii].std = 1;
        }
    }
    else if(normType==MEAN_SUB_NORM)
    {
        float mu;
        for(int ii=0; ii<nSubSeqs; ii++)
        {
            mu = subSeqPtr[ii].mean;

            for (int jj=0;jj<subSeqPtr[ii].len;jj++)
            {
                subSeqPtr[ii].pData[jj] = (subSeqPtr[ii].pData[jj]-mu);
            }
            subSeqPtr[ii].mean=0;

        }
    }
    else if(normType==MEDIAN_SUB_NORM)
    {
        float median;
        TSADATA *tempStr = (TSADATA *)malloc(sizeof(TSADATA)*subSeqPtr[0].len*2);

        for(int ii=0; ii<nSubSeqs; ii++)
        {
            memcpy(tempStr, subSeqPtr[ii].pData, sizeof(TSADATA)*subSeqPtr[ii].len);
            median = computeMedian(tempStr, subSeqPtr[ii].len);

            for (int jj=0;jj<subSeqPtr[ii].len;jj++)
            {
                subSeqPtr[ii].pData[jj] = (subSeqPtr[ii].pData[jj]-median);
            }
            subSeqPtr[ii].mean= subSeqPtr[ii].mean - median;

        }
    }
    else if(normType==QMEDIAN_SUB_NORM)
    {
        float median;
        TSADATA *tempStr = (TSADATA *)malloc(sizeof(TSADATA)*subSeqPtr[0].len*2);
        float binsPSemiTone = (procParams.repParams.binsPOct/12);

        for(int ii=0; ii<nSubSeqs; ii++)
        {
            memcpy(tempStr, subSeqPtr[ii].pData, sizeof(TSADATA)*subSeqPtr[ii].len);
            median = computeMedian(tempStr, subSeqPtr[ii].len);
            median = floor((median/binsPSemiTone)+0.5)*binsPSemiTone;

            for (int jj=0;jj<subSeqPtr[ii].len;jj++)
            {
                subSeqPtr[ii].pData[jj] = (subSeqPtr[ii].pData[jj]-median);
            }
            subSeqPtr[ii].mean= subSeqPtr[ii].mean - median;

        }
    }

    else if(normType==MAD_NORM)
    {
        float median, mad;
        TSADATA *tempStr = (TSADATA *)malloc(sizeof(TSADATA)*subSeqPtr[0].len*2);


        for(int ii=0; ii<nSubSeqs; ii++)
        {
            memcpy(tempStr, subSeqPtr[ii].pData, sizeof(TSADATA)*subSeqPtr[ii].len);
            median = computeMedian(tempStr, subSeqPtr[ii].len);
            mad = computeMAD(tempStr, subSeqPtr[ii].len, median);

            for (int jj=0;jj<subSeqPtr[ii].len;jj++)
            {
                subSeqPtr[ii].pData[jj] = (subSeqPtr[ii].pData[jj]-median)/mad;
            }
            subSeqPtr[ii].mean = (subSeqPtr[ii].mean - median)/mad;
            subSeqPtr[ii].std = subSeqPtr[ii].std/mad;
        }
    }

    t2 = clock();
    procLogPtr->tProcTS += (t2-t1)/CLOCKS_PER_SEC;

    return 1;
}


int TSAdataHandler::downSampleSubSeqs()
{
    //determining downsampling factor
    procParams.repParams.dsFactor = (int)floor((procParams.repParams.sampleRate/(procParams.repParams.pitchHop*1000))+0.5);

    int dsFac =  procParams.repParams.dsFactor;

    //this changes length of the subsequences in samples
    procParams.pattParams.subSeqLen = (int)floor(procParams.pattParams.subSeqLen/procParams.repParams.dsFactor);

    //this changes affective sampling rate of the pitch sequence stored in the subseqs
    procParams.repParams.pitchHop = procParams.repParams.pitchHop*procParams.repParams.dsFactor;

    pHop =  procParams.repParams.pitchHop;


    TSAsubSeq_t *subSeqPtr_new = (TSAsubSeq_t *)malloc(sizeof(TSAsubSeq_t)*nSubSeqs);

    for(TSAIND ii=0;ii<nSubSeqs; ii++)
    {
        subSeqPtr_new[ii].pData = (TSADATA *)malloc(sizeof(TSADATA)*procParams.pattParams.subSeqLen);
        subSeqPtr_new[ii].pTStamps = (float *)malloc(sizeof(float)*procParams.pattParams.subSeqLen);
        for(int jj=0;jj<procParams.pattParams.subSeqLen;jj++)
        {
            subSeqPtr_new[ii].pData[jj] = subSeqPtr[ii].pData[jj*dsFac];
            subSeqPtr_new[ii].pTStamps[jj] = subSeqPtr[ii].pTStamps[jj*dsFac];
        }
        subSeqPtr_new[ii].len = (int)floor(subSeqPtr[ii].len/procParams.repParams.dsFactor);
        subSeqPtr_new[ii].id = subSeqPtr[ii].id;
        subSeqPtr_new[ii].sTime = subSeqPtr[ii].sTime;
        subSeqPtr_new[ii].eTime = subSeqPtr[ii].eTime;
        free(subSeqPtr[ii].pData);
        free(subSeqPtr[ii].pTStamps);
    }
    free(subSeqPtr);

    subSeqPtr=subSeqPtr_new;

    return 1;

}

// int TSAdataHandler::storeTimeElapseInAdjSamples(){
//     
//     timeElapse = (float *)malloc(sizeof(float)*lenTS);
//     memset(timeElapse, 0, lenTS*sizeof(float));
//     
//     for (int ii=1; ii<lenTS;ii++){
//         //NOTE: we start from 1 becasue time elapse is stored in the upcoming sample. hence starting from 1
//         timeElapse[ii] = samPtr[ii].tStamp-samPtr[ii-1].tStamp;
//         //printf("Time elapse: %f\n", timeElapse[ii]);
//     }
//     return 1;
// }


int TSAdataHandler::downSampleTS()
{
    float t1,t2;
    t1 = clock();
    if (procParams.repParams.dsFactor <=0)
    {
        return 1;
    }

    TSAIND lenTS_new = ceil(lenTS/procParams.repParams.dsFactor);
    TSAsam_t * samPtr_new = (TSAsam_t *)malloc(sizeof(TSAsam_t)*lenTS_new);
    for (int ii=0;ii<lenTS_new; ii++)
    {
        samPtr_new[ii] = samPtr[ii*procParams.repParams.dsFactor];
    }

    //correcting hop size as well
    pHop = pHop*procParams.repParams.dsFactor;

    //free old memory
    free(samPtr);
    samPtr = samPtr_new;
    lenTS = lenTS_new;

    t2 = clock();
    procLogPtr->tProcTS += (t2-t1)/CLOCKS_PER_SEC;

    return 1;
}

int TSAdataHandler::readSubSeqInfo(char *fileName)
{
    FILE *fp;
    fp =fopen(fileName,"r");
    int tempi[4];
    float tempf[2];
    //TSAIND tempi;

    if (fp==NULL)
    {
        printf("Error opening file %s\n", fileName);
        return 0;
    }
    TSAIND cnt = 0;
    while (fscanf(fp, "%d\t%f\t%f\n",&tempi[0], &tempf[0], &tempf[1])!=EOF)    //read till the end of the file
    {
        subSeqPtr[cnt].id = tempi[0];
        subSeqPtr[cnt].sTime = tempf[0];
        subSeqPtr[cnt].eTime = tempf[1];
        subSeqPtr[cnt].len = (int)floor((tempf[1]-tempf[0])/procParams.repParams.pitchHop);
        cnt++;
    }
    fclose(fp);

    return 1;

}

int TSAdataHandler::setSubSeqLengthsTStamps()
{
    for(TSAIND ii=0;ii<nSubSeqs; ii++)
    {
        for(TSAIND jj=0;jj<procParams.pattParams.subSeqLen; jj++)
        {
            subSeqPtr[ii].pTStamps[jj] = subSeqPtr[ii].sTime + jj*procParams.repParams.pitchHop;
        }
    }

}

int TSAdataHandler::setSubSeqLengthsFIX(int motifLen)
{
    for(TSAIND ii=0;ii<nSubSeqs; ii++)
    {
        subSeqPtr[ii].len = motifLen;
    }
}

int TSAdataHandler::readSubSeqLengths(char *fileName)
{
    FILE *fp;
    fp =fopen(fileName,"r");
    int tempi[4];
    float tempf[2];

    if (fp==NULL)
    {
        printf("Error opening file %s\n", fileName);
        return 0;
    }
    TSAIND cnt = 0;
    while (fscanf(fp, "%f\t%f\t%d\t%d\t%d\n",&tempf[0], &tempf[1], &tempi[0], &tempi[1], &tempi[2])!=EOF)    //read till the end of the file
    {
        subSeqPtr[cnt].len = (int)floor(tempf[1]/procParams.repParams.pitchHop);
        cnt++;
    }
    fclose(fp);

    return 1;

}

int TSAdataHandler::readSubSeqData(char *fileName, TSAIND nSubs)
{
    FILE *fp;
    fp =fopen(fileName,"r");
    int nreads=0;

    subSeqPtr = (TSAsubSeq_t *)malloc(sizeof(TSAsubSeq_t)*nSubs);

    if (fp==NULL)
    {
        printf("Error opening file %s\n", fileName);
        return 0;
    }
    for (int ii=0;ii< nSubs; ii++)    //read till the end of the file
    {
        subSeqPtr[ii].pData = (TSADATA *)malloc(sizeof(TSADATA)*procParams.pattParams.subSeqLen);
        subSeqPtr[ii].pTStamps = (float *)malloc(sizeof(float)*procParams.pattParams.subSeqLen);//dummy allocation
        nreads = fread(subSeqPtr[ii].pData, sizeof(TSADATA), procParams.pattParams.subSeqLen, fp);
    }

    fclose(fp);
    nSubSeqs = nSubs;
    return 1;
}

int TSAdataHandler::readTSData(char *fileName)
{
    FILE *fp;
    float tsData, t1,t2;
    float tsTime;
    TSAIND cnt;

    t1 = clock();
    //count number of lines in the file specified
    nLinesFile = getNumLines(fHandle.getTSFileName());
    procLogPtr->lenTS += nLinesFile;
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
    while (fscanf(fp, "%f\t%f\n",&tsTime, &tsData)!=EOF)    //read till the end of the file
    {
        samPtr[cnt].value  = tsData;
        samPtr[cnt].tStamp  = tsTime;
        cnt++;
    }
    fclose(fp);
    t2 = clock();
    procLogPtr->tLoadTS += (t2-t1)/CLOCKS_PER_SEC;

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

int TSAdataHandler::dumpDiscMotifInfo(char *motifFile, TSAmotifInfo_t *priorityQDisc, int K, int verbos)
{
    FILE *fp;
    int ii;
    fp =fopen(motifFile,"w");
    for(ii=0;ii<K;ii++)
    {
        fprintf(fp, "%f\t%f\t%f\t%f\t%lf\t%lld\t%lld\n", subSeqPtr[priorityQDisc[ii].ind1].sTime, subSeqPtr[priorityQDisc[ii].ind1].eTime, subSeqPtr[priorityQDisc[ii].ind2].sTime, subSeqPtr[priorityQDisc[ii].ind2].eTime, priorityQDisc[ii].dist, priorityQDisc[ii].ind1, priorityQDisc[ii].ind2);
        if (verbos)
        {
            printf("motif pair is %f\t%f\t%f\t%lld\t%lld\n", subSeqPtr[priorityQDisc[ii].ind1].sTime,subSeqPtr[priorityQDisc[ii].ind2].sTime, priorityQDisc[ii].dist, priorityQDisc[ii].ind1%3, priorityQDisc[ii].ind2%3);
        }
    }
    fclose(fp);

    return 1;
}

int TSAdataHandler::dumpSearMotifInfo(char *motifFile, TSAmotifInfoExt_t **priorityQSear, TSAIND nQueries, int K, int verbos)
{
    FILE *fp;
    fp = fopen(motifFile, "w");

    for(TSAIND ii=0;ii<K;ii++)
    {
        for(TSAIND jj=0;jj < nQueries;jj++)
        {
            if ( priorityQSear[jj][ii].dist < INF)
            {
                fprintf(fp, "%f\t%f\t%f\t%f\t%f\t%d\t", subSeqPtr[priorityQSear[jj][ii].ind1].sTime, subSeqPtr[priorityQSear[jj][ii].ind1].eTime, priorityQSear[jj][ii].sTime, priorityQSear[jj][ii].eTime, priorityQSear[jj][ii].dist, priorityQSear[jj][ii].searchFileID);
            }
            else
            {
                fprintf(fp, "%f\t%f\t%f\t%f\t%f\t%d\t", -1.0,-1.0,-1.0,-1.0,-1.0,-1);
            }

        }
        fprintf(fp, "\n");
    }
    fclose(fp);

    return 1;
}

int TSAdataHandler::dumpPatternKNNInfo(char *motifFile, TSAmotifInfoExt_t **priorityQSear, TSAIND nQueries, int KNN, int verbos)
{
    FILE *fp;
    fp = fopen(motifFile, "w");

    for(TSAIND ii=0;ii<nQueries;ii++)
    {
        for(TSAIND jj=0;jj < KNN;jj++)
        {
            if ( priorityQSear[ii][jj].dist < INF)
            {
                fprintf(fp, "%lld\t%lld\t%f\n", subSeqPtr[priorityQSear[ii][jj].ind1].id, priorityQSear[ii][jj].patternID2, priorityQSear[ii][jj].dist);
            }
            else
            {
                fprintf(fp, "%lld\t%lld\t%f\n", -1,-1,INF);
            }

        }

    }
    fclose(fp);

    return 1;
}

int TSAdataHandler::dumpPatternDISTInfo(char *outputFile, TSAmotifInfoExt_t **priorityQSear, TSAIND nPatterns, TSAIND *pattPerQ, int verbos)
{
    FILE *fp;
    fp = fopen(outputFile, "w");

    for(TSAIND ii=0;ii<nPatterns;ii++)
    {
        for(TSAIND jj=0;jj < pattPerQ[ii];jj++)
        {
            if ( priorityQSear[ii][jj].dist < INF)
            {
                fprintf(fp, "%lld\t%lld\t%f\n", subSeqPtr[priorityQSear[ii][jj].ind1].id, priorityQSear[ii][jj].patternID2, priorityQSear[ii][jj].dist);
            }
            else
            {
                fprintf(fp, "%lld\t%lld\t%f\n", -1,-1,INF);
            }

        }

    }
    fclose(fp);

    return 1;
}


fileNameHandler::fileNameHandler()
{
    nSearchFiles=-1;
}
fileNameHandler::~fileNameHandler()
{
    if(nSearchFiles>0)
    {
        for(int ii=0; ii< nSearchFiles; ii++)
        {
            free(searchFileNames[ii]);
        }
        free(searchFileNames);
    }
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

char* fileNameHandler::getFlatNoteFileName()
{
    memset(fileName, '\0', sizeof(char)*MAX_FNAME_CHARS);
    strcat(fileName, baseName);
    strcat(fileName, fileExtPtr->FNFileExt);

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
char* fileNameHandler::getOutFileNamePostRR(int similarityMeasure)
{
    memset(fileName, '\0', sizeof(char)*MAX_FNAME_CHARS);
    strcat(fileName, baseName);
    strcat(fileName, fileExtPtr->outFileExt);
    strcat(fileName, SimMeasureNames[similarityMeasure]);

    return fileName;

}

char* fileNameHandler::getMappFileName()
{
    memset(fileName, '\0', sizeof(char)*MAX_FNAME_CHARS);
    strcat(fileName, baseName);
    strcat(fileName, fileExtPtr->mappFileExt);

    return fileName;

}

char* fileNameHandler::getSearchListFileName()
{
    memset(fileName, '\0', sizeof(char)*MAX_FNAME_CHARS);
    strcat(fileName, baseName);
    strcat(fileName, fileExtPtr->searchListExt);

    return fileName;
}
char* fileNameHandler::getQueryFileName()
{
    memset(fileName, '\0', sizeof(char)*MAX_FNAME_CHARS);
    strcat(fileName, baseName);
    strcat(fileName, fileExtPtr->queryFileExt);

    return fileName;
}

char* fileNameHandler::getSubSeqFileName()
{
    memset(fileName, '\0', sizeof(char)*MAX_FNAME_CHARS);
    strcat(fileName, baseName);
    strcat(fileName, fileExtPtr->subSeqFileExt);

    return fileName;
}

char* fileNameHandler::getSubSeqTNFileName()
{
    memset(fileName, '\0', sizeof(char)*MAX_FNAME_CHARS);
    strcat(fileName, baseName);
    strcat(fileName, fileExtPtr->subSeqTNFileExt);

    return fileName;
}

char* fileNameHandler::getSubSeqInfoFileName()
{
    memset(fileName, '\0', sizeof(char)*MAX_FNAME_CHARS);
    strcat(fileName, baseName);
    strcat(fileName, fileExtPtr->subSeqInfoFileExt);

    return fileName;
}

char* fileNameHandler::getPatternKNNFileName()
{
    memset(fileName, '\0', sizeof(char)*MAX_FNAME_CHARS);
    strcat(fileName, baseName);
    strcat(fileName, fileExtPtr->patternKNNExt);

    return fileName;
}




int fileNameHandler::loadSearchFileList()
{
    char tempFilename[MAX_FNAME_CHARS];

    FILE *fp1 = fopen(getSearchListFileName(), "r");
    if (fp1==NULL)
    {
        printf("Error opening file %s\n", getSearchListFileName());
        return 0;
    }
    int ii=0;
    searchFileNames = (char **)malloc(sizeof(char *)*MAX_NUM_SEARCHFILES);
    while(fgets(tempFilename, MAX_FNAME_CHARS, fp1))
    {
        searchFileNames[ii] = (char *)malloc(sizeof(char)*MAX_FNAME_CHARS);
        sscanf(tempFilename, "%[^\n]s\n", searchFileNames[ii]);
        ii++;
    }
    nSearchFiles  = ii;
    fclose(fp1);
    return 1;
}

