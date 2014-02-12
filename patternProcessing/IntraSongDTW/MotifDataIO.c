/******************************************************************************
*******************************************************************************
****** Sankalp Gulati                                                   *******
****** MUSIC TECHNOLOGY GROUP, UPF, BARCELONA                           *******
******                                                                  *******
*******************************************************************************
******************************************************************************/

#include "MotifDataIO.h"


/*
 * This function quickly returns number of lines in a text file
 */
INDTYPE getNumLines(const char *file)
{
    int fp;
    struct stat fs;
    char *buf;
    INDTYPE line=0, ii;

    fp = open(file,O_RDONLY);
    if (fp == -1) 
    {
                printf("Error opening file %s\n",file);
                return 0;
    }
    if (fstat(fp, &fs) == -1)
    {
                printf("Error opening file %s\n",file);
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


INDTYPE generateSubsequenceDB(DATATYPE ***d, segInfo_t **t, int *motifLen, char *baseName, fileExts_t *myFileExts, procParams_t *myProcParams, procLogs_t *myProcLogs, int verbos)
{
     FILE *fp;
    char pitchFile[200]={'\0'}, tonicFile[200]={'\0'}, segmentFile[200]={'\0'};
    
    INDTYPE numLinesInFile, ind, ii, jj, ll, lenTS, N, blacklistCNT;
    DATATYPE *pitchSamples, **data, **dataInterp;
    
    float tonic, pHop, temp1, temp[4]={0}, t1, t2, ex, pitchTemp, timeTemp;
    float *timeSamples, *mean, *stdVec, *indLow, *indHigh;
    
    int lenMotifReal,lenMotifRealM1,lenMotifInterpHM1, lenMotifInterpLM1, lenMotifInterpH, lenMotifInterpL,dsFactor, nRead, varSam; 
    int *blacklist;
    
    segInfoInterp_t *tStamps;
    segInfo_t *taniSegs, *tStampsInterp;
    
    
    //pitch file name
    strcat(pitchFile,baseName);
    strcat(pitchFile,myFileExts->pitchExt);
    //tonic file name
    strcat(tonicFile,baseName);
    strcat(tonicFile,myFileExts->tonicExt);
    //segment file name
    strcat(segmentFile,baseName);
    strcat(segmentFile,myFileExts->segExt);
    
    //########################## READING PITCH DATA ##########################
    // Reading number of lines in the pitch file
    numLinesInFile = getNumLines(pitchFile);
    myProcLogs->totalPitchSamples += numLinesInFile;
    // after downsampling we will be left with these many points
    numLinesInFile = floor(numLinesInFile/myProcParams->dsFactor) +1;
    
    //allocating memory for pitch and time samples
    pitchSamples = (DATATYPE*)malloc(sizeof(DATATYPE)*numLinesInFile);     // since we don't know silence regions, allocate maximum possible number of samples
    timeSamples = (float*)malloc(sizeof(float)*numLinesInFile);
    
    blacklist = (int*)malloc(sizeof(int)*numLinesInFile);  //subsequences which we don't have to consider. Unfortunately we know them after we already stored them
    memset(blacklist,0,sizeof(int)*numLinesInFile);
    
    data = (DATATYPE **)malloc(sizeof(DATATYPE *)*numLinesInFile);         //since we don't know valid subsequences, we allocate max possible subsequences and later discard them and free the memory
    tStamps = (segInfoInterp_t *)malloc(sizeof(segInfoInterp_t)*numLinesInFile);                
    
    mean = (float *)malloc(sizeof(float)*numLinesInFile);
    memset(mean,0,sizeof(float)*numLinesInFile);
    stdVec = (float *)malloc(sizeof(float)*numLinesInFile);
    memset(stdVec,0,sizeof(float)*numLinesInFile);
    
    //opening pitch file JUST FOR OBTAINING HOP SIZE OF THE PITCH SEQUENCE
    fp =fopen(pitchFile,"r");
    if (fp==NULL)
    {
        printf("Error opening file %s\n", pitchFile);
        return 0;
    }
    //reading just first two lines, in order to obtain hopsize//
    nRead = fscanf(fp, "%f\t%f\n",&temp[0],&temp[1]);
    nRead = fscanf(fp, "%f\t%f\n",&temp[2],&temp[3]);
    fclose(fp);
    pHop = (temp[2]-temp[0])*myProcParams->dsFactor;  //final hop size afte downsampling
    lenMotifReal = (int)round(myProcParams->durMotif/pHop);
    *motifLen = lenMotifReal;
    lenMotifInterpH = (int)ceil((myProcParams->durMotif*myProcParams->factorHigh)/pHop)+1;  //adding one because we need one extra sample for cubic interpolation
    lenMotifInterpL = (int)round((myProcParams->durMotif*myProcParams->factorLow)/pHop);
    varSam = (int)round(myProcParams->varDur/pHop);
    temp1 = ((float)myProcParams->binsPOct)/LOG2;
    lenMotifRealM1 = lenMotifReal-1;
    lenMotifInterpHM1 = lenMotifInterpH-1;
    lenMotifInterpLM1 = lenMotifInterpL-1;
    
    //Opening the tonic file
    fp =fopen(tonicFile,"r");
    if (fp==NULL)
    {
        printf("Error opening file %s\n", tonicFile);
        return 0;
    }
    nRead = fscanf(fp, "%f\n",&tonic);
    fclose(fp);
    
    
    //Finally opening the pitch file to read the data
    fp =fopen(pitchFile,"r");
    t1 = clock();
    ind=0;
    jj=0;
    while (fscanf(fp, "%f\t%f\n",&timeTemp,&pitchTemp)!=EOF)    //read till the end of the file
    {
        
        if (pitchTemp>myProcParams->minPossiblePitch) //only fill in meaningful pitch data, reject other pitch samples
        {
            pitchSamples[ind]= (DATATYPE)round(temp1*log((pitchTemp+EPS)/tonic));
            timeSamples[ind] = timeTemp;
            ind++;
            jj++;
        }
        
        for(ii=0;ii<myProcParams->dsFactor-1;ii++)
        {
            // just bypass other samples to downsample the pitch sequence
            nRead = fscanf(fp, "%f\t%f\n",&timeTemp,&pitchTemp);
            jj++;
        }
        
              
    }
    fclose(fp);
    t2 = clock();
    myProcLogs->timeDataLoad += (t2-t1)/CLOCKS_PER_SEC;
    myProcLogs->totalPitchNonSilSamples += ind;
    if (verbos)
    {printf("Time taken to load the pitch data :%f\n",(t2-t1)/CLOCKS_PER_SEC);}
    
    //########################## Subsequence generation + selection step ##########################
    // In subsequence selection our aim is to discard those subsequences which result into trivial matches, 
    // which are flat regions. For removing flat regions the obvious choice is to consider variance of a 
    // subsequence and filter using a myProcParams->threshold. The problem here is large amounts of octave errors because
    // of which the variance increases but affectively large chunk of data is still flat. 
    
    // For solving problem mentioned above we resort to short duration variance for deciding wheather a 
    // given sample belongs to a flat region or not. Later we accumulate total number of samples in a 
    // subsequence which belong to a flat region and filter the subsequence based on a myProcParams->threshold.
    t1 = clock();
    //computing local mean and variance
    for(ii=0;ii<2*varSam+1;ii++)
    {
        mean[varSam] +=  pitchSamples[ii] ;
    }    
    for(ii=varSam+1;ii<ind-varSam;ii++)
    {
        mean[ii] =  mean[ii-1] - pitchSamples[ii-varSam-1] + pitchSamples[ii+varSam] ;  //running mean
    }
    N = 2*varSam+1 ; 
    
    //computing running variance
    for(ii=0;ii<2*varSam+1;ii++)
    {
        temp[0] = (pitchSamples[ii]- mean[ii]/N);
        stdVec[varSam] += temp[0]*temp[0];
    } 
    for(ii=varSam+1;ii<ind-varSam;ii++)
    {
        temp[0] = (pitchSamples[ii-varSam-1] -mean[ii-varSam-1]/N);
        stdVec[ii] =  stdVec[ii-1] - temp[0]*temp[0] ; 
        temp[1] = (pitchSamples[ii+varSam] -mean[ii+varSam]/N);
        stdVec[ii] += temp[1]*temp[1] ;
    }
    // after computing running variance, applying a myProcParams->threshold and selecting the onces which are corresponsing to non flat regions
    for(ii=varSam;ii<ind-varSam;ii++)
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
    
    lenTS  = ind;       //number of non trivial pitch samples, they might still carry number of subsequences which have to be discarded
    //blackCnt=lenTS;
    ind=0;
    ex=0;
    
    //############################## generating subsequences ##########################################
    while(ind<lenTS)
    {
        tStamps[ind].str = timeSamples[ind];
        
        if (ind<lenTS-lenMotifInterpHM1)
        {
            data[ind] = (DATATYPE *)malloc(sizeof(DATATYPE*)*lenMotifInterpH);
            tStamps[ind].end = timeSamples[ind+lenMotifRealM1];
            tStamps[ind].endInterpH = timeSamples[ind+lenMotifInterpHM1];
            tStamps[ind].endInterpL = timeSamples[ind+lenMotifInterpLM1];
            if (fabs(tStamps[ind].str - tStamps[ind].endInterpH) > myProcParams->durMotif*myProcParams->factorHigh + myProcParams->maxPauseDur)       //allow 200 ms pauses in total not more than that
                {
                    blacklist[ind]=1;
                }
        }
        
        for(ll = min(ind,lenTS-lenMotifInterpHM1-1) ; ll >= max(0,ind-lenMotifInterpHM1) ; ll--)
        {
            data[ll][ind-ll] = pitchSamples[ind]; 
        }
        
        ex+=stdVec[ind];
        if (ind >= lenMotifInterpHM1)
        {
            if (blacklist[ind-lenMotifInterpHM1]==0)
            {
                if(ex < myProcParams->flatThreshold*(float)lenMotifInterpH)
                {
                    
                    blacklist[ind-lenMotifInterpHM1]=1;
                }
                
            }
            
            ex -= stdVec[ind-lenMotifInterpHM1];
        }
       ind++;        
    }
    lenTS = lenTS-lenMotifInterpHM1;   //we only have lenMotifInterpHM1 samples less than original number lenTS. it still counts blacklisted candidates
    myProcLogs->totalSubsGenerated += lenTS;
    t2 = clock();
    printf("Time taken to create subsequences:%f\n",(t2-t1)/CLOCKS_PER_SEC);
    
    //free memory which is not needed further
    free(mean);
    free(stdVec);
    free(pitchSamples);
    free(timeSamples);
    
    
    // Finding all the subsequences that should be blacklisted because they fall in TANI regions
    t1 = clock();
    numLinesInFile = getNumLines(segmentFile);
    taniSegs = (segInfo_t *)malloc(sizeof(segInfo_t)*numLinesInFile);
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
        taniSegs[ii].str = temp[0];
        taniSegs[ii].end = temp[1];
        ii++;
    }
    fclose(fp);
    //blacklisting the indexes whose time stamps lie between these segments
    for(ii=0;ii<lenTS;ii++)
    {
        for(jj=0;jj<numLinesInFile;jj++)
        {
            if ((tStamps[ii].end >= taniSegs[jj].str)&&(tStamps[ii].str <= taniSegs[jj].end))
                blacklist[ii]=1;
        }
        
    }
    t2 = clock();
    if (verbos)
    {printf("Time taken to blacklist tani sections :%f\n",(t2-t1)/CLOCKS_PER_SEC);}
    
    //%%%%%%%%%%%%%%%% Removing blacklisted subsequences %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    t1 = clock();
    //finding number of non blacklisted subsequences
    N=0;
    for (ii=0;ii<lenTS; ii++)
    {
        if (blacklist[ii]==0)
            N++;
    }
    blacklistCNT = lenTS-N;
    myProcLogs->totalSubsBlacklisted += blacklistCNT;
    dataInterp = (DATATYPE **)malloc(sizeof(DATATYPE *)*N*3);   //allocating memory also to have interpolated subsequences
    tStampsInterp = (segInfo_t *)malloc(sizeof(segInfo_t)*(N)*3);         //allocating memory also to have interpolated subsequences
    jj=0;
    
    //we do interpolation as well
    indLow = (float *)malloc(sizeof(float)*lenMotifReal);
    indHigh = (float *)malloc(sizeof(float)*lenMotifReal);
    for (ii=0;ii<lenMotifReal; ii++)
    {
        indLow[ii] = myProcParams->factorLow*ii;
        indHigh[ii] = myProcParams->factorHigh*ii;
    }
    
    for (ii=0;ii<lenTS;ii++)
    {
        if (blacklist[ii]==0)
        {
            
            //low stretched subsequence
            dataInterp[jj] = (DATATYPE *)malloc(sizeof(DATATYPE)*lenMotifReal);
            cubicInterpolate(data[ii], dataInterp[jj], indLow, lenMotifReal);
            tStampsInterp[jj].str = tStamps[ii].str;
            tStampsInterp[jj].end = tStamps[ii].endInterpL;
            jj++;
            
            //normal
            dataInterp[jj] = data[ii];
            tStampsInterp[jj].str = tStamps[ii].str;
            tStampsInterp[jj].end = tStamps[ii].end;
            jj++;
            
            //compacted subsequence
            dataInterp[jj] = (DATATYPE *)malloc(sizeof(DATATYPE)*lenMotifReal);
            cubicInterpolate(data[ii], dataInterp[jj], indHigh, lenMotifReal);
            tStampsInterp[jj].str = tStamps[ii].str;
            tStampsInterp[jj].end = tStamps[ii].endInterpH;
            jj++;
            
            
        }
        else
        {
            free(data[ii]);
        }
    }
    lenTS = N*3; 
    free(data);
    free(tStamps);
    free(blacklist);
    
    t2 = clock();
    myProcLogs->timeRemBlacklist += (t2-t1)/CLOCKS_PER_SEC;
    myProcLogs->totalSubsInterpolated += lenTS;
    if (verbos)
    {printf("Time taken to remove blacklist subsequences :%f\n",(t2-t1)/CLOCKS_PER_SEC);}
    
    if( verbos == 1 )
        printf("Finally number of subsequences are: %lld\nNumber of subsequences removed are: %lld\n",lenTS,blacklistCNT);
        printf("Length of Each Time Series : %d\n\n",lenMotifReal);

    *d = dataInterp;
    *t = tStampsInterp;
    return lenTS;
}

INDTYPE loadSeedMotifSequence(DATATYPE ***d, segInfo_t **t, int *motifLen, char *baseName, fileExts_t *myFileExts, procParams_t *myProcParams, int maxNMotifsPairs, int verbos)
{
     FILE *fp;
    char pitchFile[200]={'\0'}, tonicFile[200]={'\0'}, motifFile[200]={'\0'};
    
    INDTYPE numLinesInFile, ind, ii, jj, ll, lenTS, N, totalPitchNonSilSamples, *seedMotifInd;
    DATATYPE *pitchSamples, **data, **dataInterp;
    
    float tonic, pHop, temp1, temp[10]={0}, ex, pitchTemp, timeTemp;
    float *timeSamples, *indLow, *indHigh, min_val;
    
    int lenMotifReal,lenMotifRealM1,lenMotifInterpHM1, lenMotifInterpLM1, lenMotifInterpH, lenMotifInterpL,dsFactor, nRead, NMotifs; 
    
    segInfoInterp_t *tStamps;
    segInfo_t *taniSegs, *tStampsInterp, *seedMotifs;
    
    
    //pitch file name
    strcat(pitchFile,baseName);
    strcat(pitchFile,myFileExts->pitchExt);
    //tonic file name
    strcat(tonicFile,baseName);
    strcat(tonicFile,myFileExts->tonicExt);
    //motif file name
    strcat(motifFile,baseName);
    strcat(motifFile,myFileExts->seedMotifExt);
    
    //########################## READING PITCH DATA ##########################
    // Reading number of lines in the pitch file
    numLinesInFile = getNumLines(pitchFile);
    // after downsampling we will be left with these many points
    numLinesInFile = floor(numLinesInFile/myProcParams->dsFactor) +1;
    
    //allocating memory for pitch and time samples
    pitchSamples = (DATATYPE*)malloc(sizeof(DATATYPE)*numLinesInFile);     // since we don't know silence regions, allocate maximum possible number of samples
    timeSamples = (float*)malloc(sizeof(float)*numLinesInFile);
    
    //opening pitch file JUST FOR OBTAINING HOP SIZE OF THE PITCH SEQUENCE
    fp =fopen(pitchFile,"r");
    if (fp==NULL)
    {
        printf("Error opening file %s\n", pitchFile);
        return 0;
    }
    //reading just first two lines, in order to obtain hopsize//
    nRead = fscanf(fp, "%f\t%f\n",&temp[0],&temp[1]);
    nRead = fscanf(fp, "%f\t%f\n",&temp[2],&temp[3]);
    fclose(fp);
    pHop = (temp[2]-temp[0])*myProcParams->dsFactor;  //final hop size afte downsampling
    lenMotifReal = (int)round(myProcParams->durMotif/pHop);
    *motifLen = lenMotifReal;
    lenMotifInterpH = (int)ceil((myProcParams->durMotif*myProcParams->factorHigh)/pHop)+1;  //adding one because we need one extra sample for cubic interpolation
    lenMotifInterpL = (int)round((myProcParams->durMotif*myProcParams->factorLow)/pHop);
    temp1 = ((float)myProcParams->binsPOct)/LOG2;
    lenMotifRealM1 = lenMotifReal-1;
    lenMotifInterpHM1 = lenMotifInterpH-1;
    lenMotifInterpLM1 = lenMotifInterpL-1;
    
    //Opening the tonic file
    fp =fopen(tonicFile,"r");
    if (fp==NULL)
    {
        printf("Error opening file %s\n", tonicFile);
        return 0;
    }
    nRead = fscanf(fp, "%f\n",&tonic);
    fclose(fp);
    
    
    //Finally opening the pitch file to read the data
    fp =fopen(pitchFile,"r");
    ind=0;
    jj=0;
    while (fscanf(fp, "%f\t%f\n",&timeTemp,&pitchTemp)!=EOF)    //read till the end of the file
    {
        
        if (pitchTemp>myProcParams->minPossiblePitch) //only fill in meaningful pitch data, reject other pitch samples
        {
            pitchSamples[ind]= (DATATYPE)round(temp1*log((pitchTemp+EPS)/tonic));
            timeSamples[ind] = timeTemp;
            ind++;
            jj++;
        }
        
        for(ii=0;ii<myProcParams->dsFactor-1;ii++)
        {
            // just bypass other samples to downsample the pitch sequence
            nRead = fscanf(fp, "%f\t%f\n",&timeTemp,&pitchTemp);
            jj++;
        }
        
              
    }
    totalPitchNonSilSamples = ind;
    fclose(fp);
    
    //################## reading seed motif information %%%%%%%%%%%%%%%%%%%%%
    fp = fopen(motifFile,"r");
    if (fp==NULL)
    {
        printf("Error opening file %s\n", motifFile);
        return 0;
    }
    seedMotifs = (segInfo_t*)malloc(sizeof(segInfo_t)*maxNMotifsPairs*2);
    jj=0;
    ii=0;
    while(fscanf(fp,"%f\t%f\t%f\t%f\t%f\t%f\t%f\n", &temp[0], &temp[1], &temp[2], &temp[3], &temp[4], &temp[5], &temp[6])!=EOF)
    {
        if(ii>=maxNMotifsPairs)
            break;
        
        if(temp[4]<INF)
        {
            seedMotifs[jj].str=temp[0];
            seedMotifs[jj].end=temp[1];
            jj++;
            
            seedMotifs[jj].str=temp[2];
            seedMotifs[jj].end=temp[3];
            jj++;
        }
        ii++;
    }
    NMotifs=jj;
    fclose(fp);
    
    //############ finding indexes where to start filling seed motif sequence #################
    seedMotifInd = (INDTYPE*)malloc(sizeof(INDTYPE)*NMotifs);
    for(ii=0;ii<NMotifs;ii++)
    {
        min_val = INF;
        for(jj=0;jj<totalPitchNonSilSamples;jj++)
        {
            if (fabs(timeSamples[jj]-seedMotifs[ii].str)<min_val)
            {
                min_val = fabs(timeSamples[jj]-seedMotifs[ii].str);
                seedMotifInd[ii] = jj;
            }
        }
    }
    
    data = (DATATYPE **)malloc(sizeof(DATATYPE *)*NMotifs);         //since we don't know valid subsequences, we allocate max possible subsequences and later discard them and free the memory
    tStamps = (segInfoInterp_t *)malloc(sizeof(segInfoInterp_t)*NMotifs); 
    
    for(ii=0;ii<NMotifs;ii++)
    {
        data[ii] = (DATATYPE *)malloc(sizeof(DATATYPE*)*lenMotifInterpH);
        tStamps[ii].str = timeSamples[seedMotifInd[ii]];
        tStamps[ii].end = timeSamples[seedMotifInd[ii]+lenMotifRealM1];
        tStamps[ii].endInterpH = timeSamples[seedMotifInd[ii]+lenMotifInterpHM1];
        tStamps[ii].endInterpL = timeSamples[seedMotifInd[ii]+lenMotifInterpLM1];
        for(jj=0;jj<lenMotifInterpH;jj++)
        {
            data[ii][jj]=pitchSamples[seedMotifInd[ii]+jj];
        }
        
    }
    N = NMotifs;
    dataInterp = (DATATYPE **)malloc(sizeof(DATATYPE *)*N*3);   //allocating memory also to have interpolated subsequences
    tStampsInterp = (segInfo_t *)malloc(sizeof(segInfo_t)*(N)*3);         //allocating memory also to have interpolated subsequences
    jj=0;
    //we do interpolation as well
    indLow = (float *)malloc(sizeof(float)*lenMotifReal);
    indHigh = (float *)malloc(sizeof(float)*lenMotifReal);
    for (ii=0;ii<lenMotifReal; ii++)
    {
        indLow[ii] = myProcParams->factorLow*ii;
        indHigh[ii] = myProcParams->factorHigh*ii;
    }
    
    for (ii=0;ii<N;ii++)
    {
        //low stretched subsequence
        dataInterp[jj] = (DATATYPE *)malloc(sizeof(DATATYPE)*lenMotifReal);
        cubicInterpolate(data[ii], dataInterp[jj], indLow, lenMotifReal);
        tStampsInterp[jj].str = tStamps[ii].str;
        tStampsInterp[jj].end = tStamps[ii].endInterpL;
        jj++;
        
        //normal
        dataInterp[jj] = data[ii];
        tStampsInterp[jj].str = tStamps[ii].str;
        tStampsInterp[jj].end = tStamps[ii].end;
        jj++;
        
        //compacted subsequence
        dataInterp[jj] = (DATATYPE *)malloc(sizeof(DATATYPE)*lenMotifReal);
        cubicInterpolate(data[ii], dataInterp[jj], indHigh, lenMotifReal);
        tStampsInterp[jj].str = tStamps[ii].str;
        tStampsInterp[jj].end = tStamps[ii].endInterpH;
        jj++;
        

    }
    lenTS = N*3; 
    free(data);
    free(tStamps);
    free(pitchSamples);
    free(timeSamples);
    
    *d = dataInterp;
    *t = tStampsInterp;
    return lenTS;
}

void dumpSearchMotifInfo(char *motifFile, char *mappFile, char *searchFile, motifInfo** topKmotifs, segInfo_t *tStampsInterpSeed, segInfo_t *tStampsInterp, int NSeeds, INDTYPE K, mappInfo_t *mapp, int verbos)
{
    FILE *fp;
    INDTYPE ii=0, lineWritten;
    int jj=0;
    int terminate=1;
    fp = fopen(motifFile, "ab");
    lineWritten=0;
    for(ii=0;ii<K;ii++)
    {
        terminate=1;
        for(jj=0;jj<NSeeds/3;jj++)
        {
            if(topKmotifs[jj][ii].dist<INF)
            {
                fprintf(fp, "%f\t%f\t%f\t%f\t%f\t", tStampsInterpSeed[topKmotifs[jj][ii].ind1].str, tStampsInterpSeed[topKmotifs[jj][ii].ind1].end, tStampsInterp[topKmotifs[jj][ii].ind2].str, tStampsInterp[topKmotifs[jj][ii].ind2].end, topKmotifs[jj][ii].dist);
                terminate=0;
                
            }
            else
            {
                fprintf(fp, "%f\t%f\t%f\t%f\t%f\t", -1.0, -1.0, -1.0, -1.0,-1.0);
                
            }
        }
        fprintf(fp, "\n");
        lineWritten++;
        if (terminate)
            break;
    }
    fclose(fp);
    
    fp = fopen(mappFile, "ab");
    fprintf(fp, "%s\t%lld\t%lld\n", searchFile, mapp->last_line, mapp->last_line+lineWritten-1);
    mapp->last_line+=lineWritten;
    fclose(fp);
    
}

void dumpDiscoveredMotifInfo(char *motifFile, motifInfo *topKmotifs, segInfo_t *tStampsInterp, int K, int verbos)
{
    FILE *fp;
    int ii;
    fp =fopen(motifFile,"w");
    for(ii=0;ii<K;ii++)
    {
        fprintf(fp, "%f\t%f\t%f\t%f\t%f\t%lld\t%lld\n", tStampsInterp[topKmotifs[ii].ind1].str, tStampsInterp[topKmotifs[ii].ind1].end, tStampsInterp[topKmotifs[ii].ind2].str, tStampsInterp[topKmotifs[ii].ind2].end, topKmotifs[ii].dist, topKmotifs[ii].ind1, topKmotifs[ii].ind2);
        if (verbos)
        {
            printf("motif pair is %f\t%f\t%f\t%lld\t%lld\n", tStampsInterp[topKmotifs[ii].ind1].str,tStampsInterp[topKmotifs[ii].ind2].str, topKmotifs[ii].dist, topKmotifs[ii].ind1%3, topKmotifs[ii].ind2%3);
        }
    }
    fclose(fp);
}

void dumpDiscoveryLogs(char *logFile, procLogs_t myProcLogs, int verbos)
{
    FILE *fp;
    
    fp =fopen(logFile,"w");
    fprintf(fp, "\nCommit Id of the code used for processing:\t%s\n", myProcLogs.commitID);
    fprintf(fp, "\n#################### TIME RELATED STATS ####################\n");
    fprintf(fp, "Time taken to load the pitch data:\t%f\n", myProcLogs.timeDataLoad);
    fprintf(fp, "Time taken to generate subsequences:\t%f\n", myProcLogs.timeGenSubs);
    fprintf(fp, "Time taken to remove blacklisted subsequences:\t%f\n", myProcLogs.timeRemBlacklist);
    fprintf(fp, "Time taken to generate envelops:\t%f\n", myProcLogs.timeGenEnvelops);
    fprintf(fp, "Time taken to discover patterns:\t%f\n", myProcLogs.timeDiscovery);
    fprintf(fp, "Time taken to write data:\t%f\n", myProcLogs.timeWriteData);
    fprintf(fp, "Total time taken by the process:\t%f\n", myProcLogs.timeTotal);
    
    fprintf(fp, "\n#################### DATA POINTS RELATED STATS ####################\n");
    fprintf(fp, "Total number of pitch samples in the file:\t%lld\n", myProcLogs.totalPitchSamples);
    fprintf(fp, "Total number of non zero pitch samples in the file:\t%lld\n", myProcLogs.totalPitchNonSilSamples);
    fprintf(fp, "Total number of subsequences generated originally:\t%lld\n", myProcLogs.totalSubsGenerated);
    fprintf(fp, "Total number of subsequences blacklisted:\t%lld\n", myProcLogs.totalSubsBlacklisted);
    fprintf(fp, "Total number of subsequences after interpolation:\t%lld\n", myProcLogs.totalSubsInterpolated);
    
    fprintf(fp, "\n#################### FNC CALLS RELATED STATS ####################\n");
    fprintf(fp, "Number of FL lowerbound is computed:\t%lld\n", myProcLogs.totalFLDone);
    fprintf(fp, "Number of LB_Keogh_EQ computed:\t%lld\n", myProcLogs.totalLBKeoghEQ);
    fprintf(fp, "Number of LB_Keogh_EC computed:\t%lld\n", myProcLogs.totalLBKeoghEC);
    fprintf(fp, "Number of times DTW computed:\t%lld\n", myProcLogs.totalDTWComputations);
    fprintf(fp, "Number of updates of priority list:\t%lld\n", myProcLogs.totalPriorityUpdates);
    fclose(fp);
    
    
    if (verbos)
    {
        printf("\nCommit Id of the code used for processing:\t%s\n", myProcLogs.commitID);
        printf("\n#################### TIME RELATED STATS ####################\n");
        printf("Time taken to load the pitch data:\t%f\n", myProcLogs.timeDataLoad);
        printf("Time taken to generate subsequences:\t%f\n", myProcLogs.timeGenSubs);
        printf("Time taken to remove blacklisted subsequences:\t%f\n", myProcLogs.timeRemBlacklist);
        printf("Time taken to generate envelops:\t%f\n", myProcLogs.timeGenEnvelops);
        printf("Time taken to discover patterns:\t%f\n", myProcLogs.timeDiscovery);
        printf("Time taken to write data:\t%f\n", myProcLogs.timeWriteData);
        printf("Total time taken by the process:\t%f\n", myProcLogs.timeTotal);
        
        printf("\n#################### DATA POINTS RELATED STATS ####################\n");
        printf("Total number of pitch samples in the file:\t%lld\n", myProcLogs.totalPitchSamples);
        printf("Total number of non zero pitch samples in the file:\t%lld\n", myProcLogs.totalPitchNonSilSamples);
        printf("Total number of subsequences generated originally:\t%lld\n", myProcLogs.totalSubsGenerated);
        printf("Total number of subsequences blacklisted:\t%lld\n", myProcLogs.totalSubsBlacklisted);
        printf("Total number of subsequences after interpolation:\t%lld\n", myProcLogs.totalSubsInterpolated);
        
        printf("\n#################### FNC CALLS RELATED STATS ####################\n");        
        printf("Number of FL lowerbound is computed:\t%lld\n", myProcLogs.totalFLDone);
        printf("Number of LB_Keogh_EQ computed:\t%lld\n", myProcLogs.totalLBKeoghEQ);
        printf("Number of LB_Keogh_EC computed:\t%lld\n", myProcLogs.totalLBKeoghEC);
        printf("Number of times DTW computed:\t%lld\n", myProcLogs.totalDTWComputations);
        printf("Number of updates of priority list:\t%lld\n", myProcLogs.totalPriorityUpdates);
    }
}

