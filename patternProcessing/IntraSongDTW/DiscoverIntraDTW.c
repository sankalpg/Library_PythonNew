/******************************************************************************
*******************************************************************************
****** Sankalp Gulati                                                   *******
****** MUSIC TECHNOLOGY GROUP, UPF, BARCELONA                           *******
******                                                                  *******
*******************************************************************************
******************************************************************************/


#include "DiscoverIntraDTW.h"



/*
 * This function quickly returns number of lines in a text file
 */
long long getNumLines(const char *file)
{
    int fp;
    struct stat fs;
    char *buf;
    long long line=0, ii;

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


int main( int argc , char *argv[])
{
    char *baseName, *pitchExt, *tonicExt, *segExt, *motifExt, pitchFile[200]={'\0'}, tonicFile[200]={'\0'}, segmentFile[200]={'\0'}, motifFile[200]={'\0'}, logFile[200]={'\0'};
    float tonic,blackDur, durMotif, t1,t2, pHop, *timeSamples, minPossiblePitch, allowedSilDur, temp1, pitchTemp, timeTemp, *stdVec, *mean, std, varDur, threshold,ex;
    int lenMotifReal,lenMotifRealM1,lenMotifInterpHM1, lenMotifInterpLM1, lenMotifInterpH, lenMotifInterpL, verbos=0, bandDTW, numReads,dsFactor, *blacklist,allowedSilSam, binsPOct, nRead; 
    INDTYPE    lenTS, count_DTW=0, ind, blackCnt=0;
    FILE *fp, *fp_out;
    long long numLinesInFile, pp=0;
    long long         K,ii,jj,ll, varSam, N;
    DATATYPE **data,**dataInterp, **U, **L, *accLB, pitch , ex2, *pitchSamples;
    DISTTYPE LB_Keogh_EQ, realDist,LB_Keogh_EC,bsf=INF,**costMTX, LB_kim_FL;
    motifInfo *topKmotifs;
    segInfo *taniSegs, *tStampsInterp;
    segInfoInterp *tStamps;
    float temp[4]={0}, maxPauseDur, flatThreshold,factorLow, factorHigh;
    double *indNormal;
    float *indHigh, *indLow;
    procLogs myProcLogs;
    char commitID[] = "6f58f1c6eba3b863ba7945e4bc0a8cd997e84f97";
    myProcLogs.commitID = commitID;
    myProcLogs.timeDataLoad=0;
    myProcLogs.timeGenSubs=0;
    myProcLogs.timeRemBlacklist=0; 
    myProcLogs.timeGenEnvelops=0;
    myProcLogs.timeDiscovery=0;
    myProcLogs.totalPitchSamples=0;
    myProcLogs.totalPitchNonSilSamples=0;
    myProcLogs.totalSubsGenerated=0;
    myProcLogs.totalSubsBlacklisted=0;
    myProcLogs.totalSubsInterpolated=0;
    myProcLogs.totalFLDone=0;
    myProcLogs.totalLBKeoghEQ=0;
    myProcLogs.totalLBKeoghEC=0;
    myProcLogs.totalDTWComputations=0;
    myProcLogs.totalPriorityUpdates=0;
    
    if(argc < 11 || argc > 12)
    {
        printf("\nInvalid number of arguments!!!\n");
        exit(1);
    }
    
    baseName = argv[1];
    pitchExt = argv[2];
    tonicExt = argv[3];
    segExt = argv[4];
    motifExt = argv[5];
    durMotif = atof(argv[6]);
    K = atoi(argv[7]);
    blackDur = atof(argv[8]);
    if (atof(argv[9])>0)
    {
        bsf = atof(argv[9]);
    }
    dsFactor = atoi(argv[10]);
    
    if( argc == 12 ){verbos = atoi(argv[11]);}
    
    //############ CRUCIAL PARAMETERS ##################
    minPossiblePitch = 60.0;
    allowedSilDur = 0.15;
    binsPOct = 120;
    varDur = 0.1;
    threshold = 225;
    flatThreshold = 0.8;
    maxPauseDur = 0.5;
    factorLow = 0.9;
    factorHigh = 1.1;
    
    //########################## READING PITCH DATA ##########################
    //pitch file name
    strcat(pitchFile,baseName);
    strcat(pitchFile,pitchExt);
    //tonic file name
    strcat(tonicFile,baseName);
    strcat(tonicFile,tonicExt);
    //segment file name
    strcat(segmentFile,baseName);
    strcat(segmentFile,segExt);
    //motif file name
    strcat(motifFile,baseName);
    strcat(motifFile,motifExt);
    //log file name
    strcat(logFile,baseName);
    strcat(logFile,".proclog");    
    
    // Reading number of lines in the pitch file
    numLinesInFile = getNumLines(pitchFile);
    myProcLogs.totalPitchSamples = numLinesInFile;
    
    // after downsampling we will be left with these many points
    numLinesInFile = floor(numLinesInFile/dsFactor) +1;
    
    //allocating memory for pitch and time samples
    pitchSamples = (DATATYPE*)malloc(sizeof(DATATYPE)*numLinesInFile);     // since we don't know silence regions, allocate maximum possible number of samples
    timeSamples = (float*)malloc(sizeof(float)*numLinesInFile);
    
    blacklist = (int*)malloc(sizeof(int)*numLinesInFile);  //subsequences which we don't have to consider. Unfortunately we know them after we already stored them
    memset(blacklist,0,sizeof(int)*numLinesInFile);
    
    data = (DATATYPE **)malloc(sizeof(DATATYPE *)*numLinesInFile);         //since we don't know valid subsequences, we allocate max possible subsequences and later discard them and free the memory
    tStamps = (segInfoInterp *)malloc(sizeof(segInfoInterp)*numLinesInFile);                
    
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
    pHop = (temp[2]-temp[0])*dsFactor;  //final hop size afte downsampling
    lenMotifReal = (int)round(durMotif/pHop);
    lenMotifInterpH = (int)ceil((durMotif*factorHigh)/pHop)+1;  //adding one because we need one extra sample for cubic interpolation
    lenMotifInterpL = (int)round((durMotif*factorLow)/pHop);
    varSam = (int)round(varDur/pHop);
    temp1 = ((float)binsPOct)/LOG2;
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
        
        if (pitchTemp>minPossiblePitch) //only fill in meaningful pitch data, reject other pitch samples
        {
            pitchSamples[ind]= (DATATYPE)round(temp1*log((pitchTemp+EPS)/tonic));
            timeSamples[ind] = timeTemp;
            ind++;
            jj++;
        }
        
        for(ii=0;ii<dsFactor-1;ii++)
        {
            // just bypass other samples to downsample the pitch sequence
            nRead = fscanf(fp, "%f\t%f\n",&timeTemp,&pitchTemp);
            jj++;
        }
        
              
    }
    fclose(fp);
    t2 = clock();
    myProcLogs.timeDataLoad = (t2-t1)/CLOCKS_PER_SEC;
    myProcLogs.totalPitchNonSilSamples = ind;
    if (verbos)
    {printf("Time taken to load the pitch data :%f\n",myProcLogs.timeDataLoad);}
    
    //########################## Subsequence generation + selection step ##########################
    // In subsequence selection our aim is to discard those subsequences which result into trivial matches, 
    // which are flat regions. For removing flat regions the obvious choice is to consider variance of a 
    // subsequence and filter using a threshold. The problem here is large amounts of octave errors because
    // of which the variance increases but affectively large chunk of data is still flat. 
    
    // For solving problem mentioned above we resort to short duration variance for deciding wheather a 
    // given sample belongs to a flat region or not. Later we accumulate total number of samples in a 
    // subsequence which belong to a flat region and filter the subsequence based on a threshold.
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
    // after computing running variance, applying a threshold and selecting the onces which are corresponsing to non flat regions
    for(ii=varSam;ii<ind-varSam;ii++)
    {
        if (stdVec[ii]>threshold)
        {
            stdVec[ii] = 1;
        }
        else
        {
            stdVec[ii] = 0;
        }            
            
    }
    
    lenTS  = ind;       //number of non trivial pitch samples, they might still carry number of subsequences which have to be discarded
    blackCnt=lenTS;
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
            if (fabs(tStamps[ind].str - tStamps[ind].endInterpH) > durMotif*factorHigh + maxPauseDur)       //allow 200 ms pauses in total not more than that
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
                if(ex < flatThreshold*(float)lenMotifInterpH)
                {
                    
                    blacklist[ind-lenMotifInterpHM1]=1;
                }
                
            }
            
            ex -= stdVec[ind-lenMotifInterpHM1];
        }
       ind++;        
    }
    lenTS = lenTS-lenMotifInterpHM1;   //we only have lenMotifInterpHM1 samples less than original number lenTS. it still counts blacklisted candidates
    myProcLogs.totalSubsGenerated = lenTS;
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
    taniSegs = (segInfo *)malloc(sizeof(segInfo)*numLinesInFile);
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
    myProcLogs.totalSubsBlacklisted = myProcLogs.totalSubsGenerated-N;
    dataInterp = (DATATYPE **)malloc(sizeof(DATATYPE *)*N*3);   //allocating memory also to have interpolated subsequences
    tStampsInterp = (segInfo *)malloc(sizeof(segInfo)*(N)*3);         //allocating memory also to have interpolated subsequences
    jj=0;
    
    //we do interpolation as well
    
    indLow = (float *)malloc(sizeof(float)*lenMotifReal);
    indHigh = (float *)malloc(sizeof(float)*lenMotifReal);
    for (ii=0;ii<lenMotifReal; ii++)
    {
        indLow[ii] = factorLow*ii;
        indHigh[ii] = factorHigh*ii;
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
    myProcLogs.timeRemBlacklist = (t2-t1)/CLOCKS_PER_SEC;
    myProcLogs.totalSubsInterpolated = lenTS;
    if (verbos)
    {printf("Time taken to remove blacklist subsequences :%f\n", myProcLogs.timeRemBlacklist);}
    
    if( verbos == 1 )
        printf("Finally number of subsequences are: %lld\nNumber of subsequences removed are: %lld\n",lenTS,myProcLogs.totalSubsBlacklisted);
        printf("Length of Each Time Series : %d\n\n",lenMotifReal);
    
    //################# Precomputing envelope of each subsequence for the LB Keogh lower bound ###########################
    bandDTW = int(lenMotifReal*0.1);
    U = (DATATYPE **)malloc(sizeof(DATATYPE *)*lenTS);
    L= (DATATYPE **)malloc(sizeof(DATATYPE *)*lenTS);
    accLB = (DATATYPE *)malloc(sizeof(DATATYPE)*lenMotifReal);
    topKmotifs = (motifInfo *)malloc(sizeof(motifInfo)*K);
    costMTX = (DISTTYPE **)malloc(sizeof(DISTTYPE *)*lenMotifReal);

    t1 = clock();
    for (ii=0;ii<lenTS;ii++)
    {
        U[ii] = (DATATYPE *)malloc(sizeof(DATATYPE)*lenMotifReal);
        L[ii] = (DATATYPE *)malloc(sizeof(DATATYPE)*lenMotifReal);
        computeRunningMinMax(dataInterp[ii], U[ii], L[ii], lenMotifReal, bandDTW);
        
    }
    t2 = clock();
    myProcLogs.timeGenEnvelops = (t2-t1)/CLOCKS_PER_SEC;
    if (verbos)
    {printf("Time taken to generate envelopes for lower bound :%f\n",myProcLogs.timeGenEnvelops);}
    
    //initialization
    for(ii=0;ii<K;ii++)
    {
        topKmotifs[ii].dist = INF;
        topKmotifs[ii].ind1 = 0;
        topKmotifs[ii].ind2 = 0;
    }

    for (ii=0;ii<lenMotifReal;ii++)
        {
          costMTX[ii] = (DISTTYPE *)malloc(sizeof(DISTTYPE)*lenMotifReal);
          for (jj=0;jj<lenMotifReal;jj++)
          {
              costMTX[ii][jj]=FLT_MAX;
          }
          
        }
    
    //timing
    t1 = clock();
    
    //Computing all combinations to obtain top K best matches
    for(ii=0;ii<lenTS;ii++)
    {
        for(jj=ii+1;jj<lenTS;jj++)
        {
            // So we have three subs for every original sub. Its low interp, normal and high interp (if we consider 2 original subs we have now 9 combinations of 3 variants of each). Based on some intuitions we remove combinations which repeat (approximately repeat) like low1-normal2 would be really be close to normal1-high2 and so on and so forth. We removed 4 such combinations out of 9. This saves us a lot of computations
            if (((ii%3==0)&&(jj%3==0))||((ii%3==2)&&(jj%3==2))||((ii%3==0)&&(jj%3==1))||((ii%3==2)&&(jj%3==1)))
                continue;
            if (fabs(tStampsInterp[ii].str-tStampsInterp[jj].str)<blackDur)
            {
                continue;
            }
            
            LB_kim_FL = computeLBkimFL(dataInterp[ii][0], dataInterp[jj][0], dataInterp[ii][lenMotifReal-1], dataInterp[jj][lenMotifReal-1]);
            myProcLogs.totalFLDone++;
            if (LB_kim_FL< bsf) 
            {
                LB_Keogh_EQ = computeKeoghsLB(U[ii],L[ii],accLB, dataInterp[jj],lenMotifReal, bsf);
                myProcLogs.totalLBKeoghEQ++;
                if(LB_Keogh_EQ < bsf)
                {
                    LB_Keogh_EC = computeKeoghsLB(U[jj],L[jj],accLB, dataInterp[ii],lenMotifReal, bsf);
                    myProcLogs.totalLBKeoghEC++;
                    if(LB_Keogh_EC < bsf)
                    {
                        realDist = dtw1dBandConst(dataInterp[ii], dataInterp[jj], lenMotifReal, lenMotifReal, costMTX, 0, bandDTW, bsf, accLB);
                        myProcLogs.totalDTWComputations++;
                        if(realDist<bsf)
                        {
                            bsf = manageTopKMotifs(topKmotifs, tStampsInterp, K, ii, jj, realDist, blackDur);
                            myProcLogs.totalPriorityUpdates++;
                        }
                    }
                }
            }
        }
    }

    //timing
    t2 = clock();
    myProcLogs.timeDiscovery = (t2-t1)/CLOCKS_PER_SEC;
    if (verbos)
    {printf("Time taken to compute all combinations :%f\n",myProcLogs.timeDiscovery);}
    
    fp =fopen(motifFile,"w");
    for(ii=0;ii<K;ii++)
    {
        fprintf(fp, "%f\t%f\t%f\t%f\t%f\t%lld\t%lld\n", tStampsInterp[topKmotifs[ii].ind1].str, tStampsInterp[topKmotifs[ii].ind1].end, tStampsInterp[topKmotifs[ii].ind2].str, tStampsInterp[topKmotifs[ii].ind2].end, topKmotifs[ii].dist, topKmotifs[ii].ind1, topKmotifs[ii].ind2);
        printf("motif pair is %f\t%f\t%f\t%lld\t%lld\n", tStampsInterp[topKmotifs[ii].ind1].str,tStampsInterp[topKmotifs[ii].ind2].str, topKmotifs[ii].dist, topKmotifs[ii].ind1%3, topKmotifs[ii].ind2%3);
    }
    fclose(fp);
    
    
    fp =fopen(logFile,"w");
    fprintf(fp, "\nCommit Id of the code used for processing:\t%s\n", myProcLogs.commitID);
    fprintf(fp, "\n#################### TIME RELATED STATS ####################\n");
    fprintf(fp, "Time taken to load the pitch data:\t%f\n", myProcLogs.timeDataLoad);
    fprintf(fp, "Time taken to generate subsequences:\t%f\n", myProcLogs.timeGenSubs);
    fprintf(fp, "Time taken to remove blacklisted subsequences:\t%f\n", myProcLogs.timeRemBlacklist);
    fprintf(fp, "Time taken to generate envelops:\t%f\n", myProcLogs.timeGenEnvelops);
    fprintf(fp, "Time taken to discover patterns:\t%f\n", myProcLogs.timeDiscovery);
    
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
    
    //Memory clearing
    for(ii=0;ii<lenTS;ii++)
    {
        free(dataInterp[ii]);
        free(U[ii]);
        free(L[ii]);
        
    }
    free(dataInterp);
    free(accLB);
    free(tStampsInterp);
    free(U);
    free(L);
    free(topKmotifs);
    for(ii=0;ii<lenMotifReal;ii++)
    {
        free(costMTX[ii]);
    }
    free(costMTX);
    return 1;
    
}



DISTTYPE manageTopKMotifs(motifInfo *topKmotifs, segInfo *tStamps, int K, INDTYPE ind1 , INDTYPE ind2, DISTTYPE dist, float blackDur)
{
    int ii=0;
    int sortInd = -1;
    int matchInd = -1;
    
    for(ii=0;ii<K;ii++)
    {
        if ((topKmotifs[ii].dist > dist)&&(sortInd==-1))
        {
            sortInd=ii;
        }
        // searching if we already have a motif in out top K list which is near to the currently good match
        if ((fabs(tStamps[topKmotifs[ii].ind1].str-tStamps[ind1].str) < blackDur) || (fabs(tStamps[topKmotifs[ii].ind2].str-tStamps[ind1].str) < blackDur) || (fabs(tStamps[topKmotifs[ii].ind1].str-tStamps[ind2].str) < blackDur) || (fabs(tStamps[topKmotifs[ii].ind2].str-tStamps[ind2].str) < blackDur))
        {
            matchInd=ii;
            break;
        }
        
    }
    if (sortInd==-1)//we couldn't get any satisfactory replacement before we get a close neighbour
    {
        return topKmotifs[K-1].dist;
    }
    //There are three possibilities
    //1) There is no match found in the existing top motifs, simplest
    if (matchInd==-1)
    {
        memmove(&topKmotifs[sortInd+1], &topKmotifs[sortInd], sizeof(motifInfo)*(K-(sortInd+1)));
        topKmotifs[sortInd].dist = dist;
        topKmotifs[sortInd].ind1 = ind1;
        topKmotifs[sortInd].ind2 = ind2;
    }
    else if (sortInd == matchInd)
    {
        topKmotifs[sortInd].dist = dist;
        topKmotifs[sortInd].ind1 = ind1;
        topKmotifs[sortInd].ind2 = ind2;
    }
    else if (sortInd < matchInd)
    {
        memmove(&topKmotifs[sortInd+1], &topKmotifs[sortInd], sizeof(motifInfo)*(matchInd-sortInd));
        topKmotifs[sortInd].dist = dist;
        topKmotifs[sortInd].ind1 = ind1;
        topKmotifs[sortInd].ind2 = ind2;
    }
    
    return topKmotifs[K-1].dist;
    
}

void linearlyInterpolate(float *array, int size, float val1, float val2)
{
    int ii=0;
    float diff;
    
    diff = (val2-val1)/(size+1);
    
    
    for(ii=1;ii<=size;ii++)
    {
        array[ii] = val1 + ii*diff;
    }
}



