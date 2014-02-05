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
    char *baseName, *pitchExt, *tonicExt, *segExt, *motifExt, pitchFile[200]={'\0'}, tonicFile[200]={'\0'}, segmentFile[200]={'\0'}, motifFile[200]={'\0'};
    float tonic,blackDur, durMotif, t1,t2, pHop, *timeSamples, minPossiblePitch, allowedSilDur, temp1, pitchTemp, timeTemp, *stdVec, *mean, std, varDur, threshold,ex;
    int lenMotif,lenMotifM1,  verbos=0, bandDTW, numReads,dsFactor, *blacklist,allowedSilSam, binsPOct, nRead; 
    INDTYPE    lenTS, count_DTW=0, ind, blackCnt=0;
    FILE *fp, *fp_out;
    long long numLinesInFile, pp=0;
    long long         K,ii,jj,ll, varSam, N;
    DATATYPE **data,**dataInterp, **U, **L, *accLB, pitch , ex2, *pitchSamples;
    DISTTYPE LB_Keogh_EQ, realDist,LB_Keogh_EC,bsf=INF,**costMTX, LB_kim_FL;
    motifInfo *topKmotifs;
    segInfo *taniSegs, *tStamps, *tStampsInterp;
    float temp[4]={0}, maxPauseDur, flatThreshold;
    
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
    
    //
    
    
    
    //############ CRUCIAL PARAMETERS ##################
    minPossiblePitch = 60.0;
    allowedSilDur = 0.15;
    binsPOct = 120;
    varDur = 0.1;
    threshold = 225;
    flatThreshold = 0.8;
    maxPauseDur = 0.5;
    
    
    //DERIVED
   
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
    //segment file name
    strcat(motifFile,baseName);
    strcat(motifFile,motifExt);
    
    // Reading number of lines in the pitch file
    numLinesInFile = getNumLines(pitchFile);
    
    // after downsampling we will be left with these many points
    numLinesInFile = floor(numLinesInFile/dsFactor) +1;
    
    //allocating memory for pitch and time samples
    pitchSamples = (DATATYPE*)malloc(sizeof(DATATYPE)*numLinesInFile);     // since we don't know silence regions, allocate maximum possible number of samples
    timeSamples = (float*)malloc(sizeof(float)*numLinesInFile);
    
    blacklist = (int*)malloc(sizeof(int)*numLinesInFile);  //subsequences which we don't have to consider. Unfortunately we know them after we already stored them
    memset(blacklist,0,sizeof(int)*numLinesInFile);
    
    data = (DATATYPE **)malloc(sizeof(DATATYPE *)*numLinesInFile);         //since we don't know valid subsequences, we allocate max possible subsequences and later discard them and free the memory
    tStamps = (segInfo *)malloc(sizeof(segInfo)*numLinesInFile);                
    
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
    lenMotif = (int)round(durMotif/pHop);
    varSam = (int)round(varDur/pHop);
    temp1 = ((float)binsPOct)/LOG2;
    lenMotifM1 = lenMotif-1;

    
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
    
    //########################## Subsequence generation + selection step ##########################
    // In subsequence selection our aim is to discard those subsequences which result into trivial matches, 
    // which are flat regions. For removing flat regions the obvious choice is to consider variance of a 
    // subsequence and filter using a threshold. The problem here is large amounts of octave errors because
    // of which the variance increases but affectively large chunk of data is still flat. 
    
    // For solving problem mentioned above we resort to short duration variance for deciding wheather a 
    // given sample belongs to a flat region or not. Later we accumulate total number of samples in a 
    // subsequence which belong to a flat region and filter the subsequence based on a threshold.
    
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
        
        if (ind<lenTS-lenMotifM1)
        {
            data[ind] = (DATATYPE *)malloc(sizeof(DATATYPE*)*lenMotif);
            tStamps[ind].end = timeSamples[ind+lenMotifM1];
            if (fabs(tStamps[ind].str - tStamps[ind].end) > durMotif + maxPauseDur)       //allow 200 ms pauses in total not more than that
                {
                    blacklist[ind]=1;
                }
        }
        
        for(ll = min(ind,lenTS-lenMotifM1-1) ; ll >= max(0,ind-lenMotif) ; ll--)
        {
            data[ll][ind-ll] = pitchSamples[ind]; 
        }
        
        ex+=stdVec[ind];
        if (ind >= lenMotifM1)
        {
            if (blacklist[ind-lenMotifM1]==0)
            {
                if(ex < flatThreshold*(float)lenMotif)
                {
                    
                    blacklist[ind-lenMotifM1]=1;
                }
                
            }
            

            ex -= stdVec[ind-lenMotifM1];
        }
       ind++;        
    }
    lenTS = lenTS-lenMotifM1;   //we only have lenMotifM1 samples less than original number lenTS. it still counts blacklisted candidates
    t2 = clock();
    printf("Time taken to load the data and create subsequences:%f\n",(t2-t1)/CLOCKS_PER_SEC);
    
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
    printf("Time taken to blacklist tani sections :%f\n",(t2-t1)/CLOCKS_PER_SEC);
    
    //%%%%%%%%%%%%%%%% Removing blacklisted subsequences %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    t1 = clock();
    //finding number of non blacklisted subsequences
    N=0;
    for (ii=0;ii<lenTS; ii++)
    {
        if (blacklist[ii]==0)
            N++;
    }
    dataInterp = (DATATYPE **)malloc(sizeof(DATATYPE *)*N*3);   //allocating memory also to have interpolated subsequences
    tStampsInterp = (segInfo *)malloc(sizeof(segInfo)*(N)*3);         //allocating memory also to have interpolated subsequences
    jj=0;
    for (ii=0;ii<lenTS;ii++)
    {
        if (blacklist[ii]==0)
        {
            dataInterp[jj] = data[ii];
            tStampsInterp[jj].str = tStamps[ii].str;
            tStampsInterp[jj].end = tStamps[ii].end;
            jj++;
        }
        else
        {
            free(data[ii]);
        }
    }
    lenTS = N; 
    free(data);
    free(tStamps);
    free(blacklist);
    
    t2 = clock();
    printf("Time taken to remove blacklist subsequences :%f\n",(t2-t1)/CLOCKS_PER_SEC);
    
    if( verbos == 1 )
        printf("Finally number of subsequences are: %lld\nNumber of subsequences removed are: %lld\n",lenTS,blackCnt-lenTS-lenMotifM1);
        printf("Length of Each Time Series : %d\n\n",lenMotif);
    
    //################# Precomputing envelope of each subsequence for the LB Keogh lower bound ###########################
    bandDTW = int(lenMotif*0.1);
    U = (DATATYPE **)malloc(sizeof(DATATYPE *)*lenTS);
    L= (DATATYPE **)malloc(sizeof(DATATYPE *)*lenTS);
    accLB = (DATATYPE *)malloc(sizeof(DATATYPE)*lenMotif);
    topKmotifs = (motifInfo *)malloc(sizeof(motifInfo)*K);
    costMTX = (DISTTYPE **)malloc(sizeof(DISTTYPE *)*lenMotif);

    t1 = clock();
    for (ii=0;ii<lenTS;ii++)
    {
        U[ii] = (DATATYPE *)malloc(sizeof(DATATYPE)*lenMotif);
        L[ii] = (DATATYPE *)malloc(sizeof(DATATYPE)*lenMotif);
        computeRunningMinMax(dataInterp[ii], U[ii], L[ii], lenMotif, bandDTW);
        
    }
    t2 = clock();
    printf("Time taken to generate envelopes for lower bound :%f\n",(t2-t1)/CLOCKS_PER_SEC);
    
    //initialization
    for(ii=0;ii<K;ii++)
    {
        topKmotifs[ii].dist = INF;
        topKmotifs[ii].ind1 = 0;
        topKmotifs[ii].ind2 = 0;
    }

    for (ii=0;ii<lenMotif;ii++)
        {
          costMTX[ii] = (DISTTYPE *)malloc(sizeof(DISTTYPE)*lenMotif);
          for (jj=0;jj<lenMotif;jj++)
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
            if (fabs(tStampsInterp[ii].str-tStampsInterp[jj].str)<blackDur)
            {
                continue;
            }
            
            LB_kim_FL = computeLBkimFL(dataInterp[ii][0], dataInterp[jj][0], dataInterp[ii][lenMotif-1], dataInterp[jj][lenMotif-1]);
            if (LB_kim_FL< bsf) 
            {
                LB_Keogh_EQ = computeKeoghsLB(U[ii],L[ii],accLB, dataInterp[jj],lenMotif, bsf);
                if(LB_Keogh_EQ < bsf)
                {
                    LB_Keogh_EC = computeKeoghsLB(U[jj],L[jj],accLB, dataInterp[ii],lenMotif, bsf);
                    if(LB_Keogh_EC < bsf)
                    {
                        realDist = dtw1dBandConst(dataInterp[ii], dataInterp[jj], lenMotif, lenMotif, costMTX, 0, bandDTW, bsf, accLB);
                        count_DTW+=1;
                        if(realDist<bsf)
                        {
                            bsf = manageTopKMotifs(topKmotifs, tStampsInterp, K, ii, jj, realDist, blackDur);
                        }
                    }
                }
            }
        }
    }

    //timing
    t2 = clock();
    
    if (verbos)
    {
        printf("Time taken to compute all combinations :%f\n",(t2-t1)/CLOCKS_PER_SEC);
    }
    fp =fopen(motifFile,"w");
    for(ii=0;ii<K;ii++)
    {
        fprintf(fp, "%f\t%f\t%f\t%f\t%f\n", tStampsInterp[topKmotifs[ii].ind1].str, tStampsInterp[topKmotifs[ii].ind1].end, tStampsInterp[topKmotifs[ii].ind2].str, tStampsInterp[topKmotifs[ii].ind2].end, topKmotifs[ii].dist);
        printf("motif pair is %f\t%f\t%f\n", tStampsInterp[topKmotifs[ii].ind1].str,tStampsInterp[topKmotifs[ii].ind2].str, topKmotifs[ii].dist);
    }
    printf("Total dtw computations %lld\n", count_DTW);
    fclose(fp);

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
    for(ii=0;ii<lenMotif;ii++)
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

