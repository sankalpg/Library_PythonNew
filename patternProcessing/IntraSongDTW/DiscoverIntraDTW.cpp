/******************************************************************************
*******************************************************************************
****** Sankalp Gulati                                                   *******
****** MUSIC TECHNOLOGY GROUP, UPF, BARCELONA                           *******
******                                                                  *******
*******************************************************************************
******************************************************************************/


#include "DiscoverIntraDTW.h"




#ifdef DEBUG
FILE *fpdb;
fpdb = fopen("dump.txt","w");
#endif


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

long long getNumLines2(const char *file)
{
    FILE *fp;
    char ch;
    long long line=0, ii;

    fp = fopen(file,"r");
    
    ch = getc(fp);
    while(ch!=EOF)
    {
        line+=1;
        ch = getc(fp);
    }
    
    return line;
    fclose(fp);
    
}



int main( int argc , char *argv[])
{
    char *baseName, *pitchExt, *tonicExt, pitchFile[200]={'\0'}, tonicFile[200]={'\0'};
    float tonic,blackDur, durMotif, t1,t2,*tStamps, pHop, *timeSamples, minPossiblePitch, allowedSilDur, temp1, pitchTemp, timeTemp, *stdVec, *mean, std, varDur, threshold,ex;
    int lenMotif,lenMotifM1,  verbos=0, band, numReads,dsFactor, *blacklist,allowedSilSam, binsPOct; 
    INDTYPE    lenTS, count_DTW=0, ind, blackCnt=0;
    FILE *fp, *fp_out;
    long long numPitchSam, pp=0;
    int         K,ii,jj,ll, mm, varSam, N;
    DATATYPE **data, *U, *L, *U2, *L2, *accLB, pitch , ex2, *pitchSamples;
    DISTTYPE LB_Keogh_EQ, realDist,LB_Keogh_EC,bsf=INF,**cost2, LB_kim_FL;
    motifInfo *topKmotifs;
    float temp[4]={0};
    
    if(argc < 9 || argc > 10)
    {
        printf("\nInvalid number of arguments!!!\n");
        exit(1);
    }
    
    baseName = argv[1];
    pitchExt = argv[2];
    tonicExt = argv[3];
    durMotif = atof(argv[4]);
    K = atoi(argv[5]);
    blackDur = atof(argv[6]);
    if (atof(argv[7])>0)
    {
        bsf = atof(argv[7]);
    }
    dsFactor = atoi(argv[8]);
    
    if( argc == 10 ){verbos = atoi(argv[9]);}
    
    //
    
    
    
    //############ CRUCIAL PARAMETERS ##################
    minPossiblePitch = 60.0;
    allowedSilDur = 0.15;
    binsPOct = 120;
    varDur = 0.1;
    threshold = 225;
    
    
    //DERIVED
   
    //########################## READING PITCH DATA ##########################
    //pitchFile[0]="\0";
    //tonicFile[0]="\0";
    strcat(pitchFile,baseName);
    strcat(pitchFile,pitchExt);
    
    strcat(tonicFile,baseName);
    strcat(tonicFile,tonicExt);
        
    // Reading number of lines in the pitch file
    numPitchSam = getNumLines(pitchFile);
    // after downsampling we will be left with these many points
    numPitchSam = floor(numPitchSam/dsFactor) +1;
    //allocating memory for pitch and time samples
    pitchSamples = (DATATYPE*)malloc(sizeof(DATATYPE)*numPitchSam);
    timeSamples = (float*)malloc(sizeof(float)*numPitchSam);
    blacklist = (int*)malloc(sizeof(int)*numPitchSam);
    
    data = (DATATYPE **)malloc(sizeof(DATATYPE *)*numPitchSam);
    tStamps = (float*)malloc(sizeof(float)*numPitchSam);
    stdVec = (float *)malloc(sizeof(float)*numPitchSam);
    mean = (float *)malloc(sizeof(float)*numPitchSam);
    
    fp =fopen(pitchFile,"r");
    fscanf(fp, "%f\t%f\n",&temp[0],&temp[1]);
    fscanf(fp, "%f\t%f\n",&temp[2],&temp[3]);
    fclose(fp);
    
    fp =fopen(strcat(baseName, tonicExt),"r");
    fscanf(fp, "%f\n",&tonic);
    fclose(fp);
    
    ind=0;
    jj=0;
    pHop = (temp[2]-temp[0])*dsFactor;
    lenMotif = (int)round(durMotif/pHop);
    varSam = (int)round(varDur/pHop);
    temp1 = ((float)binsPOct)/LOG2;
    lenMotifM1 = lenMotif-1;
    
    fp =fopen(pitchFile,"r");
    t1 = clock();
    ind=0;
    
    while ((fscanf(fp, "%f\t%f\n",&timeTemp,&pitchTemp)!=EOF)&&ind<numPitchSam)
    {
        if (pitchTemp>minPossiblePitch)
        {
            pitchSamples[ind]= (DATATYPE)(temp1*log((pitchTemp+EPS)/tonic));
            timeSamples[ind] = timeTemp;
            ind++;
        }
        
        for(ii=0;ii<dsFactor-1;ii++)
        {
            fscanf(fp, "%f\t%f\n",&timeTemp,&pitchTemp);
            jj++;
        }
        
        jj++;        
    }
    fclose(fp);
    
    //computing local mean and variance
    memset(mean,0,sizeof(float)*numPitchSam);
    memset(stdVec,0,sizeof(float)*numPitchSam);
    for(ii=0;ii<2*varSam+1;ii++)
    {
        mean[varSam] +=  pitchSamples[ii] ;
    }    
    for(ii=varSam+1;ii<ind-varSam;ii++)
    {
        mean[ii] =  mean[ii-1] - pitchSamples[ii-varSam-1] + pitchSamples[ii+varSam] ;
    }
    N = 2*varSam+1 ; 
    
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
    for(ii=0;ii<varSam;ii++)
    {
        stdVec[ii] = 1;
    }
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
    for(ii=ind-varSam;ii<ind;ii++)
    {
        stdVec[ii] = 1;
    }
    
    lenTS  = ind;
    ind=0;
    ex=0;
    blackCnt=0;
    while(ind<lenTS)
    {
        if (ind<lenTS-lenMotifM1)
            data[ind] = (DATATYPE *)malloc(sizeof(DATATYPE*)*lenMotif);
        
        for(ll = min(ind,lenTS-lenMotifM1-1) ; ll >= max(0,ind-lenMotif) ; ll--)
        {
            data[ll][ind-ll] = pitchSamples[ind]; 
            mm++;
        }
        //printf("%lld\n",ind);
        tStamps[ind] = timeSamples[ind];
        ex+=stdVec[ind];
        
        if (ind >= lenMotifM1)
        {
            if(ex>0.8*lenMotif)
            {
                blacklist[ind-lenMotif]=0;
            }
            else
            {
                blacklist[ind-lenMotif]=1;
                blackCnt++;
            }
            ex -= stdVec[ind-lenMotif];
            
        }
        ind++;        
    }
    lenTS  = ind-lenMotifM1;
    printf("%lld\t%lld",lenTS,blackCnt);
    t2 = clock();
    printf("Time taken to load the data :%f\n",(t2-t1)/CLOCKS_PER_SEC);
    
    /*fp = fopen("std.txt","w");
    for(ii=0;ii<lenTS;ii++)
    {
        fprintf(fp,"%f\t%f\n",tStamps[ii],stdVec[ii]);
    }
    fclose(fp);*/
    

    

    
    
    
    if( verbos == 1 )
        printf("\nNumber of Time Series : %lld\nLength of Each Time Series : %d\n\n",lenTS,lenMotif);
    
    
    //Memory allocation
    //data = (DATATYPE **)malloc(sizeof(DATATYPE *)*lenTS);
    //tStamps = (double *)malloc(sizeof(double)*lenTS);
    topKmotifs = (motifInfo *)malloc(sizeof(motifInfo)*K);
    cost2 = (DISTTYPE **)malloc(sizeof(DISTTYPE *)*lenMotif);
    U = (DATATYPE *)malloc(sizeof(DATATYPE)*lenMotif);
    L = (DATATYPE *)malloc(sizeof(DATATYPE)*lenMotif);
    U2 = (DATATYPE *)malloc(sizeof(DATATYPE)*lenMotif);
    L2 = (DATATYPE *)malloc(sizeof(DATATYPE)*lenMotif);    
    accLB = (DATATYPE *)malloc(sizeof(DATATYPE)*lenMotif);
    //initialization
    for(ii=0;ii<K;ii++)
    {
        topKmotifs[ii].dist = INF;
        topKmotifs[ii].ind1 = 0;
        topKmotifs[ii].ind2 = 0;
    }
    
    for(ii=0;ii<lenMotif;ii++)
    {
        cost2[ii] = (DISTTYPE *)malloc(sizeof(DISTTYPE)*lenMotif);
    }
    

    for (ii=0;ii<lenMotif;ii++)
        {
          for (jj=0;jj<lenMotif;jj++)
          {
              cost2[ii][jj]=FLT_MAX;
          }
          
        }
    
    //timing
    t1 = clock();
    //Computing all combinations to obtain top K best matches
    band = int(lenMotif*0.1);
    
    for(ii=0;ii<lenTS;ii++)
    {
        //computing lower and uper envelope for Keogh's lower bound
        computeRunningMinMax(data[ii], U, L, lenMotif, band);
        
        for(jj=ii+1;jj<lenTS;jj++)
        {

            if (fabs(tStamps[ii]-tStamps[jj])<blackDur)
            {
                continue;
            }
            
            LB_kim_FL = computeLBkimFL(data[ii][0], data[jj][0], data[ii][lenMotif-1], data[jj][lenMotif-1]);
            if (LB_kim_FL< bsf) 
            {
                LB_Keogh_EQ = computeKeoghsLB(U,L,accLB, data[jj],lenMotif, bsf);
                if(LB_Keogh_EQ < bsf)
                {
                    realDist = dtw1dBandConst(data[ii], data[jj], lenMotif, lenMotif, cost2, 0, band, bsf, accLB);
                    count_DTW+=1;
                    if(realDist<bsf)
                    {
                        bsf = manageTopKMotifs(topKmotifs, tStamps, K, ii, jj, realDist, blackDur);

#ifdef DEBUG                        
                        for(dd=0;dd<K;dd++)
                        {
                            fprintf(fpdb,"%lld\t%lld\t",ii,jj);
                            fprintf(fpdb,"%f\t",topKmotifs[dd].dist);
                            fprintf(fpdb,"\n",);
                        }
#endif

                        
                    } 
                        
                        
                    /*computeRunningMinMax(data[jj], U2, L2, lenMotif, band);
                    LB_Keogh_EC = computeKeoghsLB(U2,L2,accLB, data[ii],lenMotif, bsf);
                    
                    if(LB_Keogh_EC < bsf)
                    {
                        realDist = dtw1dBandConst(data[ii], data[jj], lenMotif, lenMotif, cost2, 0, band, bsf, accLB);
                        count_DTW+=1;
                        if(realDist<bsf)
                        {
                            bsf = manageTopKMotifs(topKDist, topKInd, tStamps, K, ii, jj, realDist, blackDur);
                        } 
                    }*/
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
    for(ii=0;ii<K;ii++)
    {
        printf("motif pair is %f\t%f\t%f\n", tStamps[topKmotifs[ii].ind1],tStamps[topKmotifs[ii].ind2], topKmotifs[ii].dist);
    }
    printf("Total dtw computations %lld\n", count_DTW);    

#ifdef DEBUG                        
                  fclose(fpdb);
#endif
    
    return 1;
    
}



DISTTYPE manageTopKMotifs(motifInfo *topKmotifs, float *tStamps, int K, INDTYPE ind1 , INDTYPE ind2, DISTTYPE dist, float blackDur)
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
        if ((fabs(tStamps[topKmotifs[ii].ind1]-tStamps[ind1]) < blackDur) || (fabs(tStamps[topKmotifs[ii].ind2]-tStamps[ind1]) < blackDur) || (fabs(tStamps[topKmotifs[ii].ind1]-tStamps[ind2]) < blackDur) || (fabs(tStamps[topKmotifs[ii].ind2]-tStamps[ind2]) < blackDur))
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