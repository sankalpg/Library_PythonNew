/******************************************************************************
*******************************************************************************
****** Sankalp Gulati                                                   *******
****** MUSIC TECHNOLOGY GROUP, UPF, BARCELONA                           *******
******                                                                  *******
*******************************************************************************
******************************************************************************/


#include "SearchInterDTW.h"
//#define DEBUG_GENERATION

int compareMotifInfo(const void *a, const void *b)
{
    if (((motifInfo*)a)->dist > ((motifInfo*)b)->dist)
    {
        return 1;
    }
    else if (((motifInfo*)a)->dist < ((motifInfo*)b)->dist)
    {
        return -1;
    }
    return 0;
}


int main( int argc , char *argv[])
{
    FILE *fp;
    char *baseName, motifFile[400]={'\0'}, logFile[400]={'\0'}, searchFileList[400]={'\0'}, searchFile[400] = {'\0'}, mappFile[400] = {'\0'}, paramOutFile[400]={'\0'} ;
    float t1,t2, t3,t4;
    int lenMotifReal, verbos=0, bandDTW, maxNMotifsPairs, nInterFact, **combMTX; 
    INDTYPE    NSeed, lenTS, K,ii,jj;
    bool sameFile;
    
    DATATYPE **dataInterpSeed, **dataInterp, **U, **L, **USeed, **LSeed, *accLB;
    DISTTYPE LB_Keogh_EQ, realDist,LB_Keogh_EC,bsf=INF,**costMTX, LB_kim_FL;
    motifInfo **topKmotifs;
    segInfo_t *tStampsInterp, *tStampsInterpSeed;
    procLogs_t myProcLogs;
    procParams_t myProcParams;
    fileExts_t myFileExts;
    mappInfo_t mapp;
    mapp.last_line =1;
    
    t3=clock();
    
    char commitID[] = "6f58f1c6eba3b863ba7945e4bc0a8cd997e84f97";
    myProcLogs.commitID = commitID;
    myProcLogs.timeDataLoad=0;
    myProcLogs.timeGenSubs=0;
    myProcLogs.timeRemBlacklist=0; 
    myProcLogs.timeGenEnvelops=0;
    myProcLogs.timeDiscovery=0;
    myProcLogs.timeWriteData=0;
    myProcLogs.timeTotal=0;
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
    
    
    if(argc < 18 || argc > 19)
    {
        printf("\nInvalid number of arguments!!!\n");
        exit(1);
    }
    
    baseName = argv[1];
    myFileExts.pitchExt = argv[2];
    myFileExts.tonicExt = argv[3];
    myFileExts.segExt = argv[4];
    myFileExts.searchExt = argv[5];
    myFileExts.seedMotifExt = argv[6];
    myFileExts.motifExt = argv[7];
    myFileExts.mappExt = argv[8];
    myFileExts.logExt = argv[9];
    myFileExts.paramOutExt = argv[10];
    myProcParams.durMotif = atof(argv[11]);
    K = atoi(argv[12]);
    myProcParams.blackDur = atof(argv[13]);
    if (atof(argv[14])>0)
    {
        bsf = atof(argv[14]);
    }
     myProcParams.dsFactor = atoi(argv[15]);
     maxNMotifsPairs = atoi(argv[16]);
     myProcParams.nInterpFac=atoi(argv[17]);
    
    if( argc == 19 ){verbos = atoi(argv[18]);}
    
    //############ CRUCIAL PARAMETERS ##################
    myProcParams.minPossiblePitch = 60.0;
    myProcParams.binsPOct = 1200;
    myProcParams.varDur = 0.1;
    myProcParams.threshold = 45.2267016866645;
    myProcParams.flatThreshold = 0.8;
    myProcParams.maxPauseDur = 0.5;
    myProcParams.DTWBand = 0.1;
    myProcParams.removeTaniSegs=1;
    
    
    if (myProcParams.nInterpFac==1)
    {
        myProcParams.interpFac[0]=1.0;
        
        myProcParams.combMTX = (int **)malloc(sizeof(int*)*myProcParams.nInterpFac);
        for(ii=0;ii<myProcParams.nInterpFac;ii++)
        {
            myProcParams.combMTX[ii] =  (int *)malloc(sizeof(int)*myProcParams.nInterpFac);
            for(jj=0;jj<myProcParams.nInterpFac;jj++)
            {
                myProcParams.combMTX[ii][jj] = 1;
            }
        }        
    }
    else if (myProcParams.nInterpFac==3)
    {
        myProcParams.interpFac[0]=0.9;
        myProcParams.interpFac[1]=1.0;
        myProcParams.interpFac[2]=1.1;
        
        myProcParams.combMTX = (int **)malloc(sizeof(int*)*myProcParams.nInterpFac);
        for(ii=0;ii<myProcParams.nInterpFac;ii++)
        {
            myProcParams.combMTX[ii] =  (int *)malloc(sizeof(int)*myProcParams.nInterpFac);
            for(jj=0;jj<myProcParams.nInterpFac;jj++)
            {
                myProcParams.combMTX[ii][jj] = combAllwd_3[ii][jj];
            }
        }        
    }
    else if (myProcParams.nInterpFac==5)
    {
        myProcParams.interpFac[0]=0.9;
        myProcParams.interpFac[1]=0.95;
        myProcParams.interpFac[2]=1.0;
        myProcParams.interpFac[3]=1.05;
        myProcParams.interpFac[4]=1.1;
        
        myProcParams.combMTX = (int **)malloc(sizeof(int*)*myProcParams.nInterpFac);
        for(ii=0;ii<myProcParams.nInterpFac;ii++)
        {
            myProcParams.combMTX[ii] =  (int *)malloc(sizeof(int)*myProcParams.nInterpFac);
            for(jj=0;jj<myProcParams.nInterpFac;jj++)
            {
                myProcParams.combMTX[ii][jj] = combAllwd_5[ii][jj];
            }
        }        
    }
    nInterFact = myProcParams.nInterpFac;
    combMTX = myProcParams.combMTX;    
    
    //####################################################
    //motif file name
    strcat(motifFile,baseName);
    strcat(motifFile,myFileExts.motifExt);
    //motif file name
    strcat(mappFile,baseName);
    strcat(mappFile,myFileExts.mappExt);
    //log file name
    strcat(logFile,baseName);
    strcat(logFile,myFileExts.logExt);
    //search file name
    strcat(searchFileList,baseName);
    strcat(searchFileList,myFileExts.searchExt); 
    //paramOut file name
    strcat(paramOutFile,baseName);
    strcat(paramOutFile,myFileExts.paramOutExt); 
    
    //erasing motif and mapp file
    fp = fopen(motifFile,"w");
    fclose(fp);
    fp = fopen(mappFile,"w");
    fclose(fp);
    
    
    // loading sequence corresponding to seed motifs
    NSeed = loadSeedMotifSequence(&dataInterpSeed, &tStampsInterpSeed, &lenMotifReal, baseName, &myFileExts, &myProcParams, &myProcLogs, maxNMotifsPairs, verbos);
    
#ifdef DEBUG_GENERATION
    fp = fopen("subsequences.bin","wb");
    for(ii=0;ii<NSeed;ii++)
    {
        for(jj=0;jj<lenMotifReal;jj++)
        {
            fprintf(fp, "%f\t", dataInterpSeed[ii][jj]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
    return 1;
    
#endif      
    
    // generating envelops for the seed motifs
    bandDTW = (int)floor(lenMotifReal*myProcParams.DTWBand);
    USeed = (DATATYPE **)malloc(sizeof(DATATYPE *)*NSeed);
    LSeed= (DATATYPE **)malloc(sizeof(DATATYPE *)*NSeed);
    accLB = (DATATYPE *)malloc(sizeof(DATATYPE)*lenMotifReal);
    topKmotifs = (motifInfo **)malloc(sizeof(motifInfo*)*NSeed/nInterFact);
    
    for (ii=0;ii<NSeed;ii++)
    {
        USeed[ii] = (DATATYPE *)malloc(sizeof(DATATYPE)*lenMotifReal);
        LSeed[ii] = (DATATYPE *)malloc(sizeof(DATATYPE)*lenMotifReal);
        computeRunningMinMax(dataInterpSeed[ii], USeed[ii], LSeed[ii], lenMotifReal, bandDTW);
        
    }
    
    
    costMTX = (DISTTYPE **)malloc(sizeof(DISTTYPE *)*lenMotifReal);
    for (ii=0;ii<lenMotifReal;ii++)
    {
        costMTX[ii] = (DISTTYPE *)malloc(sizeof(DISTTYPE)*lenMotifReal);
        for (jj=0;jj<lenMotifReal;jj++)
        {
            costMTX[ii][jj]=FLT_MAX;
        }
        
    }
    
    
    // iterating over all files 
    fp = fopen(searchFileList, "r");
    while(fscanf(fp, "%s\n",searchFile)!=EOF)
    {
        //generating subsequence database for file to be searched
        lenTS = readPreProcessGenDB(&dataInterp, &tStampsInterp, &lenMotifReal, searchFile, &myFileExts, &myProcParams, &myProcLogs, verbos);
        
        t1=clock();
        //computing envelops for the file to be searched
        U = (DATATYPE **)malloc(sizeof(DATATYPE *)*lenTS);
        L= (DATATYPE **)malloc(sizeof(DATATYPE *)*lenTS);
        
        for (ii=0;ii<lenTS;ii++)
        {
            U[ii] = (DATATYPE *)malloc(sizeof(DATATYPE)*lenMotifReal);
            L[ii] = (DATATYPE *)malloc(sizeof(DATATYPE)*lenMotifReal);
            computeRunningMinMax(dataInterp[ii], U[ii], L[ii], lenMotifReal, bandDTW);
            
        }
        t2=clock();
        myProcLogs.timeGenEnvelops += (t2-t1)/CLOCKS_PER_SEC;
        for(jj=0;jj<NSeed/nInterFact;jj++)
        {
            topKmotifs[jj] = (motifInfo *)malloc(sizeof(motifInfo)*K);
            for(ii=0;ii<K;ii++)
            {
                topKmotifs[jj][ii].dist = INF;
                topKmotifs[jj][ii].ind1 = 0;
                topKmotifs[jj][ii].ind2 = 0;
            }
        }
        
        //############## Performing a search using basic DTW ########################
        t1=clock();
        sameFile = strcmp(baseName, searchFile);
        for(ii=0;ii<NSeed;ii++)
        {
            for(jj=0;jj<lenTS;jj++)
            {
                if (myProcParams.combMTX[ii%nInterFact][jj%nInterFact]==0)
                {
                    continue;
                }
                if ((sameFile)&&(fabs(tStampsInterpSeed[ii].str-tStampsInterp[jj].str)< myProcParams.blackDur))
                {
                    continue;
                }
                
                LB_kim_FL = computeLBkimFL(dataInterpSeed[ii][0], dataInterp[jj][0], dataInterpSeed[ii][lenMotifReal-1], dataInterp[jj][lenMotifReal-1], SqEuclidean);
                myProcLogs.totalFLDone++;
                if (LB_kim_FL< bsf)
                {
                    LB_Keogh_EQ = computeKeoghsLB(USeed[ii],LSeed[ii], accLB, dataInterp[jj],lenMotifReal, bsf, SqEuclidean);
                    myProcLogs.totalLBKeoghEQ++;
                    if(LB_Keogh_EQ < bsf)
                    {
                        LB_Keogh_EC = computeKeoghsLB(U[jj],L[jj],accLB, dataInterpSeed[ii],lenMotifReal, bsf, SqEuclidean);
                        myProcLogs.totalLBKeoghEC++;
                        if(LB_Keogh_EC < bsf)
                        {
                            realDist = dtw1dBandConst(dataInterpSeed[ii], dataInterp[jj], lenMotifReal, lenMotifReal, costMTX, SqEuclidean, bandDTW, bsf, accLB);
                            myProcLogs.totalDTWComputations++;
                            if(realDist<bsf)
                            {
                                bsf = manageTopKMotifs(topKmotifs[ii/nInterFact], tStampsInterpSeed, tStampsInterp, K, ii, jj, realDist, myProcParams.blackDur);
                                myProcLogs.totalPriorityUpdates++;
                            }
                        }
                    }
                }
            }
        } 

        //############## Rank Refinement using sophisticated DTW ########################
        
        //recomputing the distance
        for(ii=0;ii<NSeed/nInterFact;ii++)
        {
            for(jj=0;jj<K;jj++)
            {
                topKmotifs[ii][jj].dist = dtw1dBandConst_localConst(dataInterpSeed[topKmotifs[ii][jj].ind1], dataInterp[topKmotifs[ii][jj].ind2], lenMotifReal, lenMotifReal, costMTX, SqEuclidean, bandDTW, INF, accLB);
            }
            //sorting the priority list
            qsort (topKmotifs[ii], K, sizeof(motifInfo), compareMotifInfo);
        }   
            
        t2=clock();
        myProcLogs.timeDiscovery += (t2-t1)/CLOCKS_PER_SEC;
        t1=clock();
        dumpSearchMotifInfo(motifFile, mappFile, searchFile, topKmotifs, tStampsInterpSeed, tStampsInterp, NSeed, K, &mapp, nInterFact, verbos);
        t2=clock();
        myProcLogs.timeWriteData += (t2-t1)/CLOCKS_PER_SEC;
        
        for(ii=0;ii<lenTS;ii++)
        {
            free(dataInterp[ii]);
            free(U[ii]);
            free(L[ii]);
        }
        for(jj=0;jj<NSeed/nInterFact;jj++)
            {
                free(topKmotifs[jj]);
            }
        free(U);
        free(L);
        free(dataInterp);
        free(tStampsInterp);
            
        
    }
    
    fclose(fp);
    
    //Memory clearing
    for(ii=0;ii<NSeed;ii++)
    {
        free(dataInterpSeed[ii]);
        free(USeed[ii]);
        free(LSeed[ii]);
        
    }
    free(dataInterpSeed);
    free(USeed);
    free(LSeed);
    free(accLB);
    free(tStampsInterpSeed);
    free(topKmotifs);
    for(ii=0;ii<lenMotifReal;ii++)
    {
        free(costMTX[ii]);
    }
    free(costMTX);
    
    t4=clock();
    myProcLogs.timeTotal += (t4-t3)/CLOCKS_PER_SEC;
    
    dumpDiscoveryLogs(logFile, myProcLogs, verbos);
    dumpParameterValuesUsed(paramOutFile, &myProcParams);
    
    return 1;
    
}

DISTTYPE manageTopKMotifs(motifInfo *topKmotifs, segInfo_t *tStamps1, segInfo_t *tStamps2, int K, INDTYPE ind1 , INDTYPE ind2, DISTTYPE dist, float blackDur)
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
         if ((tStamps1[topKmotifs[ii].ind1].str == tStamps1[ind1].str) && (fabs(tStamps2[topKmotifs[ii].ind2].str-tStamps2[ind2].str)<blackDur))
        {
            matchInd=ii;
            break;
        }
        
    }
    if (sortInd==-1)//we couldn't get any satisfactory replacement before we get a close neighbour
    {
        return topKmotifs[K-1].dist;
    }
    //There are threbsfe possibilities
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



