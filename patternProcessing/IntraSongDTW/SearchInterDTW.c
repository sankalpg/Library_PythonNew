/******************************************************************************
*******************************************************************************
****** Sankalp Gulati                                                   *******
****** MUSIC TECHNOLOGY GROUP, UPF, BARCELONA                           *******
******                                                                  *******
*******************************************************************************
******************************************************************************/


#include "SearchInterDTW.h"
#define DEBUG_GENERATION



int main( int argc , char *argv[])
{
    FILE *fp, *fp_out;
    char *baseName, motifFile[400]={'\0'}, logFile[400]={'\0'}, searchFileList[400]={'\0'}, searchFile[400] = {'\0'}, mappFile[400] = {'\0'} ;
    float t1,t2, t3,t4;
    int lenMotifReal, verbos=0, bandDTW, maxNMotifsPairs; 
    INDTYPE    NSeed, lenTS, count_DTW=0, numLinesInFile, K,ii,jj;
    
    DATATYPE **dataInterpSeed, **dataInterp, **U, **L, **USeed, **LSeed, *accLB;
    DISTTYPE LB_Keogh_EQ, realDist,LB_Keogh_EC,bsf=INF,**costMTX, LB_kim_FL;
    motifInfo **topKmotifs;
    segInfo_t *taniSegs, *tStampsInterp, *tStampsInterpSeed;
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
    
    
    if(argc < 16 || argc > 17)
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
    myProcParams.durMotif = atof(argv[10]);
    K = atoi(argv[11]);
    myProcParams.blackDur = atof(argv[12]);
    if (atof(argv[13])>0)
    {
        bsf = atof(argv[13]);
    }
     myProcParams.dsFactor = atoi(argv[14]);
     maxNMotifsPairs = atoi(argv[15]);
    
    if( argc == 17 ){verbos = atoi(argv[16]);}
    
    //############ CRUCIAL PARAMETERS ##################
    myProcParams.minPossiblePitch = 60.0;
    myProcParams.allowedSilDur = 0.15;
    myProcParams.binsPOct = 120;
    myProcParams.varDur = 0.1;
    myProcParams.threshold = 225;
    myProcParams.flatThreshold = 0.8;
    myProcParams.maxPauseDur = 0.5;
    myProcParams.factorLow = 0.9;
    myProcParams.factorHigh = 1.1;
    myProcParams.DTWBand = 0.1;
    myProcParams.removeTaniSegs=1;
    
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
    //return 1;
    
#endif      
    
    // generating envelops for the seed motifs
    bandDTW = (int)floor(lenMotifReal*myProcParams.DTWBand);
    USeed = (DATATYPE **)malloc(sizeof(DATATYPE *)*NSeed);
    LSeed= (DATATYPE **)malloc(sizeof(DATATYPE *)*NSeed);
    accLB = (DATATYPE *)malloc(sizeof(DATATYPE)*lenMotifReal);
    topKmotifs = (motifInfo **)malloc(sizeof(motifInfo*)*NSeed/3);
    
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
        for(jj=0;jj<NSeed/3;jj++)
        {
            topKmotifs[jj] = (motifInfo *)malloc(sizeof(motifInfo)*lenTS);
            for(ii=0;ii<lenTS;ii++)
            {
                topKmotifs[jj][ii].dist = INF;
                topKmotifs[jj][ii].ind1 = 0;
                topKmotifs[jj][ii].ind2 = 0;
            }
        }
        t1=clock();
        if (strcmp(baseName, searchFile))
        {
            for(ii=0;ii<NSeed;ii++)
            {
                for(jj=0;jj<lenTS;jj++)
                {
                    if (((ii%3==0)&&(jj%3==0))||((ii%3==2)&&(jj%3==2))||((ii%3==0)&&(jj%3==1))||((ii%3==2)&&(jj%3==1)))
                        continue;
                    if (fabs(tStampsInterpSeed[ii].str-tStampsInterp[jj].str)< myProcParams.blackDur)
                    {
                        continue;
                    }
                    
                    
                    LB_kim_FL = computeLBkimFL(dataInterpSeed[ii][0], dataInterp[jj][0], dataInterpSeed[ii][lenMotifReal-1], dataInterp[jj][lenMotifReal-1]);
                    myProcLogs.totalFLDone++;
                    if (LB_kim_FL< bsf)
                    {
                        LB_Keogh_EQ = computeKeoghsLB(USeed[ii],LSeed[ii], accLB, dataInterp[jj],lenMotifReal, bsf);
                        myProcLogs.totalLBKeoghEQ++;
                        if(LB_Keogh_EQ < bsf)
                        {
                            LB_Keogh_EC = computeKeoghsLB(U[jj],L[jj],accLB, dataInterpSeed[ii],lenMotifReal, bsf);
                            myProcLogs.totalLBKeoghEC++;
                            if(LB_Keogh_EC < bsf)
                            {
                                realDist = dtw1dBandConst(dataInterpSeed[ii], dataInterp[jj], lenMotifReal, lenMotifReal, costMTX, 0, bandDTW, bsf, accLB);
                                myProcLogs.totalDTWComputations++;
                                if(realDist<bsf)
                                {
                                    realDist = dtw1dBandConst_localConst(dataInterpSeed[ii], dataInterp[jj], lenMotifReal, lenMotifReal, costMTX, 0, bandDTW, bsf, accLB);
                                    manageTopKMotifs(topKmotifs[ii/3], tStampsInterpSeed, tStampsInterp, lenTS, ii, jj, realDist, myProcParams.blackDur);
                                    myProcLogs.totalPriorityUpdates++;
                                }
                            }
                        }
                    }
                }
            } 
            
        }
        else
        {
            for(ii=0;ii<NSeed;ii++)
            {
                for(jj=0;jj<lenTS;jj++)
                {
                    if (((ii%3==0)&&(jj%3==0))||((ii%3==2)&&(jj%3==2))||((ii%3==0)&&(jj%3==1))||((ii%3==2)&&(jj%3==1)))
                        continue;
                
                    LB_kim_FL = computeLBkimFL(dataInterpSeed[ii][0], dataInterp[jj][0], dataInterpSeed[ii][lenMotifReal-1], dataInterp[jj][lenMotifReal-1]);
                    myProcLogs.totalFLDone++;
                    if (LB_kim_FL< bsf)
                    {
                        LB_Keogh_EQ = computeKeoghsLB(USeed[ii],LSeed[ii], accLB, dataInterp[jj],lenMotifReal, bsf);
                        myProcLogs.totalLBKeoghEQ++;
                        if(LB_Keogh_EQ < bsf)
                        {
                            LB_Keogh_EC = computeKeoghsLB(U[jj],L[jj],accLB, dataInterpSeed[ii],lenMotifReal, bsf);
                            myProcLogs.totalLBKeoghEC++;
                            if(LB_Keogh_EC < bsf)
                            {
                                realDist = dtw1dBandConst(dataInterpSeed[ii], dataInterp[jj], lenMotifReal, lenMotifReal, costMTX, 0, bandDTW, bsf, accLB);
                                myProcLogs.totalDTWComputations++;
                                if(realDist<bsf)
                                {
                                    realDist = dtw1dBandConst_localConst(dataInterpSeed[ii], dataInterp[jj], lenMotifReal, lenMotifReal, costMTX, 0, bandDTW, bsf, accLB);
                                    manageTopKMotifs(topKmotifs[ii/3], tStampsInterpSeed, tStampsInterp, lenTS, ii, jj, realDist, myProcParams.blackDur);
                                    myProcLogs.totalPriorityUpdates++;
                                }
                            }
                        }
                    }
                }
            }            
        }
    t2=clock();
    myProcLogs.timeDiscovery += (t2-t1)/CLOCKS_PER_SEC;
    t1=clock();
    dumpSearchMotifInfo(motifFile, mappFile, searchFile, topKmotifs, tStampsInterpSeed, tStampsInterp, NSeed, lenTS, &mapp, verbos);
    t2=clock();
    myProcLogs.timeWriteData += (t2-t1)/CLOCKS_PER_SEC;
    
    for(ii=0;ii<lenTS;ii++)
    {
        free(dataInterp[ii]);
        free(U[ii]);
        free(L[ii]);
    }
    for(jj=0;jj<NSeed/3;jj++)
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



