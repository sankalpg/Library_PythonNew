/******************************************************************************
*******************************************************************************
****** Sankalp Gulati                                                   *******
****** MUSIC TECHNOLOGY GROUP, UPF, BARCELONA                           *******
******                                                                  *******
*******************************************************************************
******************************************************************************/


#include "DiscoverIntraDTW.h"
#define DEBUG_GENERATION



int main( int argc , char *argv[])
{
    FILE *fp;
    char *baseName, motifFile[400]={'\0'}, logFile[400]={'\0'};
    float t1,t2;
    int lenMotifReal, verbos=0, bandDTW; 
    INDTYPE    lenTS, count_DTW=0, numLinesInFile, K,ii,jj;
    
    DATATYPE **dataInterp, **U, **L, *accLB;
    DISTTYPE LB_Keogh_EQ, realDist,LB_Keogh_EC,bsf=INF,**costMTX, LB_kim_FL;
    motifInfo *topKmotifs;
    segInfo_t *taniSegs, *tStampsInterp;
    procLogs_t myProcLogs;
    procParams_t myProcParams;
    fileExts_t myFileExts;
    
    
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
    
    if(argc < 12 || argc > 13)
    {
        printf("\nInvalid number of arguments!!!\n");
        exit(1);
    }
    
    baseName = argv[1];
    myFileExts.pitchExt = argv[2];
    myFileExts.tonicExt = argv[3];
    myFileExts.segExt = argv[4];
    myFileExts.motifExt = argv[5];
    myFileExts.logExt = argv[6];
    myProcParams.durMotif = atof(argv[7]);
    K = atoi(argv[8]);
    myProcParams.blackDur = atof(argv[9]);
    if (atof(argv[10])>0)
    {
        bsf = atof(argv[10]);
    }
     myProcParams.dsFactor = atoi(argv[11]);
    
    if( argc == 13 ){verbos = atoi(argv[12]);}
    
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
    //log file name
    strcat(logFile,baseName);
    strcat(logFile,myFileExts.logExt);  
    
    
    lenTS = readPreProcessGenDB(&dataInterp, &tStampsInterp, &lenMotifReal, baseName, &myFileExts, &myProcParams, &myProcLogs, verbos);
    
#ifdef DEBUG_GENERATION
    fp = fopen("subsequences.bin","wb");
    for(ii=0;ii<lenTS;ii++)
    {
        for(jj=0;jj<lenMotifReal;jj++)
        {
            fprintf(fp, "%f\t", dataInterp[ii][jj]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
    return 1;
    
#endif    
    //################# Precomputing envelope of each subsequence for the LB Keogh lower bound ###########################
    bandDTW = (int)floor(lenMotifReal*myProcParams.DTWBand);
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
            if (fabs(tStampsInterp[ii].str-tStampsInterp[jj].str)< myProcParams.blackDur)
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
                            bsf = manageTopKMotifs(topKmotifs, tStampsInterp, K, ii, jj, realDist, myProcParams.blackDur);
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
    

    dumpDiscoveredMotifInfo(motifFile, topKmotifs, tStampsInterp, K, verbos);
    dumpDiscoveryLogs(logFile, myProcLogs, verbos);
    
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



DISTTYPE manageTopKMotifs(motifInfo *topKmotifs, segInfo_t *tStamps, int K, INDTYPE ind1 , INDTYPE ind2, DISTTYPE dist, float blackDur)
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


