/******************************************************************************
*******************************************************************************
****** Sankalp Gulati                                                   *******
****** MUSIC TECHNOLOGY GROUP, UPF, BARCELONA                           *******
******                                                                  *******
*******************************************************************************
******************************************************************************/


#include "ComputePatternKNN.h"
//#define DEBUG_GENERATION



int main( int argc , char *argv[])
{
    FILE *fp, *fp1, *fp2;
    char *baseName, *patternInfoExt, *patternDataExt, *flistExt, motifFile[N_SIM_MEASURES][400]={'\0'}, *knnExt, filelistFilename[400]={'\0'}, searchFileNames[2000][400] = {'\0'}, tempFilename[400]= {'\0'}, pitchFile[400]={'\0'}, patternInfoFile[400]={'\0'}, patternDataFile[400]={'\0'};
    float t1,t2, t3,t4, temp[10], pHop, max_factor;
    int lenMotifReal, verbos=0, bandDTW, maxNMotifsPairs, nInterFact, **combMTX, searchFileID, *emptySpaceInd, emptySpaceCnt, priorityListInd, nPriorityList, *emptySpacePtr, match_found, NFilesSearch, nRead, NPatternsFile1, NPatternsFile2; 
    INDTYPE    NSeed, K,ii,jj, pp,mm, ss, ind1, ind2, lenTS1, lenTS2;
    bool sameFile;
    patternInfo_t *patternInfo1, *patternInfo2;
    
    DATATYPE **data1Interp, **data2Interp, *dataStr1, **dataPtr1, *dataStr2, **dataPtr2, **U1, **L1, **U2, **L2, *accLB;
    DISTTYPE LB_Keogh_EQ, realDist,LB_Keogh_EC,bsf=INF,**costMTX, LB_kim_FL, *bsfArray, bsf_local;
    patternDist_t **topKmotifs;
    segInfo_t *tStampsInterpDummy1, *tStampsInterpDummy2;
    procLogs_t myProcLogs;
    procParams_t myProcParams;
    fileExts_t myFileExts;
    longTermDataStorage_t **longTermDataStorage;
    INDTYPE patternID;
    segInfoInterp_t *tStampsDummy1, *tStampsDummy2;
    
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
    
    
    if(argc < 13 || argc > 14)
    {
        printf("\nInvalid number of arguments!!!\n");
        exit(1);
    }
    
    baseName = argv[1];
    myFileExts.pitchExt = argv[2];
    myFileExts.tonicExt = argv[3];
    patternInfoExt = argv[4];
    patternDataExt = argv[5];
    flistExt = argv[6];
    knnExt = argv[7];
    myProcParams.durMotif = atof(argv[8]);
    K = atoi(argv[9]);
    myProcParams.nInterpFac=atoi(argv[10]);
    myProcParams.dsFactor = atoi(argv[11]);
    if (atoi(argv[12])>0)
    {
        myProcParams.simMeasureRankRefinement = (int*)malloc(sizeof(int)*1);
        myProcParams.simMeasureRankRefinement[0] = atoi(argv[18]);
        myProcParams.nSimMeasuresUsed =1;
    }
    else    //In such case use 4 different similarity measures needed for experimentation
    {
        myProcParams.simMeasureRankRefinement = (int*)malloc(sizeof(int)*4);
        myProcParams.simMeasureRankRefinement[0] = SqEuclidean;
        myProcParams.simMeasureRankRefinement[1] = CityBlock;
        myProcParams.simMeasureRankRefinement[2] = ShiftCityBlock;
        myProcParams.simMeasureRankRefinement[3] = ShiftLinExp;
        myProcParams.nSimMeasuresUsed =4;
    }
     
    if( argc == 14 ){verbos = atoi(argv[13]);}
    
    //############ CRUCIAL PARAMETERS ##################
    myProcParams.minPossiblePitch = 60.0;
    myProcParams.binsPOct = 1200;
    myProcParams.varDur = 0.1;
    myProcParams.threshold = 45.0;
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
    for(ii=0;ii<myProcParams.nSimMeasuresUsed;ii++)
    {
        strcat(motifFile[ii],baseName);
        strcat(motifFile[ii],knnExt);
        strcat(motifFile[ii],SimMeasureNames[myProcParams.simMeasureRankRefinement[ii]]);
    }

    //search file name
    strcat(filelistFilename,baseName);
    strcat(filelistFilename,flistExt); 
    
    
    
    fp1 = fopen(filelistFilename, "r");
    ii=0;
    while(fgets(tempFilename, 400, fp1))
    {
        sscanf(tempFilename, "%[^\n]s\n", searchFileNames[ii]);
        ii++;
        
    }
    fclose(fp1);
    
    NFilesSearch = ii;
    
    //since we need to know hop size and compute motif lengths in terms of samples and since we have stored only the processed patterns (downsampled etc). We need to read one pitch file with the same processing parameters to fetch these values.
    // Opening pitch file (JUST TO OBTAIN HOP SIZE)
    //pitch file name
    strcat(pitchFile,searchFileNames[0]);
    strcat(pitchFile,myFileExts.pitchExt);
    fp1 =fopen(pitchFile,"r");
    if (fp1==NULL)
    {
        printf("Error opening file %s\n", pitchFile);
        return 0;
    }
    //reading just first two lines, in order to obtain hopsize//
    nRead = fscanf(fp1, "%f\t%f\n",&temp[0],&temp[1]);
    nRead = fscanf(fp1, "%f\t%f\n",&temp[2],&temp[3]);
    pHop = (temp[2]-temp[0])*myProcParams.dsFactor;  //final hop size afte downsampling
    fclose(fp1);
    
    // calculating lengths of the patterns
    for (ii=0;ii<myProcParams.nInterpFac; ii++)
    {
        myProcParams.motifLengths[ii] = (int)ceil((myProcParams.durMotif*myProcParams.interpFac[ii])/pHop);
        if (myProcParams.interpFac[ii]==1)
        {
            myProcParams.indexMotifLenReal = ii;
        }
        if (myProcParams.interpFac[ii]>max_factor)
        {
            max_factor = myProcParams.interpFac[ii];
            myProcParams.indexMotifLenLongest = ii;
        }
    }
    
    //CRUCIAL POINT !!! since cubic interpolation needs 4 points (2 ahead) just store 
    myProcParams.motifLengths[myProcParams.indexMotifLenLongest]+=1;    
 
    strcat(patternInfoFile,baseName);
    strcat(patternInfoFile,patternInfoExt);
    readPatternDump(patternInfoFile, &patternInfo1, &NPatternsFile1);
    
    //since we know the length of the data by now. Lets assign memory for it
    dataStr1 = (DATATYPE *)malloc(sizeof(DATATYPE)*(NPatternsFile1*myProcParams.motifLengths[myProcParams.indexMotifLenLongest]));
    dataPtr1 = (DATATYPE **)malloc(sizeof(DATATYPE *)*NPatternsFile1);
    for (ss = 0 ; ss < NPatternsFile1 ; ss++)
    {
        dataPtr1[ss] = &dataStr1[ss*myProcParams.motifLengths[myProcParams.indexMotifLenLongest]];
    }
    
    strcat(patternDataFile,baseName);
    strcat(patternDataFile,patternDataExt);
    fp1 = fopen(patternDataFile, "rb");
    fread ( dataStr1, sizeof(DATATYPE), NPatternsFile1, fp1);
    fclose(fp1);
    
    //generate multiple interpolated versions
    tStampsDummy1 = (segInfoInterp_t*)malloc(sizeof(segInfoInterp_t)*NPatternsFile1);
    generateInterpolatedSequences(dataPtr1, tStampsDummy1, &data1Interp,  &tStampsInterpDummy1, (INDTYPE)NPatternsFile1, &myProcParams);

    lenTS1 = NPatternsFile1*myProcParams.nInterpFac;
    nPriorityList = NPatternsFile1;

    // generating envelops for the seed motifs
    lenMotifReal = myProcParams.motifLengths[myProcParams.indexMotifLenReal];
    bandDTW = (int)floor(lenMotifReal*myProcParams.DTWBand);
    U1 = (DATATYPE **)malloc(sizeof(DATATYPE *)*lenTS1);
    L1= (DATATYPE **)malloc(sizeof(DATATYPE *)*lenTS1);
    accLB = (DATATYPE *)malloc(sizeof(DATATYPE)*lenMotifReal);
    topKmotifs = (patternDist_t **)malloc(sizeof(patternDist_t*)*nPriorityList);
    
    for (ii=0;ii<lenTS1;ii++)
    {
        U1[ii] = (DATATYPE *)malloc(sizeof(DATATYPE)*lenMotifReal);
        L1[ii] = (DATATYPE *)malloc(sizeof(DATATYPE)*lenMotifReal);
        computeRunningMinMax(data1Interp[ii], U1[ii], L1[ii], lenMotifReal, bandDTW);   
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

    for(jj=0;jj<nPriorityList;jj++)
    {
        topKmotifs[jj] = (patternDist_t *)malloc(sizeof(patternDist_t)*K);
        for(ii=0;ii<K;ii++)
        {
            topKmotifs[jj][ii].dist = INF;
            topKmotifs[jj][ii].patternID1 = PID_DEFAULT1;
            topKmotifs[jj][ii].patternID1 = PID_DEFAULT2;
        }
    }
    bsfArray = (DISTTYPE*)malloc(sizeof(DISTTYPE)*nPriorityList);
    for(ii=0;ii<nPriorityList;ii++)
    {
        bsfArray[ii] = bsf;
    }

    memset(tempFilename, '\0', sizeof(char)*400);
    for (ss=0;ss<NFilesSearch; ss++)
    {
        //read the data 
        strcat(tempFilename,searchFileNames[ss]);
        strcat(tempFilename,patternInfoExt);
        readPatternDump(tempFilename, &patternInfo2, &NPatternsFile2);
        memset(tempFilename, '\0', sizeof(char)*400);
        
        //since we know the length of the data by now. Lets assign memory for it
        dataStr2 = (DATATYPE *)malloc(sizeof(DATATYPE)*(NPatternsFile2*myProcParams.motifLengths[myProcParams.indexMotifLenLongest]));
        dataPtr2 = (DATATYPE **)malloc(sizeof(DATATYPE *)*NPatternsFile2);
        for (ii = 0 ; ii < NPatternsFile2 ; ii++)
        {
            dataPtr2[ii] = &dataStr2[ii*myProcParams.motifLengths[myProcParams.indexMotifLenLongest]];
        }
        
        strcat(tempFilename,searchFileNames[ss]);
        strcat(tempFilename,patternDataExt);
        fp1 = fopen(tempFilename, "rb");
        fread ( dataStr2, sizeof(DATATYPE), NPatternsFile2, fp1);
        fclose(fp1);
        memset(tempFilename, '\0', sizeof(char)*400);
        
        //generate multiple interpolated versions
        tStampsDummy2 = (segInfoInterp_t*)malloc(sizeof(segInfoInterp_t)*NPatternsFile2);
        generateInterpolatedSequences(dataPtr2, tStampsDummy2, &data2Interp,  &tStampsInterpDummy2, (INDTYPE) NPatternsFile2, &myProcParams);
        
        lenTS2 = NPatternsFile2*myProcParams.nInterpFac;

        t1=clock();
        //computing envelops for the file to be searched
        U2 = (DATATYPE **)malloc(sizeof(DATATYPE *)*lenTS2);
        L2 = (DATATYPE **)malloc(sizeof(DATATYPE *)*lenTS2);
        
        for (ii=0;ii<lenTS2;ii++)
        {
            U2[ii] = (DATATYPE *)malloc(sizeof(DATATYPE)*lenMotifReal);
            L2[ii] = (DATATYPE *)malloc(sizeof(DATATYPE)*lenMotifReal);
            computeRunningMinMax(data2Interp[ii], U2[ii], L2[ii], lenMotifReal, bandDTW);
            
        }
        t2=clock();
        myProcLogs.timeGenEnvelops += (t2-t1)/CLOCKS_PER_SEC;

        for (ii=0; ii < NPatternsFile1; ii++)
        {
            for (jj=0; jj < NPatternsFile2; jj++)  
            {
                ind1 = ii*myProcParams.nInterpFac;
                ind2 = jj*myProcParams.nInterpFac;

                bsf_local = bsfArray[ii];
                for (pp = 0; pp < myProcParams.nInterpFac; pp++)
                {
                    for (mm = 0; mm < myProcParams.nInterpFac; mm++)
                    {
                        // skip if its not a valid combinations of the interpolation factor (for optimizing and doing only 9 out of 25 combinations)
                        if (myProcParams.combMTX[pp][mm]==0)
                        {
                            continue;
                        }

                        LB_kim_FL = computeLBkimFL(data1Interp[ind1+pp][0], data2Interp[ind2+mm][0], data1Interp[ind1+pp][lenMotifReal-1], data2Interp[ind2+mm][lenMotifReal-1], SqEuclidean);
                        myProcLogs.totalFLDone++;
                        if (LB_kim_FL< bsf_local)
                        {
                            LB_Keogh_EQ = computeKeoghsLB(U1[ind1+pp],L1[ind1+pp], accLB, data2Interp[ind2+mm],lenMotifReal, bsf_local, SqEuclidean);
                            myProcLogs.totalLBKeoghEQ++;
                            if(LB_Keogh_EQ < bsf_local)
                            {
                                LB_Keogh_EC = computeKeoghsLB(U2[ind2+mm],L2[ind2+mm],accLB, data1Interp[ind1+pp],lenMotifReal, bsf_local, SqEuclidean);
                                myProcLogs.totalLBKeoghEC++;
                                if(LB_Keogh_EC < bsf_local)
                                {
                                    realDist = dtw1dBandConst(data1Interp[ind1+pp], data2Interp[ind2+mm], lenMotifReal, lenMotifReal, costMTX, SqEuclidean, bandDTW, bsf_local, accLB);
                                    myProcLogs.totalDTWComputations++;
                                    if(realDist < bsf_local)
                                    {
                                        bsf_local = realDist;
                                    }
                                }
                            }
                        }

                    }   
                }

                if (bsf_local < bsfArray[ii])
                {
                    bsfArray[ii] = manageTopKMotifs(topKmotifs[ii], nPriorityList, patternInfo1[ii].id, patternInfo2[jj].id, bsf_local);
                    myProcLogs.totalPriorityUpdates++;
                }


            }
        }

        for(ii=0;ii<(lenTS2);ii++)
        {   
            free(U2[ii]);
            free(L2[ii]);
        }
        for(ii=0;ii<(NPatternsFile2);ii++)
        {   
            for (jj=0; jj < myProcParams.nInterpFac; jj++)
            {
                if (myProcParams.interpFac[jj] !=1)
                {
                    free(data2Interp[ii]);
                }
                
            }
        }


        free(U2);
        free(L2);
        free(data2Interp);
        free(dataStr2);
        free(dataPtr2);
        free(tStampsDummy2);
        free(tStampsInterpDummy2);
    
    }

    for(ii=0;ii<(lenTS1);ii++)
    {   
        free(U1[ii]);
        free(L1[ii]);
    }
    for(ii=0;ii<(NPatternsFile1);ii++)
    {   
        for (jj=0; jj < myProcParams.nInterpFac; jj++)
        {
            if (myProcParams.interpFac[jj] ==1)
            {
                continue;
            }
            free(data1Interp[ii]);    
        }
    }


    free(U1);
    free(L1);
    free(data1Interp);
    free(dataStr1);
    free(dataPtr1);
    free(tStampsDummy1);
    free(tStampsInterpDummy1);

    
        


}




DISTTYPE manageTopKMotifs(patternDist_t *topKmotifs, int K, INDTYPE id1, INDTYPE id2, DISTTYPE dist)
{
    int ii=0;
    int sortInd = -1;
    int matchInd = -1;
    int a=0;
    
    
    for(ii=0;ii<K;ii++)
    {
        if ((topKmotifs[ii].dist > dist)&&(sortInd==-1))
        {
            sortInd=ii;
        }

        // searching if we already have a motif in out top K list which is near to the currently good match
        if ((topKmotifs[ii].patternID1 == id1) && (topKmotifs[ii].patternID2 == id2))
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
        memmove(&topKmotifs[sortInd+1], &topKmotifs[sortInd], sizeof(patternDist_t)*(K-(sortInd+1)));
        topKmotifs[sortInd].dist = dist;
        topKmotifs[sortInd].patternID1 = id1;
        topKmotifs[sortInd].patternID2 = id2;
    }
    else if (sortInd == matchInd)
    {
        topKmotifs[sortInd].dist = dist;
        topKmotifs[sortInd].patternID1 = id1;
        topKmotifs[sortInd].patternID2 = id2;
    }
    else if (sortInd < matchInd)
    {
        memmove(&topKmotifs[sortInd+1], &topKmotifs[sortInd], sizeof(patternDist_t)*(matchInd-sortInd));
        topKmotifs[sortInd].dist = dist;
        topKmotifs[sortInd].patternID1 = id1;
        topKmotifs[sortInd].patternID2 = id2;
    }
    
    return topKmotifs[K-1].dist;
    
}
    
    
    
    
    
    
    
    
    
    
    
    
    
    
/*    
     
        t2=clock();
        myProcLogs.timeDiscovery += (t2-t1)/CLOCKS_PER_SEC;

        for(ii=0;ii<lenTS;ii++)
        {
            free(dataInterp[ii]);
            free(U[ii]);
            free(L[ii]);
        }

        free(U);
        free(L);
        free(dataInterp);
        free(tStampsInterp);
        
        //in the mapp file make a mapping of searchFileID and its name
        fprintf(fp2,"%d\t%s\n", searchFileID, searchFile);
        searchFileID++;
            
        
    }
    fclose(fp);
    fclose(fp2);

    
    //############## Rank Refinement using sophisticated DTW ########################
    for (mm=0;mm<myProcParams.nSimMeasuresUsed;mm++)
    {
            // Since there can be multiple similarity measure used for rank refinement (mainly during experiment phase) this rank refinement step should be in loop, no need to loop rest of the steps
            
            //recomputing the distance
            for(ii=0;ii<nPriorityList;ii++)
            {
                for(jj=0;jj<K;jj++)
                {
                    if (topKmotifs[ii][jj].dist <INF)   //do refinement only for a valid top entry, leave the infinites!!
                    {
                        topKmotifs[ii][jj].dist = dtw1dBandConst_localConst(dataInterpSeed[topKmotifs[ii][jj].ind1], topKmotifs[ii][jj].storagePtr->data, lenMotifReal, lenMotifReal, costMTX, myProcParams.simMeasureRankRefinement[mm], bandDTW, INF, accLB);
                    }
                    else
                    {
                        if (verbos)
                        {
                            printf("There is some serious problem in rank refinement step %lld,%lld",jj,ii);
                        }
                    }
                    
                }
                //sorting the priority list
                qsort (topKmotifs[ii], K, sizeof(motifInfo), compareMotifInfo);
            }   
                
            
            t1=clock();
            dumpSearchMotifInfo(motifFile[mm], topKmotifs, tStampsInterpSeed, NSeed, K, nInterFact, verbos);
            t2=clock();
            myProcLogs.timeWriteData += (t2-t1)/CLOCKS_PER_SEC;
    }

    
    //Memory clearing
    for(ii=0;ii<NSeed;ii++)
    {
        free(dataInterpSeed[ii]);
        free(USeed[ii]);
        free(LSeed[ii]);
        
    }
    for(jj=0;jj<nPriorityList;jj++)
    {
        free(topKmotifs[jj]);
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
    
    for(ii=0;ii<nInterFact;ii++)
    {
        free(combMTX[ii]);
    }
    for(jj=0;jj<nPriorityList;jj++)
    {
        for(ii=0;ii<K;ii++)
        {
            free(longTermDataStorage[jj][ii].data);
        }
        
        free(longTermDataStorage[jj]);
        
    }
    free(longTermDataStorage);
    free(emptySpaceInd);
    free(combMTX);
    free(myProcParams.simMeasureRankRefinement);
    free(bsfArray);
    
    t4=clock();
    myProcLogs.timeTotal += (t4-t3)/CLOCKS_PER_SEC;
    
    dumpDiscoveryLogs(logFile, myProcLogs, verbos);
    dumpParameterValuesUsed(paramOutFile, &myProcParams);

    */






