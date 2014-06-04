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
    FILE *fp, *fp2;
    char *baseName, *patternInfoExt, *patternDataExt, *flistExt, motifFile[N_SIM_MEASURES][400]={'\0'}, *knnExt, filelistFilename[400]={'\0'}, searchFileNames[2000][400] = {'\0'}, tempFilename[400]= {'\0'};
    float t1,t2, t3,t4, temp[10], pHop;
    int lenMotifReal, verbos=0, bandDTW, maxNMotifsPairs, nInterFact, **combMTX, mm, searchFileID, *emptySpaceInd, emptySpaceCnt, priorityListInd, nPriorityList, *emptySpacePtr, match_found, NFilesSearch; 
    INDTYPE    NSeed, lenTS, K,ii,jj, pp, ss, NPatternsFile1, NPatternsFile2;
    bool sameFile;
    patternInfo_t *patternInfo1, *patternInfo2;
    
    DATATYPE **data1Interp, **data2Interp, *dataStr1, **dataPtr1, *dataStr2, **dataPtr2, **U1, **L1, **U2, **L2, *accLB;
    DISTTYPE LB_Keogh_EQ, realDist,LB_Keogh_EC,bsf=INF,**costMTX, LB_kim_FL, *bsfArray;
    motifInfo **topKmotifs;
    segInfo_t *tStampsInterpDummy;
    procLogs_t myProcLogs;
    procParams_t myProcParams;
    fileExts_t myFileExts;
    longTermDataStorage_t **longTermDataStorage;
    INDTYPE patternID;
    segInfoInterp_t *tStampsDummy;
    
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
    strcat(pitchFile,myFileExts->pitchExt);
    fp1 =fopen(pitchFile,"r");
    if (fp1==NULL)
    {
        printf("Error opening file %s\n", pitchFile);
        return 0;(myProcParams->motifLengths[myProcParams->indexMotifLenLongest]
    }
    //reading just first two lines, in order to obtain hopsize//
    nRead = fscanf(fp1, "%f\t%f\n",&temp[0],&temp[1]);
    nRead = fscanf(fp1, "%f\t%f\n",&temp[2],&temp[3]);
    pHop = (temp[2]-temp[0])*myProcParams->dsFactor;  //final hop size afte downsampling
    fclose(fp1);
    
    // calculating lengths of the patterns
    for (ii=0;ii<myProcParams->nInterpFac; ii++)
    {
        myProcParams->motifLengths[ii] = (int)ceil((myProcParams->durMotif*myProcParams->interpFac[ii])/pHop);
        if (myProcParams->interpFac[ii]==1)
        {
            myProcParams->indexMotifLenReal = ii;
        }
        if (myProcParams->interpFac[ii]>max_factor)
        {
            max_factor = myProcParams->interpFac[ii];
            myProcParams->indexMotifLenLongest = ii;
        }
    }
    
    //CRUCIAL POINT !!! since cubic interpolation needs 4 points (2 ahead) just store 
    myProcParams->motifLengths[myProcParams->indexMotifLenLongest]+=1;    
    
    
    for (ii=0;ii<NFilesSearch;ii++)
    {
        //read the data 
        strcat(tempFilename,searchFileNames[ii]);
        strcat(tempFilename,patternInfoExt);
        readPatternDump(tempFilename, &patternInfo1, &NPatternsFile1);
        
        //since we know the length of the data by now. Lets assign memory for it
        dataStr1 = (DATATYPE *)malloc(sizeof(DATATYPE)*(NPatternsFile1*myProcParams->motifLengths[myProcParams->indexMotifLenLongest]))
        dataPtr1 = (DATATYPE **)malloc(sizeof(DATATYPE *)*NPatternsFile1);
        for (ss = 0 ; ss < NPatternsFile1 ; ss++)
        {
            dataPtr1[ss] = &dataStr1[ss*myProcParams->motifLengths[myProcParams->indexMotifLenLongest]];
        }
        
        strcat(tempFilename,searchFileNames[ii]);
        strcat(tempFilename,patternDataExt);
        fp1 = fopen(tempFilename, "rb");
        fread ( dataStr1, sizeof(DATATYPE), NPatternsFile1, fp1);
        fclose(fp1);
        
        //generate multiple interpolated versions
        tStampsDummy = (segInfoInterp_t*)malloc(sizeof(segInfoInterp_t)*NPatternsFile1);
        generateInterpolatedSequences(dataPtr1, tStampsDummy, &data1Interp,  &tStampsInterpDummy, NPatternsFile1, myProcParams);
        
        for (jj=ii;jj<NFilesSearch; jj++)
        {
            //read the data 
            strcat(tempFilename,searchFileNames[jj]);
            strcat(tempFilename,patternInfoExt);
            readPatternDump(tempFilename, &patternInfo1, &NPatternsFile1);
            
            //since we know the length of the data by now. Lets assign memory for it
            dataStr1 = (DATATYPE *)malloc(sizeof(DATATYPE)*(NPatternsFile1*myProcParams->motifLengths[myProcParams->indexMotifLenLongest]))
            dataPtr1 = (DATATYPE **)malloc(sizeof(DATATYPE *)*NPatternsFile1);
            for (ss = 0 ; ss < NPatternsFile1 ; ss++)
            {
                dataPtr1[ss] = &dataStr1[ss*myProcParams->motifLengths[myProcParams->indexMotifLenLongest]];
            }
            
            strcat(tempFilename,searchFileNames[jj]);
            strcat(tempFilename,patternDataExt);
            fp1 = fopen(tempFilename, "rb");
            fread ( dataStr1, sizeof(DATATYPE), NPatternsFile1, fp1);
            fclose(fp1);
            
            //generate multiple interpolated versions
            tStampsDummy = (segInfoInterp_t*)malloc(sizeof(segInfoInterp_t)*NPatternsFile1);
            generateInterpolatedSequences(dataPtr1, tStampsDummy, &data1Interp,  &tStampsInterpDummy, NPatternsFile1, myProcParams);
        
        }
    }


    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

    nPriorityList = (int)floor(NSeed/nInterFact);
    
#ifdef DEBUG_GENERATION
    fp = fopen("subsequences.bin","wb");
    for(ii=0;ii<NSeed;ii++)
    {
        for(jj=0;jj<lenMotifReal;jj++)
        {
            fprintf(fp, "%f\t", dataInterpSeed[ii][jj]);
        }
        fprintf(fp, "\n");
    }read
    fclose(fp);
    return 1;
    
#endif      
    
    // generating envelops for the seed motifs
    bandDTW = (int)floor(lenMotifReal*myProcParams.DTWBand);
    USeed = (DATATYPE **)malloc(sizeof(DATATYPE *)*NSeed);
    LSeed= (DATATYPE **)malloc(sizeof(DATATYPE *)*NSeed);
    accLB = (DATATYPE *)malloc(sizeof(DATATYPE)*lenMotifReal);
    topKmotifs = (motifInfo **)malloc(sizeof(motifInfo*)*nPriorityList);
    
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

    for(jj=0;jj<nPriorityList;jj++)
    {
        topKmotifs[jj] = (motifInfo *)malloc(sizeof(motifInfo)*K);
        for(ii=0;ii<K;ii++)
        {
            topKmotifs[jj][ii].dist = INF;
            topKmotifs[jj][ii].ind1 = 0;
            topKmotifs[jj][ii].ind2 = 0;
            topKmotifs[jj][ii].patternID = PID_DEFAULT2;
            topKmotifs[jj][ii].searchFileID = FID_DEFAULT1;
            
            
        }
    }
    bsfArray = (DISTTYPE*)malloc(sizeof(DISTTYPE)*nPriorityList);
    for(ii=0;ii<nPriorityList;ii++)
    {
        bsfArray[ii] = bsf;
    }
    
    longTermDataStorage = (longTermDataStorage_t **)malloc(sizeof(longTermDataStorage_t*)*nPriorityList);
    for(jj=0;jj<nPriorityList;jj++)
    {
        longTermDataStorage[jj] = (longTermDataStorage_t *)malloc(sizeof(longTermDataStorage_t)*K);
        
        for(ii=0;ii<K;ii++)
        {
            longTermDataStorage[jj][ii].data = (DATATYPE *)malloc(sizeof(DATATYPE)*lenMotifReal);
            longTermDataStorage[jj][ii].patternID = PID_DEFAULT1;
        }
        
    }
    
    emptySpaceInd = (int *)malloc(sizeof(int*)*K);
    
    // iterating over all files 
    fp2 = fopen(mappFile,"w");
    searchFileID=0;
    
    fp = fopen(searchFileList, "r");
    patternID = 0;
    while(fgets(tempFilename, 400, fp))
    {
        sscanf(tempFilename, "%[^\n]s\n", searchFile);
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
        
        //############## Performing a search using basic DTW ########################
        t1=clock();
        for(jj=0;jj<lenTS;jj++)
        {
            for(ii=0;ii<NSeed;ii++)
            {
                priorityListInd = (int)floor(ii/nInterFact);
        
                if (myProcParams.combMTX[ii%nInterFact][jj%nInterFact]==0)
                {
                    continue;
                }
                if ((strcmp(baseName, searchFile)==0)&&(fabs(tStampsInterpSeed[ii].str-tStampsInterp[jj].str)< myProcParams.blackDur))
                    //beware that basename and searchFile name should both have either full path or relative path.
                {
                    continue;
                }
                
                LB_kim_FL = computeLBkimFL(dataInterpSeed[ii][0], dataInterp[jj][0], dataInterpSeed[ii][lenMotifReal-1], dataInterp[jj][lenMotifReal-1], SqEuclidean);
                myProcLogs.totalFLDone++;
                if (LB_kim_FL< bsfArray[priorityListInd])
                {
                    LB_Keogh_EQ = computeKeoghsLB(USeed[ii],LSeed[ii], accLB, dataInterp[jj],lenMotifReal, bsfArray[priorityListInd], SqEuclidean);
                    myProcLogs.totalLBKeoghEQ++;
                    if(LB_Keogh_EQ < bsfArray[priorityListInd])
                    {
                        LB_Keogh_EC = computeKeoghsLB(U[jj],L[jj],accLB, dataInterpSeed[ii],lenMotifReal, bsfArray[priorityListInd], SqEuclidean);
                        myProcLogs.totalLBKeoghEC++;
                        if(LB_Keogh_EC < bsfArray[priorityListInd])
                        {
                            realDist = dtw1dBandConst(dataInterpSeed[ii], dataInterp[jj], lenMotifReal, lenMotifReal, costMTX, SqEuclidean, bandDTW, bsfArray[priorityListInd], accLB);
                            myProcLogs.totalDTWComputations++;
                            if(realDist<bsfArray[priorityListInd])
                            {
                                bsfArray[priorityListInd] = manageTopKMotifs(topKmotifs[priorityListInd], tStampsInterpSeed, tStampsInterp, K, ii, jj, realDist, myProcParams.blackDur, searchFileID);
                                myProcLogs.totalPriorityUpdates++;
                            }
                        }
                    }
                }
            }
            
        }
        for(ii=0;ii<NSeed;ii++)
        {
            priorityListInd = (int)floor(ii/nInterFact);
            manageTopKMotifsData(topKmotifs[priorityListInd], longTermDataStorage[priorityListInd], dataInterp, tStampsInterp, &patternID, emptySpaceInd, lenMotifReal, K, searchFileID);
            
        }
        
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
    
    return 1;
    
}

/*
 * This function manages the data for top K motifs so that we can perform rank refinement at the end. The logic is simple but steps may seesm a bit complex.
 */
int manageTopKMotifsData(motifInfo *topKmotifs, longTermDataStorage_t *longTermDataStorage, DATATYPE** dataInterp, segInfo_t *tStampsInterp, INDTYPE *patternID, int *emptySpaceInd, int lenMotifReal, int K, int searchFileID)
{
    int match_found, *emptySpacePtr;
    int emptySpaceCnt=0;
    int pp,ss;
    
    //Lets find out what patterns are removed from the priority list after processing one seed over entire search file. if patterns are removed, store there location is emptySapceInd buffer which we later use to copy new data (data which is newly added to priority list)
    match_found=0;
    for(pp=0;pp<K;pp++)
    {
        match_found=0;
        for (ss=0;ss<K;ss++)
        {
            if (longTermDataStorage[pp].patternID == topKmotifs[ss].patternID)
            {
                match_found=1;
                break;
            }
        }
        if(match_found==0)
        {
            emptySpaceInd[emptySpaceCnt] = pp;
            emptySpaceCnt++;
        }
    }
    
    emptySpacePtr = emptySpaceInd;
    //After processing every seed motif over a file, check what are the new patterns added in the priority list and store data corresponding to them in the longTermStorage structure and assign them pattern ID so that we can track these patterns in the future
    for(pp=0;pp<K;pp++)
    {
        //newly added patterns will be the ones added from the current search file so just search only in that domain
        if ((topKmotifs[pp].searchFileID == searchFileID)&&(topKmotifs[pp].patternID==PID_DEFAULT3))
        {
            topKmotifs[pp].patternID = *patternID;
            (*patternID)++;
            
            //Also in such case add data to the long term storage;
            if(emptySpaceCnt>0)
            {
                memcpy(longTermDataStorage[emptySpacePtr[0]].data, dataInterp[topKmotifs[pp].ind2], sizeof(DATATYPE)*lenMotifReal);
                longTermDataStorage[emptySpacePtr[0]].patternID = topKmotifs[pp].patternID;
                longTermDataStorage[emptySpacePtr[0]].strTime = tStampsInterp[topKmotifs[pp].ind2].str;
                longTermDataStorage[emptySpacePtr[0]].endTime = tStampsInterp[topKmotifs[pp].ind2].end;
                topKmotifs[pp].storagePtr = &longTermDataStorage[emptySpacePtr[0]];
                emptySpacePtr++;
                emptySpaceCnt--;
            }
            else
            {
                printf("SOMETHING TERRIBLE HAPPENED");
                return -1;
            }
        }
    }
    
    return 1;
}



DISTTYPE manageTopKMotifs(motifInfo *topKmotifs, segInfo_t *tStamps1, segInfo_t *tStamps2, int K, INDTYPE ind1 , INDTYPE ind2, DISTTYPE dist, float blackDur, int searchFileID)
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
        if ((topKmotifs[ii].searchFileID==searchFileID) && (fabs(tStamps2[topKmotifs[ii].ind2].str-tStamps2[ind2].str)<blackDur))
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
        topKmotifs[sortInd].searchFileID = searchFileID;
        topKmotifs[sortInd].patternID = PID_DEFAULT3;
    }
    else if (sortInd == matchInd)
    {
        topKmotifs[sortInd].dist = dist;
        topKmotifs[sortInd].ind1 = ind1;
        topKmotifs[sortInd].ind2 = ind2;
        topKmotifs[sortInd].searchFileID = searchFileID;
        topKmotifs[sortInd].patternID = PID_DEFAULT3;
    }
    else if (sortInd < matchInd)
    {
        memmove(&topKmotifs[sortInd+1], &topKmotifs[sortInd], sizeof(motifInfo)*(matchInd-sortInd));
        topKmotifs[sortInd].dist = dist;
        topKmotifs[sortInd].ind1 = ind1;
        topKmotifs[sortInd].ind2 = ind2;
        topKmotifs[sortInd].searchFileID = searchFileID;
        topKmotifs[sortInd].patternID = PID_DEFAULT3;
    }
    
    return topKmotifs[K-1].dist;
    
}




