
/*******************************************************************************
****** Sankalp Gulati                                                   *******
****** MUSIC TECHNOLOGY GROUP, UPF, BARCELONA                           *******
******                                                                  *******
*******************************************************************************
******************************************************************************/


#include "ComputePatternDistances.h"
//#define DEBUG_GENERATION

#define MOTIFDBSIZE     500

int main( int argc , char *argv[])
{
    FILE *fp, *fp1, *fp2;
    char *baseName, *patternInfoExt, *patternDataExt, *flistExt, motifFile[N_SIM_MEASURES][400]={'\0'}, *knnExt, filelistFilename[400]={'\0'}, searchFileNames[2000][400] = {'\0'}, tempFilename[400]= {'\0'}, pitchFile[400]={'\0'}, patternInfoFile[400]={'\0'}, patternDataFile[400]={'\0'}, kNNOutFile[400]={'\0'}, *blackListExt;
    float t1,t2, t3,t4, temp[10], pHop, max_factor;
    int err, lenMotifReal, verbos=0, bandDTW, maxNMotifsPairs, nInterFact, **combMTX, searchFileID, *emptySpaceInd, emptySpaceCnt, priorityListInd, nPriorityList, *emptySpacePtr, match_found, NFilesSearch, nRead, NPatternsFile1, NPatternsFile2, overWrite, *isBlackListed1, *isBlackListed2; 
    INDTYPE    NSeed, K,ii,jj, pp,mm, ss, ind1, ind2, lenTS1, lenTS2;
    bool sameFile;
    patternInfo_t *patternInfo1, *patternInfo2;
    pattCntManager_t  *pattCndMan;
    
    DATATYPE **data1Interp, **data2Interp, **data1, **data2, **U1, **L1, **U2, **L2, *accLB;
    DISTTYPE LB_Keogh_EQ, realDist,LB_Keogh_EC,bsf=INF,**costMTX, LB_kim_FL, *bsfArray, bsf_local, distThshld;
    patternDist_t **patternDists;
    segInfo_t *tStampsInterpDummy1, *tStampsInterpDummy2;
    procLogs_t myProcLogs;
    procParams_t myProcParams;
    fileExts_t myFileExts;
  
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
    
    
    if(argc < 15 || argc > 16)
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
    blackListExt = argv[8];
    myProcParams.durMotif = atof(argv[9]);
    distThshld = atoi(argv[10]);
    myProcParams.nInterpFac=atoi(argv[11]);
    myProcParams.dsFactor = atoi(argv[12]);
    if (atoi(argv[13])>0)
    {
        myProcParams.simMeasureRankRefinement = (int*)malloc(sizeof(int)*1);
        myProcParams.simMeasureRankRefinement[0] = atoi(argv[13]);
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
    overWrite = atoi(argv[14]);
     
    if( argc == 16 ){verbos = atoi(argv[15]);}
    
    //############ CRUCIAL PARAMETERS ##################
    myProcParams.minPossiblePitch = 60.0;
    myProcParams.binsPOct = 1200;
    myProcParams.varDur = 0.1;
    myProcParams.threshold = 45.0;
    myProcParams.flatThreshold = 0.8;
    myProcParams.maxPauseDur = 0.5;
    myProcParams.DTWBand = 0.1;
    myProcParams.removeTaniSegs=1;
    
    strcat(kNNOutFile, baseName);
    strcat(kNNOutFile, knnExt);
    
    if (access(kNNOutFile, F_OK) == 0)
    {
        if (verbos ==1)
        {
            printf("Files does exists");
        }
        
        if (overWrite == 0)
        {
            if (verbos ==1)
            {
                printf("File does exists, I am returning");
            }
            
            return 1;
        }
    }
    else
    {
        if (verbos ==1)
        {
            printf("Output file doesn't exist, need to compute the output");
        }
        
    }
    
    
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
    
    
    
    //search file name
    strcat(filelistFilename,baseName);
    strcat(filelistFilename,flistExt); 

    
    fp1 = fopen(filelistFilename, "r");
    if (fp1==NULL)
    {
        printf("Error opening file %s\n", filelistFilename);
        return 0;
    }
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
    err = readPatternDump(patternInfoFile, &patternInfo1, &NPatternsFile1);
    if (err==0)
    {
        return 0;
    }
    lenTS1 = NPatternsFile1*myProcParams.nInterpFac;
    
    memset(tempFilename, '\0', sizeof(char)*400);
    strcat(tempFilename,baseName);
    strcat(tempFilename,blackListExt);         
    isBlackListed1 = (int*)malloc(sizeof(int)*NPatternsFile1);
    readBlackListDump(tempFilename, isBlackListed1);
    memset(tempFilename, '\0', sizeof(char)*400);
    
	

    //since we know the length of the data by now. Lets assign memory for it
    data1 = (DATATYPE **)malloc(sizeof(DATATYPE *)*NPatternsFile1);
    
    strcat(patternDataFile,baseName);
    strcat(patternDataFile,patternDataExt);
    fp1 = fopen(patternDataFile, "rb");

    for (ss = 0 ; ss < NPatternsFile1 ; ss++)
    {
        data1[ss] = (DATATYPE *)malloc(sizeof(DATATYPE)*myProcParams.motifLengths[myProcParams.indexMotifLenLongest]);
        fread(data1[ss], sizeof(DATATYPE), myProcParams.motifLengths[myProcParams.indexMotifLenLongest], fp1);

    }
    
    fclose(fp1);
    
    //generate multiple interpolatedpattCndMan versions
    tStampsDummy1 = (segInfoInterp_t*)malloc(sizeof(segInfoInterp_t)*NPatternsFile1);
    generateInterpolatedSequences(data1, tStampsDummy1, &data1Interp,  &tStampsInterpDummy1, (INDTYPE)NPatternsFile1, &myProcParams);

    
    nPriorityList = NPatternsFile1;

    // generating envelops for the seed motifs
    lenMotifReal = myProcParams.motifLengths[myProcParams.indexMotifLenReal];
    bandDTW = (int)floor(lenMotifReal*myProcParams.DTWBand);
    U1 = (DATATYPE **)malloc(sizeof(DATATYPE *)*lenTS1);
    L1= (DATATYPE **)malloc(sizeof(DATATYPE *)*lenTS1);
    accLB = (DATATYPE *)malloc(sizeof(DATATYPE)*lenMotifReal);
    patternDists = (patternDist_t **)malloc(sizeof(patternDist_t*)*nPriorityList);
    pattCndMan = (pattCntManager_t *)malloc(sizeof(pattCntManager_t)*nPriorityList);
    
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
	if (isBlackListed1[jj]==1)
	{
	  pattCndMan[jj].pattCnt =0;
	  continue;
	}
        patternDists[jj] = (patternDist_t *)malloc(sizeof(patternDist_t)*MOTIFDBSIZE);
	pattCndMan[jj].allocSpace =MOTIFDBSIZE;
	pattCndMan[jj].pattCnt =0;

    }


    for (ss=0;ss<NFilesSearch; ss++)
    {
        //read the data 
    	memset(tempFilename, '\0', sizeof(char)*400);
        strcat(tempFilename,searchFileNames[ss]);
        strcat(tempFilename,patternInfoExt);
        err =  readPatternDump(tempFilename, &patternInfo2, &NPatternsFile2);
	
        if (err==0)
        {
            continue;
        }
        memset(tempFilename, '\0', sizeof(char)*400);
        strcat(tempFilename,baseName);
        strcat(tempFilename,blackListExt);         
        isBlackListed2 = (int*)malloc(sizeof(int)*NPatternsFile2);
        readBlackListDump(tempFilename, isBlackListed2);
        memset(tempFilename, '\0', sizeof(char)*400);
        
        
        printf("File number %d\t%s\n",(int)ss, searchFileNames[ss]);
        
    	if (ss==1019){printf("Checkpoint -3\n");}

        lenTS2 = NPatternsFile2*myProcParams.nInterpFac;

        //since we know the length of the data by now. Lets assign memory for it
        data2 = (DATATYPE **)malloc(sizeof(DATATYPE *)*NPatternsFile2);
        strcat(tempFilename,searchFileNames[ss]);
        strcat(tempFilename,patternDataExt);
        fp1 = fopen(tempFilename, "rb");
        if (fp1==NULL)
        {
            printf("Error opening file %s\n", tempFilename);
            return 0;
        }
        if (ss==1019){printf("Checkpoint -2\n");}
        for (ii = 0 ; ii < NPatternsFile2 ; ii++)
        {
            if (ss==1019){printf("index %d\n",ii);}
            data2[ii] = (DATATYPE *)malloc(sizeof(DATATYPE)*myProcParams.motifLengths[myProcParams.indexMotifLenLongest]);
            fread(data2[ii], sizeof(DATATYPE), myProcParams.motifLengths[myProcParams.indexMotifLenLongest], fp1);
        }
        fclose(fp1);
        memset(tempFilename, '\0', sizeof(char)*400);
        
        if (ss==1019){printf("Checkpoint 0\n");}
        
        //generate multiple interpolated versions
        tStampsDummy2 = (segInfoInterp_t*)malloc(sizeof(segInfoInterp_t)*NPatternsFile2);
        generateInterpolatedSequences(data2, tStampsDummy2, &data2Interp,  &tStampsInterpDummy2, (INDTYPE) NPatternsFile2, &myProcParams);

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

        if (ss==1019){printf("Checkpoint 1\n");}
        for (ii=0; ii < NPatternsFile1; ii++)
        {
            if (isBlackListed1[ii]==1)
            {
                continue;
            }
            for (jj=0; jj < NPatternsFile2; jj++)  
            {
                
                if(isBlackListed2[jj]==1)
                {
                    continue;
                }
                
                ind1 = ii*myProcParams.nInterpFac;
                ind2 = jj*myProcParams.nInterpFac;

                bsf_local = distThshld;
                for (pp = 0; pp < myProcParams.nInterpFac; pp++)
                {
                    for (mm = 0; mm < myProcParams.nInterpFac; mm++)
                    {
                        // skip if its not a valid combinations of the interpolation factor (for optimizing and doing only 9 out of 25 combinations)
                        if (myProcParams.combMTX[pp][mm]==0)
                        {
                            continue;
                        }
                        if (ss==1019){printf("Checkpoint 2\n");}
                        LB_kim_FL = computeLBkimFL(data1Interp[ind1+pp][0], data2Interp[ind2+mm][0], data1Interp[ind1+pp][lenMotifReal-1], data2Interp[ind2+mm][lenMotifReal-1], myProcParams.simMeasureRankRefinement[0]);
                        myProcLogs.totalFLDone++;
                        if (LB_kim_FL< bsf_local)
                        {
                            LB_Keogh_EQ = computeKeoghsLB(U1[ind1+pp],L1[ind1+pp], accLB, data2Interp[ind2+mm],lenMotifReal, bsf_local, myProcParams.simMeasureRankRefinement[0]);
                            myProcLogs.totalLBKeoghEQ++;
                            if(LB_Keogh_EQ < bsf_local)
                            {
                                LB_Keogh_EC = computeKeoghsLB(U2[ind2+mm],L2[ind2+mm],accLB, data1Interp[ind1+pp],lenMotifReal, bsf_local, myProcParams.simMeasureRankRefinement[0]);
                                myProcLogs.totalLBKeoghEC++;
                                if(LB_Keogh_EC < bsf_local)
                                {
                                    if (ss==1019){printf("Checkpoint 3\n");}
                                    realDist = dtw1dBandConst(data1Interp[ind1+pp], data2Interp[ind2+mm], lenMotifReal, lenMotifReal, costMTX, myProcParams.simMeasureRankRefinement[0], bandDTW, bsf_local, accLB);
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
                if (ss==1019){printf("Checkpoint 4\n");}

                if (bsf_local < distThshld)
                {
                    if (bsf_local > 0) // update only when the patterns are different
                    {
                        manageMotifStorage(&patternDists[jj], pattCndMan, patternInfo1[ii].id, patternInfo2[jj].id, bsf_local);
                        myProcLogs.totalPriorityUpdates++;
                    }
                }


            }
            if (ss==1019){printf("Checkpoint 5\n");}
        }
        
        if (ss==1019){printf("Checkpoint 6\n");}

        for(ii=0;ii<lenTS2;ii++)
        {   
            free(U2[ii]);
            free(L2[ii]);
            free(data2Interp[ii]);
        }
        if (ss==1019){printf("Checkpoint 7\n");}
        
        free(U2);
        free(L2);
        free(data2Interp);
        free(data2);
        free(tStampsDummy2);
        free(tStampsInterpDummy2);
        free(isBlackListed2);
    
    }
    if (ss==1019){printf("Checkpoint 8\n");}

    for(ii=0;ii<(lenTS1);ii++)
    {   
        free(U1[ii]);
        free(L1[ii]);
        free(data1Interp[ii]);
    }
    
    free(U1);
    free(L1);
    free(data1Interp);
    free(data1);
    free(tStampsDummy1);
    free(tStampsInterpDummy1);
    free(isBlackListed1);

    dumpKNNPatterns(kNNOutFile, patternDists, pattCndMan, nPriorityList); 
    return 1;


}

int updatePatternStorage(patternDist_t **patternDist, pattCntManager_t *patCntMan)
{
    patternDist_t *newStorage;
    
    newStorage = (patternDist_t *)malloc(sizeof(patternDist_t)*(patCntMan->allocSpace + MOTIFDBSIZE));
    patCntMan->allocSpace += MOTIFDBSIZE;
    
    memcpy(newStorage, *patternDist, sizeof(patternDist_t)*patCntMan->pattCnt);
    
    free(*patternDist);
    *patternDist = newStorage;
  
}

int manageMotifStorage(patternDist_t **patternDist, pattCntManager_t *patCntMan, INDTYPE id1, INDTYPE id2, DISTTYPE dist)
{
    if (patCntMan->allocSpace == patCntMan->pattCnt)
    {
	updatePatternStorage(patternDist, patCntMan);
    }
  
    (*patternDist)[patCntMan->pattCnt].dist = dist ; 
    (*patternDist)[patCntMan->pattCnt].patternID1 = id1 ; 
    (*patternDist)[patCntMan->pattCnt].patternID1 = id2 ;
    patCntMan->pattCnt++;
    

}

int dumpKNNPatterns(char *filename, patternDist_t **patternDists, pattCntManager_t *patCntMan, int nPriorityList) 
{	
	int ii, jj;
	FILE *fp;

	fp = fopen(filename, "w"); 
	for(ii=0;ii<nPriorityList;ii++)
	{
		for(jj=0;jj<patCntMan[ii].pattCnt;jj++)
		{
			fprintf(fp, "%lld\t%lld\t%f\n", patternDists[ii][jj].patternID1, patternDists[ii][jj].patternID2, patternDists[ii][jj].dist);
		}
	}
	fclose(fp);

}

