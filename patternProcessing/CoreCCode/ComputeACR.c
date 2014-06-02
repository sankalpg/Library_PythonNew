/******************************************************************************
*******************************************************************************
****** Sankalp Gulati                                                   *******
****** MUSIC TECHNOLOGY GROUP, UPF, BARCELONA                           *******
******                                                                  *******
*******************************************************************************
******************************************************************************/



#include "MotifDataIO.h"


void computeACR(DATATYPE *data, int lenMotifReal, double *acrArray, int nPoints)
{
    int ii,jj;
    double mean=0, std=0, temp, sum, max_val=0;
    //mean subtraction
    for(ii=0;ii<lenMotifReal;ii++)
    {
        mean+=data[ii];
    }
    mean = mean/lenMotifReal;
    
    for(ii=0;ii<lenMotifReal;ii++)
    {
        temp = (data[ii]-mean);
        std+=temp*temp;
    }
    std = sqrt(std/lenMotifReal);
    
    //normalizing data
    for(ii=0;ii<lenMotifReal;ii++)
    {
        data[ii] = (data[ii]-mean)/std;
    }
    
    for(ii=0;ii<nPoints;ii++)
    {
        sum =0;
        for (jj=0;jj<lenMotifReal-ii;jj++)
        {
            sum+= data[jj]*data[jj+ii];
        }
        acrArray[ii] =sum*lenMotifReal/(lenMotifReal-ii);
        if (acrArray[ii] > max_val)
        {
            max_val = acrArray[ii];
        }
    }
    
    for(ii=0;ii<nPoints;ii++)
    {
        acrArray[ii] = acrArray[ii]/max_val;
    }
    
    
}



int main( int argc , char *argv[])
{
    FILE *fp, *fp_out;
    char *baseName, dumpFile[400]={'\0'}, logFile[400]={'\0'}, *dumpExt;
    float t1,t2;
    int lenMotifReal, verbos=0, bandDTW, nPoints, **combMTX, nInterFact; 
    INDTYPE    lenTS, count_DTW=0, numLinesInFile, K,ii,jj, cntr1=0;;
    double **acrArray;
    
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
    
    if(argc < 8 || argc > 9)
    {
        printf("\nInvalid number of arguments!!!\n");
        exit(1);
    }
    
    baseName = argv[1];
    myFileExts.pitchExt = argv[2];
    myFileExts.tonicExt = argv[3];
    myFileExts.segExt = argv[4];
    dumpExt = argv[5];
    myProcParams.durMotif = atof(argv[6]);
    nPoints = atoi(argv[7]);
    
    if( argc == 9 ){verbos = atoi(argv[8]);}
    
    
    
    //############ CRUCIAL PARAMETERS ##################
    myProcParams.minPossiblePitch = 60.0;
    myProcParams.binsPOct = 1200;
    myProcParams.varDur = 0.1;
    myProcParams.threshold = 45.0;
    myProcParams.flatThreshold = 0.8;
    myProcParams.maxPauseDur = 0.5;
    myProcParams.DTWBand = 0.1;
    myProcParams.removeTaniSegs=1;
    
    myProcParams.nInterpFac = 1;
    myProcParams.dsFactor=1;
    
    
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
    strcat(dumpFile,baseName);
    strcat(dumpFile,dumpExt);
    
    
    lenTS = readPreProcessGenDB(&dataInterp, &tStampsInterp, &lenMotifReal, baseName, &myFileExts, &myProcParams, &myProcLogs, verbos);
    
   
    acrArray = (double **)malloc(sizeof(double*)*lenTS);
    for(ii=0;ii<lenTS;ii++)
    {
        acrArray[ii] = (double *)malloc(sizeof(double)*nPoints);
    }
    //timing
    t1 = clock();
    
    //Computing all combinations to obtain top K best matches
    for(ii=0;ii<lenTS;ii++)
    {
        computeACR(dataInterp[ii],lenMotifReal, acrArray[ii], nPoints);
    }

    //timing
    t2 = clock();
    myProcLogs.timeDiscovery = (t2-t1)/CLOCKS_PER_SEC;
    if (verbos)
    {printf("Time taken to compute all combinations :%f\n",myProcLogs.timeDiscovery);}
    
    
    fp = fopen(dumpFile,"w");
    for(ii=0;ii<lenTS;ii++)
    {
        for(jj=0;jj<nPoints;jj++)
        {
                fprintf(fp,"%f\t", acrArray[ii][jj]);
        }
        
        fprintf(fp,"\n");
    }
    fclose(fp);
    
    //Memory clearing
    for(ii=0;ii<lenTS;ii++)
    {
        free(dataInterp[ii]);
        free(acrArray[ii]);
    }
    free(dataInterp);
    free(acrArray);
    free(tStampsInterp);
    return 1;
    
}


