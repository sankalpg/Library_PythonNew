/******************************************************************************
*******************************************************************************
****** Sankalp Gulati                                                   *******
****** MUSIC TECHNOLOGY GROUP, UPF, BARCELONA                           *******
******                                                                  *******
*******************************************************************************
******************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <signal.h>
#include <stdint.h>
#include <string.h>
#include "../../similarityMeasures/dtw/dtw.h"
#include "../../basicDSPFuncs/basicDSPCFuncs.h"
#include <float.h>

//preprocessor defines!!
//#define FIXPDATA
#define FLOATDATA


#ifdef FIXPDATA
#define DATATYPE       uint32_t
#else
#define DATATYPE       double
#endif

#define DISTTYPE        double
#define INDTYPE        long long

#define INF 999999999999.0
# define computeLBkimFL(a,b,c,d) ((a-b)*(a-b)) + ((c-d)*(c-d))


DATATYPE computedist(DATATYPE * data1, DATATYPE * data2, int len)
{
    DISTTYPE sum=0;
    DATATYPE diff;
    int ii=0;
    
    for(ii=0;ii<len;ii++)
    {
        diff = data1[ii] - data2[ii];
        sum+=(DISTTYPE)(diff*diff);
    }
    
    return sum;
    
}

DISTTYPE manageTopKMotifs(DISTTYPE *topKDist, INDTYPE *topKInd, double *tStamps, int K, INDTYPE ind1 , INDTYPE ind2, DISTTYPE dist, double blackDur)
{
    int ii=0;
    int sortInd = -1;
    int matchInd = -1;
    
    for(ii=0;ii<K;ii++)
    {
        if ((topKDist[ii] > dist)&&(sortInd==-1))
        {
            sortInd=ii;
        }
        // searching if we already have a motif in out top K list which is near to the currently good match
        if ((fabs(tStamps[topKInd[2*ii]]-tStamps[ind1]) < blackDur) || (fabs(tStamps[topKInd[2*ii+1]]-tStamps[ind1]) < blackDur) || (fabs(tStamps[topKInd[2*ii]]-tStamps[ind2]) < blackDur) || (fabs(tStamps[topKInd[2*ii+1]]-tStamps[ind2]) < blackDur))
        {
            matchInd=ii;
            break;
        }
        
    }
    if (sortInd==-1)//we couldn't get any satisfactory replacement before we get a close neighbour
    {
        return topKDist[K-1];
    }
    //There are three possibilities
    //1) There is no match found in the existing top motifs, simplest
    if (matchInd==-1)
    {
        memcpy(&topKDist[sortInd+1], &topKDist[sortInd], sizeof(DISTTYPE)*(K-(sortInd+1)));
        memcpy(&topKInd[2*(sortInd+1)], &topKInd[2*sortInd], sizeof(INDTYPE)*2*(K-(sortInd+1)));
        topKDist[sortInd]=dist;
        topKInd[2*sortInd]=ind1;
        topKInd[2*sortInd+1]=ind2;
    }
    else if (sortInd == matchInd)
    {
        topKDist[sortInd]=dist;
        topKInd[2*sortInd]=ind1;
        topKInd[2*sortInd+1]=ind2;
    }
    else if (sortInd < matchInd)
    {
        memcpy(&topKDist[sortInd+1], &topKDist[sortInd], sizeof(DISTTYPE)*(matchInd-sortInd));
        memcpy(&topKInd[2*(sortInd+1)], &topKInd[2*sortInd], sizeof(INDTYPE)*2*(matchInd-sortInd));
        topKDist[sortInd]=dist;
        topKInd[2*sortInd]=ind1;
        topKInd[2*sortInd+1]=ind2;
    }
    
    return topKDist[K-1];
    
}

FILE *fp_DB, *fp_TS, *fp_out;

int main( int argc , char *argv[])
{
    INDTYPE    lenTS, count_DTW=0,offset;
    int         lenMotif;
    double      blackDur;
    int         K,ii,jj,ll, numReads;
    int verbos = 0;
    double t1,t2;
    double *tStamps;
    DATATYPE **data;
    DISTTYPE LB_Keogh_EQ, realDist,LB_Keogh_EC;
    DISTTYPE bsf=INF;
    DISTTYPE *topKDist;
    INDTYPE *topKInd;
    DISTTYPE *cost, **cost2;
    DATATYPE *U,*L,*U2,*L2, *accLB, LB_kim_FL;
    int band;
    
    
    
    if(argc < 8 || argc > 9)
    {
        printf("Invalid number of arguments!!!");
        exit(1);
    }
    
    fp_DB = fopen(argv[1],"rb");
    fp_TS = fopen(argv[2],"r");
    lenTS  = atol(argv[3]);
    lenMotif = atoi(argv[4]);
    K = atoi(argv[5]);
    blackDur = atof(argv[6]);
    fp_out = fopen(argv[7],"r");
    
    if( argc == 9 ){verbos = atoi(argv[8]);}
    
    if( verbos == 1 )
        printf("\nNumber of Time Series : %lld\nLength of Each Time Series : %d\n\n",lenTS,lenMotif);
    
    
    //Memory allocation
    data = (DATATYPE **)malloc(sizeof(DATATYPE *)*lenTS);
    tStamps = (double *)malloc(sizeof(double)*lenTS);
    topKDist = (DISTTYPE *)malloc(sizeof(DISTTYPE)*K);
    topKInd = (INDTYPE *)malloc(sizeof(INDTYPE)*K*2);
    cost = (DISTTYPE *)malloc(sizeof(DISTTYPE)*lenMotif*lenMotif);
    cost2 = (DISTTYPE **)malloc(sizeof(DISTTYPE *)*lenMotif);
    U = (DATATYPE *)malloc(sizeof(DATATYPE)*lenMotif);
    L = (DATATYPE *)malloc(sizeof(DATATYPE)*lenMotif);
    U2 = (DATATYPE *)malloc(sizeof(DATATYPE)*lenMotif);
    L2 = (DATATYPE *)malloc(sizeof(DATATYPE)*lenMotif);    
    accLB = (DATATYPE *)malloc(sizeof(DATATYPE)*lenMotif);
    //initialization
    for(ii=0;ii<K;ii++)
    {
        topKInd[2*ii]=0;
        topKInd[2*ii+1]=0;
        topKDist[ii]=(DISTTYPE)INF;
    }
    
    for(ii=0;ii<lenMotif;ii++)
    {
        cost2[ii] = (DISTTYPE *)malloc(sizeof(DISTTYPE)*lenMotif);
    }
    
    //timing
    t1 = clock();
    
    //reading data, loading into memory
    for (ii=0;ii<lenTS;ii++)
    {
        data[ii] = (DATATYPE *)malloc(sizeof(DATATYPE)*lenTS);
        numReads = fread( data[ii], sizeof(DATATYPE), lenMotif, fp_DB);
        //printf("%d\t%d\n",ii,numReads);
        if (numReads<lenMotif)
        {
            printf("Binary data size was wrongly specified");
            exit(1);
        }
    }
    fclose(fp_DB);
    
    //timing
    t2 = clock();
    
    if (verbos)
    {
        printf("Time taken to load the data :%f\n",(t2-t1)/CLOCKS_PER_SEC);
    }
    
    //timing
    t1 = clock();
    
    //reading time stamps, loading into memory
    numReads = fread(tStamps, sizeof(double), lenTS, fp_TS);
    if (numReads<lenTS)
    {
        printf("Error in reading time stamps, size doesn't match specific length of the time series");
        exit(1);
    }
    fclose(fp_TS);
    
    //timing
    t2 = clock();
    
    if (verbos)
    {
        printf("Time taken to load time stamps :%f\n",(t2-t1)/CLOCKS_PER_SEC);
    }
    
    for (ii=1;ii<lenMotif;ii++)
        {
          for (jj=1;jj<lenMotif;jj++)
          {
              cost[(ii*lenMotif)+ jj] = FLT_MAX;
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
        
        for(jj=ii;jj<lenTS;jj++)
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
                        if (realDist ==1482.0)
                        {
                            realDist = 1482.0;
                        }
                        bsf = manageTopKMotifs(topKDist, topKInd, tStamps, K, ii, jj, realDist, blackDur);
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
        printf("motif pair is %f\t%f\t%f\n", tStamps[topKInd[2*ii]],tStamps[topKInd[2*ii+1]], topKDist[ii]);
    }
    printf("Total dtw computations %lld\n", count_DTW);    
    
    return 1;
    
}