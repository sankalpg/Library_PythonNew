/******************************************************************************
*******************************************************************************
******                                                                  *******
******     This code is written by Abdullah Al Mueen at the department  *******
******     of Computer Science and Engineering of University of         *******
******     California - RIverside.                                      *******
******                                                                  *******
*******************************************************************************
******************************************************************************/

/*#############################################################################
######                                                                  #######
######     This code is open to use, reuse and distribute at the user's #######
######     own risk and as long as the above credit is ON. The complete #######
######     description of the algorithm and methods applied can be      #######
######     found in the paper - EXACT DISCOVERY OF TIME SERIES MOTIFS   #######
######     by Abdullah Mueen, Eamonn Keogh and Qiang Zhu.               #######
######                                                                  #######
#############################################################################*/







#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <signal.h>
#include <stdint.h>
#include <string.h>



FILE *fp, *fp_out, *fp_black, *fp_timeStamps; 
int ref;

#define FIXPOINT


#ifdef FIXPOINT

#define DATATYPE uint32_t
#define DISTTYPE double
#define INDTYPE long long 

#define MAX_UINT32 4294967296 -1 
#define INF MAX_UINT32


#else   //floatingpoint

#define DATATYPE double
#define DISTTYPE double
#define INDTYPE long long 

#define INF 999999999999.0
#endif //FIXPOINT

DATATYPE **data;
DATATYPE **dref;
DISTTYPE **indices, *topKDist;
INDTYPE *ind, *blacklist, *topKInd;
INDTYPE loc1 , loc2;
double *timeStamps;

double *stdRef;
INDTYPE TIMESERIES;
int LENGTH;
int MAXREF;
DISTTYPE bsf,bsf_old;
    
/* Calculates the distance between two time series x and y. If the distance is
larger than the best so far (bsf) it stops computing and returns the approximate
distance. To get exact distance the bsf argument should be omitted.*/

DISTTYPE distance(DATATYPE *x, DATATYPE *y, int length , DISTTYPE best_so_far = INF )
{
    int i;
    DATATYPE sum = 0;
    DISTTYPE bsf2 = best_so_far*best_so_far;
    for ( i = 0 ; i < length && sum < bsf2 ; i++ )
        sum += (x[i]-y[i])*(x[i]-y[i]);
    
#ifdef FIXPOINT 
    if (sum>MAX_UINT32)
    {
        return sqrt(MAX_UINT32);
    }
    else
    {
        return sqrt(sum);
    }
#else
    return sqrt(sum);
#endif  //FIXPOINT
    }

/*Comparison function for qsort function. Compares two time series by using their
distances from the reference time series. */

int comp1(const void *a,const void *b)
{
    INDTYPE *x=(INDTYPE *)a;
    INDTYPE *y=(INDTYPE *)b;

    if (indices[ref][*x]>indices[ref][*y])
        return 1;
    else if (indices[ref][*x]<indices[ref][*y])
        return -1;
    else
        return 0;
    }


int comp2(const void *a,const void *b)
{
    int *x=(int *)a;
    int *y=(int *)b;

    if (stdRef[*x]<stdRef[*y])
        return 1;
    else if (stdRef[*x]>stdRef[*y])
        return -1;
    else
        return 0;
    }
    

void error(int id)
{
    if(id==1)
        printf("ERROR : Memory can't be allocated!!!\n\n");
    else if ( id == 2 )
        printf("ERROR : File not Found!!!\n\n");
    else if ( id == 3 )
        printf("ERROR : Can't create Output File!!!\n\n");
    else if ( id == 4 )
        printf("ERROR : Invalid Number of Arguments!!!\n\n");

    exit(1);

    }

void stop_exec(int sig)
{

    printf("Current Motif is (%lld",loc1);
    printf(",%lld)",loc2);
    printf(" and the distance is %lf\n",bsf);
    exit(1);
    }

DISTTYPE maintainTopKMotifs(DISTTYPE*motifBuffer, int K, DISTTYPE dist, INDTYPE *topKInd, INDTYPE ind1, INDTYPE ind2)
{
    int ii=0;
    
    for (ii=0;ii<K;ii++)
    {
        if (motifBuffer[ii]>dist)
        {
            memcpy(&motifBuffer[ii+1], &motifBuffer[ii], sizeof(DISTTYPE)*(K-ii));
            motifBuffer[ii] = dist;
            break;
        }
    }
    memcpy(&topKInd[2*(ii+1)], &topKInd[2*ii], 2*sizeof(DISTTYPE)*(K-ii));
    topKInd[2*ii] = ind1;
    topKInd[2*ii+1] = ind2;
    return motifBuffer[K];
}

void updateBlackList(INDTYPE*blacklist, double*timeStamps, double blackDur, INDTYPE loc)
{
    INDTYPE ii=1, loc_curr;
    while (1)
    {
        loc_curr = loc + ii;
        if ((loc_curr>=TIMESERIES) || (fabs(timeStamps[loc]-timeStamps[loc_curr])>blackDur))
        {
            ii=ii;
            break;
        }
        else
        {
            blacklist[loc_curr]=1;
            ii++;
        }
    }
    ii=-1;
    while (1)
    {
        loc_curr = loc + ii;
        if ((loc_curr<=0) || (fabs(timeStamps[loc]-timeStamps[loc_curr])>blackDur))
        {
            ii=ii;
            break;
        }
        else
        {
            blacklist[loc_curr]=1;
            ii--;
        }
    }
    
}
int main(  int argc , char *argv[] )
{

    DISTTYPE d;
    DATATYPE dd;
    INDTYPE i , j , blackInd, black, ii,tsInd=0;
    INDTYPE offset = 0;
    int abandon = 0 ,  r = 0;
    DISTTYPE ex , ex2 , mean, std;
    int *rInd;
    INDTYPE length;
    int clear = 0;
    double *cntr;
    double t1,t2;
    int verbos = 0;
    double count = 0;
    int topK = 1;
    double ts, blackDur;
    
    /* taking inpput time series from the file in the data array and Normalizing
    them as well. */

    bsf = INF;
    i = 0;
    j = 0;
    blackInd=0;
    ex = ex2 = 0;

    t1 = clock();
    signal(SIGINT,stop_exec);


    if(argc < 11 || argc > 12)
    {
        printf("Invalid number of arguments!!!");
        exit(1);
    }


    fp = fopen(argv[1],"r");
    TIMESERIES = atol(argv[2]);
    LENGTH = atoi(argv[3]);
    MAXREF = atoi(argv[4]);
    clear = LENGTH;
    bsf_old = atof(argv[5])+0.000001;
    fp_out = fopen(argv[6],"w");
    fp_black = fopen(argv[7],"r");
    topK = atoi(argv[8]);
    fp_timeStamps = fopen(argv[9],"r");
    blackDur = atof(argv[10]);
    
    if( argc == 12 )
        verbos = atoi(argv[11]);

    if( verbos == 1 )
        printf("\nNumber of Time Series : %lld\nLength of Each Time Series : %d\n\n",TIMESERIES,LENGTH);

    data = (DATATYPE **)malloc(sizeof(DATATYPE *)*TIMESERIES);
    ind = (INDTYPE *)malloc(sizeof(INDTYPE)*TIMESERIES);
    blacklist = (INDTYPE *)calloc(TIMESERIES, sizeof(INDTYPE));
    topKDist = (DISTTYPE *)malloc(sizeof(DISTTYPE)*topK);
    topKInd =  (INDTYPE *)malloc(sizeof(INDTYPE)*topK*2);
    timeStamps = (double*)malloc(sizeof(double)*TIMESERIES);
    
    for (ii=0;ii<topK;ii++)
    {
        topKDist[ii]=INF;
    }
    
    //filling blacklist array
    while(fscanf(fp_black,"%lld",&black) != EOF && blackInd < TIMESERIES)
    {
        blacklist[black]=1;
        blackInd++;
    }
    
    //filling timeStamps array
    while(fscanf(fp_timeStamps,"%lf",&ts) != EOF && tsInd < TIMESERIES)
    {
        timeStamps[tsInd]=ts;
        tsInd++;
    }
    
    
    if( data == NULL || ind == NULL )
    {
        error(1);
    }
    data[0] = (DATATYPE *)malloc(sizeof(DATATYPE)*LENGTH);
    
   // printf("%ld %ld\n\n",data,data[0]);
    if( data[0] == NULL )
    {   
        error(1);
    }
    //printf("Hello2\n");
#ifdef FIXPOINT    
    while(fscanf(fp,"%d",&dd) != EOF && i < TIMESERIES)
#else
    while(fscanf(fp,"%lf",&dd) != EOF && i < TIMESERIES)
#endif         
    {
        data[i][j] = dd;
        //ex += d;
        //ex2 += d*d;
        if( j == LENGTH - 1 )
        {
            //mean = ex = ex/LENGTH;
            //ex2 = ex2/LENGTH;
            //std = sqrt(ex2-ex*ex);
            /*for( int k = 0 ; k < LENGTH ; k++ )
            {
                data[i][k] = (data[i][k]-mean)/std;
                
            }*/
            //ex = ex2 = 0;
            ind[i] = i;
            i++;
            if( i <= TIMESERIES )
            {
                data[i] = (DATATYPE *)malloc(sizeof(DATATYPE)*LENGTH);
                
            }
            //printf("Hello2.5\n");
            if( data[i] == NULL )
            { 
                error(1);
                
            }

            j = 0;
        }
        else
            j++;
    }

    fclose(fp);
    //printf("Hello3\n");
    
    if(verbos == 1)
    {
        printf("Data Have been Read and Normalized\n\n");
    }

    dref = (DATATYPE **)malloc(MAXREF*sizeof(DATATYPE *));
    indices = (DISTTYPE **)malloc(MAXREF*sizeof(DISTTYPE *));
    stdRef = (double *)malloc(MAXREF*sizeof(double));
    cntr = (double *)malloc(MAXREF*sizeof(double));
    rInd = (int *)malloc(MAXREF*sizeof(int));
    
    if ( dref == NULL || indices == NULL || stdRef == NULL || cntr == NULL || rInd == NULL )
    {
        error(1);
    }

    //////////////////////////////////////////////////////////////////////////

    /*Generating the reference time series. Here it is a random time series*/

    srand ( time(NULL) );

    for( r = 0 ; r < MAXREF ; r++ )
    {

        dref[r] = (DATATYPE *)malloc(sizeof(DATATYPE)*LENGTH);
        indices[r] = (DISTTYPE *)malloc(sizeof(DISTTYPE)*TIMESERIES);
        if( dref[r] == NULL || indices[r] == NULL )
        {
            error(1);
        }

        INDTYPE random_pick = rand() % TIMESERIES;
        for( i = 0 ; i < LENGTH ; i++ )
        {
            dref[r][i] = data[random_pick][i];
        }

        if( verbos == 1 )
        {
            printf("\nThe %lldth Time Series is chosen as %dth reference\n",random_pick,r);
        }

            /*Creating the indices array which is a 2 by TIMESERIES
            sized vector having the indices (to the data array) and distances (from the
            reference time series) in each row.*/


            ex = 0;
            ex2 = 0;
            
            for( i = 0 ; i < TIMESERIES ; i++ )
            {
                if ((i == random_pick) || (blacklist[i]==1))
                { indices[r][i] = INF; continue; }
                d = indices[r][i] = distance(data[i],dref[r],LENGTH);
                count = count + 1;
                ex += d;
                ex2 += d*d;
                if ( abs(i - random_pick ) <= clear )  continue;
                if (( d < bsf )&&(d > bsf_old))
                {
                    loc1 = i; loc2 = random_pick;
                    bsf = maintainTopKMotifs(topKDist, topK, d, topKInd, loc1, loc2);
                    //updateBlackList(blacklist, timeStamps, blackDur, loc1);
                    //updateBlackList(blacklist, timeStamps, blackDur, loc2);
                    
                    if(verbos == 1)
                        printf("New best-so-far is %f and (%lld , %lld) are the new motif pair\n",bsf,loc1,loc2);
                }

            }

            ex = ex/(TIMESERIES-1);
            ex2 = ex2/(TIMESERIES-1);
            std = sqrt(ex2-ex*ex);
            


            rInd[r] = r;
            stdRef[r] = std;
            cntr[r] = 0;
            ////////////////////////////////////////////////////////////////////
  }
  
        if(verbos == 1)
        {
            printf("\nReferences are picked and Dist has been Computed\n\n");
        }

        /*Sort the standard Deviations*/
        qsort(rInd,MAXREF,sizeof(int),comp2);

        ref = rInd[0];

        INDTYPE remaining = TIMESERIES;

        /*Sort indices using the distances*/
        qsort(ind,TIMESERIES,sizeof(INDTYPE),comp1);
        ///////////////////////////////////


         /*Motif Search loop of the algorithm that finds the motif. The algorithm
        computes the distances between a pair of time series only when it thinks
        them as a potential motif'*/

        if(verbos == 1)
        {
            printf("Orderings have been Computed and Search has begun\n\n");
        }
        offset = 0;
        abandon = 0;
        
        while (!abandon && offset < remaining)
        {
            abandon = 1;
            offset++;

            for(i = 0 ; i < remaining - offset ; i++ )
            {
                INDTYPE left = ind[i];
                INDTYPE right = ind[i + offset];
                if (blacklist[left] + blacklist[right] >=1)
                {
                    continue;
                }
                if( abs(left-right) <= clear )
                {
                    continue;
                }

                //According to triangular inequality distances between left and right
                //is obviously greater than lower_bound.
                double lower_bound = 0;
                r = 0;
                do
                {
                    lower_bound = fabs(indices[rInd[r]][right] - indices[rInd[r]][left]);
                    r++;
                }while( r < MAXREF && lower_bound < bsf );


                if (r >= MAXREF && lower_bound < bsf)
                {

                    abandon = 0;
                    count =  count + 1;
                    d = distance( data[left] , data[right] , LENGTH , bsf );
                    signal(SIGINT,SIG_IGN);
                    if (( d < bsf )&&(d > bsf_old))
                    {
                        t2 = clock();
                        loc1 = left;
                        loc2 = right;
                        bsf = maintainTopKMotifs(topKDist, topK, d, topKInd, loc1, loc2);
                        //updateBlackList(blacklist, timeStamps, blackDur, loc1);
                        //updateBlackList(blacklist, timeStamps, blackDur, loc2);
                        if(verbos == 1)
                        {
                            printf("New best-so-far is %lf and (%lld , %lld) are the new motif pair\n",bsf,loc1,loc2);
                    
                        }

                        }
                    signal(SIGINT,stop_exec);
                    }
                }
            }
            

        t2 = clock();
        if(verbos == 1)
            printf("\nExecution Time was : %lf seconds\n",(t2-t1)/CLOCKS_PER_SEC);
        printf("\n\nFinal Motif is the pair ( %lld",loc1);
        printf(", %lld ) and the Motif distance is %lf\n",loc2,bsf);
        for (ii=0;ii<topK;ii++)
        {
            fprintf(fp_out, "%lld\t%lld\t%f\n",topKInd[2*ii],topKInd[2*ii+1],topKDist[ii]);
        }
        fclose(fp_out);
        }
