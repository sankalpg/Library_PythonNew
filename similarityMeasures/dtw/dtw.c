/*  
Author: Sankalp gulati
email: sankalp.gulati@gmai.com
Affiliation: Universitat Pompeu Fabra

License: to be decided !!!   
*/

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include "dtw.h"
#include <float.h>
#include "tables.h"

#define BINPOMAKAM 53.0



#define N_PER_CENT 10   //number of rows for one cent
#define CENT_PER_BIN 1
#define MULT_FACTOR_TABLE 10   //N_PER_CENT*CENT_PER_BIN



//#define ENABLE_EA


//########################################## Similarity Measures #####################################
double distEuclidean(double a, double b)
{
    double diff;
    diff = a-b;
    return sqrt(diff*diff);
}

double distSqEuclidean(double a, double b)
{
    double diff;
    diff = a-b;
    return (diff*diff);
}

double distCityBlock(double a, double b)
{
    return fabs(a-b);
}

double distShiftCityBlock(double a, double b)
{
    int index;
    index = (int)(fabs(a-b)*MULT_FACTOR_TABLE);
    return simShifCB[min(SIM_TB_SIZE, index)];
}

double distTMM_CityBlock(double a, double b)
{
    double diff1, diff2;
    diff1 = fabs(a-b);
    diff2 = fmod(diff1, BINPOMAKAM);
    return min(diff2, BINPOMAKAM-diff2);
}

double distShiftLinExp(double a, double b)
{
    int index;
    index = (int)(fabs(a-b)*MULT_FACTOR_TABLE);
    return simShifCBExp[min(SIM_TB_SIZE, index)];
}

simMeasure mySimMeasure[N_SIM_MEASURES] = {distEuclidean, distSqEuclidean, distCityBlock, distShiftCityBlock, distTMM_CityBlock, distShiftLinExp};
const char*SimMeasureNames[N_SIM_MEASURES] = {"Euclidean", "SqEuclidean", "CityBlock", "ShiftCityBlock", "TMM_CityBlock", "ShiftLinExp"};

distMeasures myDistMeasures[N_SIM_MEASURES] = {distEuclidean, distSqEuclidean, distCityBlock, distShiftCityBlock, distTMM_CityBlock, distShiftLinExp};






//////////////////////////////// TO BE REMOVE AFTER SOME TIME TEMP ONES \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
// euclidean distance
double EucDist(double a, double b)
{
    double diff;
    diff = a-b;
        return (diff*diff);
}

double octBy2WrappedCitiblock(double a, double b)
{
    double diff1, diff2;
    diff1 = fabs(a-b);
    diff2 = fmod(diff1, BINPOMAKAM);
    return min(diff2, BINPOMAKAM-diff2);
    
}


double LocalDist(double *a, double*b, int dim, DistMethods fptr, double *wghts)
{
        switch (dim)
        {
        case 1:
            return ((*fptr)(a[0],b[0]))*wghts[0];
        case 2:
            return ((*fptr)(a[0],b[0]))*wghts[0] + ((*fptr)(a[1],b[1]))*wghts[1];
        case 3:
            return ((*fptr)(a[0],b[0]))*wghts[0] + ((*fptr)(a[1],b[1]))*wghts[1] + ((*fptr)(a[2],b[2]))*wghts[2];
        }
        return 0.0;
}

double min3(double a, double b, double c)
{

if (b < a)
    a = b;
if (c < a)
    return c;
return a;
}

/*Euclidean distance between two sequences*/

double euclideanSeq(double *x, double*y, int x_len, int y_len, double**cost, int dist_type, int bandwidth, double bsf, double *accLB)
{
    int ii;
    double dist=0;
    
    for (ii=0;ii<min(x_len, y_len);ii++)
    {
        dist+=mySimMeasure[dist_type](x[ii], y[ii]);
    }
    
    return dist;
}

int path_Euclidean(double **cost, int n, int m)
{
    return n;
}



/*
 * All the DTW variants for one dimensional time series
 * 
 * NOTE: this first category of dtw functions are quite generic and can be used in C/Cython/Python/Mex/MATLAB. They do not require any specific consideration from the calling function.
 * Whereas there are some specific dtw versions I wrote specifically for my usage which incorporate lower bounding, early abandoning and do not initialize cost matrix for faster processing. (find them after this category)
 * 
 * We should Ideally be keeping input param structure same for all variants so that they can be later indexed by a function pointer
*   # ==    Global  Local   Subsequence
*           No      No      No
*           Yes     No      No
*           No      Yes     No
*           Yes     Yes     No
*           No      Yes     Yes (with sub version it always make sense to have local alignment)
*           Yes     Yes     Yes (with sub version it always make sense to have local alignment)
 */

double dtw_GLS(double *x, double*y, int x_len, int y_len, double*cost, dtwParams_t params)
{
        // declarations of variables
        int i,j, bandwidth;    
        float min_vals, factor, max_del, bound; 
        
        max_del = max(params.delStep, params.moveStep);
        bandwidth = params.bandwidth;
        
        if (params.globalType==0)
        {
            /*This is along 45 degree, the usual way to do in dtw*/
            factor=1;
            /*Since this is subsequence dtw, even if bandwidth is smaller than abs(x_len-y_len) its not a problem*/
        }
        else if (params.globalType==1)
        {
            factor = (float)y_len/(float)x_len;
        }
        
        if (params.initCostMtx==1)
        {
            //putting infi in all cost mtx
            for (i=0;i<x_len;i++)
            {
                for (j=0;j<y_len;j++)
                {
                    cost[(i*y_len)+ j] = FLT_MAX;
                }
            
            }            
        }
        if (params.hasGlobalConst==0)
        {
            bandwidth = max(x_len, y_len);  
        }
        
        if (params.reuseCostMtx==0)
        {
            

            //Initializing the row and columns of cost matrix
            cost[0]= (*myDistMeasures[params.distType])(x[0],y[0]);
            if (params.initFirstCol==1)
            {
                for (i=1;i<min(bandwidth+1, x_len);i++)
                {
                    cost[i*y_len]=(*myDistMeasures[params.distType])(x[i],y[0]) + cost[((i-1)*y_len)];
                }
            }
            for (j=1;j< min(bandwidth+1, y_len);j++)
            {
                cost[j]=(*myDistMeasures[params.distType])(x[0],y[j]);
            }
            
            /*Initializing other rows and colms till we can support out local constraints to be applied*/
            for (i=1;i<=min(bandwidth+1, x_len-1);i++)
            {
                for (j=1;j<=max_del;j++)
                {
                    cost[(i*y_len)+ j] = (*myDistMeasures[params.distType])(x[i],y[j]) + min3(((i-1)*y_len)+ j, cost[((i-1)*y_len)+(j-1)], cost[((i)*y_len)+(j-1)]);
                }
            }
            for (j=1;j<= min(bandwidth+1,y_len-1);j++)
            {
                for (i=1;i<=max_del;i++)
                {
                    cost[(i*y_len)+ j] = (*myDistMeasures[params.distType])(x[i],y[j]) + min3(cost[((i-1)*y_len)+(j)], cost[((i-1)*y_len)+(j-1)], cost[((i)*y_len)+(j-1)]);
                }
            }
            
            //filling in all the cumulative cost matrix
            for (i=max_del;i<x_len;i++)
            {
                for (j=max(max_del, i-bandwidth);j<=min(y_len-1, i+bandwidth);j++)
                {
                    cost[(i*y_len)+ j] = (*myDistMeasures[params.distType])(x[i],y[j]) + 
                                            min3(cost[(i-params.moveStep)*y_len+(j-params.delStep)], cost[((i-params.diagStep)*y_len)+(j-params.diagStep)], cost[((i-params.delStep)*y_len)+(j-params.moveStep)]);
                }
            
            }
        }
        
        if (params.isSubsequence==1)
        {   
            min_vals = FLT_MAX;
            for (i=x_len-1;i>=max(x_len-1-bandwidth, 0);i--)
            {
                j = y_len -1;
                if(cost[(i*y_len)+ j] < min_vals)
                {
                    min_vals = cost[(i*y_len)+ j];
                }
            }
            for (j=y_len-1;j>=max(y_len-1-bandwidth,0);j--)
            {
                i = x_len -1;
                if(cost[(i*y_len)+ j] < min_vals)
                {
                    min_vals = cost[(i*y_len)+ j];
                }
            }
        }
        else
        {
            min_vals = cost[x_len*y_len -1];
        }
        return min_vals;

}


// This is the standard dtw implementation for multidimensional time series. In case the dimension is > 1 the point to point distance is computed as weighted point to point 
//distance between each dimension. The weight function is input to the function.
double dtwNd_std(double *x, double*y, MatrixSize*size_x, MatrixSize*size_y, int NFeatDim, double*cost, Config*  myConfig)
{
        // declarations of variables
        int i,j, NRows, NCols;	
        DistMethods myDistMethods[5]={NULL};
        
        //some initialization
        NRows=size_x->Nrow;
        NCols=size_y->Nrow;
        myDistMethods[Euclidean]=&EucDist;
        
        //printf("%f\t%f\n",myConfig->DistWghts[0],myConfig->DistWghts[1]);
        //Initializing the row and columns of cost matrix
        cost[0]=LocalDist(&x[0],&y[0],NFeatDim, myDistMethods[myConfig->DistMthd],myConfig->DistWghts);
        for (i=1;i<NRows;i++)
        {
        cost[i*NCols]=LocalDist(&x[i*size_x->Ncol],&y[0], NFeatDim, myDistMethods[myConfig->DistMthd], myConfig->DistWghts) + cost[(i-1)*NCols];
        }
        for (j=1;j<NCols;j++)
        {
        cost[j]=LocalDist(&x[0],&y[j*size_y->Ncol], NFeatDim, myDistMethods[myConfig->DistMthd], myConfig->DistWghts) + cost[(j-1)];
        }
        //filling in all the cumulative cost matrix
        for (i=1;i<NRows;i++)
        {
        for (j=1;j<NCols;j++)
        {
            cost[(i*NCols)+ j] = LocalDist(&x[i*size_x->Ncol],&y[j*size_y->Ncol], NFeatDim, myDistMethods[myConfig->DistMthd], myConfig->DistWghts) + 
                                    min3(cost[(i-1)*NCols+j], cost[(i-1)*NCols+(j-1)], cost[i*NCols+(j-1)]);
        }
        
        }
        return cost[(NRows*NCols)-1];
}

// Compute the warp path starting at cost[startx, starty]
// If startx = -1 -> startx = n-1; if starty = -1 -> starty = m-1
int pathLocal(double *cost, int n, int m, int startx, int starty, DTW_path *p, dtwParams_t params)
{
//params.delStep, params.moveStep and params.diagStep
int i, j, k, z1, z2, a;
int *px;
int *py;
double min_cost;

if ((startx >= n) || (starty >= m))
    return 0;

if (startx < 0)
    startx = n - 1;

if (starty < 0)
    starty = m - 1;
    
i = startx;
j = starty;
k = 1;

// allocate path for the worst case
px = (int *) malloc ((startx+1) * (starty+1) * sizeof(int));
py = (int *) malloc ((startx+1) * (starty+1) * sizeof(int));

px[0] = i;
py[0] = j;

while ((i > 0) || (j > 0))
    {
    if (i == 0)
        j--;
    else if (j == 0)
        i--;
    else if ((i==params.moveStep) && (j==params.moveStep)) // cannot do del
        {
        for (a=1;a<=params.diagStep;a++) // step diagonally
        {
            i--;
            j--;
        }
        }
    else if ((i==params.diagStep) && (j==params.diagStep))
        {
        for (a=1;a<=params.diagStep;a++)
        {
            i--;
            j--;
        }
        }
    else if (i==params.moveStep)
        {
        min_cost = min(cost[(i-params.diagStep)*m+(j-params.diagStep)], 
                       cost[(i-params.moveStep)*m+(j-params.delStep)]);
        
        if (cost[(i-params.diagStep)*m+(j-params.diagStep)] == min_cost)
            {
            for (a=1;a<=params.diagStep;a++)
            {
                i--;
                j--;
            }
            }
        else
            {
            for (a=1;a<=params.moveStep;a++)
            {
                i--;
            }
            for (a=1;a<=params.delStep;a++)
            {
                j--;
            }
        }
        }
    else if (j==params.moveStep)
        {
        min_cost = min(cost[(i-params.diagStep)*m+(j-params.diagStep)], 
                        cost[(i-params.delStep)*m+(j-params.moveStep)]);
        
        if (cost[(i-params.diagStep)*m+(j-params.diagStep)] == min_cost)
            {
            for (a=1;a<=params.diagStep;a++)
            {
                i--;
                j--;
            }
            }
        else
            {
            for (a=1;a<=params.delStep;a++)
            {
                i--;
            }
            for (a=1;a<=params.moveStep;a++)
            {
                j--;
            }
        }
        }
    else
        {
        min_cost = min3(cost[(i-params.delStep)*m+(j-params.moveStep)],
                        cost[(i-params.diagStep)*m+(j-params.diagStep)], 
                        cost[(i-params.moveStep)*m+(j-params.delStep)]);
        
        if (cost[(i-params.diagStep)*m+(j-params.diagStep)] == min_cost)
            {
            for (a=1;a<=params.diagStep;a++)
            {
                i--;
                j--;
            }
            }
        else if (cost[(i-params.moveStep)*m+(j-params.delStep)] == min_cost)
            {
            for (a=1;a<=params.moveStep;a++)
            {
                i--;
            }
            for (a=1;a<=params.delStep;a++)
            {
                j--;
            }
            }
        else
            {
            for (a=1;a<=params.delStep;a++)
            {
                i--;
            }
            for (a=1;a<=params.moveStep;a++)
            {
                j--;
            }
            }
        }
    
    px[k] = i;
    py[k] = j;
    k++;      
    }

p->px = (int *) malloc (k * sizeof(int));
p->py = (int *) malloc (k * sizeof(int));
for (z1=0, z2=k-1; z1<k; z1++, z2--)
    {
    p->px[z1] = px[z2];
    p->py[z1] = py[z2];
    }
p->plen = k;

free(px);
free(py);

return 1;
}

// Compute the warp path starting at cost[startx, starty]
// If startx = -1 -> startx = n-1; if starty = -1 -> starty = m-1
int path(double *cost, int n, int m, int startx, int starty, DTW_path *p)
{
int i, j, k, z1, z2;
int *px;
int *py;
double min_cost;

if ((startx >= n) || (starty >= m))
    return 0;

if (startx < 0)
    startx = n - 1;

if (starty < 0)
    starty = m - 1;
    
i = startx;
j = starty;
k = 1;

// allocate path for the worst case
px = (int *) malloc ((startx+1) * (starty+1) * sizeof(int));
py = (int *) malloc ((startx+1) * (starty+1) * sizeof(int));

px[0] = i;
py[0] = j;

while ((i > 0) || (j > 0))
    {
    if (i == 0)
        j--;
    else if (j == 0)
        i--;
    else
        {
        min_cost = min3(cost[(i-1)*m+j],
                        cost[(i-1)*m+(j-1)], 
                        cost[i*m+(j-1)]);
        
        if (cost[(i-1)*m+(j-1)] == min_cost)
            {
            i--;
            j--;
            }
        else if (cost[i*m+(j-1)] == min_cost)
            j--;
        else
            i--;
        }
    
    px[k] = i;
    py[k] = j;
    k++;      
    }

p->px = (int *) malloc (k * sizeof(int));
p->py = (int *) malloc (k * sizeof(int));
for (z1=0, z2=k-1; z1<k; z1++, z2--)
    {
    p->px[z1] = px[z2];
    p->py[z1] = py[z2];
    }
p->plen = k;

free(px);
free(py);

return 1;
}


int path_11(double **cost, int n, int m)
{
int i, j, k, z1, z2;
int *px;
int *py;
double min_cost;

    
i =  n - 1;
j = m - 1;
k = 1;

while ((i > 0) || (j > 0))
    {
        if (i == 0)
        {
            j--;
        }
        else if (j == 0)
        {
            i--;
        }
        else
        {
            min_cost = min3(cost[i-1][j],
                            cost[i-1][j-1], 
                            cost[i][j-1]);
            
            if (cost[i-1][j-1] == min_cost)
            {
                i--;
                j--;
            }
            else if (cost[i][j-1] == min_cost)
            {
                j--;
            }
            else
            {
                i--;
            }
        }
        k++;      
    }

return k;
}

int path_12(double **cost, int n, int m)
{
int i, j, k, z1, z2;
int *px;
int *py;
double min_cost;

    
i =  n - 1;
j = m - 1;
k = 1;

while ((i > 0) || (j > 0))
    {
        if (i == 0)
        {
            j--;
        }
        else if (j == 0)
        {
            i--;
        }
        else if (i == 1)
            {
                min_cost = min(cost[i-1][j-1], cost[i-1][j-2]);
                
                if (cost[i-1][j-1] == min_cost)
                    {
                        i--;
                        j--;
                    }
                else
                {
                    j-=2;
                    i--;
                }
            }
        else if (j ==1)
            {
                min_cost = min(cost[i-2][j-1], cost[i-1][j-1]);
                
                if (cost[i-1][j-1] == min_cost)
                {
                    i--;
                    j--;
                }
                else
                {
                    i-=2;
                    j--;
                }
            }        
        
        else
            {
                min_cost = min3(cost[i-2][j-1],
                                cost[i-1][j-1], 
                                cost[i-1][j-2]);
                
                if (cost[i-1][j-1] == min_cost)
                {
                    i--;
                    j--;
                }
                else if (cost[i-1][j-2] == min_cost)
                {
                    j-=2;
                    i--;
                }
                else
                {
                    i-=2;
                    j--;
                }
            }
        k++;      
    }

return k;
}


double dist4Path(double *x, double*y, MatrixSize*size_x, MatrixSize*size_y, int NFeatDim, DTW_path* path_t, int path_len, Config* myConfig)
{
        // declarations of variables
        int i,j, NRows, NCols, ind_x, ind_y;	
        DistMethods myDistMethods[5]={NULL};
        double dist=0;
        float *a;
        a = (float*)x;
        
        //some initialization
        NRows=size_x->Nrow;
        NCols=size_y->Nrow;
        myDistMethods[Euclidean]=&EucDist;
        
        
        for (i=0;i<path_len;i++)
        {
        ind_x = path_t->px[i];
        ind_y = path_t->py[i];
        dist+=LocalDist(&x[ind_x*size_x->Ncol],&y[ind_y*size_y->Ncol], NFeatDim, myDistMethods[myConfig->DistMthd], myConfig->DistWghts);
        }
        return dist;
}

//This is a basic (standard) DTW measure without any constraint or fansy stuff. This is for one dimensional time series.
double dtw1d_std(double *x, double*y, int x_len, int y_len, double*cost, int dist_type)
{
        // declarations of variables
        int i,j;	
        DistMethods myDistMethods[5]={NULL};
        
        //setting up types of methods availale for measuring point to point distance
        myDistMethods[Euclidean]=&EucDist;
        myDistMethods[TMM_CityBlock]=&octBy2WrappedCitiblock;
        
        //Initializing the row and columns of cost matrix
        cost[0]= (*myDistMethods[dist_type])(x[0],y[0]);
        for (i=1;i<x_len;i++)
        {
        cost[i*y_len]=(*myDistMethods[dist_type])(x[i],y[0]) + cost[(i-1)*y_len];
        }
        for (j=1;j<y_len;j++)
        {
        cost[j]=(*myDistMethods[dist_type])(x[0],y[j]) + cost[(j-1)];
        }
        
        //filling in all the cumulative cost matrix
        for (i=1;i<x_len;i++)
        {
        for (j=1;j<y_len;j++)
        {
            cost[(i*y_len)+ j] = (*myDistMethods[dist_type])(x[i],y[j]) + 
                                    min3(cost[(i-1)*y_len+j], cost[((i-1)*y_len)+(j-1)], cost[(i*y_len)+(j-1)]);
        }
        
        }
        
        return cost[(x_len*y_len)-1];
}

//This is a constrained DTW, constraint is a band of certain width along the 45 degree line. (quite and standard constraint talked about by Sakoe, H. and Chiba)
double dtw1d_BandConstraint45(double *x, double*y, int x_len, int y_len, double*cost, int dist_type, int bandwidth)
{
        // declarations of variables
        int i,j;	
        DistMethods myDistMethods[5]={NULL};
        
        
        //CHANGES DUE TO CONSTRAINTS
        //the bandwidth of the constraint can't go beyong the abs(y_len-x_len), so
        bandwidth = max(bandwidth, abs(y_len-x_len)); // adapt constraint width
        //putting infi in all cost mtx
        for (i=0;i<x_len;i++)
        {
            for (j=0;j<y_len;j++)
            {
                cost[(i*y_len)+ j] = FLT_MAX;
            }
        
        }

        //setting up types of methods availale for measuring point to point distance
        myDistMethods[Euclidean]=&EucDist;
        
        //Initializing the row and columns of cost matrix
        cost[0]= (*myDistMethods[dist_type])(x[0],y[0]);
        for (i=1;i<x_len;i++)
        {
        cost[i*y_len]=(*myDistMethods[dist_type])(x[i],y[0]) + cost[(i-1)*y_len];
        }
        for (j=1;j<y_len;j++)
        {
        cost[j]=(*myDistMethods[dist_type])(x[0],y[j]) + cost[(j-1)];
        }
        
        //filling in all the cumulative cost matrix
        for (i=1;i<x_len;i++)
        {
        for (j=max(1, i-bandwidth);j<=min(y_len-1, i+bandwidth);j++)
        {
            cost[(i*y_len)+ j] = (*myDistMethods[dist_type])(x[i],y[j]) + 
                                    min3(cost[(i-1)*y_len+j], cost[((i-1)*y_len)+(j-1)], cost[(i*y_len)+(j-1)]);
        }
        
        }
        
        return cost[(x_len*y_len)-1];
}


double dtw1d_BandConst_LocalConst(double *x, double*y, int x_len, int y_len, double*cost, int dist_type, int bandwidth)
{
        // declarations of variables
        int i,j;    
        float min_vals; 
        DistMethods myDistMethods[5]={NULL};
        
        
        //CHANGES DUE TO CONSTRAINTS
        //the bandwidth of the constraint can't go beyong the abs(y_len-x_len), so
        bandwidth = max(bandwidth, abs(y_len-x_len)); // adapt constraint width
        //putting infi in all cost mtx
        for (i=0;i<x_len;i++)
        {
            for (j=0;j<y_len;j++)
            {
                cost[(i*y_len)+ j] = FLT_MAX;
            }
        
        }
        
        //setting up types of methods availale for measuring point to point distance
        myDistMethods[Euclidean]=&EucDist;
        myDistMethods[TMM_CityBlock]=&octBy2WrappedCitiblock;
        
        //Initializing the row and columns of cost matrix
        cost[0]= (*myDistMethods[dist_type])(x[0],y[0]);
        for (i=1;i<bandwidth+1;i++)
        {
            cost[i*y_len]=(*myDistMethods[dist_type])(x[i],y[0])  + cost[(i-1)*y_len];
        }
        for (j=1;j<bandwidth+1;j++)
        {
            cost[j]=(*myDistMethods[dist_type])(x[0],y[j]) + cost[j-1];
        }
        for (i=1;i<=bandwidth+1;i++)
        {
            j=1;
            cost[(i*y_len)+ j] = (*myDistMethods[dist_type])(x[i],y[j]) + min3(((i-1)*y_len)+ j, cost[((i-1)*y_len)+(j-1)], cost[((i)*y_len)+(j-1)]);
        }
        for (j=1;j<=bandwidth+1;j++)
        {
            i=1;
            cost[(i*y_len)+ j] = (*myDistMethods[dist_type])(x[i],y[j]) + min3(cost[((i-1)*y_len)+(j)], cost[((i-1)*y_len)+(j-1)], cost[((i)*y_len)+(j-1)]);
        }
        
        //filling in all the cumulative cost matrix
        for (i=2;i<x_len;i++)
        {
        for (j=max(2, i-bandwidth);j<=min(y_len-1, i+bandwidth);j++)
        {
            cost[(i*y_len)+ j] = (*myDistMethods[dist_type])(x[i],y[j]) + 
                                    min3(cost[(i-1)*y_len+(j-2)], cost[((i-1)*y_len)+(j-1)], cost[((i-2)*y_len)+(j-1)]);
        }
        
        }

        return cost[(x_len*y_len)-1];
        

}


double dtw1d_BandConst_LocalConst_Subsequence(double *x, double*y, int x_len, int y_len, double*cost, int dist_type, int bandwidth)
{
        // declarations of variables
        int i,j;    
        float min_vals; 
        DistMethods myDistMethods[5]={NULL};
        
        
        //CHANGES DUE TO CONSTRAINTS
        //the bandwidth of the constraint can't go beyong the abs(y_len-x_len), so
        bandwidth = max(bandwidth, abs(y_len-x_len)); // adapt constraint width
        //putting infi in all cost mtx
        for (i=0;i<x_len;i++)
        {
            for (j=0;j<y_len;j++)
            {
                cost[(i*y_len)+ j] = FLT_MAX;
            }
        
        }
        
        //setting up types of methods availale for measuring point to point distance
        myDistMethods[Euclidean]=&distEuclidean;
        myDistMethods[SqEuclidean]=&distSqEuclidean;
        myDistMethods[CityBlock]=&distCityBlock;
        myDistMethods[ShiftCityBlock]=&distShiftCityBlock;
        myDistMethods[TMM_CityBlock]=&distTMM_CityBlock;
        
        
        
        //Initializing the row and columns of cost matrix
        cost[0]= (*myDistMethods[dist_type])(x[0],y[0]);
        for (i=1;i<min(bandwidth+1, x_len);i++)
        {
            cost[i*y_len]=(*myDistMethods[dist_type])(x[i],y[0]);
        }
        for (j=1;j< min(bandwidth+1, y_len);j++)
        {
            cost[j]=(*myDistMethods[dist_type])(x[0],y[j]);
        }
        for (i=1;i<=min(bandwidth+1, x_len-1);i++)
        {
            j=1;
            cost[(i*y_len)+ j] = (*myDistMethods[dist_type])(x[i],y[j]) + min3(((i-1)*y_len)+ j, cost[((i-1)*y_len)+(j-1)], cost[((i)*y_len)+(j-1)]);
        }
        for (j=1;j<= min(bandwidth+1,y_len-1);j++)
        {
            i=1;
            cost[(i*y_len)+ j] = (*myDistMethods[dist_type])(x[i],y[j]) + min3(cost[((i-1)*y_len)+(j)], cost[((i-1)*y_len)+(j-1)], cost[((i)*y_len)+(j-1)]);
        }
        
        //filling in all the cumulative cost matrix
        for (i=2;i<x_len;i++)
        {
        for (j=max(2, i-bandwidth);j<=min(y_len-1, i+bandwidth);j++)
        {
            cost[(i*y_len)+ j] = (*myDistMethods[dist_type])(x[i],y[j]) + 
                                    min3(cost[(i-1)*y_len+(j-2)], cost[((i-1)*y_len)+(j-1)], cost[((i-2)*y_len)+(j-1)]);
        }
        
        }
        min_vals = FLT_MAX;
        for (i=x_len-1;i>=max(x_len-1-bandwidth, 0);i--)
        {
            j = y_len -1;
            if(cost[(i*y_len)+ j] < min_vals)
            {
                min_vals = cost[(i*y_len)+ j];
            }
        }
        for (j=y_len-1;j>=max(y_len-1-bandwidth,0);j--)
        {
            i = x_len -1;
            if(cost[(i*y_len)+ j] < min_vals)
            {
                min_vals = cost[(i*y_len)+ j];
            }
        } 
        return min_vals;
        

}



/*
 * This is an optimized DTW code with ability to band contraint the path. Additionally it has early abandoning included which is enabled if #define ENABLE_EA is there.
 * NOTE that to use this code pre initialization of cost to INF is required,  otherwise this code wont work
 * Inputs
 * x - time series (subsequence)
 * y - time series (subsequence)
 * x_len - length of x subsequence
 * y_len - length of y subsequence
 * cost - cost matrix. Memory should be initilized beforehand and should be initialized to INF
 * dist_type - type of distance to be used,  Euclidean is hardcoded for now. Also this euclidean doesn't have a square root calculation
 * bandwidth - band of the shakochiba constraint
 * bsf - best so far distance needed in the case of early abandoning
 * accLB - accumulated lower bound distance,  needed for the case of early abandoning. Something that improves early abandoning. 
  */
double dtw1dBandConst_old(double *x, double*y, int x_len, int y_len, double*cost, int dist_type, int bandwidth, double bsf, double *accLB)
{
        // declarations of variables
        int i,j, ind, overflow; 
        double min_vals, leftLB, temp;

#ifdef ENABLE_EA
        temp = bsf - accLB[y_len-1] ;
#endif        
        
        //Initializing the row and columns of cost matrix
        cost[0]= EucDist(x[0],y[0]);

        for (i=1;i<=bandwidth;i++)
        {
        cost[i*y_len]=EucDist(x[i],y[0]) + cost[(i-1)*y_len];
        }
        for (j=1;j<=bandwidth;j++)
        {
        cost[j]=EucDist(x[0],y[j]) + cost[(j-1)];
        }
        
        //filling in all the cumulative cost matrix
        for (i=1;i<x_len;i++)
        {
            
#ifdef ENABLE_EA
            leftLB = temp + accLB[i] ;
            overflow=1;
#endif                
            ind = i*y_len;
            
            for (j=max(1, i-bandwidth);j<min(y_len, i+bandwidth);j++)
            {
                min_vals = min3(cost[ind-y_len+j], cost[ind-y_len+(j-1)], cost[ind+(j-1)]);
                
#ifdef ENABLE_EA  
                if (min_vals > leftLB)
                {
                    cost[ind+ j]  =  FLT_MAX;
                }
                else
                {
                    cost[ind+ j] = EucDist(x[i],y[j]) + min_vals;
                }
                
                if (cost[ind+ j] < leftLB)
                {
                overflow=0;
                }

#else
                cost[ind+ j] = EucDist(x[i],y[j]) + min_vals;
#endif 

            }
        
#ifdef ENABLE_EA         
        if(overflow==1)
        {
            return FLT_MAX;
        }
#endif        
        
        }
        
        return cost[(x_len*y_len)-1];
}


//######################################## NEW SET OF FUNCTIONS TO BE CALLED IN C FUNCTIONS #####################################

double dtw1dBandConst(double *x, double*y, int x_len, int y_len, double**cost, int dist_type, int bandwidth, double bsf, double *accLB)
{
        // declarations of variables
        int i,j, ind, overflow; 
        double min_vals, leftLB, temp;

#ifdef ENABLE_EA
        temp = bsf - accLB[y_len-1] ;
#endif        
        
        //Initializing the row and columns of cost matrix
        cost[0][0]= mySimMeasure[dist_type](x[0],y[0]);

        for (i=1;i<=bandwidth;i++)
        {
        cost[i][0]=mySimMeasure[dist_type](x[i],y[0]) + cost[i-1][0];
        }
        for (j=1;j<=bandwidth;j++)
        {
        cost[0][j]=mySimMeasure[dist_type](x[0],y[j]) + cost[0][j-1];
        }
        
        //filling in all the cumulative cost matrix
        for (i=1;i<x_len;i++)
        {
            
#ifdef ENABLE_EA
            leftLB = temp + accLB[i] ;
            overflow=1;
#endif                
           
            for (j=max(1, i-bandwidth);j<=min(y_len-1, i+bandwidth);j++)
            {
                min_vals = min3(cost[i-1][j], cost[i-1][j-1], cost[i][j-1]);
                
#ifdef ENABLE_EA  
                if (min_vals > leftLB)
                {
                    cost[i][j]  =  FLT_MAX;
                }
                else
                {
                    cost[i][j] = mySimMeasure[dist_type](x[i],y[j]) + min_vals;
                }
                
                if (cost[i][j] < leftLB)
                {
                overflow=0;
                }

#else
                cost[i][j] = mySimMeasure[dist_type](x[i],y[j]) + min_vals;
#endif 

            }
        
#ifdef ENABLE_EA         
        if(overflow==1)
        {
            return FLT_MAX;
        }
#endif        
        
        }
        return cost[x_len-1][y_len-1];
}

double dtw1dBandConst_localConst(double *x, double*y, int x_len, int y_len, double**cost, int dist_type, int bandwidth, double bsf, double *accLB)
{
        // declarations of variables
        int i,j, ind, overflow; 
        double min_vals, leftLB, temp;
        
#ifdef ENABLE_EA
        temp = bsf - accLB[y_len-1] ;
#endif        
        
        //Initializing the row and columns of cost matrix
        cost[0][0]= mySimMeasure[dist_type](x[0],y[0]);

        for (i=1;i<=bandwidth;i++)
        {
        cost[i][0]=mySimMeasure[dist_type](x[i],y[0]) + cost[i-1][0];
        }
        for (j=1;j<=bandwidth;j++)
        {
        cost[0][j]=mySimMeasure[dist_type](x[0],y[j]) + cost[0][j-1];
        }
        
        for (i=1;i<=bandwidth+1;i++)
        {
            j=1;
            min_vals = min3(cost[i-1][j], cost[i-1][j-1], cost[i][j-1]);
            cost[i][j] = mySimMeasure[dist_type](x[i],y[j]) + min_vals;
        }
        for (j=1;j<=bandwidth+1;j++)
        {
            i=1;
            min_vals = min3(cost[i-1][j], cost[i-1][j-1], cost[i][j-1]);
            cost[i][j] = mySimMeasure[dist_type](x[i],y[j]) + min_vals;
        }
        
        //filling in all the cumulative cost matrix
        for (i=2;i<x_len;i++)
        {
            
#ifdef ENABLE_EA
            leftLB = temp + accLB[i] ;
            overflow=1;
#endif                
           
            for (j=max(2, i-bandwidth);j<=min(y_len-1, i+bandwidth);j++)
            {
                min_vals = min3(cost[i-2][j-1], cost[i-1][j-1], cost[i-1][j-2]);
                
#ifdef ENABLE_EA  
                if (min_vals > leftLB)
                {
                    cost[i][j]  =  FLT_MAX;
                }
                else
                {
                    cost[i][j] = mySimMeasure[dist_type](x[i],y[j]) + min_vals;
                }
                
                if (cost[i][j] < leftLB)
                {
                overflow=0;
                }

#else
                cost[i][j] = mySimMeasure[dist_type](x[i],y[j]) + min_vals;
#endif 

            }
        
#ifdef ENABLE_EA         
        if(overflow==1)
        {
            return FLT_MAX;
        }
#endif        
        
        }
        return cost[x_len-1][y_len-1];
}


double dtw1dBandConst_subsequence(double *x, double*y, int x_len, int y_len, double**cost, int dist_type, int bandwidth, double bsf, double *accLB)
{
        // declarations of variables
        int i,j, ind, overflow; 
        double min_vals, leftLB, temp;

#ifdef ENABLE_EA
        temp = bsf - accLB[y_len-1] ;
#endif        
        
        //Initializing the row and columns of cost matrix
        cost[0][0]= mySimMeasure[dist_type](x[0],y[0]);

        for (i=1;i<=bandwidth;i++)
        {
            cost[i][0]=mySimMeasure[dist_type](x[i],y[0]);
        }
        for (j=1;j<=bandwidth;j++)
        {
            cost[0][j]=mySimMeasure[dist_type](x[0],y[j]);
        }
        
        //filling in all the cumulative cost matrix
        for (i=1;i<x_len;i++)
        {
            
#ifdef ENABLE_EA
            leftLB = temp + accLB[i] ;
            overflow=1;
#endif                
           
            for (j=max(1, i-bandwidth);j<=min(y_len-1, i+bandwidth);j++)
            {
                min_vals = min3(cost[i-1][j], cost[i-1][j-1], cost[i][j-1]);
                
#ifdef ENABLE_EA  
                if (min_vals > leftLB)
                {
                    cost[i][j]  =  FLT_MAX;
                }
                else
                {
                    cost[i][j] = mySimMeasure[dist_type](x[i],y[j]) + min_vals;
                }
                
                if (cost[i][j] < leftLB)
                {
                overflow=0;
                }

#else
                cost[i][j] = mySimMeasure[dist_type](x[i],y[j]) + min_vals;
#endif 

            }
        
#ifdef ENABLE_EA         
        if(overflow==1)
        {
            return FLT_MAX;
        }
#endif        
        
        }
        min_vals = FLT_MAX;
        for (i=x_len-1;i>=x_len-1-bandwidth;i--)
        {
            if(cost[i][y_len-1] < min_vals)
            {
                min_vals = cost[i][y_len-1];
            }
        }
        for (i=y_len-1;i>=y_len-1-bandwidth;i--)
        {
            if(cost[x_len-1][i] < min_vals)
            {
                min_vals = cost[x_len-1][i];
            }
        }        

        
        
        return min_vals;
}


double dtw1dBandConst_subsequence_localConst(double *x, double*y, int x_len, int y_len, double**cost, int dist_type, int bandwidth, double bsf, double *accLB)
{
        // declarations of variables
        int i,j, ind, overflow; 
        double min_vals, leftLB, temp;

#ifdef ENABLE_EA
        temp = bsf - accLB[y_len-1] ;
#endif        
        
        //Initializing the row and columns of cost matrix
        cost[0][0]= mySimMeasure[dist_type](x[0],y[0]);

        for (i=1;i<=bandwidth;i++)
        {
            cost[i][0]=mySimMeasure[dist_type](x[i],y[0]);
        }
        for (j=1;j<=bandwidth;j++)
        {
            cost[0][j]=mySimMeasure[dist_type](x[0],y[j]);
        }
        
        for (i=1;i<=bandwidth+1;i++)
        {
            j=1;
            min_vals = min3(cost[i-1][j], cost[i-1][j-1], cost[i][j-1]);
            cost[i][j] = mySimMeasure[dist_type](x[i],y[j]) + min_vals;
        }
        for (j=1;j<=bandwidth+1;j++)
        {
            i=1;
            min_vals = min3(cost[i-1][j], cost[i-1][j-1], cost[i][j-1]);
            cost[i][j] = mySimMeasure[dist_type](x[i],y[j]) + min_vals;
        }
        
        //filling in all the cumulative cost matrix
        for (i=2;i<x_len;i++)
        {
            
#ifdef ENABLE_EA
            leftLB = temp + accLB[i] ;
            overflow=1;
#endif                
           
            for (j=max(2, i-bandwidth);j<=min(y_len-1, i+bandwidth);j++)
            {
                min_vals = min3(cost[i-2][j-1], cost[i-1][j-1], cost[i-1][j-2]);
                
#ifdef ENABLE_EA  
                if (min_vals > leftLB)
                {
                    cost[i][j]  =  FLT_MAX;
                }
                else
                {
                    cost[i][j] = mySimMeasure[dist_type](x[i],y[j]) + min_vals;
                }
                
                if (cost[i][j] < leftLB)
                {
                overflow=0;
                }

#else
                cost[i][j] = mySimMeasure[dist_type](x[i],y[j]) + min_vals;
#endif 

            }
        
#ifdef ENABLE_EA         
        if(overflow==1)
        {
            return FLT_MAX;
        }
#endif        
        
        }
        min_vals = FLT_MAX;
        for (i=x_len-1;i>=x_len-1-bandwidth;i--)
        {
            if(cost[i][y_len-1] < min_vals)
            {
                min_vals = cost[i][y_len-1];
            }
        }
        for (i=y_len-1;i>=y_len-1-bandwidth;i--)
        {
            if(cost[x_len-1][i] < min_vals)
            {
                min_vals = cost[x_len-1][i];
            }
        }        

        
        
        return min_vals;
}





/*#####################################################################################################################
################################################ LOWER BOUND FUNCTIONS ##############################################
#####################################################################################################################*/
double computeLBkimFL(double a1, double a2, double b1, double b2, int dist_type)
{
    double distF, distL;
    distF = mySimMeasure[dist_type](a1, a2);
    distL = mySimMeasure[dist_type](b1, b2);
    return (distF + distL);
}

/*
 * This function computes Keogh lower bound of DTW distance,  as an input you have to precompute
 * U - running max of one of the subsequence
 * L - running min of one of the subsequence
 * accLB - accumulated lower bound at each index,  incrementally computed
 * data - subsequence to be used for the computation of lower bound
 * lenMotif - length of the subsequence
 * bsf - best so far distance,  which is kind of the threshold for every computation.
 * 
 * Additionally this function also has a early abandoning step,  which is dependent on the preprocessor #define. Use 'ENABLE_EA' to enable it 
*/
double computeKeoghsLB(double *U, double *L, double * accLB, double *data,int lenMotif, double bsf, int dist_type)
{
    int ii;
    double sum=0;
    for(ii=0;ii<lenMotif;ii++)
    {
        if (data[ii]>U[ii])
        {
            sum+=mySimMeasure[dist_type](data[ii],U[ii]);
        }
        else if (data[ii]<L[ii])
        {
            sum+=mySimMeasure[dist_type](data[ii],L[ii]);
        }
        accLB[ii] = sum;
#ifdef ENABLE_EA             
        if (sum>bsf)
        {
            return FLT_MAX;
        }
#endif        
    }
    return sum;
    
}


//Probably not functinally correct functions but kept to just have as a backup


/*
double computeLBkimFL_extended(double *U1, double *L1, double *data1, double *U2, double *L2, double *data2, int lenMotif)
{
    double sum1=0, sum2=0;
    int ind=0;
    
    // data1 versus envelope of data2
    ind = 0;
    if (data1[ind]>U2[ind])
    {
        sum1+=EucDist(data1[ind],U2[ind]);
    }
    else if (data1[ind]<L2[ind])
    {
        sum1+=EucDist(data1[ind],L2[ind]);
    }
    
    ind = lenMotif-1;
    if (data1[ind]>U2[ind])
    {
        sum1+=EucDist(data1[ind],U2[ind]);
    }
    else if (data1[ind]<L2[ind])
    {
        sum1+=EucDist(data1[ind],L2[ind]);
    }
    
    // data2 versus envelope of data1
    ind = 0;
    if (data2[ind]>U1[ind])
    {
        sum1+=EucDist(data2[ind],U1[ind]);
    }
    else if (data2[ind]<L1[ind])
    {
        sum1+=EucDist(data2[ind],L1[ind]);
    }
    
    ind = lenMotif-1;
    if (data2[ind]>U1[ind])
    {
        sum1+=EucDist(data2[ind],U1[ind]);
    }
    else if (data2[ind]<L1[ind])
    {
        sum1+=EucDist(data2[ind],L1[ind]);
    }
    
    return max(sum1,sum2);
}

*/



