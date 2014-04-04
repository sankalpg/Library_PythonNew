/*
 * Author: Sankalp gulati
 * email: sankalp.gulati@upf.edu.com
 * Affiliation: Universitat Pompeu Fabra
 *
 *THIS IS A MEX WRAPPER FOR SEVERAL DTW FUNCTIONS
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include "mex.h"
#include "matrix.h"
#include "../dtw/dtw.h"

/*Enable this if you are compiling on old mac*/
//#define char16_t uint16_t

/*#defs*/
#define  N_DTW_VARIANTS 6

/*typedefs*/
typedef double (*dtwVariants)(double*, double*, int, int, double*, int, int);



dtwVariants myDtwVariants[N_DTW_VARIANTS] = {   dtw1d_BandConst_LocalConst_Subsequence,
dtw1d_BandConst_LocalConst_Subsequence,
dtw1d_BandConst_LocalConst_Subsequence,
dtw1d_BandConst_LocalConst_Subsequence,
dtw1d_BandConst_LocalConst_Subsequence,
dtw1d_BandConst_LocalConst_Subsequence};

void mexFunction ( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    /*
     *  USAGE:
     *  <filename>(x, y, distType, dtwType, bandWidth, cost[optional])
     *
     *  x = one dimensional time series of length N (can be treated as query)
     *  y = one dimensional time series of length M (can be treated as reference)

     *  Params are: 
     *  
     *  dtwParams.distType (default 0)-->Distance type used in dtw
     *              OPTIONS
     *              # == distance type
     *              0 ==Euclidean;
     *              1 = SqEuclidean; (faster than Euclidean)
     *              2 = CityBlock
     *              3 = ShiftCityBlock
     *              4 = TMM_CityBlock (Customized city block for TMM)
     *
     *  dtwParams.hasGlobalConst (default 0)-->Flag selecting whether or not to apply global constraint
     *  
     *  dtwParams.globalType (default 0)--> type of global constraint
     *  
     *  dtwParams.bandwidth (default is 10% of query)-->bandwidth of the global constraint in samples
     *  
     *  dtwParams.initCostMtx (default 1) --> whether or not to initialize cost matrix inside dtw. 
     *          This can be very handy tool to apply constraint from outside, just make a design 
     *          and pass it to the function. Also one can do iterative stuff with this.
     *
     *  dtwParams.reuseCostMtx (default 0) --> if enabled will not recompute the cost matrix, 
     *          instead use the one provided directly
     *
     *  dtwParams.delStep (default 1) --> step size for local consraint (horizontal)
     *
     *  dtwParams.moveStep (default 1) --> step size for local consraint (verticle)
     *
     *  dtwParams.diagStep (default 1) --> step size for local consraint (diagonal value)
     *
     *  dtwParams.initFirstCol (default 1) -->whether to initialize the fist column of cost or not!!.
     *
     *  dtwParams.isSubsequence (default 0) --> whether or not its a subsequence dtw variant

     *  cost [optional] = 2d arary of NxM. if provided used as cost matrix for DTW computation.
     *                    Providing this matrix avoids allocation of memory. So if you are using DTW inside a loop with
     *                    same duration time series (NxM) across iterations, allocate it once and just pass it.
     *
     *
     */
    
    
    double *inp1, *inp2, *costPtr, dist , *pathOutPtr, min_val;
    mxArray *cost, *tmp;
    int distType, dtwType, bandWidth, ii, jj, min_x, min_y;
    mwSize len1, len2;
    DTW_path path_t;
    mwSize     NStructElems;
    int nfields, ifield;
    const char *fieldName;       /* pointers to field names */
    double fieldVal; 
    mwSize     sizebuf;
    dtwParams_t params;
    
    
    /*fetching pointers for input data*/
    inp1 = mxGetPr(prhs[0]);
    inp2 = mxGetPr(prhs[1]);
    
    /*Computing lengths of the input vectors*/
    len1 = mxGetDimensions(prhs[0])[0];
    len2 = mxGetDimensions(prhs[1])[0];
    
    if ((nrhs<3)|(nrhs>4))
    {
        printf("Please provide right number of parameters\n");
        return ;
    }
    else if (nrhs==3)
    {
        //allocate memory for cost matrix
        /*NOTE that the mateix is (len2xlen1) because matlab interprets colxrow as 1D unline cython*/
        cost = mxCreateDoubleMatrix(len2,len1, mxREAL);
        costPtr = mxGetPr(cost);
    }
    else if (nrhs==4)
    {
        costPtr = mxGetPr(prhs[3]);
    }

    nfields = mxGetNumberOfFields(prhs[2]);
    NStructElems = mxGetNumberOfElements(prhs[2]);
    
    if(NStructElems>1)
    {
        printf("something is wrong in the way you have passed dtw param structure");
    }
    
    /*setting default values of the dtw parameters*/
    params.distType = 0;
    params.hasGlobalConst = 0;
    params.globalType = 0;
    params.bandwidth = (int)(0.1*len1);
    params.initCostMtx = 1;
    params.reuseCostMtx = 0;
    params.delStep = 1;
    params.moveStep = 1;
    params.diagStep = 1;
    params.initFirstCol = 1;
    params.isSubsequence = 0;
    
    
    /* Fetching dtw parameters */
    for (ifield=0; ifield< nfields; ifield++)
    {
        fieldName = mxGetFieldNameByNumber(prhs[2],ifield);
        tmp =  mxGetFieldByNumber(prhs[2],0,ifield);
        sizebuf = mxGetElementSize(tmp);
        memcpy(&fieldVal, mxGetData(tmp), sizebuf);
        
        if (strcmp(fieldName, "distType")==0)
        {
            //printf("%s\t%f\n",fieldName, fieldVal);
            params.distType = (int)fieldVal;
        }
        else if (strcmp(fieldName, "hasGlobalConst")==0)
        {
            //printf("%s\t%f\n",fieldName, fieldVal);
            params.hasGlobalConst = (int)fieldVal;
        }
        else if (strcmp(fieldName, "globalType")==0)
        {
            //printf("%s\t%f\n",fieldName, fieldVal);
            params.globalType = (int)fieldVal;
        }
        else if (strcmp(fieldName, "bandwidth")==0)
        {
            //printf("%s\t%f\n",fieldName, fieldVal);
            params.bandwidth =(float)fieldVal;
        }
        else if (strcmp(fieldName, "initCostMtx")==0)
        {
            //printf("%s\t%f\n",fieldName, fieldVal);
            params.initCostMtx = (int)fieldVal;
        }
        else if (strcmp(fieldName, "reuseCostMtx")==0)
        {
            //printf("%s\t%f\n",fieldName, fieldVal);
            params.reuseCostMtx = (int)fieldVal;
        }
        else if (strcmp(fieldName, "delStep")==0)
        {
            //printf("%s\t%f\n",fieldName, fieldVal);
            params.delStep = (int)fieldVal;
        }
        else if (strcmp(fieldName, "moveStep")==0)
        {
            //printf("%s\t%f\n",fieldName, fieldVal);
            params.moveStep = (int)fieldVal;
        }
        else if (strcmp(fieldName, "diagStep")==0)
        {
            //printf("%s\t%f\n",fieldName, fieldVal);
            params.diagStep = (int)fieldVal;
        }
        else if (strcmp(fieldName, "initFirstCol")==0)
        {
            //printf("%s\t%f\n",fieldName, fieldVal);
            params.initFirstCol = (int)fieldVal;
        }
        else if (strcmp(fieldName, "isSubsequence")==0)
        {
            //printf("%s\t%f\n",fieldName, fieldVal);
            params.isSubsequence = (int)fieldVal;
        }
    }


    //dist = myDtwVariants[dtwType](inp1, inp2, len1, len2, costPtr, distType, bandWidth);
    /*printf(" params.distType %d\n", params.distType);
    printf(" params.hasGlobalConst %d\n", params.hasGlobalConst);
    printf(" params.globalType %d\n", params.globalType);
    printf(" params.bandwidth %d\n", params.bandwidth);
    printf(" params.initCostMtx %d\n", params.initCostMtx);
    printf(" params.reuseCostMtx %d\n", params.reuseCostMtx);
    printf(" params.delStep %d\n", params.delStep);
    printf(" params.moveStep %d\n", params.moveStep);
    printf(" params.diagStep %d\n", params.diagStep);
    printf(" params.initFirstCol %d\n", params.initFirstCol);
    printf(" params.isSubsequence %d\n", params.isSubsequence);*/
    dist = dtw_GLS(inp1, inp2, len1, len2, costPtr, params);
    
    
    if (nlhs == 1)
    {
        //user has only asked for distance
        plhs[0] =  mxCreateDoubleScalar(dist);
        return;
    }
    else
    {
        //user is greedy he askss for more stuff than just the distance, well, lets make him happy :|
        plhs[0] =  mxCreateDoubleScalar(dist);
        min_x=-1;
        min_y=-1;
        if (params.isSubsequence)
        {
            /*In subsequence dtw we have to search for the path end you remember right?*/
            min_val=100000000000000000000.0;
            jj = len2-1;
            for (ii=0;ii<len1;ii++)
            {
                if(costPtr[ii*len2+jj]< min_val)
                {
                    min_val = costPtr[ii*len2+jj];
                    min_x = ii;
                    min_y = jj;
                }
            }
            ii = len1-1;
            for (jj=0;jj<len2;jj++)
            {
                if(costPtr[ii*len2+jj] <min_val)
                {
                    min_val = costPtr[ii*len2+jj];
                    min_x = ii;
                    min_y = jj;
                }
            }
        }
        
        //This means that dtw is not a subsequence version, in which case we start backtracking from the n,m of cost
        
        path(costPtr, len1, len2, min_x, min_y, &path_t);

        plhs[1] =  mxCreateDoubleScalar(path_t.plen);
        
        if (nlhs==2)
        {
            return;
        }
        else
        {
            plhs[2] = mxCreateDoubleMatrix(path_t.plen,2, mxREAL);
            pathOutPtr = mxGetPr(plhs[2]);
            
            for(ii=0;ii<path_t.plen;ii++)
            {
                pathOutPtr[ii] = (double)path_t.px[ii];
                pathOutPtr[ii+path_t.plen] = (double)path_t.py[ii];
            
            }
            /*free the memory for path, that was allocated inside the path function*/
            free(path_t.px);
            free(path_t.py);

            if(nlhs==3)
            {
                
                return;
            }
            else
            {
                if (nrhs==4)
                {
                    plhs[3]=prhs[3];
                }
                else
                {
                    plhs[3]=cost;
                    return;
                }
                
                
                
            }
            
            
        }
        
    }
    

}
