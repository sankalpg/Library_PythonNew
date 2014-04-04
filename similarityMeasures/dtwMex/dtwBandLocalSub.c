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
     *  distType =  # == distance type
     *              0 ==Euclidean;
     *              1 = SqEuclidean; (faster than Euclidean)
     *              2 = CityBlock
     *              3 = ShiftCityBlock
     *              4 = TMM_CityBlock (Customized city block for TMM)
     *
     *  dtwType =   # ==    Global  Local   Subsequence
     *              0 ==    No      No      No
     *              1 ==    Yes     No      No
     *              2 ==    No      Yes     No
     *              3 ==    Yes     Yes     No
     *              4 ==    No      Yes     Yes (with sub version it always make sense to have local alignment)
     *              5 ==    Yes     Yes     Yes (with sub version it always make sense to have local alignment)
     *
     *  bandWidth = -1 whenever Global constraint is No, else band-width of the Global constraint
     *
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
    dtwParams_t params;
    mwSize     NStructElems;
    int nfields, ifield;
    mwIndex    jstruct;
    mxClassID  *classIDflags;
    const char **fnames;       /* pointers to field names */
    double *fvals; 
    
    
    /*fetching pointers for input data*/
    inp1 = mxGetPr(prhs[0]);
    inp2 = mxGetPr(prhs[1]);
    distType = mxGetScalar(prhs[2]);
    dtwType = mxGetScalar(prhs[3]);
    bandWidth = mxGetScalar(prhs[4]);
    
    /*Computing lengths of the input vectors*/
    len1 = mxGetDimensions(prhs[0])[0];
    len2 = mxGetDimensions(prhs[1])[0];
    
    if (nrhs==6)
    {
        //Cost matrix is also passed into the function
        costPtr = mxGetPr(prhs[5]);
    }
    else if (nrhs < 6)
    {
        //allocate memory for cost matrix
        /*NOTE that the mateix is (len2xlen1) because matlab interprets colxrow as 1D unline cython*/
        cost = mxCreateDoubleMatrix(len2,len1, mxREAL);
        costPtr = mxGetPr(cost);
    }
    else if (nrhs >6)
    {
        
        nfields = mxGetNumberOfFields(prhs[6]);
        NStructElems = mxGetNumberOfElements(prhs[6]);
        
        
        
        fnames = mxCalloc(nfields, sizeof(*fnames));
        fvals = mxCalloc(nfields, sizeof(double));
        mwSize     sizebuf;
        double temp;
                
                
        /* get field name pointers */
        for (ifield=0; ifield< nfields; ifield++)
        {
            fnames[ifield] = mxGetFieldNameByNumber(prhs[6],ifield);
            tmp =  mxGetFieldByNumber(prhs[6],0,ifield);
            sizebuf = mxGetElementSize(tmp);
            memcpy(&temp, mxGetData(tmp), sizebuf);
            printf("%f\n",temp);
        }
        
        for (ifield=0; ifield< nfields; ifield++)
        {
            //printf("%s\t%f\n",fnames[ifield], fvals[ifield]);
        }
        plhs[0] =  mxCreateDoubleScalar(-1);
        return ;
    }
    
    /*Other useful check*/
    if (dtwType%2 ==0)
    {
        //This means Global constraint has to be applied, which means user better pass bandWidth as sensible number
        if (bandWidth==-1)
        {
            plhs[0] =  mxCreateDoubleScalar(-1);
            return ;
        }
    }
    
    /*Allocating cost matrix for dtw computation*/
    
    /*Computing distance using DTW function*/
    /*params.distType=distType;
    params.globalType=0;
    params.bandwidth=bandWidth;
    params.initCostMtx=1;
    params.reuseCostMtx=0;
    params.delStep=2;
    params.moveStep=1;
    params.diagStep=1;
    params.initFirstCol=1;
    params.hasGlobalConst=1;
    if (dtwType>=4)
    {
        params.isSubsequence=1;
    }*/
    dist = myDtwVariants[dtwType](inp1, inp2, len1, len2, costPtr, distType, bandWidth);
    /*dist = dtw_GLS(inp1, inp2, len1, len2, costPtr, params);*/
    
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
        if (dtwType>=4)
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
                if (nrhs==6)
                {
                    plhs[3]=prhs[5];
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
