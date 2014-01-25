/*
 * This implementation of running min and max is proposed by "Daniel Lemire" in 
 * STREAMING MAXIMUM-MINIMUM FILTER USING NO
MORE THAN THREE COMPARISONS PER ELEMENT
 * This is supposed to be quite a fast implementation of running min and max
 *
 * Author: Sankalp Gulati
 * No part of this code is taken from anywhere and is solely based on the pseudocode in above paper
 * NOTE: this code can be made memory efficient by using push pop kind of things rather than allocating big memory chunk.
 * 
 * I made it computationally efficient without worrying much about memory inefficiency
 * Note that min and max is returned for a range of 2*winLen+1 indices
 */

#include <stdio.h>
#include <stdlib.h>

int computeRunningMinMax(double *data, double *U, double *L, int lenData, int winLen)
{
    
    int UstartInd, UnumPoints, LstartInd, LnumPoints, ii;
    int *utemp, *ltemp;
    
    
    //without worrying much about memory just allocating memory for temporary monotonic wedges used for this implementation
    utemp = (int*)malloc(sizeof(int)*lenData);
    ltemp = (int*)malloc(sizeof(int)*lenData);
    
    //Initializing
    UstartInd=0;
    utemp[UstartInd] = 0;
    UnumPoints=1;
    
    LstartInd=0;
    ltemp[LstartInd] = 0;
    LnumPoints=1;
    
    for(ii=1;ii<lenData;ii++)
    {
        if (ii>winLen)
        {
            U[ii-winLen-1] = data[utemp[UstartInd]];
            L[ii-winLen-1] = data[ltemp[LstartInd]];
        }
        
        if(data[ii]>data[ii-1])
        {
            UnumPoints-=1;
            
            while ((UnumPoints>0)&&(data[ii]>data[utemp[UstartInd+UnumPoints-1]]))
            {
                UnumPoints-=1;
            }
        }
        else
        {
            LnumPoints-=1;
            
            while ((LnumPoints>0)&&(data[ii]< data[ltemp[LstartInd+LnumPoints-1]]))
            {
                LnumPoints-=1;
            }
            
        }
        UnumPoints+=1;
        LnumPoints+=1;
        utemp[UstartInd+UnumPoints-1]=ii;
        ltemp[LstartInd+LnumPoints-1]=ii;
        
        if (ii == 2*winLen+1 + utemp[UstartInd])
        {
            UstartInd+=1;
            UnumPoints-=1;
        }        
        else if (ii == 2*winLen+1 + ltemp[LstartInd])
        {
            LstartInd+=1;
            LnumPoints-=1;
        }
            
        
    }
    
    for (ii = lenData; ii < lenData+winLen+1; ii++) 
    {

        U[ii-winLen-1] = data[utemp[UstartInd]];
        L[ii-winLen-1] = data[ltemp[LstartInd]];

        if (ii-utemp[UstartInd] >= 2*winLen+1)           
        {
            UstartInd+=1;
            UnumPoints-=1;
        }

        if (ii-ltemp[LstartInd] >= 2*winLen+1)           
        {
            LstartInd+=1;
            LnumPoints-=1;
        }
    }
    
    free(utemp);
    free(ltemp);
    
    return 1;
    
}