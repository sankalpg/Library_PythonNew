
/*******************************************************************************
****** Sankalp Gulati                                                   *******
****** MUSIC TECHNOLOGY GROUP, UPF, BARCELONA                           *******
******                                                                  *******
*******************************************************************************
******************************************************************************/


#include "MotifDataIO.h"

int comparePairs(const void *a, const void *b)
{
    if (((patternDist_t*)a)->dist > ((patternDist_t*)b)->dist)
    {
        return 1;
    }
    else if (((patternDist_t*)a)->dist < ((patternDist_t*)b)->dist)
    {
        return -1;
    }
    return 0;
}


int main( int argc , char *argv[])
{
    FILE *fp1;
    char *patternInfoExt, *patternKNNExt, *filelistFilename, *blackListPatternExt, searchFileNames[2000][400] = {'\0'}, tempFilename[400]= {'\0'},  blackListFile[400]={'\0'};
    int verbos =1, ii, jj, mm, nPatterns, NFilesSearch, NPairs;
    patternInfo_t *patternInfo;
    patternDist_t *patternDist;
    int *blackListArray;
    float timeMax=0;
    float resolution = 0.01;
    int nTimeSamples, str, end, str_old, end_old, hasOverlap, OverlapIndex;
    int *overlapArray, err;
    


   if(argc < 5 || argc > 6)
    {
        printf("\nInvalid number of arguments!!!\n");
        exit(1);
    }
    
    patternInfoExt = argv[1];
    patternKNNExt = argv[2];
    filelistFilename = argv[3];
    blackListPatternExt = argv[4];

    if( argc == 6 ){verbos = atoi(argv[5]);}


    fp1 = fopen(filelistFilename, "r");
    ii=0;
    while(fgets(tempFilename, 400, fp1))
    {
        sscanf(tempFilename, "%[^\n]s\n", searchFileNames[ii]);
        ii++;
        
    }
    fclose(fp1);
    NFilesSearch = ii;
    memset(tempFilename, '\0', sizeof(char)*400);

    for (ii=0;ii<NFilesSearch;ii++)
    {
        strcat(tempFilename,searchFileNames[ii]);
        strcat(tempFilename,patternInfoExt);
        err =  readPatternDump(tempFilename, &patternInfo, &nPatterns);
        memset(tempFilename, '\0', sizeof(char)*400);

        strcat(tempFilename,searchFileNames[ii]);
        strcat(tempFilename,patternKNNExt);
        err = readKNNDump(tempFilename, &patternDist, &NPairs);
        memset(tempFilename, '\0', sizeof(char)*400);

        if ( nPatterns != NPairs)
        {
            if (verbos ==1)
            {
                printf("There is some problem with the info and knn file for the file %s\n",searchFileNames[ii]);
            }
        }

        blackListArray = (int *)malloc(sizeof(int)*nPatterns);
        for(jj=0; jj< nPatterns; jj++)
        {
            blackListArray[jj]=1; //by default all of them are blacklisted
        }
        
        //sorting the pairs according to ascending order of the distance values
        //qsort (patternDist, nPatterns, sizeof(patternDist_t), comparePairs);

        // finding out the longest time stamp in the current file
        for(jj=0; jj< nPatterns; jj++)
        {
            if (patternInfo[jj].end > timeMax)
            {
                timeMax = patternInfo[jj].end;
            }
        }
        
        //allocating array for storing ovrlapping history
        nTimeSamples = (int)ceil(timeMax/resolution);
        
        overlapArray = (int*)malloc(sizeof(int)*nTimeSamples);
        for(jj=0; jj< nTimeSamples; jj++)
        {
            overlapArray[jj]=-1;
        }
        
        for (jj=0;jj<nPatterns;jj++)
        {
            hasOverlap=0;
            OverlapIndex=-1;
            if (patternInfo[jj].id!=patternDist[jj].patternID1)
            {
                if (verbos==1)
                {
                    printf("There is a mismatch in Info and Dist file for filename %s\n",searchFileNames[ii]);
                }
                break;
            }
                
                
            str = (int)floor(patternInfo[jj].str/resolution);
            end = (int)floor(patternInfo[jj].end/resolution);
            
            // first check if there is an overlap
            for(mm=str;mm<=end; mm++)
            {
                if (overlapArray[mm]!=-1)
                {
                    hasOverlap =1;
                    OverlapIndex = overlapArray[mm];
                    break;
                }
            }
            
            if (hasOverlap==1)
            {
                if (patternDist[OverlapIndex].dist > patternDist[jj].dist)
                {
                    //this means that the current pattern is better and so the older one should be removed
                    str_old = (int)floor(patternInfo[OverlapIndex].str/resolution);
                    end_old = (int)floor(patternInfo[OverlapIndex].end/resolution);
                    for(mm=str_old;mm<=end_old; mm++)
                    {
                        overlapArray[mm] = -1;
                    }
                    hasOverlap=0;
                    blackListArray[OverlapIndex]=1;
                }
                //the other case where the older pattern is the better one leave thigns as it is. Since hasOverlap flag is 1 there won't be any action
                
            }
            if (hasOverlap==0)
            {
                for(mm=str;mm<=end; mm++)
                {
                    overlapArray[mm] = jj;
                }
                blackListArray[jj]=0;
                
            }
          
        }
        
        strcat(blackListFile,searchFileNames[ii]);
        strcat(blackListFile,blackListPatternExt);
        fp1 = fopen(blackListFile, "w");
        for(jj=0;jj<nPatterns;jj++)
        {
            fprintf(fp1, "%d\n",blackListArray[jj]);
        }
        fclose(fp1);
        memset(blackListFile, '\0', sizeof(char)*400);
     
    }
}
