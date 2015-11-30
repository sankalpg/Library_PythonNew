#include "TSAlogs.h"

TSAlogs::TSAlogs()
{
    procLogs.tLoadTS = 0;
    procLogs.tProcTS = 0;
    procLogs.tGenEnv = 0; 
    procLogs.tGenUniScal = 0;
    procLogs.tTotal = 0;
    procLogs.tDump = 0;
    
    procLogs.lenTS = 0;
    procLogs.lenProcTS = 0;
    procLogs.nSubSeqs = 0;
    procLogs.nSubSeqsBL = 0;
    procLogs.nProcSubSeq = 0;
    
    procLogs.nLB_KIM_FL = 0;
    procLogs.nLB_Keogh_EQ = 0;
    procLogs.nLB_Keogh_EC = 0;
    procLogs.nDTW_EA = 0;
    procLogs.nPriorityUpdates = 0;
}

int  TSAlogs::dumpProcLogs(char *logFile, int verbos)
{
    FILE *fp;
    
    fp =fopen(logFile,"w");
    fprintf(fp, "\n#################### TIME RELATED STATS ####################\n");
    fprintf(fp, "Time taken to load the time series :\t%f\n", procLogs.tLoadTS);
    fprintf(fp, "Time taken to process time series:\t%f\n", procLogs.tProcTS);
    fprintf(fp, "Time taken to generate envelops:\t%f\n", procLogs.tGenEnv);
    fprintf(fp, "Time taken to uniformly scale subsequences:\t%f\n", procLogs.tGenUniScal);
    fprintf(fp, "Time taken to write data:\t%f\n", procLogs.tDump);
    fprintf(fp, "Total time taken by the process:\t%f\n", procLogs.tTotal);
    
    fprintf(fp, "\n#################### DATA POINTS RELATED STATS ####################\n");
    fprintf(fp, "Total number of time series samples:\t%lld\n", procLogs.lenTS);
    fprintf(fp, "Total number of time series samples after processing:\t%lld\n", procLogs.lenProcTS);
    fprintf(fp, "Total number of generated subsequences:\t%lld\n", procLogs.nSubSeqs);
    fprintf(fp, "Total number of subsequences blacklisted:\t%lld\n", procLogs.nSubSeqsBL);
    fprintf(fp, "Total number of subsequences after scaling:\t%lld\n", procLogs.nProcSubSeq);
    
    fprintf(fp, "\n#################### FNC CALLS RELATED STATS ####################\n");
    fprintf(fp, "Number of FL lowerbound is computed:\t%lld\n", procLogs.nLB_KIM_FL);
    fprintf(fp, "Number of LB_Keogh_EQ computed:\t%lld\n", procLogs.nLB_Keogh_EQ);
    fprintf(fp, "Number of LB_Keogh_EC computed:\t%lld\n", procLogs.nLB_Keogh_EC);
    fprintf(fp, "Number of times DTW computed:\t%lld\n", procLogs.nDTW_EA);
    fprintf(fp, "Number of updates of priority list:\t%lld\n", procLogs.nPriorityUpdates);
    fclose(fp);
    
    
    if (verbos)
    {
        printf("\n#################### TIME RELATED STATS ####################\n");
        printf("Time taken to load the time series :\t%f\n", procLogs.tLoadTS);
        printf("Time taken to process time series:\t%f\n", procLogs.tProcTS);
        printf("Time taken to generate envelops:\t%f\n", procLogs.tGenEnv);
        printf("Time taken to uniformly scale subsequences:\t%f\n", procLogs.tGenUniScal);
        printf("Time taken to write data:\t%f\n", procLogs.tDump);
        printf("Total time taken by the process:\t%f\n", procLogs.tTotal);
        
        printf("\n#################### DATA POINTS RELATED STATS ####################\n");
        printf("Total number of time series samples:\t%lld\n", procLogs.lenTS);
        printf("Total number of time series samples after processing:\t%lld\n", procLogs.lenProcTS);
        printf("Total number of generated subsequences:\t%lld\n", procLogs.nSubSeqs);
        printf("Total number of subsequences blacklisted:\t%lld\n", procLogs.nSubSeqsBL);
        printf("Total number of subsequences after scaling:\t%lld\n", procLogs.nProcSubSeq);
        
        printf("\n#################### FNC CALLS RELATED STATS ####################\n");
        printf("Number of FL lowerbound is computed:\t%lld\n", procLogs.nLB_KIM_FL);
        printf("Number of LB_Keogh_EQ computed:\t%lld\n", procLogs.nLB_Keogh_EQ);
        printf("Number of LB_Keogh_EC computed:\t%lld\n", procLogs.nLB_Keogh_EC);
        printf("Number of times DTW computed:\t%lld\n", procLogs.nDTW_EA);
        printf("Number of updates of priority list:\t%lld\n", procLogs.nPriorityUpdates);
    }
}


