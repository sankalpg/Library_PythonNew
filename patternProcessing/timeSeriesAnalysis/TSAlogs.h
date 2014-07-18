#ifndef TSALOG_H

#define TSALOG_H


#include "TSAhashDefs.h"
#include "TSAdataStructs.h"

class TSAlogs
{
public:
    procLogs_t procLogs;
    
    TSAlogs();
    int dumpProcLogs(char *logFile, int verbos);
    
};

#endif //TSALOG_H