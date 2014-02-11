

#include "DiscoverIntraDTW.h"

INDTYPE getNumLines(const char *file);

INDTYPE generateSubsequenceDB(DATATYPE ***d, segInfo_t **t, int *motifLen, char *baseName, fileExts_t *myFileExts, procParams_t *myProcParams, procLogs_t *myProcLogs, int verbos);

void dumpLogs(char *logFile, procLogs_t myProcLogs, int verbos);
void dumpMotifInfo(char *motifFile, motifInfo *topKmotifs, segInfo_t *tStampsInterp, int K, int verbos);