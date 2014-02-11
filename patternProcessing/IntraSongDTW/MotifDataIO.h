
#ifndef MOTIFDATAIO_H

#define MOTIFDATAIO_H

#include "SearchInterDTW.h"
#include <string.h>

typedef struct mappInfo
{
    INDTYPE last_line;
    
}mappInfo_t;

INDTYPE getNumLines(const char *file);

INDTYPE generateSubsequenceDB(DATATYPE ***d, segInfo_t **t, int *motifLen, char *baseName, fileExts_t *myFileExts, procParams_t *myProcParams, procLogs_t *myProcLogs, int verbos);
INDTYPE loadSeedMotifSequence(DATATYPE ***d, segInfo_t **t, int *motifLen, char *baseName, fileExts_t *myFileExts, procParams_t *myProcParams, int maxNMotifs, int verbos);

void dumpSearchMotifInfo(char *motifFile, char *mappFile, char *searchFile, motifInfo** topKmotifs, segInfo_t *tStampsInterpSeed, segInfo_t *tStampsInterp, int NSeeds, INDTYPE K, mappInfo_t *mapp, int verbos);
void dumpDiscoveryLogs(char *logFile, procLogs_t myProcLogs, int verbos);
void dumpDiscoveredMotifInfo(char *motifFile, motifInfo *topKmotifs, segInfo_t *tStampsInterp, int K, int verbos);

#endif //MOTIFDATAIO_H