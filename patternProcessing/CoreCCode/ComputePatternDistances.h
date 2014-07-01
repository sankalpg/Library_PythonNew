
#ifndef DiscoverIntraDTW_H

#define DiscoverIntraDTW_H


#include "MotifDataIO.h"
#include <unistd.h>

//functions
int updatePatternStorage(patternDist_t **patternDist, pattCntManager_t *patCntMan);
int manageMotifStorage(patternDist_t **patternDist, pattCntManager_t *patCntMan, INDTYPE id1, INDTYPE id2, DISTTYPE dist);
int dumpKNNPatterns(char *filename, patternDist_t **patternDists, pattCntManager_t *patCntMan, int nPriorityList);
#endif //#ifndef DiscoverIntraDTW_H




