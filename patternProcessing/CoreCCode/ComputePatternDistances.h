
#ifndef DiscoverIntraDTW_H

#define DiscoverIntraDTW_H


#include "MotifDataIO.h"
#include <unistd.h>

//functions
DISTTYPE manageTopKMotifs(patternDist_t *topKmotifs, int K, INDTYPE id1, INDTYPE id2, DISTTYPE dist);

int dumpKNNPatterns(char *filename, patternDist_t **topKmotifs, int K, int nPriorityList); 
#endif //#ifndef DiscoverIntraDTW_H




