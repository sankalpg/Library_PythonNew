
#ifndef TSA_HASHDEFS_H

#define TSA_HASHDEFS_H

#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <fcntl.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <signal.h>
#include <stdint.h>
#include <string.h>
#include <float.h>
#include <string.h>
#include <algorithm>
#include "../../similarityMeasures/dtw/dtw.h"
#include "../../basicDSPFuncs/basicDSPCFuncs.h"


//typdefs
#define TSADATA      double
#define TSADIST        double
#define TSAIND        long long int


// other Constants
#define INF FLT_MAX
#define LOG2  0.693147180559945
#define EPS 0.0000000000000001
#define MAXNTEMPOFACTORS 10
#define MAX_FNAME_CHARS 400
#define MAX_FEXT_CHARS 100
#define MAX_NUM_SEARCHFILES 2000

#define PID_DEFAULT1    -1
#define PID_DEFAULT2    -2
#define PID_DEFAULT3    -3

#define FID_DEFAULT1    -1

enum 
{
MOTIFPAIR_DUMP_FORMAT=0,
VIGNESH_MOTIF_ANNOT_FORMAT
};


//Extern declarations
extern const int combAllwd_5[5][5];
extern const int combAllwd_3[3][3];

#endif //TSA_HASHDEFS_H

