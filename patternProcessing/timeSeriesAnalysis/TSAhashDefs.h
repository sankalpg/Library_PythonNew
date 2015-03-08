
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
#define Q_DIST_BLOCKSIZE 500

#define PID_DEFAULT1    -1
#define PID_DEFAULT2    -2
#define PID_DEFAULT3    -3

#define FID_DEFAULT1    -1

enum 
{
MOTIFPAIR_DUMP_FORMAT=0,
VIGNESH_MOTIF_ANNOT_FORMAT,
MY_MOTIF_ANNOT_FORMAT,//which is basically edited vbersion of vignesh's annotation with start and end and id instead of start duration and id
PATTERNS_PER_FILE_DUMP,
MOTIFID1_MOTIFID2_DIST,
};

enum 
{
NO_NORM=0,
TONIC_NORM,
Z_NORM,
MEAN_SUB_NORM,
MEDIAN_SUB_NORM,
MAD_NORM,
QMEDIAN_SUB_NORM, 
TONIC_NORM_PASAPA
};

enum
{
NORM_1=0,
PATT_LEN,
PATH_LEN,
MAXLEN_NO_NORM,
MAXLEN_PATH_LEN
};

enum 
{
NO_QUANT=0,
NEAREST_NOTE_QUANT,
};

enum
{
	Var1=1,
	Var2
};//These are mthod variants for supervised analysis. Var1: ground truth segmentation is considered and patterns are linearly stretched to the same length. Var2: candiadte patterns are taken to be of the length of the query.


//Extern declarations
extern const int combAllwd_5[5][5];
extern const int combAllwd_3[3][3];

#endif //TSA_HASHDEFS_H

