
#ifndef  TSA_DATAIO_H

#define TSA_DATAIO_H

#include "TSAhashDefs.h"
#include "TSAdataStructs.h"

class TSAdataIO
{

private:
    char*       fileNameTS;
    char*       fileNameMotifs;    
public:
    void*       readTSSeq(char* fileName, void* seq);
    void*       readTSSubSeq(char* fileName, void *subSeq, int len, int sizeSample);
    int         dumpMotifInfo();
    int         countNumberLines();
}

class TSAparamIO
{
    
public:
    procParams_t procParams;
    fileExts_t fileExts;
    
    int readCommandLineIO(int argc , char *argv[]);
    int readFileExtsInfoFile(char* infoFileExts);
    int readProcessingParams(char* paramFileName);
}


#endif //TSA_DATAIO_H

