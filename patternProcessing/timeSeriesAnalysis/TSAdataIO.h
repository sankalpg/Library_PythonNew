
#ifndef  TSA_DATAIO_H

#define TSA_DATAIO_H

#include "TSAhashDefs.h"
#include "TSAdataStructs.h"


class TSAparamHandle
{
    
public:
    procParams_t procParams;
    fileExts_t fileExts;
    
    TSAparamHandle();
    ~TSAparamHandle();
    
    int readParamsFromFile(char *paramFile);
    int readFileExtsInfoFile(char* infoFileExts);
    
    procParams_t* getParamPtr();
    fileExts_t* getExtPtr();
    
};

class fileNameHandler
{
    
public:
    char fileName[MAX_FNAME_CHARS];
    fileExts_t *fileExtPtr;
    char *baseName;
    char **searchFileNames;
    int nSearchFiles;
    
    fileNameHandler();
    ~fileNameHandler();
    int initialize(char *bName, fileExts_t *fExtPtr);
    char *getTSFileName();
    char *getTonicFileName();
    char *getBlackListSegFileName();
    char *getLogFileName();
    char *getParamDumpFileName();
    char *getOutFileName();
    char *getMappFileName();
    char *getSearchListFileName();
    char* getQueryFileName();
    int  loadSearchFileList();
    char *getOutFileNamePostRR(int similarityMeasure);
    char *getSubSeqFileName();
    char *getSubSeqTNFileName();
    char *getSubSeqInfoFileName();
};


class TSAdataHandler
{
    
public:
    
    procLogs_t *procLogPtr;
    fileExts_t *fileExtPtr;
    procParams_t procParams;
    char *baseName;
    fileNameHandler fHandle;
    
    
    TSAIND lenTS;
    TSAIND nSubSeqs;
    TSAIND nLinesFile;
    
    
    int *blacklist;
    int isBlackListAlloc;
    
    TSAseg_t *queryTStamps;
    TSAIND nQueries;
    
    
    
    TSAsam_t *samPtr;
    TSAsubSeq_t *subSeqPtr;
    
    float pHop;
    
    int         readTSData(char *fileName);
    int         readSubSeqData(char *fileName, TSAIND nSubs);
    int         readSubSeqLengths(char *fileName);
    
    //void*       readTSSubSeq(char* fileName, void *subSeq, int len, int sizeSample);
    int       readHopSizeTS(char *fileName);
    int         dumpMotifInfo();
    int         countNumberLines();
    int         readSubSeqInfo(char *fileName);
    int         genTemplate1SubSeqs();
    TSAIND      getNumLines(const char *file);
    int         setSubSeqLengthsFIX(int motifLen);
    int         downSampleTS();
    int         downSampleSubSeqs();
    int         quantizeSampleTS(int quantizationType);
    int         quantizeSampleSubSeqs(int quantizationType);
    int         filterSamplesTS();
    int         convertHz2Cents(char *tonicFileName);
    int         initializeBlackList();
    int         updateBLDurThsld();
    int         updateBLStdThsld();
    int         updateBLInvalidSegment(char *fileName);
    int         calculateDiffMotifLengths();
    int         genSlidingWindowSubSeqs();
    int         computeStdTSLocal(float **std, int varSam);
    int         filterBlackListedSubSeqs();
    int         genUniScaledSubSeqs();
    int         genUniScaledSubSeqsVarLen();
    int         dumpDiscMotifInfo(char *motifFile, TSAmotifInfo_t *priorityQDisc, int K, int verbos);
    int         dumpSearMotifInfo(char *motifFile, TSAmotifInfoExt_t **priorityQSear, TSAIND nQueries, int K, int verbos);
    int         readQueryTimeStamps(char *queryFileName, int format);
    int         genSubSeqsWithTStarts(TSAseg_t *queryTStamps, TSAIND nQueries);
    int         genSubSeqsWithTStamps(TSAseg_t *qTStamps, TSAIND nQueries);
    int         loadMotifDataTemplate1();
    int         freeSubSeqsMem();
    int         normalizeSubSeqs(int normType);
    
    TSAdataHandler(char *bName, procLogs_t *procLogs, fileExts_t *fileExts, procParams_t *pParams);
    ~TSAdataHandler();
    
};


#endif //TSA_DATAIO_H
