g++ -I. ../../similarityMeasures/dtw/dtw.c ../../basicDSPFuncs/basicDSPCFuncs.c ../../similarityMeasures/dtw/tables.c patternProTable.c ComputePatternDistances.c MotifDataIO.c -O3 -o ComputePatternDistances_O3

