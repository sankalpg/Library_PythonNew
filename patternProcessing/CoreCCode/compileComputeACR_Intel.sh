g++ -I. ../../similarityMeasures/dtw/dtw.c ../../basicDSPFuncs/basicDSPCFuncs.c ../../similarityMeasures/dtw/tables.c patternProTable.c ComputeACR.c MotifDataIO.c -O3 -o ComputeACR_O3
