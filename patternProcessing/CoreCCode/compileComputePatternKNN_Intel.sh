g++ -I. ../../similarityMeasures/dtw/dtw.c ../../basicDSPFuncs/basicDSPCFuncs.c ../../similarityMeasures/dtw/tables.c patternProTable.c ComputePatternKNN.c MotifDataIO.c -O3 -o ComputePatternKNN_O3

