g++ -I. ../../similarityMeasures/dtw/dtw.c ../../basicDSPFuncs/basicDSPCFuncs.c ../../similarityMeasures/dtw/tables.c patternProTable.c DumpTimeSeries.c MotifDataIO.c -O3 -o DumpTimeSeries_O3
