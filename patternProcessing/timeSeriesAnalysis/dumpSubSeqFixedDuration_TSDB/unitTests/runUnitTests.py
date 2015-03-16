import numpy as np
import os, sys
import filecmp


if __name__ == "__main__":
    
    baseDir = "/home/sankalp/Work/Work_PhD/library_pythonnew/patternProcessing/timeSeriesAnalysis/dumpSubSeqFixedDuration_TSDB/unitTests/unitTestFiles/"
    outExt = '.patternData'
    outExt_REF = '.patternData_REF'
    #running reference and current code on allthe configurations
    
    print "################################"
    print "testing for dumpSubSeqFixedDuration"
    print "################################"
    
    ########## test1 #########
    testInd = 1
    fileName = '01_Inta_Chala'
    
    cmd1 = './dumpSubSeqFixedDuration_O3_REF ' +'"'+ baseDir +fileName +'"'+ ' '+ 'procParamsDump.txt fileExtensionsdumpSubSeqFixedDuration_REF.txt 0'
    cmd2 = './../dumpSubSeqFixedDuration_O3 ' +'"'+ baseDir +fileName +'"'+ ' '+ 'procParamsDump.txt fileExtensionsdumpSubSeqFixedDuration.txt 0'
    
    os.system(cmd1)
    os.system(cmd2)
    if filecmp.cmp(baseDir +fileName + outExt_REF, baseDir +fileName + outExt) and os.path.getsize(baseDir +fileName + outExt_REF) > 0:
        print "################"
        print "Test %d passed."%(testInd)
        print "################"
    else:
        print "################"
        print "Test %d failed."%(testInd)
        print "################"


    ########## test2 #########
    testInd = 2
    fileName = '02_Ragaratnamalika'
    
    cmd1 = './dumpSubSeqFixedDuration_O3_REF ' +'"'+ baseDir +fileName +'"'+ ' '+ 'procParamsDump.txt fileExtensionsdumpSubSeqFixedDuration_REF.txt 0'
    cmd2 = './../dumpSubSeqFixedDuration_O3 ' +'"'+ baseDir +fileName +'"'+ ' '+ 'procParamsDump.txt fileExtensionsdumpSubSeqFixedDuration.txt 0'
    
    os.system(cmd1)
    os.system(cmd2)
    if filecmp.cmp(baseDir +fileName + outExt_REF, baseDir +fileName + outExt) and os.path.getsize(baseDir +fileName + outExt_REF) > 0:
        print "################"
        print "Test %d passed."%(testInd)
        print "################"
    else:
        print "################"
        print "Test %d failed."%(testInd)
        print "################"