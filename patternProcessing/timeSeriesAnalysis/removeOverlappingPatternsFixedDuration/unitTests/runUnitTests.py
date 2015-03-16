import numpy as np
import os, sys
import filecmp


if __name__ == "__main__":
    
    baseDir = "/home/sankalp/Work/Work_PhD/library_pythonnew/patternProcessing/timeSeriesAnalysis/removeOverlappingPatternsFixedDuration/unitTests/unitTestFiles/"
    outExt = '.blackPatterns'
    outExt_REF = '.blackPatterns_REF'
    #running reference and current code on allthe configurations
    
    print "################################"
    print "testing for removeOverlappingPatternsFixedDuration"
    print "################################"
    
    ########## test1 #########
    testInd = 1
    fileName = '01_Inta_Chala'
    
    cmd1 = './removeOverlappingPatternsFixedDuration_O3_REF ' +'"'+ baseDir +fileName +'"'+ ' '+ 'procParamsRemoveOverlap.txt fileExtensionsRemoveOverlap_REF.txt'
    cmd2 = './../removeOverlappingPatternsFixedDuration_O3 ' +'"'+ baseDir +fileName +'"'+ ' '+ 'procParamsRemoveOverlap.txt fileExtensionsRemoveOverlap.txt'
    
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
    
    cmd1 = './removeOverlappingPatternsFixedDuration_O3_REF ' +'"'+ baseDir +fileName +'"'+ ' '+ 'procParamsRemoveOverlap.txt fileExtensionsRemoveOverlap_REF.txt'
    cmd2 = './../removeOverlappingPatternsFixedDuration_O3 ' +'"'+ baseDir +fileName +'"'+ ' '+ 'procParamsRemoveOverlap.txt fileExtensionsRemoveOverlap.txt'
    
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
