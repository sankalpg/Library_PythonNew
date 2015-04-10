import numpy as np
import os, sys
import filecmp


if __name__ == "__main__":
    
    baseDir = "/home/sankalp/Work/Work_PhD/library_pythonnew/patternProcessing/timeSeriesAnalysis/searchVariableDuration_TSDB/unitTests/unitTestFiles/"
    outExt = '.motifSearch'
    outExt_REF = '.motifSearch_REF'
    #running reference and current code on allthe configurations
    
    print "################################"
    print "testing for searchPatternsVarLen"
    print "################################"
    
    ########## test1 #########
    testInd = 1
    fileName = 'ALB_Bhairavi1_70albmgpnganesan'
    
    cmd1 = './searchPatternsVarLen_O3_REF ' +'"'+ baseDir +fileName +'"'+ ' '+ 'procParamsSearch.txt fileExtensionsSearchVarLen_REF.txt 10 -1 0'
    cmd2 = './../searchPatternsVarLen_O3 ' +'"'+ baseDir +fileName +'"'+ ' '+ 'procParamsSearch.txt fileExtensionsSearchVarLen.txt 10 -1 0'
    
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
        
        
    ########## test1 #########
    testInd = 2
    fileName = 'Sanjay_Bhairavi1_FASfeb2001'
    
    cmd1 = './searchPatternsVarLen_O3_REF ' +'"'+ baseDir +fileName +'"'+ ' '+ 'procParamsSearch.txt fileExtensionsSearchVarLen_REF.txt 10 -1 0'
    cmd2 = './../searchPatternsVarLen_O3 ' +'"'+ baseDir +fileName +'"'+ ' '+ 'procParamsSearch.txt fileExtensionsSearchVarLen.txt 10 -1 0'
    
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
        
        
        
        
    print "####################################"
    print "testing for searchPatternsAllDBQuery"
    print "####################################"  
    outExt = '.motifSearch_ALL'
    outExt_REF = '.motifSearch_REF_ALL'
    ########## test1 #########
    testInd = 1
    fileName = 'ALB_Bhairavi1_70albmgpnganesan'
    
    cmd1 = './searchPatternsAllDBQuery_O3_REF ' +'"'+ baseDir +fileName +'"'+ ' '+ 'procParamsSearch.txt fileExtensionsSearchVarLen_REF_ALL.txt 10 -1 0'
    cmd2 = './../searchPatternsAllDBQuery_O3 ' +'"'+ baseDir +fileName +'"'+ ' '+ 'procParamsSearch.txt fileExtensionsSearchVarLen_ALL.txt 10 -1 0'
    
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
        
        
    ########## test1 #########
    testInd = 2
    fileName = 'Sanjay_Bhairavi1_FASfeb2001'
    
    cmd1 = './searchPatternsAllDBQuery_O3_REF ' +'"'+ baseDir +fileName +'"'+ ' '+ 'procParamsSearch.txt fileExtensionsSearchVarLen_REF_ALL.txt 10 -1 0'
    cmd2 = './../searchPatternsAllDBQuery_O3 ' +'"'+ baseDir +fileName +'"'+ ' '+ 'procParamsSearch.txt fileExtensionsSearchVarLen_ALL.txt 10 -1 0'
    
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
            
        
 