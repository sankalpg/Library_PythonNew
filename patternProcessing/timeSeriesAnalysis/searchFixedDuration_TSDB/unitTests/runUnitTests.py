import numpy as np
import os, sys
import filecmp

def createFileList(filename, baseDir, flistExt, flistExtLocal):

    lines = open(os.path.join(baseDir, filename + flistExt),'r').readlines()
    fid = open(os.path.join(baseDir, filename + flistExtLocal),'w')
    for line in lines:
        fname = os.path.basename(line.strip())
        fid.write("%s\n"%(os.path.join(baseDir, fname)))
    fid.close()

if __name__ == "__main__":
    
    baseDir = sys.argv[1]
    fileList = '.flist'
    fileListLocal = '.flistLocal'
    #baseDir = "/home/sankalp/Work/Work_PhD/library_pythonnew/patternProcessing/timeSeriesAnalysis/searchFixedDuration_TSDB/unitTests/unitTestFiles/"
    outExt = '.motifSearchSqEuclidean'
    outExt_REF = '.motifSearch_REFSqEuclidean'
    #running reference and current code on allthe configurations
    
    ########## test1 #########
    testInd = 1
    fileName = '03_Teliyaleru_Rama'
    createFileList(fileName, baseDir, fileList, fileListLocal)
 
    cmd1 = './searchPatterns_O3_REF ' +'"'+ baseDir +fileName +'"'+ ' '+ 'procParamsSearch.txt fileExtensionsSearch_REF.txt 5 -1 0'
    cmd2 = './../searchPatterns_O3 ' +'"'+ baseDir +fileName +'"'+ ' '+ 'procParamsSearch.txt fileExtensionsSearch.txt 5 -1 0'
    
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
    fileName = '3-02_Aduvum_Sholluval_(Padam)'
    createFileList(fileName, baseDir, fileList, fileListLocal)

    cmd1 = './searchPatterns_O3_REF ' +'"'+ baseDir +fileName +'"'+ ' '+ 'procParamsSearch.txt fileExtensionsSearch_REF.txt 5 -1 0'
    cmd2 = './../searchPatterns_O3 ' +'"'+ baseDir +fileName +'"'+ ' '+ 'procParamsSearch.txt fileExtensionsSearch.txt 5 -1 0'
    
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
        
        
    ########## test3 #########
    testInd = 3
    fileName = '01_Inta_Chala'
    createFileList(fileName, baseDir, fileList, fileListLocal)

    cmd1 = './searchPatterns_O3_REF ' +'"'+ baseDir +fileName +'"'+ ' '+ 'procParamsSearch.txt fileExtensionsSearch_REF.txt 5 -1 0'
    cmd2 = './../searchPatterns_O3 ' +'"'+ baseDir +fileName +'"'+ ' '+ 'procParamsSearch.txt fileExtensionsSearch.txt 5 -1 0'
    
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
