import numpy as np
import os, sys
import filecmp


if __name__ == "__main__":
    
    baseDir = "unitTestFiles/"
    outExt = '.motif'
    outExt_REF = '.motif_REF'
    #running reference and current code on allthe configurations
    
    ########## test1 #########
    testInd = 1
    fileName = '03_Teliyaleru_Rama'
    
    cmd1 = './discoverPatterns_O3_REF ' +'"'+ baseDir +fileName +'"'+ ' '+ 'procParamsDiscovery.txt fileExtensionsDiscovery_REF.txt 5 -1 0'
    cmd2 = './../discoverPatterns_O3 ' +'"'+ baseDir +fileName +'"'+' '+ 'procParamsDiscovery.txt fileExtensionsDiscovery.txt 5 -1 0'
    
    os.system(cmd1)
    os.system(cmd2)
    if filecmp.cmp(baseDir +fileName + outExt_REF, baseDir +fileName + outExt)  and os.path.getsize(baseDir +fileName + outExt_REF) > 0:
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
    
    cmd1 = './discoverPatterns_O3_REF ' +'"'+ baseDir +fileName +'"'+ ' '+ 'procParamsDiscovery.txt fileExtensionsDiscovery_REF.txt 5 -1 0'
    cmd2 = './../discoverPatterns_O3 ' +'"'+ baseDir +fileName +'"'+' '+ 'procParamsDiscovery.txt fileExtensionsDiscovery.txt 5 -1 0'
    
    os.system(cmd1)
    os.system(cmd2)
    if filecmp.cmp(baseDir +fileName + outExt_REF, baseDir +fileName + outExt)  and os.path.getsize(baseDir +fileName + outExt_REF) > 0:
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
    
    cmd1 = './discoverPatterns_O3_REF ' +'"'+ baseDir +fileName +'"'+ ' '+ 'procParamsDiscovery.txt fileExtensionsDiscovery_REF.txt 5 -1 0'
    cmd2 = './../discoverPatterns_O3 ' +'"'+ baseDir +fileName +'"'+' '+ 'procParamsDiscovery.txt fileExtensionsDiscovery.txt 5 -1 0'
    
    os.system(cmd1)
    os.system(cmd2)
    if filecmp.cmp(baseDir +fileName + outExt_REF, baseDir +fileName + outExt)  and os.path.getsize(baseDir +fileName + outExt_REF) > 0:
        print "################"
        print "Test %d passed."%(testInd)
        print "################"
    else:
        print "################"
        print "Test %d failed."%(testInd)
        print "################"        