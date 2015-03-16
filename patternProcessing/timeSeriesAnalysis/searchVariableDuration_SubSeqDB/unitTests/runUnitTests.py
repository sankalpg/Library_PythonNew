import numpy as np
import os, sys
import filecmp


if __name__ == "__main__":
    
    searchExt = '.motifSearch.motifSearch'
    referenceOut = 'output_reference/'
    currentOut = 'output_current/'
    
    #running reference and current code on allthe configurations
    
    
    ########## test1 #########
    testInd = 1
    configFile = 'configFiles_533_1'
    cmd1 = './ICASSP2015_Experiment_O3_REF ' + referenceOut +configFile + ' '+ referenceOut +configFile +   '.txt fileExtensions.txt 198 -1 1'
    cmd2 = './../searchVariableDurationSubSeqDB_O3 ' + currentOut +configFile + ' '  + currentOut +configFile   +    '.txt fileExtensions.txt 198 -1 1'
    os.system(cmd1)
    os.system(cmd2)
    if filecmp.cmp(referenceOut +configFile + searchExt, currentOut +configFile+ searchExt) and os.path.getsize(referenceOut +configFile + searchExt) > 0:
        print "################"
        print "Test %d passed."%(testInd)
        print "################"
    else:
        print "################"
        print "Test %d failed."%(testInd)
        print "################"
        

    ########## test2 #########
    testInd = 2
    configFile = 'configFiles_533_2'
    cmd1 = './ICASSP2015_Experiment_O3_REF ' + referenceOut +configFile + ' '+ referenceOut +configFile +   '.txt fileExtensions.txt 198 -1 1'
    cmd2 = './../searchVariableDurationSubSeqDB_O3 ' + currentOut +configFile + ' '  + currentOut +configFile   +    '.txt fileExtensions.txt 198 -1 1'
    os.system(cmd1)
    os.system(cmd2)
    if filecmp.cmp(referenceOut +configFile + searchExt, currentOut +configFile+ searchExt) and os.path.getsize(referenceOut +configFile + searchExt) > 0:
        print "################"
        print "Test %d passed."%(testInd)
        print "################"
    else:
        print "################"
        print "Test %d failed."%(testInd)
        print "################"        
    
    
    
    ########## test3 #########
    testInd = 3
    configFile = 'configFiles_533_3'
    cmd1 = './ICASSP2015_Experiment_O3_REF ' + referenceOut +configFile + ' '+ referenceOut +configFile +   '.txt fileExtensions.txt 198 -1 1'
    cmd2 = './../searchVariableDurationSubSeqDB_O3 ' + currentOut +configFile + ' '  + currentOut +configFile   +    '.txt fileExtensions.txt 198 -1 1'
    os.system(cmd1)
    os.system(cmd2)
    if filecmp.cmp(referenceOut +configFile + searchExt, currentOut +configFile+ searchExt) and os.path.getsize(referenceOut +configFile + searchExt) > 0:
        print "################"
        print "Test %d passed."%(testInd)
        print "################"
    else:
        print "################"
        print "Test %d failed."%(testInd)
        print "################"        
    
    
    
    
    ########## test4 #########
    testInd = 4
    configFile = 'configFiles_533_4'
    cmd1 = './ICASSP2015_Experiment_O3_REF ' + referenceOut +configFile + ' '+ referenceOut +configFile +   '.txt fileExtensions.txt 198 -1 1'
    cmd2 = './../searchVariableDurationSubSeqDB_O3 ' + currentOut +configFile + ' '  + currentOut +configFile   +    '.txt fileExtensions.txt 198 -1 1'
    os.system(cmd1)
    os.system(cmd2)
    if filecmp.cmp(referenceOut +configFile + searchExt, currentOut +configFile+ searchExt) and os.path.getsize(referenceOut +configFile + searchExt) > 0:
        print "################"
        print "Test %d passed."%(testInd)
        print "################"
    else:
        print "################"
        print "Test %d failed."%(testInd)
        print "################"      
        
        
        
    ########## test5 #########
    testInd = 5
    configFile = 'configFiles_533_5'
    cmd1 = './ICASSP2015_Experiment_O3_REF ' + referenceOut +configFile + ' '+ referenceOut +configFile +   '.txt fileExtensions.txt 198 -1 1'
    cmd2 = './../searchVariableDurationSubSeqDB_O3 ' + currentOut +configFile + ' '  + currentOut +configFile   +    '.txt fileExtensions.txt 198 -1 1'
    os.system(cmd1)
    os.system(cmd2)
    if filecmp.cmp(referenceOut +configFile + searchExt, currentOut +configFile+ searchExt) and os.path.getsize(referenceOut +configFile + searchExt) > 0:
        print "################"
        print "Test %d passed."%(testInd)
        print "################"
    else:
        print "################"
        print "Test %d failed."%(testInd)
        print "################"            
        
        
        
    ########## test6 #########
    testInd = 6
    configFile = 'configFiles_533_6'
    cmd1 = './ICASSP2015_Experiment_O3_REF ' + referenceOut +configFile + ' '+ referenceOut +configFile +   '.txt fileExtensions.txt 198 -1 1'
    cmd2 = './../searchVariableDurationSubSeqDB_O3 ' + currentOut +configFile + ' '  + currentOut +configFile   +    '.txt fileExtensions.txt 198 -1 1'
    os.system(cmd1)
    os.system(cmd2)
    if filecmp.cmp(referenceOut +configFile + searchExt, currentOut +configFile+ searchExt) and os.path.getsize(referenceOut +configFile + searchExt) > 0:
        print "################"
        print "Test %d passed."%(testInd)
        print "################"
    else:
        print "################"
        print "Test %d failed."%(testInd)
        print "################"            
        
        
        
    ########## test7 #########
    testInd = 7
    configFile = 'configFiles_533_7'
    cmd1 = './ICASSP2015_Experiment_O3_REF ' + referenceOut +configFile + ' '+ referenceOut +configFile +   '.txt fileExtensions.txt 198 -1 1'
    cmd2 = './../searchVariableDurationSubSeqDB_O3 ' + currentOut +configFile + ' '  + currentOut +configFile   +    '.txt fileExtensions.txt 198 -1 1'
    os.system(cmd1)
    os.system(cmd2)
    if filecmp.cmp(referenceOut +configFile + searchExt, currentOut +configFile+ searchExt) and os.path.getsize(referenceOut +configFile + searchExt) > 0:
        print "################"
        print "Test %d passed."%(testInd)
        print "################"
    else:
        print "################"
        print "Test %d failed."%(testInd)
        print "################"            
        
        
        
        
    ########## test8 #########
    testInd = 8
    configFile = 'configFiles_533_8'
    cmd1 = './ICASSP2015_Experiment_O3_REF ' + referenceOut +configFile + ' '+ referenceOut +configFile +   '.txt fileExtensions.txt 198 -1 1'
    cmd2 = './../searchVariableDurationSubSeqDB_O3 ' + currentOut +configFile + ' '  + currentOut +configFile   +    '.txt fileExtensions.txt 198 -1 1'
    os.system(cmd1)
    os.system(cmd2)
    if filecmp.cmp(referenceOut +configFile + searchExt, currentOut +configFile+ searchExt) and os.path.getsize(referenceOut +configFile + searchExt) > 0:
        print "################"
        print "Test %d passed."%(testInd)
        print "################"
    else:
        print "################"
        print "Test %d failed."%(testInd)
        print "################"            
        
        
        
    ########## test9 #########
    testInd = 9
    configFile = 'configFiles_533_9'
    cmd1 = './ICASSP2015_Experiment_O3_REF ' + referenceOut +configFile + ' '+ referenceOut +configFile +   '.txt fileExtensions.txt 198 -1 1'
    cmd2 = './../searchVariableDurationSubSeqDB_O3 ' + currentOut +configFile + ' '  + currentOut +configFile   +    '.txt fileExtensions.txt 198 -1 1'
    os.system(cmd1)
    os.system(cmd2)
    if filecmp.cmp(referenceOut +configFile + searchExt, currentOut +configFile+ searchExt) and os.path.getsize(referenceOut +configFile + searchExt) > 0:
        print "################"
        print "Test %d passed."%(testInd)
        print "################"
    else:
        print "################"
        print "Test %d failed."%(testInd)
        print "################"            
        
    ########## test10 #########
    testInd = 10
    configFile = 'configFiles_533_10'
    cmd1 = './ICASSP2015_Experiment_O3_REF ' + referenceOut +configFile + ' '+ referenceOut +configFile +   '.txt fileExtensions.txt 198 -1 1'
    cmd2 = './../searchVariableDurationSubSeqDB_O3 ' + currentOut +configFile + ' '  + currentOut +configFile   +    '.txt fileExtensions.txt 198 -1 1'
    os.system(cmd1)
    os.system(cmd2)
    if filecmp.cmp(referenceOut +configFile + searchExt, currentOut +configFile+ searchExt) and os.path.getsize(referenceOut +configFile + searchExt) > 0:
        print "################"
        print "Test %d passed."%(testInd)
        print "################"
    else:
        print "################"
        print "Test %d failed."%(testInd)
        print "################"                    