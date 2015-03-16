import numpy as np
import os, sys

sys.path.append(os.path.join(os.path.dirname(__file__), '../../../batchProcessing/'))
import batchProcessing as BP


def batchProcessDumpSubSeqFixedDuration(root_dir, paramFile, fileExtFile, searchExt = '.mp3'):
    
    filenames = BP.GetFileNamesInDir(root_dir, searchExt)
    
    for ii,filename in enumerate(filenames):
        fname,ext = os.path.splitext(filename)
        print "Processing file %d of %d"%(ii+1, len(filenames))
        cmd1 = './dumpSubSeqFixedDuration_O3 ' +'"'+ fname +'"'+ ' '+ paramFile+ ' '+ fileExtFile+' 1'
        os.system(cmd1)
        
        
        