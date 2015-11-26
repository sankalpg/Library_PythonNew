import numpy as np
import scipy.signal as sig
from scipy.io import loadmat
import string
eps = np.finfo(np.float).eps
import sys, os
import matplotlib.pyplot as plt
sys.path.append(os.path.join(os.path.dirname(__file__), '../batchProcessing'))

import basicOperations as BO
import batchProcessing as BP
import pitchHistogram as PH
import segmentation as seg
import transcription as ts



def batchProc(root_dir, audioExt = '.mp3', pitchExt = '.pitchSilIntrpPP', tonicExt = '.tonicFine'):
    
    filenames = BP.GetFileNamesInDir(root_dir, '.mp3')
    segObj = seg.melodySegmentation()
    #phObj = PH.PitchHistogram()
    
    
    for filename in filenames[:1]:
        print "Processing file %s"%filename
        
        fname, ext = os.path.splitext(filename)
        pitch, time, Hop = BO.readPitchFile(fname + pitchExt)
        tonic = np.loadtxt(fname  + tonicExt)
        pcents = BO.PitchHz2Cents(pitch, tonic)
        pdata = (time,pcents,Hop)
        
        breathPheases = segObj.ExtractBreathPhrases(pcents, Hop, 0.5)
        breathPheases_s= []
        print "Hop size is %f"%Hop
        
        for b in breathPheases:
            breathPheases_s.append([b[0]*Hop, b[1]*Hop, b[2]])
            
        #tsObj = ts.Data(pdata) 
        #tsObj.Data(pdata, thres=8, width=35, verbose=False)
        #return breathPheases_s
    
       # def findValidSvaras(pitch, tonic):
            
        phObj = PH.PitchHistogram(pitch, float(tonic))
        phObj.ComputePitchHistogram(Oct_fold=1)
        phObj.ValidSwarLocEstimation()
        phObj.SwarLoc2Cents()
        svara = phObj.swarCents
        print svara
        phObj.PlotHistogram()