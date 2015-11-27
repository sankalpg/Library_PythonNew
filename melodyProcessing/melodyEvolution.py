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
        
        breathPheases = np.array(segObj.ExtractBreathPhrases(pcents, Hop, 0.5)).astype(np.float)
        validTime = breathPheases
        validTime[:,0] = validTime[:,0]*Hop
        validTime[:,1] = validTime[:,1]*Hop
        
        print validTime
        np.savetxt((''.join([filename,'.bphrases'])),validTime, delimiter = "\t", fmt="%f\t%f\t%d")
        #print "Hop size is %f"%Hop
        
        #for b in breathPheases:
            #np.appendbreathPheases_s.append([b[0]*Hop, b[1]*Hop, b[2]])
            #validTime = breathPheases_s[-1]
            #print np.array(validTime).shape
            ##np.savetxt('filename.bphrases',np.array(validTime, delimiter = " ", fmt="%f")
            ##validStart = breathPheases_s[-1][0]
            ##validEnd = breathPheases_s[-1][1]
            
            
        tsObj = ts.Data(pdata, thres=8, width=35, verbose=False) 
        song_notes = tsObj.transcribed
        song_durs = [e-s+1 for s,e in zip(tsObj.st, tsObj.en)]
        song_tuple = zip(song_notes, song_durs)
        
        print len(tsObj.times), len(tsObj.levels)
        
        #plt.plot(tsObj.times[:len(tsObj.levels)], tsObj.ts[:len(tsObj.levels)])
        #plt.plot(tsObj.times[:len(tsObj.levels)], tsObj.levels)
        #plt.show()
        
        #print validStart,validEnd
        #validTime = breathPheases_s[-1][0]:breathPheases_s[-1][1]
        #print validTime
        
        #return breathPheases_s
    
       # def findValidSvaras(pitch, tonic):
            
        #phObj = PH.PitchHistogram(pitch, float(tonic))
        #phObj.ComputePitchHistogram(Oct_fold=1)
        #phObj.ValidSwarLocEstimation()
        #phObj.SwarLoc2Cents()
        #svara = phObj.swarCents
        #print svara
        #phObj.PlotHistogram()