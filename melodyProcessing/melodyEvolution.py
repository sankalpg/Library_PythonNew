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
        print len(pdata[1])*Hop
        
        breathPheases = np.array(segObj.ExtractBreathPhrases(pcents, Hop, 0.5)).astype(np.float)
        validTime = breathPheases
        validTime[:,0] = validTime[:,0]*Hop
        validTime[:,1] = validTime[:,1]*Hop
        
        #print validTime
        bphraseFilename = (''.join([fname,'.bphrases']))
        np.savetxt(bphraseFilename,validTime, delimiter = "\t", fmt="%f\t%f\t%d")
        #print "Hop size is %f"%Hop
        
        #for b in breathPheases:
            #np.appendbreathPheases_s.append([b[0]*Hop, b[1]*Hop, b[2]])
            #validTime = breathPheases_s[-1]
            #print np.array(validTime).shape
            ##np.savetxt('filename.bphrases',np.array(validTime, delimiter = " ", fmt="%f")
            ##validStart = breathPheases_s[-1][0]
            ##validEnd = breathPheases_s[-1][1]
            
        svaraSemitone, ignoreNotes = findValidSvaras(pitch,tonic)
        print svaraSemitone, ignoreNotes
        
        
        tsObj = ts.Data(pdata, ignoreNotes, thres=60, width=35, verbose=False) 
        song_notes = tsObj.transcribed
        song_durs = [e-s+1 for s,e in zip(tsObj.st, tsObj.en)]
        song_tuple = zip(song_notes, song_durs)
        #print song_tuple
        
        plt.plot(tsObj.times[:len(tsObj.levels)], tsObj.ts[:len(tsObj.levels)],'g--', alpha=0.7)
        
        #fig = plt.figure(figsize=(18,10))
        for i in range(validTime.shape[0]):
	    #ax = plt.subplot(2,3,i+1)
	    print validTime[i,0], validTime[i,1]
	    #ax.set_title(tsObj.get_notes(tsObj.get_string_part(validTime[i,0], validTime[i,1])[0]))
	    beg = max(m for m, n in enumerate(tsObj.times) if n < validTime[i,0])
	    fin = max(m for m, n in enumerate(tsObj.times) if n < validTime[i,1])
	    print beg, fin
	    #plt.plot(tsObj.times[beg:fin], tsObj.ts[beg:fin],'g--', alpha=0.7)
	    #plt.plot(tsObj.times[beg:fin], tsObj.levels[beg:fin],'r-', linewidth=2)
	    #ax.set_ylim([-600,1800])
	    temp = tsObj.get_notes(tsObj.get_string_part(validTime[i,0], validTime[i,1])[0])
	    print '{0}\t {1:.2f}\t {2:.2f}\t {3}'.format(len(temp), (fin-beg)/100., (len(temp)*100./(fin-beg)), temp)
	#plt.show()
        
        #print len(tsObj.times), len(tsObj.levels)
        
        #plt.plot(tsObj.times[:len(tsObj.levels)], tsObj.ts[:len(tsObj.levels)])
        plt.plot(tsObj.times[:len(tsObj.levels)], tsObj.levels, 'r-', linewidth=2)
        plt.show()

        
        #return breathPheases_s
    
def findValidSvaras(pitch,tonic):
      
    phObj = PH.PitchHistogram(pitch, float(tonic))
    phObj.ComputePitchHistogram(Oct_fold=1)
    phObj.ValidSwarLocEstimation()
    phObj.SwarLoc2Cents()
    svara = phObj.swarCents
    print svara
    phObj.PlotHistogram()
    svaraSemitone = map(lambda svara:round(svara,-2),svara)
    #print svaraSemitone
    octaveNotes = range(0,1101,100)
    ignoreNoteLevels = list(set(octaveNotes).difference(svaraSemitone))
    noteNames = ['S','r','R','g','G','m','M','P','d','D','n','N']
    ignoreNotes = []
    for i in range(len(ignoreNoteLevels)):
        ignoreNotes.append(noteNames[octaveNotes.index(ignoreNoteLevels[i])])        
    return svaraSemitone, ignoreNotes