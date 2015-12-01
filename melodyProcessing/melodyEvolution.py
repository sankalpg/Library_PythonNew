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
    
    for filename in filenames[:]:
        print "Processing file %s"%filename
        
        #======================
        ## This is done for all
        #======================
        '''
        fname, ext = os.path.splitext(filename)
        pitch, time, Hop = BO.readPitchFile(fname + pitchExt)
        tonic = np.loadtxt(fname  + tonicExt)
        pcents = BO.PitchHz2Cents(pitch, tonic)
        pdata = (time,pcents,Hop)
        
        ## Extract Breath Phrases
        #------------------------
        breathPhrases = findBreathPhrases(segObj,fname,pcents,Hop)
        
        ## Histogram processing to extract note locations
        #------------------------------------------------
        svaraSemitone, ignoreNotes = findValidSvaras(pitch,tonic)
        #print svaraSemitone, ignoreNotes
        print "Notes being ignored are: %s" % ignoreNotes
        
        ## Svara transcription
        #---------------------
        transcription = transcribePitch(fname,pdata,ignoreNotes)
        #print transcription
        
        print "-------\nDone !!\n-------"
        '''
        
        getSvaraDistInBPs(filename)
        plotSvaraDistInBPs(filename)
        
        
        
        
def findBreathPhrases(segObj,fname,pcents,Hop):
  
    breathPhrases = np.array(segObj.ExtractBreathPhrases(pcents, Hop, 0.5)).astype(np.float)
    print "Extracting breath phrases..."
    validTime = breathPhrases
    validTime[:,0] = validTime[:,0]*Hop
    validTime[:,1] = validTime[:,1]*Hop

    bphraseFilename = (''.join([fname,'.bphrases']))
    np.savetxt(bphraseFilename,validTime, delimiter = "\t", fmt="%f\t%f\t%d")
    return breathPhrases
    
    
def findValidSvaras(pitch,tonic):
      
    phObj = PH.PitchHistogram(pitch, float(tonic))
    phObj.ComputePitchHistogram(Oct_fold=1)
    phObj.ValidSwarLocEstimation()
    phObj.SwarLoc2Cents()
    svara = phObj.swarCents
    print "Computing svara from histogram..."
    #print svara
    #phObj.PlotHistogram()
    svaraSemitone = map(lambda svara:round(svara,-2),svara)
    #print svaraSemitone
    octaveNotes = range(0,1101,100)
    ignoreNoteLevels = list(set(octaveNotes).difference(svaraSemitone))
    noteNames = ['S','r','R','g','G','m','M','P','d','D','n','N']
    ignoreNotes = []
    for i in range(len(ignoreNoteLevels)):
        ignoreNotes.append(noteNames[octaveNotes.index(ignoreNoteLevels[i])])        
    return svaraSemitone, ignoreNotes
  
  
def transcribePitch(fname,pdata,ignoreNotes):  
  
    tsObj = ts.SvarTranscription(pdata)
    transcription = np.array(tsObj.perform_transcription(ignoreNotes, thres=60, width=35, verbose=False)) 
    #print transcription
    transcriptionFilename = (''.join([fname,'.transcription']))
    np.savetxt(transcriptionFilename,transcription, delimiter = "\t", fmt="%f\t%f\t%d")
    return transcription
    



def getSvaraDistInBPs(filename, bphraseExt = '.bphrases', transExt = '.transcription'):
    """
    This function computes svara distribution for all the breath phrases
    """
    svarsBP = getSvarasInBreathPhrases(filename, bphraseExt = bphraseExt, transExt = transExt)
    
    svar2binMap = {}
    for ii, val in enumerate(range(-12, 36)):
        svar2binMap[int(val*100)] = ii
    
    dist_mtx = np.zeros((len(svar2binMap.keys()), len(svarsBP)))
    track_max = np.zeros(len(svarsBP))
    #print dist_mtx.shape
    for ii, bp in enumerate(svarsBP):
        for svr in bp:
            #print svar2binMap[svr['svar']], svr['svar'], ii
            dist_mtx[svar2binMap[svr['svar']], ii] =  svr['duration']
        if np.max(dist_mtx[:,ii]) != 0:
            dist_mtx[:,ii] = dist_mtx[:,ii]/np.max(dist_mtx[:,ii])
            track_max[ii] = np.argmax(dist_mtx[:,ii])
    #print track_max
        
    return dist_mtx, track_max
            
    
def plotSvaraDistInBPs(filename, bphraseExt = '.bphrases', transExt = '.transcription', plotName = -1):
    """
    This function plots blah blah
    """
    
    fig = plt.figure(figsize=(15,60), dpi=80)
    mtx, track_max = getSvaraDistInBPs(filename, bphraseExt = bphraseExt, transExt = transExt)
    plt.imshow(mtx, interpolation = 'nearest', origin = 'lower', cmap = plt.get_cmap('OrRd'))
    #plt.imshow(mtx, interpolation = 'nearest', origin = 'lower', extent = [0,mtx.shape[1],6,30])
    #plt.plot(track_max, 'ws', linewidth=1.0)
    #plt.xlim([0,mtx.shape[1]])
    #plt.ylim([0,48])
    
    fname, ext = os.path.splitext(filename)
    svaraDistFilename = (''.join([fname,'_svaraDist','.pdf']))
    plt.savefig(svaraDistFilename, bbox_inches='tight')
    #plt.show()
    
    
def getSvarasInBreathPhrases(filename, bphraseExt = '.bphrases', transExt = '.transcription'):
    """
    This function parses all the svaras within breath phrases
    """
    fname, ext = os.path.splitext(filename)
    bpFile = fname + bphraseExt
    transFile = fname + transExt
    
    bphrases = np.loadtxt(bpFile)
    svaras = np.loadtxt(transFile)
    output = []
    for bp in bphrases:
        ind_start = np.where(svaras[:,0] >= bp[0])[0]
        ind_end = np.where(svaras[:,1] <= bp[1])[0]
        ind_svars = np.intersect1d(ind_start, ind_end)
        bpSvaras = []
        for ii in ind_svars:
            bpSvaras.append({'start': svaras[ii][0], 'end': svaras[ii][1], 'duration': svaras[ii][1]-svaras[ii][0], 'svar': svaras[ii][2]})
        output.append(bpSvaras)
    
    return output
    