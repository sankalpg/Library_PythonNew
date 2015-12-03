import numpy as np
import scipy.signal as sig
from scipy.io import loadmat
import string
import math
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
    
    for filename in filenames[:1]:
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

    for ii, bp in enumerate(svarsBP):
        for svr in bp:
            #print svar2binMap[svr['svar']], svr['svar'], ii
            dist_mtx[svar2binMap[svr['svar']], ii] =  svr['duration']
        if np.max(dist_mtx[:,ii]) != 0:
            dist_mtx[:,ii] = dist_mtx[:,ii]/np.max(dist_mtx[:,ii])
            
    return dist_mtx



def getSvaraDistInBPsTemp(filename, bphraseExt = '.bphrases', transExt = '.transcription', windowSize = 5, hopSize = 1):
    """
    This function computes svara distribution for all the breath phrases
    """
    svarsBP = getSvarasInBreathPhrases(filename, bphraseExt = bphraseExt, transExt = transExt)
    
    svar2binMap = {}
    for ii, val in enumerate(range(-12, 36)):
        svar2binMap[int(val*100)] = ii
    
    nBps = len(svarsBP)
    nFrames = int(math.floor((nBps- windowSize)/hopSize) + 1)
    dist_mtx_cum = np.zeros((len(svar2binMap.keys()), nFrames))
    
    print nFrames, nFrames*hopSize
    
    for ii in range(0, nFrames):
        for jj in range(windowSize):
            for svr in svarsBP[ii*hopSize+jj]:
                dist_mtx_cum[svar2binMap[svr['svar']], ii] =  svr['duration']
        
        if np.max(dist_mtx_cum[:,ii]) != 0:
            dist_mtx_cum[:,ii] = dist_mtx_cum[:,ii]/np.max(dist_mtx_cum[:,ii])
            
    return dist_mtx_cum


def plotSvaraDistInBPsTemp(filename, bphraseExt = '.bphrases', transExt = '.transcription', plotName = -1, windowSize = 5, hopSize = 1):
    """
    This function plots blah blah
    """
    
    fig = plt.figure(figsize=(15,60), dpi=80)
    mtx = getSvaraDistInBPsTemp(filename, bphraseExt = bphraseExt, transExt = transExt, windowSize = windowSize, hopSize = hopSize)
    plt.imshow(mtx, interpolation = 'nearest', origin = 'lower', cmap = plt.get_cmap('OrRd'))

    
    #fname, ext = os.path.splitext(filename)
    #svaraDistFilename = (''.join([fname,'_svaraDist','.pdf']))
    #plt.savefig(svaraDistFilename, bbox_inches='tight')
    plt.show()
    
            
    
def plotSvaraDistInBPs(filename, bphraseExt = '.bphrases', transExt = '.transcription', plotName = -1):
    """
    This function plots blah blah
    """
    
    fig = plt.figure(figsize=(15,60), dpi=80)
    mtx = getSvaraDistInBPs(filename, bphraseExt = bphraseExt, transExt = transExt)
    plt.imshow(mtx, interpolation = 'nearest', origin = 'lower', cmap = plt.get_cmap('OrRd'))

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
    
    bp_duration = []
    bp_num_notes = []
    bp_num_notes_norm = []
    
    for bp in bphrases:
        ind_start = np.where(svaras[:,0] >= bp[0])[0]
        ind_end = np.where(svaras[:,1] <= bp[1])[0]
        ind_svars = np.intersect1d(ind_start, ind_end)
        
        bpSvaras = []
        for ii in ind_svars:
            bpSvaras.append({'start': svaras[ii][0], 'end': svaras[ii][1], 'duration': svaras[ii][1]-svaras[ii][0], 'svar': svaras[ii][2]})
        output.append(bpSvaras)
    
    return output
    
    
def getBreathPhraseStatistics(filename, bphraseExt = '.bphrases', transExt = '.transcription'):  
    """
    This function computes feature statistics of all the svaras within breath phrases
    """
    fname, ext = os.path.splitext(filename)
    bpFile = fname + bphraseExt
    transFile = fname + transExt
    
    bphrases = np.loadtxt(bpFile)
    svaras = np.loadtxt(transFile)
    svarsBP = []
    
    bp_duration = []
    bp_num_notes = []
    bp_num_notes_norm = []
    bp_onset = []
    
    for bp in bphrases:
        ind_start = np.where(svaras[:,0] >= bp[0])[0]
        ind_end = np.where(svaras[:,1] <= bp[1])[0]
        ind_svars = np.intersect1d(ind_start, ind_end)
        
        bp_duration.append(bp[1] - bp[0])
        bp_num_notes.append(len(ind_svars))
        bp_num_notes_norm.append(len(ind_svars)/(bp[1] - bp[0]))
        bp_onset.append(bp[0])

        
    #hist, bins = plotHist(bp_duration)
    #hist, bins = plotIOIHist(bp_onset)
    
        bpSvaras = []
        for ii in ind_svars:
            bpSvaras.append({'start': svaras[ii][0], 'end': svaras[ii][1], 'duration': svaras[ii][1]-svaras[ii][0], 'svar': svaras[ii][2]})
        svarsBP.append(bpSvaras)
        

    bp_max_note_dur = []
    note_onset_all = []
    for ii, bp in enumerate(svarsBP):
        max_val = -1
        for svr in bp:
            print svr
            note_onset_all.append(svr['start'])
            if svr['duration'] > max_val:
                max_val = svr['duration']
        bp_max_note_dur.append(max_val)
        
    
    #hist, bins = plotIOIHist(note_onset_all)
    
    
    #print len(bp_duration), len(bp_num_notes), len(bp_num_notes_norm), len(bp_max_note_dur)
    
    plt.plot(np.array(bp_duration)+50, 'o--', label="BP duration")
    plt.plot(np.array(bp_num_notes)+30, 's--', label="BP # notes")
    plt.plot(bp_num_notes_norm, 'd--', label="BP # notes norm")
    plt.plot(np.array(bp_max_note_dur)+10, '^--', label="BP max note duration")
    plt.legend(loc='upper left')
    plt.show()
    

def plotHist(parameter, bins = 100, normed=1):
    
    hist, bins = np.histogram(parameter, bins = 100, normed=1)
    width = 0.7 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2
    plt.bar(center, hist, align='center', width=width)
    plt.show()
    return hist, bins


def plotIOIHist(parameter):
    '''
    This function takes aparameter and plots IOI histogram
    '''
    ioi = np.ediff1d(parameter)
    #print ioi 
    
    hist, bins = plotHist(ioi, bins = 100, normed=1)
    return hist, bins
    