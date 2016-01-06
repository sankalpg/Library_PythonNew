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
    
    for filename in filenames[:]:
        print "Processing file %s" %filename
        
        #======================
        ## This is done for all
        #======================
        
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
        
        ## Read valid region for evolution
        #---------------------------------
        endTime = readValidVistarRegion(fname)
        
        ## Svara transcription
        #---------------------
        transcription = transcribePitch(fname,pdata,ignoreNotes, endTime)
        #print transcription
        
        #print "-------\nDone !!\n-------"
        
        
        plotSvaraDistInBPs(filename, bphraseExt = '.bphrases', transExt = '.transcription', plotName = -1)
        for windowSize, hopSize in [[5,1],[10,1]]:
	    try:
                plotSvaraDistInBPsCum(filename, bphraseExt = '.bphrases', transExt = '.transcription', plotName = -1, windowSize = windowSize, hopSize = hopSize)
            except:
	        print "Problem occured in file: %s" %filename
	        pass
        getBreathPhraseStatistics(filename, endTime, bphraseExt = '.bphrases', transExt = '.transcription')
        
        print "-------\nDone !!\n-------"
        
        
        
        
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


def readValidVistarRegion(fname):
    
    validTimeFilename = (''.join([fname,'.endSec']))
    endTime = np.loadtxt(validTimeFilename)
    return endTime
  
  
def transcribePitch(fname,pdata,ignoreNotes,endTime):  
  
    tsObj = ts.SvarTranscription(pdata, endTime)
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
    
    tx_mtx = np.zeros((len(svar2binMap.keys()),len(svar2binMap.keys())))
    notes = []

    for ii, bp in enumerate(svarsBP):
        for svr in bp:
            #print svar2binMap[svr['svar']], svr['svar'], ii
            dist_mtx[svar2binMap[svr['svar']], ii] +=  svr['duration']
            
            notes.append(svar2binMap[svr['svar']])            
            
        if np.max(dist_mtx[:,ii]) != 0:
            dist_mtx[:,ii] = dist_mtx[:,ii]/np.max(dist_mtx[:,ii])
    
    for i in range(len(notes)-1):
        tx_mtx[notes[i], notes[i+1]] += 1 
    
    steadiness_feat = np.sum(tx_mtx)/np.trace(tx_mtx)
    print "Steadiness ratio across bp is: %f" %steadiness_feat
    
    #print tx_mtx
    
    return dist_mtx, tx_mtx



def getSvaraDistInBPsCum(filename, bphraseExt = '.bphrases', transExt = '.transcription', windowSize = 5, hopSize = 1):
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
    
    #print nFrames, nFrames*hopSize
    
    for ii in range(0, nFrames):
        for jj in range(windowSize):
            for svr in svarsBP[ii*hopSize+jj]:
                dist_mtx_cum[svar2binMap[svr['svar']], ii] +=  svr['duration']
        
        if np.max(dist_mtx_cum[:,ii]) != 0:
            dist_mtx_cum[:,ii] = dist_mtx_cum[:,ii]/np.max(dist_mtx_cum[:,ii])
            
    return dist_mtx_cum


def plotSvaraDistInBPsCum(filename, bphraseExt = '.bphrases', transExt = '.transcription', plotName = -1, windowSize = 5, hopSize = 1):
    """
    This function plots blah blah
    """
    
    fig = plt.figure(figsize=(15,10), dpi=80)
    mtx = getSvaraDistInBPsCum(filename, bphraseExt = bphraseExt, transExt = transExt, windowSize = windowSize, hopSize = hopSize)
    plt.imshow(mtx, interpolation = 'nearest', origin = 'lower', cmap = plt.get_cmap('OrRd'))

    
    fname, ext = os.path.splitext(filename)
    #svaraDistCumFilename = (''.join([fname,'_svaraDistCum_', str(windowSize), '_', str(hopSize), '.pdf']))
    #plt.savefig(svaraDistCumFilename, bbox_inches='tight')
    print "Plotting svara distribution (salience) for %d bp's and %d step(s)..."%(windowSize, hopSize)
    #plt.show()
    saveFigure(fig, fname, featureName = (''.join(['_svaraDistCum_', str(windowSize), '_', str(hopSize)])))
    
            
    
def plotSvaraDistInBPs(filename, bphraseExt = '.bphrases', transExt = '.transcription', plotName = -1):
    """
    This function plots blah blah
    """
    
    fig = plt.figure(figsize=(15,10), dpi=80)
    mtx, tx_mtx = getSvaraDistInBPs(filename, bphraseExt = bphraseExt, transExt = transExt)
    plt.imshow(mtx, interpolation = 'nearest', origin = 'lower', cmap = plt.get_cmap('OrRd'))

    fname, ext = os.path.splitext(filename)
    print "Plotting svara distribution (salience) for each bp..."
    #plt.show()
    saveFigure(fig, fname, featureName = '_svaraDist')
    
    fig = plt.figure(figsize=(15,10), dpi=80)
    plt.imshow(tx_mtx, interpolation = 'nearest', origin = 'lower', cmap = plt.get_cmap('OrRd'))
    #plt.show()
    print "Plotting svara transition matrix for each bp..."
    saveFigure(fig, fname, featureName = '_svaraTxMtx')
    
    
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



def plotSvarRatioDistribution(filename, bphraseExt = '.bphrases', transExt = '.transcription', within_bp = 1, plotName = -1):
    """
    This function computes ratio of NC2 possibilities of N notes in each breath phrase. Ratio is taken as length_long/length_short. The ratio is appended for all the BPs in the file
    """
    
    #obtaining svaras in bps 
    bps = getSvarasInBreathPhrases(filename, bphraseExt = bphraseExt, transExt = transExt)
    
    durations = []
    for bp in bps:
        for note in bp:
            durations.append(note['duration'])
    print durations
    plotIOIHist(durations)    
    
    ratio = []
    note_cnt = 0
    if within_bp ==1:
        for bp in bps:
            durations = []
            for note in bp:
                note_cnt+=1
                durations.append(note['duration'])
            for ii in range(len(durations)):
                for jj in range(ii+1, len(durations)):
                    ratio.append(np.max([durations[ii], durations[jj]])/np.min([durations[ii], durations[jj]]))
    else:
        durations = []
        for bp in bps:
            for note in bp:
                note_cnt+=1
                durations.append(note['duration'])
        for ii in range(len(durations)):
            for jj in range(ii+1, len(durations)):
                ratio.append(np.max([durations[ii], durations[jj]])/np.min([durations[ii], durations[jj]]))
                    
    print len(ratio), len(bps), note_cnt
    if within_bp == 0:
        bins = np.linspace(1, 16, 16*64)
    else:
        bins = np.linspace(1, 16, 16*8)
    hist = np.histogram(ratio, bins = bins)
    #plt.bar(hist[0], hist[1][:-1])
    plt.plot(hist[1][:-1], hist[0])
    plt.show()
    
    
    
def getBreathPhraseStatistics(filename, endTime, bphraseExt = '.bphrases', transExt = '.transcription'):  
    """
    This function computes feature statistics of all the svaras within breath phrases
    """
    fname, ext = os.path.splitext(filename)
    bpFile = fname + bphraseExt
    transFile = fname + transExt
    
    bphrases = np.loadtxt(bpFile)
    bphrases = bphrases[:(np.where(bphrases[:,0] >= endTime)[0][0]),:]

    svaras = np.loadtxt(transFile)
    svarsBP = []
    bpSvaraDurDist = []
    bpSvaraDurDist_sorted = []
    
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
			
    
        bpSvaras = []
        bpSvaraDur = []
        for ii in ind_svars:
            bpSvaras.append({'start': svaras[ii][0], 'end': svaras[ii][1], 'duration': svaras[ii][1]-svaras[ii][0], 'svar': svaras[ii][2]})
            bpSvaraDur.append(svaras[ii][1]-svaras[ii][0])
        svarsBP.append(bpSvaras)
        bpSvaraDurDist.append(bpSvaraDur)
        bpSvaraDurDist_sorted.append(np.fliplr(np.atleast_2d(np.sort(np.array(bpSvaraDur))))[0])
        
        
    temp1 = np.array(bpSvaraDurDist_sorted)  
    bpSvaraDurDist_sorted = np.array(temp1.tolist())
    #print bpSvaraDurDist_sorted
    bp_num_notes_max = np.max(bp_num_notes)
    plotBarStacked(fname, bpSvaraDurDist_sorted, bp_num_notes_max)
    
    
    fig = plt.figure(figsize=(15,10), dpi=80)
    hist, bins = plotHist(bp_duration)
    print "Plotting bp duration histogram..."
    #plt.show()
    saveFigure(fig, fname, featureName = '_bpDurHist')
    
    
    fig = plt.figure(figsize=(15,10), dpi=80)
    hist, bins = plotIOIHist(bp_onset)
    print "Plotting bp onset IOIH..."
    #plt.show()
    saveFigure(fig, fname, featureName = '_bpOnsetIOIH')
    
    
    bp_max_note_dur = []
    note_onset_all = []
    for ii, bp in enumerate(svarsBP):
        max_val = -1
        for svr in bp:
            note_onset_all.append(svr['start'])
            if svr['duration'] > max_val:
                max_val = svr['duration']
        bp_max_note_dur.append(max_val)
        
    
    fig = plt.figure(figsize=(15,10), dpi=80)
    hist, bins = plotIOIHist(note_onset_all)
    print "Plotting note onset (overall) IOIH..."
    #plt.show()
    saveFigure(fig, fname, featureName = '_noteOnsetIOIH')
    
    
    #print len(bp_duration), len(bp_num_notes), len(bp_num_notes_norm), len(bp_max_note_dur)
    
    fig = plt.figure(figsize=(15,10), dpi=80)
    plt.plot(np.array(bp_duration)+50, 'o--', label="BP duration")
    plt.plot(np.array(bp_num_notes)+30, 's--', label="BP # notes")
    plt.plot(bp_num_notes_norm, 'd--', label="BP # notes norm")
    plt.plot(np.array(bp_max_note_dur)+10, '^--', label="BP max note duration")
    plt.legend(loc='upper left')
    print "Plotting bp features..."
    #plt.show()
    saveFigure(fig, fname, featureName = '_bpFeature')

    
    
    
def saveFigure(fig, fname, featureName):
    '''
    '''
    #Filename_pdf = (''.join([fname, featureName,'.pdf']))
    Filename_png = (''.join([fname, featureName,'.png']))
    plt.title(featureName)
    #plt.savefig(Filename_pdf, bbox_inches='tight')
    plt.savefig(Filename_png, bbox_inches='tight')    
    

def plotHist(parameter, bins = 100, normed=0):
    
    hist, bins = np.histogram(parameter, bins = 100, normed=normed)
    width = 0.7 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2
    plt.bar(center, hist, align='center', width=width)
    #plt.show()
    return hist, bins


def plotIOIHist(parameter):
    '''
    This function takes aparameter and plots IOI histogram
    '''
    ioi = np.ediff1d(parameter)
    #print ioi 
    
    hist, bins = plotHist(ioi, bins = 100, normed=0)
    return hist, bins
    
    
def plotBarStacked(fname, bpSvaraDurDist_sorted, bp_num_notes_max):
    '''
    This function is to plot stacked duration of notes in each breath phrase
    '''
    N = len(bpSvaraDurDist_sorted)
    dur_mtx = np.zeros(shape=(N, bp_num_notes_max))
	
    for rowind, row in enumerate(bpSvaraDurDist_sorted):
        try:
	    dur_mtx[rowind, :len(row)] = np.array(row)
	except:
	    pass	
	    
    dur_mtx_tr = np.zeros(shape=(dur_mtx.shape[1],dur_mtx.shape[0]))
    for i in range(bp_num_notes_max):
        dur_mtx_tr[i,:] = dur_mtx[:,i]
    
    dur_mtx_tr[:,:] = dur_mtx_tr[::-1,:]    
    
    fig = plt.figure(figsize=(15,10), dpi=80)
    stackBarChart(dur_mtx_tr)
    print "Plotting bp note distribution stacked bar graph..."
    saveFigure(fig, fname, featureName = '_noteDistBarStacked')

    
def stackBarChart(dur_mtx_tr):
    '''
    Stack bar with offset of the cumulative sum till last row
    '''
    data = dur_mtx_tr
    width = 0.45
    color = iter(plt.cm.Paired(np.linspace(0,1,data.shape[0])))
    c = [color.next() for ii in range(data.shape[0])]
    for ii in range(data.shape[0]-1,-1,-1):
        pl = plt.bar(range(data.shape[1]), data[ii, :], width, color=c[ii], bottom=np.sum(data[ii+1:,:],axis=0))
    #plt.show()
    