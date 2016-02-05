import numpy as np
import scipy.signal as sig
from scipy.io import loadmat
import string
import math
eps = np.finfo(np.float).eps
import sys, os
import matplotlib.pyplot as plt
sys.path.append(os.path.join(os.path.dirname(__file__), '../batchProcessing'))
from scipy.interpolate import interp1d
import pickle
import basicOperations as BO
import batchProcessing as BP
import pitchHistogram as PH
import segmentation as seg
import transcription as ts
from mutagen import easyid3
from mutagen.mp3 import MP3
import copy

svar2binMap = {}
bin2svarMap = {}
for ii, val in enumerate(range(-12, 36)):
    svar2binMap[int(val*100)] = ii
    bin2svarMap[ii] = int(val*100)

svar2binMap_of = {}
bin2svarMap_of = {}    
for ii, val in enumerate(range(0, 12)):
    svar2binMap_of[int(val*100)] = ii
    bin2svarMap_of[ii] = int(val*100)
    
def get_mbid_from_mp3(mp3_file):
   """
   fetch MBID form an mp3 file
   """
   try:
       mbid = easyid3.ID3(mp3_file)['UFID:http://musicbrainz.org'].data
   except:
       print "problem reading mbid for file %s\n" % mp3_file
       raise
   return mbid
 
def wrapperSalientContourFeatures(root_dir, audioExt = '.mp3'):
  """
  """
  filenames = BP.GetFileNamesInDir(root_dir, '.mp3')
  mbid_to_ragaId_map = get_raga_id('/home/kaustuv/Documents/ISMIR-2016/MelodyEvolution/Dataset/mbid_ragaid_Hindustani_300_selected.tsv')
  raga_info = get_vadi_samvadi('/home/kaustuv/Documents/ISMIR-2016/MelodyEvolution/Dataset/vadiSamvadiMapping.tsv')
  '''
  for filename in filenames[:]:
    print "Processing file %s" %filename
    
    fname, ext = os.path.splitext(filename)
    
    endTime = readValidVistarRegion(fname)
    dist_mtx, track_max = getSvaraDistInBPsCum(filename, endTime, bphraseExt = '.bphrases', transExt = '.transcription', windowSize = 10, hopSize = 1)
    track_max_pklfile = fname + '_svaraContour.pkl'
    pickle.dump(track_max, open(track_max_pklfile,'w'))
    
    mbid = get_mbid_from_mp3(filename)    
    raga_id = mbid_to_ragaId_map[mbid]
    
    feature_pklfile = fname + '_contourNormFeatures.pkl'
    extractFeaturesNoteEvolutionContour(track_max_pklfile, feature_pklfile, raga_info[raga_id])
  '''      
  cnt = 0     
  for filename in filenames[:]:
    print "Processing file %s" %filename
    
    fname, ext = os.path.splitext(filename)
    feature_pklfile = fname + '_contourNormFeatures.pkl'
    
    output = pickle.load(open(feature_pklfile, 'r'))
    mbid = get_mbid_from_mp3(filename)    
    raga_id = mbid_to_ragaId_map[mbid]
    
    feature_order = [output['F1'], output['F2'], output['F3'], output['F4'], output['F5'], output['F6'], output['F7'], output['F8'], output['F9'], output['F10'], output['F11'], output['F12'], output['F13'], output['F14'], output['F15'], output['F16'], output['F17'], output['F18'], output['F19'], output['F20'], raga_id]
     
    if cnt == 0:
      feature_overall = np.array([feature_order])
    else:
      feature_overall = np.vstack((feature_overall, feature_order))
    cnt += 1
   
  print feature_overall
  feature_filename = os.path.join(root_dir, 'feature_for_weka.csv')
  fid = open(feature_filename,'w')
  for row in feature_overall:
    row = [str(r) for r in row]
    fid.write("%s,"*len(row)%tuple(row))
    fid.write('\n')
  fid.close()
  
        
  
def get_raga_id(mbid_raga_mapping):
  """
  """
  lines = open(mbid_raga_mapping,'r').readlines()
  
  mbid_to_ragaId_map = {}
  
  for ii, line in enumerate(lines):
    sline = line.split(',')
    sline = [s.strip() for s in sline]
    if not mbid_to_ragaId_map.has_key(sline[0]):
      mbid_to_ragaId_map[sline[0]] = sline[1]
  
  return mbid_to_ragaId_map


def get_vadi_samvadi(raga_vadi_samvadi_mapping):
  """
  """
  lines = open(raga_vadi_samvadi_mapping,'r').readlines()
  
  raga_to_vadi_map = {}
  
  for ii, line in enumerate(lines):
    sline = line.split(',')
    sline = [s.strip() for s in sline]
    if not raga_to_vadi_map.has_key(sline[0]):
      raga_to_vadi_map[sline[0]] = {'vadi':sline[1], 'samvadi':sline[2]}
  
  return raga_to_vadi_map
 

def batchProc(root_dir, audioExt = '.mp3', pitchExt = '.pitchSilIntrpPP', tonicExt = '.tonicFine'):
    
    filenames = BP.GetFileNamesInDir(root_dir, '.mp3')
    segObj = seg.melodySegmentation()
    
    #fig = plt.figure(figsize=(15,10), dpi=80)
    
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
        #endTime = readValidVistarRegion(fname)
        
        ## Svara transcription
        #---------------------
        transcription = transcribePitch(fname,pdata,ignoreNotes)
        #print transcription
        
        print "-------\nDone !!\n-------"
        
        '''
        plotSvaraDistInBPs(filename, endTime, bphraseExt = '.bphrases', transExt = '.transcription', plotName = -1)
        for windowSize, hopSize in [[10,1]]:
	    try:
                plotSvaraDistInBPsCum(filename, endTime, bphraseExt = '.bphrases', transExt = '.transcription', plotName = -1, windowSize = windowSize, hopSize = hopSize)
            except:
	        print "*** Problem occured in file: %s" %filename
	        pass
        getBreathPhraseStatistics(filename, endTime, bphraseExt = '.bphrases', transExt = '.transcription')
        
        try:
	    plotSalientSvaraContour(filename, endTime, bphraseExt = '.bphrases', transExt = '.transcription', plotName = -1, windowSize = 10, hopSize = 1)
	except:
	    print "*** Problem occured in file: %s" %filename
	    pass
        
        print "-------\nDone !!\n-------"
        '''
    #plt.show()
	       
        
        
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
  
  
def transcribePitch(fname,pdata,ignoreNotes):  
  
    tsObj = ts.SvarTranscription(pdata)
    transcription = np.array(tsObj.perform_transcription(ignoreNotes, thres=60, width=35, verbose=False)) 
    #print transcription
    transcriptionFilename = (''.join([fname,'.transcription']))
    np.savetxt(transcriptionFilename,transcription, delimiter = "\t", fmt="%f\t%f\t%d")
    return transcription
    



def getSvaraDistInBPs(filename, endTime, bphraseExt = '.bphrases', transExt = '.transcription'):
    """
    This function computes svara distribution for all the breath phrases
    """
    svarsBP = getSvarasInBreathPhrases(filename, endTime, bphraseExt = bphraseExt, transExt = transExt)
    
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



def getSvaraDistInBPsCum(filename, endTime, bphraseExt = '.bphrases', transExt = '.transcription', windowSize = 5, hopSize = 1):
    """
    This function computes svara distribution for all the breath phrases
    """
    svarsBP = getSvarasInBreathPhrases(filename, endTime, bphraseExt = bphraseExt, transExt = transExt)
    
    nBps = len(svarsBP)
    nFrames = int(math.floor((nBps- windowSize)/hopSize) + 1)
    dist_mtx_cum = np.zeros((len(svar2binMap.keys()), nFrames))
    track_max = np.zeros(nFrames)
    
    #print nFrames, nFrames*hopSize
    
    for ii in range(0, nFrames):
        for jj in range(windowSize):
            for svr in svarsBP[ii*hopSize+jj]:
                dist_mtx_cum[svar2binMap[svr['svar']], ii] +=  svr['duration']
        
        if np.max(dist_mtx_cum[:,ii]) != 0:
            dist_mtx_cum[:,ii] = dist_mtx_cum[:,ii]/np.max(dist_mtx_cum[:,ii])
            track_max[ii] = np.argmax(dist_mtx_cum[:,ii])
          
    track_max = sig.medfilt(track_max,5)
    track_max = np.trim_zeros(track_max, trim='b')  
            
    return dist_mtx_cum, track_max


def plotSvaraDistInBPsCum(filename, endTime, bphraseExt = '.bphrases', transExt = '.transcription', plotName = -1, windowSize = 5, hopSize = 1):
    """
    This function plots blah blah
    """
    
    fig = plt.figure(figsize=(15,10), dpi=80)
    mtx, track_max = getSvaraDistInBPsCum(filename, endTime, bphraseExt = bphraseExt, transExt = transExt, windowSize = windowSize, hopSize = hopSize)
    plt.imshow(mtx, interpolation = 'nearest', origin = 'lower', cmap = plt.get_cmap('OrRd'))
    plt.plot(track_max, 'c', linewidth = 3.0)

    
    fname, ext = os.path.splitext(filename)
    print "Plotting svara distribution (salience) for %d bp's and %d step(s)..."%(windowSize, hopSize)
    #plt.show()
    saveFigure(fig, fname, featureName = (''.join(['_svaraDistCum_', str(windowSize), '_', str(hopSize)])))
    
            
    
def plotSvaraDistInBPs(filename, endTime, bphraseExt = '.bphrases', transExt = '.transcription', plotName = -1):
    """
    This function plots blah blah
    """
    
    fig = plt.figure(figsize=(15,10), dpi=80)
    mtx, tx_mtx = getSvaraDistInBPs(filename, endTime, bphraseExt = bphraseExt, transExt = transExt)
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
    
    
def normalizeMinMaxUnity(array):
  min_val = np.min(array)
  y_range = np.ptp(array, axis=0)
  return (array - min_val)/y_range, min_val, y_range 



def plotSalientSvaraContour(filename, endTime, bphraseExt = '.bphrases', transExt = '.transcription', plotName = -1, windowSize = 10, hopSize = 1):
  """
  This function plots the median filtered contour of the most salient note per bp within a [0,1] square matrix
  """
  
  #fig = plt.figure(figsize=(15,10), dpi=80)
  mtx, track_max = getSvaraDistInBPsCum(filename, endTime, bphraseExt = bphraseExt, transExt = transExt, windowSize = windowSize, hopSize = hopSize)
  
  track_max = (track_max - np.min(track_max)) / np.ptp(track_max, axis=0)
  hor_axis = np.linspace(0.0, 1.0, num = len(track_max), retstep = True)
  plt.plot(hor_axis[0], track_max, linewidth = 2.0)
  
  fname, ext = os.path.splitext(filename)
  print "Plotting svara contour"
  #plt.show()
  #saveFigure(fig, fname, featureName = '_salientSvaraContour')
    
    
    
def getSvarasInBreathPhrases(filename, endTime, bphraseExt = '.bphrases', transExt = '.transcription'):
    """
    This function parses all the svaras within breath phrases
    """
    fname, ext = os.path.splitext(filename)
    bpFile = fname + bphraseExt
    transFile = fname + transExt
    
    bphrases = np.loadtxt(bpFile)
    bphrases = bphrases[:(np.where(bphrases[:,0] >= endTime)[0][0]),:]
    
    svaras = np.loadtxt(transFile)
    svaras = svaras[:(np.where(svaras[:,0] >= endTime)[0][0]),:]
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
    bps = getSvarasInBreathPhrases(filename, endTime, bphraseExt = bphraseExt, transExt = transExt)
    
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
    #plt.show()
    
    
    
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
    svaras = svaras[:(np.where(svaras[:,0] >= endTime)[0][0]),:]
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

def groupIndices(indexes):
    """
    This function groups indexes. This is often needed to produce segments given indices.
    """
    segments = []
    segStart = indexes[0]
    N = len(indexes)
    for ii in range(len(indexes)):
        if ii == N-1:
            segments.append([segStart, indexes[ii]])
            return np.array(segments)
        if indexes[ii]+1 != indexes[ii+1]:
            segments.append([segStart, indexes[ii]])
            segStart = indexes[ii+1]

    return np.array(segments)    

def extractFeaturesNoteEvolutionContour(contour_file, output_file, raga_info = {}, num_samples_x_axis = 100, use_Ynorm = 1):
  """
  Slope-
    F1: start time till start of highest note
    F2: end of first note till start of highest note
    F3: end of highest note to end of file
  Duration-
    F4: Duration of the longest note
    F5: Duration of the second longest note
    F6: Shortest duration of the note
    F7: duration of Sa
    F8: duration of Sa'
    F9: proportion of vadi to total duration
    F10: proportion of samvadi to total duration
    F11: centroid of vadi bps
    F12: centroid of samvadi bps
  Jump-
    F13: #change / #same (from tx-mtx) 
    F14: highest jump
    F15: average jump length
  Levels-
    F16: #levels explored
    F17: starting note
    F18: ending note
  Misc-
    F19: proportion of longest note to total duration
    F20: name of the longest note
    ####F--: if there are arbitrary jumps before settling the next nyas

  """
  num_samples_x_axis = float(num_samples_x_axis)
  contour = pickle.load(open(contour_file, 'r'))
  feature = {}
  
  #before feature extraction normalizing the max contour (unity range, min val =0 and equal number of samples)
  contour_norm = contour
  min_cnt = 0
  y_range = 1
  if use_Ynorm == 1:
    contour_unity, min_cnt, y_range = normalizeMinMaxUnity(contour)
    n_sam_orig = len(contour)
    f = interp1d(np.linspace(0,1,n_sam_orig), contour_unity, kind= 'nearest')
    contour_norm = f(np.linspace(0,1,num_samples_x_axis))
    contour_norm = contour_norm.round(decimals=2)
  
  #slope related features
  max_val = np.max(contour_norm)
  ind_max_start = np.min(np.where(contour_norm == max_val)[0])
  ind_max_end = np.max(np.where(contour_norm == max_val)[0])
  ind_first_note = np.where(contour_norm == contour_norm[0])[0]
  ind_end_first_note = groupIndices(ind_first_note)[0][1]
  feature['F1'] = 100.0*(1-contour_norm[0])/ind_max_start
  feature['F2'] = 100.0*(1-contour_norm[0])/(ind_max_start-ind_end_first_note)
  if ind_max_end == num_samples_x_axis -1:
    feature['F3'] = 0
  else:
    feature['F3'] = 100.0*(contour_norm[-1] - max_val)/(num_samples_x_axis-1-ind_max_end)
  
  
  #Duration related features
  u_vals = np.unique(contour_norm)
  salient_bps = []
  for ii, v in enumerate(u_vals):
    inds = np.where(contour_norm == v)[0]
    salient_bps.append([v, len(inds)])
    segs = groupIndices(inds)
    if ii == 0:
      segments = segs
    else:
      segments = np.vstack((segments, segs))
  salient_bps = np.array(salient_bps)
  durations = segments[:,1]- segments[:,0]
  ind_single_point = np.where(durations ==0)[0]
  durations = np.delete(durations, ind_single_point)
  ind_sort = np.argsort(durations)
  feature['F4'] = durations[ind_sort[-1]]
  feature['F5'] = durations[ind_sort[-2]]
  feature['F6'] = durations[ind_sort[0]]
  
  ind_sort = np.argsort(salient_bps[:,1])
  feature['F19'] = salient_bps[ind_sort[-1],1]/num_samples_x_axis
  feature['F20'] = bin2svarMap[round((salient_bps[ind_sort[-1],0]*y_range) + min_cnt)]%1200
  
  sa_norm = (svar2binMap[0]-min_cnt)/y_range
  sa_norm = round(sa_norm, 2)
  
  sa_high_norm = (svar2binMap[1200]-min_cnt)/y_range
  sa_high_norm = round(sa_high_norm, 2)
  
  ind_sa = np.where(contour_norm == sa_norm)[0]
  ind_sa_high = np.where(contour_norm == sa_high_norm)[0]
  
  feature['F7'] = len(ind_sa)/num_samples_x_axis
  feature['F8'] = len(ind_sa_high)/num_samples_x_axis
  
  
  for ii, val in enumerate([-1200.0, 0.0, 1200.0]):
    vadi_norm = (float(svar2binMap[val+float(raga_info['vadi'])])-min_cnt)/y_range
    vadi_norm = round(vadi_norm, 2)
    if ii == 0:
      ind_vadi = np.where(contour_norm == vadi_norm)[0]
    else:
      ind_vadi = np.append(ind_vadi, np.where(contour_norm == vadi_norm)[0])
  
  feature['F9'] = len(ind_vadi)/num_samples_x_axis
  feature['F11'] = np.sum(ind_vadi)/len(ind_vadi)
  
  for ii, val in enumerate([-1200.0, 0.0, 1200.0]):
    sam_vadi_norm = (float(svar2binMap[val+float(raga_info['samvadi'])])-min_cnt)/y_range
    sam_vadi_norm = round(sam_vadi_norm ,2)
    if ii == 0:
      ind_sam_vadi = np.where(contour_norm == sam_vadi_norm)[0]  
    else:
      ind_sam_vadi = np.append(ind_sam_vadi, np.where(contour_norm == sam_vadi_norm)[0])
  
  feature['F10'] = len(ind_sam_vadi)/num_samples_x_axis
  feature['F12'] = np.sum(ind_sam_vadi)/len(ind_sam_vadi)
  
  
  #Jump related features
  diff = abs(contour_norm[1:] - contour_norm[:-1])
  inds = np.where(diff >0)[0]
  feature['F13'] = len(inds)/(num_samples_x_axis - len(inds))
  
  feature['F14'] = np.max(diff)
  feature['F15'] = np.mean(diff[inds])
  
  
  #Levels
  feature['F16'] = len(u_vals)
  feature['F17'] = bin2svarMap[contour[0]]
  feature['F18'] = bin2svarMap[contour[-1]]
  
  
  pickle.dump(feature, open(output_file, 'w'))
  
  
def dumpLongestSvaraInfo(filename, endTime, bphraseExt = '.bphrases', transExt = '.transcription'):
  """
  """
  out = getSvarasInBreathPhrases(filename, endTime, bphraseExt = '.bphrases', transExt = '.transcription')
  print len(out)
  
  st = []
  for ii, val in enumerate(out):
    durs = []
    for jj in range(len(out[ii])):
      durs.append(out[ii][jj]['duration'])
    if len(durs) != 0:
      ind = np.argmax(durs)
      pos = len(durs)
      length = np.max(durs)
      start = out[ii][ind]['start']
      end = out[ii][ind]['end']
      
      st.append(start)
  
  st = np.array(st)
  diff = np.diff(st)
  
  resolution = 50.0
  gcds_fuzzy = np.zeros((1000.0/resolution)*np.max(diff)+10)
  for d in diff:
    
    sub_inds = np.round((1000.0/resolution)*d*np.array([1, 1/2.0, 1/3.0, 1/4.0, 1/5.0, 1/6.0, 1/7.0, 1/8.0])).astype(np.int)
    print d, sub_inds
    gcds_fuzzy[sub_inds] += 1
    gcds_fuzzy[sub_inds+1] += 0.2
    gcds_fuzzy[sub_inds-1] += 0.2
    
    gcds_fuzzy[sub_inds+2] += 0.05
    gcds_fuzzy[sub_inds-2] += 0.05
  
  print gcds_fuzzy
  plt.plot((resolution/1000.0)*np.arange(len(gcds_fuzzy)), gcds_fuzzy)
  plt.show()
  
  
 
def getDurDistPerNotePerFile(filename, endTime, bphraseExt = '.bphrases', transExt = '.transcription'):
  """
  """
  #out = getSvarasInBreathPhrases(filename, endTime, bphraseExt = '.bphrases', transExt = '.transcription')
  
  fname, ext = os.path.splitext(filename)
  #bpFile = fname + bphraseExt
  transFile = fname + transExt
  
  #bphrases = np.loadtxt(bpFile)
  #bphrases = bphrases[:(np.where(bphrases[:,0] >= endTime)[0][0]),:]
  
  svaras = np.loadtxt(transFile)
  svaras = svaras[:(np.where(svaras[:,0] >= endTime)[0][0]),:]
  
  validSvaras = np.unique(svaras[:,2])
  svaraHist = np.zeros(((len(svar2binMap.keys())),50))
  
  for ii in svar2binMap.keys():
    temp = []
    for svar in svaras:
      if svar[2] == ii:
        #svaraHist[svar2binMap[ii], cnt] = (svar[1] - svar[0])
        temp.append((svar[1] - svar[0]))
       
    svaraHist[svar2binMap[ii]] = np.histogram(temp, bins = (np.arange(51))*0.2)[0]
    #print svaraHist
  #plt.hold(True)
  plt.imshow(svaraHist, interpolation = 'nearest', origin = 'lower', cmap = plt.get_cmap('OrRd'))
  plt.show()
  

def getNoteDurDistPerRaga(root_dir, audioExt = '.mp3', pitchExt = '.pitchSilIntrpPP', tonicExt = '.tonicFine'):
  """
  """
  filenames = BP.GetFileNamesInDir(root_dir, '.mp3')
  
  for filename in filenames[:]:
    print "Processing file %s" %filename
  
    fname, ext = os.path.splitext(filename)
    endTime = readValidVistarRegion(fname)
    print endTime
    
    getDurDistPerNotePerFile(filename, endTime, bphraseExt = '.bphrases', transExt = '.transcription')  
  plt.show()
  
  
def getSvarFreqHistPerFile(filename, endTime, transExt = '.transcription', Oct_fold = 1):
  
  fname, ext = os.path.splitext(filename)
  
  trans = np.loadtxt(fname + transExt)
  trans = trans[np.where(trans[:,0]<= endTime)[0]]
  
  if Oct_fold == 0:
    hist = np.zeros(len(svar2binMap.keys()))
  elif Oct_fold == 1:
    hist = np.zeros(len(svar2binMap_of.keys()))
  
  for tt in trans:    
    if Oct_fold == 0:
      hist[svar2binMap[tt[2]]]+=1
    elif Oct_fold == 1:
      hist[svar2binMap_of[tt[2]%1200]]+=1
  
  return hist
  
  
  
def getLongestNoteDurDist(filename, endTime,  bphraseExt = '.bphrases', transExt = '.transcription', group1_svars = [], group2_svars = []):
    """
    This function computes distribution of long note durations within a breath phrase based on two groups inputted
    
    NOTE: swars should be input in octave folded way in cents
    """
    
    fname, ext = os.path.splitext(filename)
  
    out = getSvarasInBreathPhrases(filename, endTime, bphraseExt = bphraseExt, transExt = transExt)
    
    long_svars = {}
    for kk in svar2binMap_of.keys():
      long_svars[kk] = []    
    
    for bp in out:
      durs = []
      for svar in bp:
	durs.append(svar['duration'])
      if len(durs) >0:
	ind = np.argmax(durs)
	long_svars[bp[ind]['svar']%1200].append(bp[ind]['duration'])
    
    durs_grp1 = []
    durs_grp2 = []
    
    for g1 in group1_svars:
      durs_grp1.extend(long_svars[g1])
    
    for g2 in group2_svars:
      durs_grp2.extend(long_svars[g2])
    
    durs_grp1 = np.array(durs_grp1)
    durs_grp2 = np.array(durs_grp2)
    
    cm1 = np.sum(durs_grp1*np.arange(1,len(durs_grp1)+1))/np.sum(durs_grp1)
    cm2 = np.sum(durs_grp2*np.arange(1,len(durs_grp2)+1))/np.sum(durs_grp2)
    
    
    return durs_grp1, durs_grp2, cm1, cm2
  
def plotLongestNoteDurDist(filename, endTime,  bphraseExt = '.bphrases', transExt = '.transcription', group1_svars = [], group2_svars = []):
  
  g1, g2 = getLongestNoteDurDist(filename, endTime,  bphraseExt = bphraseExt, transExt = transExt, group1_svars = group1_svars, group2_svars = group2_svars)
  
  bins = np.arange(50)*.1
  hist1 = np.histogram(g1, bins = bins )
  hist2 = np.histogram(g2, bins = bins )
  
  plt.hold(True)
  plt.plot(hist1[0], 'r')
  plt.plot(hist2[0], 'b')
  
  #plt.show()
    
    
def getDurStdSvarDurationPerLongestSvarOfBP(filename, endTime,  bphraseExt = '.bphrases', transExt = '.transcription'):
  
    svar_dict = {}
    for kk in svar2binMap_of.keys():
      svar_dict[kk] = [] 
    
    long_svar_dict = {}
    for kk in svar2binMap_of.keys():
      long_svar_dict[kk] = copy.deepcopy(svar_dict)
    
    out = getSvarasInBreathPhrases(filename, endTime, bphraseExt = bphraseExt, transExt = transExt)
    
    for bp in out:
      durs = []
      svar_class = []
      for svar in bp:
	durs.append(svar['duration'])
	svar_class.append(svar['svar'])
      if len(durs) >0:
	ind = np.argmax(durs)
	longest_note = bp[ind]['svar']%1200
	for ii, sc in enumerate(svar_class):
	  long_svar_dict[longest_note][sc%1200].append(durs[ii])
    return long_svar_dict
	  
      
  