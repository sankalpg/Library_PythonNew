import numpy as np
import scipy.signal as sig
eps = np.finfo(np.float).eps
import sys, os
import matplotlib.pyplot as plt
sys.path.append(os.path.join(os.path.dirname(__file__), '../batchProcessing'))

import batchProcessing as BP
import segmentation as seg


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
  

########################### Pitch sequence resampling ##################################

def resamplePitchSequence(pitch, upSampleFactor, silVal, tonic=-1, hopSize=-1):
    """
    This function resamples the pitch sequence 
    Input: 
        pitch: numpy array of pitch sequence
        hopSize: hopSize in terms of samples for the current pitch sequence
        upSampleFactor: factor by which the user wants to upsample the pitch sequence, if <1 its a downsampling factor in essense
        silVal: the value of the silence in the pitch sequence, needed to take care of such regions in interpolation
    output: 
        pitchOut: the upSampled (this could be a fraction as well for downsampling) pitch sequence as a numpy array
    """
    
    if type(pitch)==str:#passed a file name
      timePitch = np.loadtxt(pitch)
      pitch = timePitch[:,1]
      hopSize = timePitch[1,0]-timePitch[0,0]
    else:
      if type(hopSize)==int and hopSize ==-1:
        print "Please provide a valid hopsize if you want to input pithc as a ndarray"
        return -1
        
    if type(tonic)==str:
      tonic = float(np.loadtxt(tonic))
    
    indSil = np.where(pitch<=silVal)[0]  
    nSamples = pitch.size
    timeArrIn = hopSize*np.arange(nSamples)
    timeArrOut = hopSize*upSampleFactor*np.arange(int(np.ceil(nSamples/upSampleFactor)))
    
    pitchOut = np.interp(timeArrOut, timeArrIn, pitch)
    
    interpFactor = float(1/upSampleFactor)
    
    indSilNew1 = np.floor(indSil*interpFactor).astype(np.int)
    indSilNew2 = np.ceil(indSil*interpFactor).astype(np.int)
    
    indSil = np.intersect1d(indSilNew1, indSilNew2)
    
    if upSampleFactor <1:
      #note that when we do upsampling, there will always be samples at the onset and offset of a valid pitched segment which we dont want. There is an obvious reason for why will they exist, think?
      #so group the silence indices and subtract 1 from the start and add 1 to the end
      silSegs = groupIndices(indSil)
      silSegs[:,0] = silSegs[:,0]-1
      silSegs[:,1] = silSegs[:,1]+1
      
      #converting them again to indices
      indSil = []
      for seg in silSegs:
        indSil.append(range(seg[0], seg[1]+1))
      indSil = np.array(indSil)
      indNeg = np.where(indSil<0)[0]
      indSat = np.where(indSil>pitchOut.size-1)[0]
      indSil = np.delete(indSil, indNeg)
      indSil = np.delete(indSil, indSat)     
    
    pitchOut[indSil] = silVal  
    
    return pitchOut
    

###########################1) Median filtering pitch contour ###########################
#Motivation
#1) Remove spurious jumps in pitch due to small pitch errors
#2) Supress the extent of ornaments like Kan swar which arise challenges to segmentation process or in identification of stable note regions
#3) In general to smooth out pitch content supressing fine movements which might not be relevant for melodic similarity.

def medianFilterPitchContour(pitch, tonic=-1, hopSize=-1, filtLen = 0.1):
  """
  This function median filters the pitch contour to remove spurious jumps and to supress some of the transient like ornamentations which pose challenge for melody processing
  
  Input: 
    pitch: numpy array of the pitch values
    tonic: tonic of the lead singer, if provided the processing is done in cents domain
    hopSize: hop size of the pitch sequence (seconds)
    filtLen: length of the median filter to be used (seconds)
    
  Output:
    pitchOut: median filtered pitch sequence
  """

  if type(pitch)==str:#passed a file name
    timePitch = np.loadtxt(pitch)
    pitch = timePitch[:,1]
    hopSize = timePitch[1,0]-timePitch[0,0]

  if type(tonic)==str:
    tonic = float(np.loadtxt(tonic))
  
  indSil = np.where(pitch<=0)[0]
  
  filtSize = np.round((float(filtLen)/float(hopSize))).astype(np.int)
  filtSize = filtSize + 1 - np.mod(filtSize,2)  # to make the filter length odd number of samples
  
  if type(tonic)!=int:
    pitchCents  = 1200*np.log2((pitch+eps)/tonic)
  else:
    pitchCents = pitch
  
  pitchOut = sig.medfilt(pitchCents, filtSize)

  if type(tonic)!=int:
    pitchOut = tonic*np.power(2,pitchOut/1200)
  
  pitchOut[indSil] = 0
  
  return pitchOut

def batchProcessMedianFilterPitchContour(root_dir, searchExt = '.wav', pitchExt= '.tpe', tonicExt = '.tonic', outExt = '.tpeMedFilt', filtDur = .05):
  """
  Simple wrapper for batch processing of median filtering
  """
  filenames = BP.GetFileNamesInDir(root_dir, searchExt)
  for filename in filenames:
    fname, ext = os.path.splitext(filename)
    timePitch = np.loadtxt(fname+pitchExt)
    hopSize = timePitch[1,0]-timePitch[0,0]
    tonic = float(np.loadtxt(fname + tonicExt))

    pitchOut = medianFilterPitchContour(timePitch[:,1], tonic= tonic, hopSize = hopSize, filtLen = filtDur)

  timePitch[:,1] = pitchOut

  np.savetxt(fname + outExt, timePitch, delimiter = "\t")
  
  
###########################1) Low-pass filtering pitch contour ###########################
#Motivation
#1) Remove spurious jumps in pitch due to octave errors or other errors
#2) Supress the extent of ornaments like Kan swar which arise challenges to segmentation process or in identification of stable note regions
#3) In general to smooth out pitch content supressing fine movements which might not be relevant for melodic similarity.

def gaussianFilterPitchContour(pitch, tonic=-1, hopSize=-1, filtLen=0.1, sigma=0.05):
  """
  This function low-pass filters the pitch contour using gaussian shape filter to remove spurious jumps and to supress some of the transient like ornamentations which pose challenge for melody processing
  
  Input: 
    pitch: numpy array of the pitch values
    hopSize: hop size of the pitch sequence (seconds)
    filtLen: length of the lowpass filter to be used (seconds)
    sigma: sigma of the gaussian shape to be used
    
  Output:
    pitchOut: low pass filtered pitch sequence
  """
  
  if type(pitch)==str:#passed a file name
    timePitch = np.loadtxt(pitch)
    pitch = timePitch[:,1]
    hopSize = timePitch[1,0]-timePitch[0,0]

  if type(tonic)==str:
    tonic = float(np.loadtxt(tonic))
  
  indSil = np.where(pitch<=0)[0]
  N = len(pitch)
  
  filtSize = np.round((float(filtLen)/float(hopSize))).astype(np.int)
  filtSize = filtSize + 1 - np.mod(filtSize,2)  # to make the filter length odd number of samples
  sigmaSamples = np.round((float(sigma)/float(hopSize))).astype(np.int)
  f = sig.gaussian(filtSize, sigmaSamples)
  f = f/np.sum(f)
  
  centConverted = 0
  if type(tonic)!=int:
    pitchCents  = 1200*np.log2((pitch+eps)/tonic)
    centConverted = 1
  else:
    pitchCents = pitch
  
  pitchOut = sig.filtfilt(f, [1], pitchCents)

  #because of the low pass filtering many points close to silence regions are pulled down. Lets restore these points by substituting them with the original pitch sequence
  silSegs = seg.groupIndices(indSil)
  # if silSegs[0,0]==0:#starting with zero
  #   silSegs = np.delete(silSegs,0,axis=0)
  # if silSegs[-1,1]>=N-1:
  #   silSegs = np.delete(silSegs,silSegs.shape[0]-1,axis=0)
  for silSeg in silSegs:
    if silSeg[0]==0:#starting with zero
      indRight = np.arange(silSeg[1]+1, silSeg[1]+(filtSize))
      pitchOut[indRight] = pitchCents[indRight]
    elif silSeg[1]>=N-1:
      indLeft = np.arange(silSeg[0]-(filtSize-1),silSeg[0])
      pitchOut[indLeft] = pitchCents[indLeft]  
    else:
      indLeft = np.arange(silSeg[0]-(filtSize-1),silSeg[0])
      indRight = np.arange(silSeg[1]+1, silSeg[1]+(filtSize))
      pitchOut[indLeft] = pitchCents[indLeft]
      pitchOut[indRight] = pitchCents[indRight]
  
  #if centConverted==1:
  #  indHighdiff = np.where(abs(pitchOut-pitchCents)>600)[0]
  #  pitchOut[indHighdiff] = pitchCents[indHighdiff]
  #else:
  #  indHighdiff = np.where(abs(pitchOut-pitchCents)>25)[0]
  #  pitchOut[indHighdiff] = pitchCents[indHighdiff]

  if type(tonic)!=int:
    pitchOut = tonic*np.power(2,pitchOut/1200)
  
  pitchOut[indSil] = 0
  
  return pitchOut

def batchProcessGaussianFilterPitchContour(root_dir, searchExt = '.wav', pitchExt= '.tpe', tonicExt = '.tonic', outExt = '.tpeMedFilt', filtDur = .05, gaussSigma = 0.03):
  """
  Simple wrapper for batch processing of gaussian filtering
  """
  filenames = BP.GetFileNamesInDir(root_dir, searchExt)
  for filename in filenames:
    fname, ext = os.path.splitext(filename)
    timePitch = np.loadtxt(fname+pitchExt)
    hopSize = timePitch[1,0]-timePitch[0,0]
    tonic = float(np.loadtxt(fname + tonicExt))

    pitchOut = gaussianFilterPitchContour(timePitch[:,1], tonic= tonic, hopSize = hopSize, filtLen = filtDur, sigma = gaussSigma)

  timePitch[:,1] = pitchOut

  np.savetxt(fname + outExt, timePitch, delimiter = "\t")

def removeSpuriousPitchJumps(pitch, tonic=-1, hopSize=-1, filtLen=0.3):
  """
  This function removes spurious pitch jumps. Median filtering can't be applied over long window length becasue that would be too much of smoothening. But definitely 
  spurious pitch jumps which span humdreds of ms of time duration can be corrected following an heuristic based non linear filtering. This function
  attemps to perform such a filtering. This should be applied before any pitch interpolation.
  """
  #some params (fixed ones)
  octBins = 1200.0
  halfOctBins = 600.0
  semitoneBins = 100.0
  silCentVal = -5000.0

  #There are 2 type of cases
  #1) Pitch jump in the middle of a breath phrase, typically lasting over 100-300 ms and then comes back to a normal pitch track
  #2) Pitch jump at the start or at the end of a breath phrase. Such jumps are preceeded or succeeded by silences.
  
  if type(pitch)==str:#passed a file name
    timePitch = np.loadtxt(pitch)
    pitch = timePitch[:,1]
    hopSize = timePitch[1,0]-timePitch[0,0]

  if type(tonic)==str:
    tonic = float(np.loadtxt(tonic))
  
  indSil = np.where(pitch<=0)[0]
  
  filtSize = np.round((float(filtLen)/float(hopSize))).astype(np.int)
  filtSize = filtSize + 1 - np.mod(filtSize,2)  # to make the filter length odd number of samples
  
  centConverted = 0
  if type(tonic)!=int:
    pitchCents  = octBins*np.log2((pitch+eps)/tonic)
    centConverted = 1
  else:
    pitchCents = pitch
  
  # lets say there are 4 types of jumps
  #1) silence->valid pitch
  #2) valid pitch -> silence
  #3) pitch->pitch octave
  #4) pitch octave -> pitch
  #5) pitch -> pitch other jumps
  #6) pitch other jumps -> pitch
  #7) one point jump
  # Annotating at each instance what kind of jump is happening
  jumpType = np.zeros(len(pitchCents))
  diff = pitchCents[1:]-pitchCents[:-1]
  indChange = np.where(abs(diff)>=halfOctBins)[0] +1 
  indChange = np.concatenate((indChange-1, indChange, indChange+1))
  indChange = np.sort(np.unique(indChange))

  UpInd = [1,3,5]
  DnInd = [2,4,6]

  for ii in indChange:
    if ii>=pitchCents.size-2:
      continue
    diff = pitchCents[ii+1]-pitchCents[ii]
    if np.min(abs(diff-octBins*np.array([-3,-2,-1,1,2,3,4])))<=3*semitoneBins and pitchCents[ii]>silCentVal: #octave jump, -5000 is regarded as best case silence cent value, in reality it will be further less
      # there has been an octave jump and that the current sample is not silence
      if diff>0:
        jumpType[ii+1]=3
      elif diff<0:
        if jumpType[ii]==0:
          jumpType[ii]=4
        else:
          jumpType[ii]=7 # one point case

    elif abs(diff)>=halfOctBins:
      
      if pitchCents[ii]<silCentVal:
        jumpType[ii+1]=1
      elif pitchCents[ii+1]<silCentVal:
        if jumpType[ii]==0:
          jumpType[ii]=2
        else:
          jumpType[ii]=7
      elif diff>0:
        jumpType[ii+1]=5
      else:
        if jumpType[ii]==0:
          jumpType[ii]=6
        else:
          jumpType[ii]=7

  for ii in range(len(pitchCents)-filtSize):
    indChange = np.where(jumpType[ii:ii+filtSize]>0)[0]
    indChange = np.sort(indChange)
    if len(indChange)==0:
      continue
    if len(indChange)==2: # if there are exactly two jumps
      i1 = indChange[0]+ii
      i2 = indChange[1]+ii
      # if both these pairs are to and fro silence to pitch, these shouldn't be considered. So
      if jumpType[i1] + jumpType[i2] == 3: # this ensures to and fro silence one has to be 1 and another 2
        jumpType[i1] =0
        continue
      # only possibility that is allowed in the pitch error is that starting change goes from low to high and second change from high to low
      if UpInd.count(jumpType[i1]) >0 and DnInd.count(jumpType[i2]) >0:
        if jumpType[i1]==3 and jumpType[i2]!=6:
          shift = octBins*np.round((pitchCents[i1] - pitchCents[i1-1])/octBins)
          pitchCents[i1:i2+1] = pitchCents[i1:i2+1]-shift
        elif jumpType[i2]==4 and jumpType[i1]!=5:
          shift = octBins*np.round((pitchCents[i2] - pitchCents[i2+1])/octBins)
          pitchCents[i1:i2+1] = pitchCents[i1:i2+1]-shift
        elif jumpType[i1]==1:
          shift = semitoneBins*np.round((pitchCents[i2] - pitchCents[i2+1])/semitoneBins)
          pitchCents[i1:i2+1] = pitchCents[i1:i2+1]-shift
        elif jumpType[i2]==2:
          shift = semitoneBins*np.round((pitchCents[i1] - pitchCents[i1-1])/semitoneBins)
          pitchCents[i1:i2+1] = pitchCents[i1:i2+1]-shift
        else:
          if (i2-i1)*hopSize>0.5:
            err1 = np.mod(abs(pitchCents[i1]-pitchCents[i1-1]),octBins)
            err2 = np.mod(abs(pitchCents[i2+1]-pitchCents[i2]),octBins)
            if err1 < err2:
              shift = octBins*np.round((pitchCents[i1] - pitchCents[i1-1])/octBins)
            else:
              shift = octBins*np.round((pitchCents[i2] - pitchCents[i2+1])/octBins)
            pitchCents[i1:i2+1] = pitchCents[i1:i2+1]-shift
          else:
            pitchCents[i1-1:i2+2] = np.linspace(pitchCents[i1-1], pitchCents[i2+2], i2-i1+3)
        jumpType[i1]=0
        jumpType[i2]=0
      else:
        jumpType[i1]=0
      #### The only exception to the above rule of going up and down is that if there is double octave error.
      ###if jumpType[i1]==3 and ( jumpType[i2]==3 or jumpType[i2]==5):
        ###shift = octBins*np.round((pitchCents[i1] - pitchCents[i1-1])/octBins)
        ###pitchCents[i1:i2] = pitchCents[i1:i2]-shift
        ###jumpType[i1]=0
        ####jumpType[i2]=0#this is not resetted purposefully
    
    # special case of only one change point    
    if jumpType[ii+indChange[0]]==7 and len(indChange)==1:
      i1 = indChange[0]+ii
      #This is a tricky point, this situation can arise between two valid pitch samples or near the silence.
      if pitchCents[i1+1] <silCentVal:
        pitchCents[i1] = pitchCents[i1-1]
      elif pitchCents[i1-1] <silCentVal:
        pitchCents[i1] = pitchCents[i1+1]        
      else:
        pitchCents[i1] = np.mean((pitchCents[i1+1],pitchCents[i1-1]))
      jumpType[i1]=0

  if type(tonic)!=int:
    pitchOut = tonic*np.power(2,pitchCents/1200)
  
  pitchOut[indSil] = 0

  return pitchOut

def batchProcessRemoveSpuriousPitchJumps(root_dir, searchExt = '.wav', pitchExt= '.tpe', tonicExt = '.tonic', outExt = '.tpeOctCorr', filtDur = .3):
  """
  Simple wrapper for batch processing of removing spurious pitch errors
  """
  filenames = BP.GetFileNamesInDir(root_dir, searchExt)
  for filename in filenames:
    fname, ext = os.path.splitext(filename)
    timePitch = np.loadtxt(fname+pitchExt)
    hopSize = timePitch[1,0]-timePitch[0,0]
    tonic = float(np.loadtxt(fname + tonicExt))

    pitchOut = removeSpuriousPitchJumps(timePitch[:,1], tonic= tonic, hopSize = hopSize, filtLen = filtDur)

  timePitch[:,1] = pitchOut

  np.savetxt(fname + outExt, timePitch, delimiter = "\t")


def InterpolateSilence(array, silence_val, hopSize, maxSilDurIntp=0.25, interpAllSil = False):
  """

  """
    
  if type(array) ==list:
      array = np.array(array)    
  
  # interpolating zeros in middle of the time series
  array = array.astype('float')
  sil_ind = np.where(array<=silence_val)[0]    
  last_sil_ind=sil_ind[0]
  sil_ind = np.append(sil_ind,0)
  length= len(array)
  for ii in range(0,len(sil_ind)-1):
      if sil_ind[ii] +1 != sil_ind[ii+1]:
          if last_sil_ind!=0 and sil_ind[ii]!=length-1:
              inter_data = np.linspace(array[last_sil_ind-1],array[sil_ind[ii]+1] , sil_ind[ii]-last_sil_ind+3)
              if (sil_ind[ii]-last_sil_ind+1)*hopSize < maxSilDurIntp or interpAllSil:
                array[last_sil_ind-1:sil_ind[ii]+2] = inter_data
              last_sil_ind=sil_ind[ii+1]
          else:
              last_sil_ind=sil_ind[ii+1]
  

  ## zeros at the beginning and at the end of the time series are replaced by the mean of the time series
  #if interpAllSil:
  #  sil_ind = np.where(array<=silence_val)[0]
  #  nonSil_ind = np.where(array>silence_val)[0]
  #  array[sil_ind] = np.mean(array[nonSil_ind])
  
  #filling silence at the starting and at the end by random values
  if interpAllSil:
    sil_ind = np.where(array<=silence_val)[0]
    nonSil_ind = np.where(array>silence_val)[0]
    array[sil_ind] = np.random.random_integers(np.min(array[nonSil_ind]), np.max(array[nonSil_ind]), len(sil_ind))
  return array
        
def BatchProcessInterpPitchSilence(RootDir, FileExt2Proc = ".tpe", NewExt = "", PitchCol=2, SilVal=0, maxSilDurIntp=0.25):
    
  audiofilenames = GetFileNamesInDir(RootDir, FileExt2Proc)
  
  if len(NewExt)==0:
      NewExt = FileExt2Proc + "Intrp"
  
  for audiofile in audiofilenames:
        time_pitch = np.loadtxt(open(audiofile,"r"))
        hopSize = time_pitch[1,0]-time_pitch[0,0]
        new_pitch = InterpolateSilence(time_pitch[:,PitchCol-1],SilVal,hopSize, maxSilDurIntp)
        time_pitch[:,PitchCol-1]=new_pitch
        file,ext = os.path.splitext(audiofile)
        np.savetxt(file + NewExt, time_pitch, delimiter = "\t")



def postProcessPitchSequence(pitch, tonic=-1, hopSize=-1, filtDurMed=0.05, filtDurGaus=0.05, winDurOctCorr=0.3, sigmaGauss=0.025, fillSilDur= 0.25, interpAllSil=False, upSampleFactor=1):
  """
  This is a function which performs all the selected post processing steps on a pitch sequence. 
  Inputs: 
    Pitch (string or ndarray): filename or numpy array of the pitch sequence to process
    tonic (string or float): file name or float value of the tonic to be used
    hopSize (float): hopSize of the pitch sequence, needed when pitch is inputted as nd.array
    filtDurMed (float): filter duration to be used for median filtering in ms. [Negative means no median filtering]
    filtDurGaus (float): filtering duration for the gaussian smoothening to be applied in ms. [Negative means no median filtering]
    winDurOctCorr (float): window length to be used for the octave/other error correction (spurious pitch jumps). [Negative means no correction to be applied]
    sigmaGauss (float): sigma (Std) of the gaussian to be used for the smoothening.
    fillSilDur (float): maximum silence duration that has to be interpolated [Negative means no interpolation]
    interpAllSil (bool): boolean to indicate whether to interpolate all silence regions or not, (needed many times, specially when generating random noise candidates)
  Output:
    pitchOut (ndarray): output numpy array of the pitch sequence after post processing

  Processing: 
    pitch -> pitch octave correction -> low dur pitch interp -> median filt -> gaussian filt -> interp all silences -> pitch out  
  """

  #some params (fixed ones)
  octBins = 1200.0
  halfOctBins = 600.0
  semitoneBins = 100.0
  silCentVal = -5000
  silVal = 0

  ### Some important checks to do 
  if type(pitch)==str:#passed a file name
    timePitch = np.loadtxt(pitch)
    pitch = timePitch[:,1]
    hopSize = timePitch[1,0]-timePitch[0,0]

  if type(tonic)==str:
    tonic = float(np.loadtxt(tonic))
  else:
    tonic = float(tonic)
  
  if tonic >400 or tonic < 80:
    print "You should provide a valid tonic value for this processing"
  
  if upSampleFactor > 0 and upSampleFactor !=1:
    pitchResampled = resamplePitchSequence(pitch,  upSampleFactor, silVal, tonic = tonic, hopSize=hopSize)
  else:
      pitchResampled = pitch
  
  if winDurOctCorr > 0:
    pitchOctCorr = removeSpuriousPitchJumps(pitchResampled, tonic=tonic, hopSize=hopSize, filtLen=winDurOctCorr)
  else:
    pitchOctCorr = pitchResampled

  if fillSilDur > 0:
    pitchSinSilShort = InterpolateSilence(pitchOctCorr, silVal, hopSize, maxSilDurIntp=fillSilDur, interpAllSil = False)
  else:
    pitchSinSilShort = pitchOctCorr

  if filtDurMed > 0:
    pitchMedFilt = medianFilterPitchContour(pitchSinSilShort, tonic=tonic, hopSize=hopSize, filtLen = filtDurMed)
  else:
    pitchMedFilt = pitchSinSilShort

  if filtDurGaus > 0:
    pitchGausFilt = gaussianFilterPitchContour(pitchMedFilt, tonic=tonic, hopSize=hopSize, filtLen=filtDurGaus, sigma=sigmaGauss)
  else:
    pitchGausFilt = pitchMedFilt

  if interpAllSil:
    pitchIntrpAll = InterpolateSilence(pitchGausFilt, silVal, hopSize, maxSilDurIntp=fillSilDur, interpAllSil = interpAllSil)
  else:
    pitchIntrpAll = pitchGausFilt

  return pitchIntrpAll

def batchProcessPitchPostProcess(root_dir, searchExt = '.wav', pitchExt= '.tpe', tonicExt = '.tonic', outExt = '.tpeOctCorr', filtDurMed=0.05, filtDurGaus=0.05, winDurOctCorr=0.3, sigmaGauss=0.025, fillSilDur= 0.25, interpAllSil=False, upSampleFactor = 1, over_write = 1):
  """
  Wrapper to batch process pitch post processing
  """
  filenames = BP.GetFileNamesInDir(root_dir, searchExt)
  for filename in filenames:
    print "Processing file %s"%filename
    fname, ext = os.path.splitext(filename)
    try:
        timePitch = np.loadtxt(fname+pitchExt)
    except:
        continue
    hopSize = timePitch[1,0]-timePitch[0,0]
    try:
        tonic = float(np.loadtxt(fname + tonicExt))
    except:
        continue
    if os.path.isfile(fname+outExt) and over_write == 0:
        continue
    pitchOut = postProcessPitchSequence(timePitch[:,1], tonic= tonic, hopSize = hopSize, filtDurMed=filtDurMed, filtDurGaus=filtDurGaus, winDurOctCorr=winDurOctCorr, sigmaGauss=sigmaGauss, fillSilDur= fillSilDur, interpAllSil=interpAllSil, upSampleFactor = upSampleFactor)
    TStamps = float(hopSize)*float(upSampleFactor)*np.arange(pitchOut.size)
    timePitch = np.array([TStamps, pitchOut]).transpose()
    np.savetxt(fname + outExt, timePitch, delimiter = "\t")


