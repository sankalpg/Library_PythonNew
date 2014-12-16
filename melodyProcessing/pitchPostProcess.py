import numpy as np
import scipy.signal as sig



###########################1) Median filtering pitch contour ###########################
#Motivation
#1) Remove spurious jumps in pitch due to octave errors or other errors
#2) Supress the extent of ornaments like Kan swar which arise challenges to segmentation process or in identification of stable note regions
#3) In general to smooth out pitch content supressing fine movements which might not be relevant for melodic similarity.

def medianFilterPitchContour(pitch, hopSize, filtLen):
  """
  This function median filters the pitch contour to remove spurious jumps and to supress some of the transient like ornamentations which pose challenge for melody processing
  
  Input: 
    pitch: numpy array of the pitch values
    hopSize: hop size of the pitch sequence (seconds)
    filtLen: length of the median filter to be used (seconds)
    
  Output:
    pitchOut: median filtered pitch sequence
  """
  
  filtSize = np.round((float(filtLen)/float(hopSize))).astype(np.int)
  filtSize = filtSize + 1 - np.mod(filtSize,2)  # to make the filter length odd number of samples
  
  pitchOut = sig.medfilt(pitch, filtSize)
  
  return pitchOut
  
  
  
###########################1) Low-pass filtering pitch contour ###########################
#Motivation
#1) Remove spurious jumps in pitch due to octave errors or other errors
#2) Supress the extent of ornaments like Kan swar which arise challenges to segmentation process or in identification of stable note regions
#3) In general to smooth out pitch content supressing fine movements which might not be relevant for melodic similarity.

def gaussianFilterPitchContour(pitch, hopSize, filtLen, sigma):
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
  
  filtSize = np.round((float(filtLen)/float(hopSize))).astype(np.int)
  filtSize = filtSize + 1 - np.mod(filtSize,2)  # to make the filter length odd number of samples
  sigmaSamples = np.round((float(sigma)/float(hopSize))).astype(np.int)
  f = sig.gaussian(filtSize, sigmaSamples)
  f = f/np.sum(f)
  pitchOut = sig.filtfilt(f, [1], pitch)
  
  return pitchOut