#! /usr/bin/python
import numpy as np


eps = np.finfo(np.float).eps


def readPitchFile(pitchfile):
    """
    Function to read the pitch file to obtain time stamps, pitch sequence, and hopesize.
    Format of the pitch file:
    <time stamps> <pitch sequence>
    returns (time,pitch,hopsize)
    """
    timepitch = np.loadtxt(pitchfile)
    timeData = timepitch[:,0]
    pitchData = timepitch[:,1]    
    phop = timeData[1]-timeData[0]
    return pitchData,timeData,phop
    
    
def PitchHz2Cents(pitch, tonic, binPOctave=1200):
    """
    Function to convert the pitch values from Hz to cent scale using provided tonic value
    """
    #print "tonic used for normalization is "+ str(tonic) + " cents"
    #ind_zero_pitch = np.where(pitch<=0)[0]  ###TODO remove this line
    pCents=binPOctave*np.log2((eps+pitch)/tonic)
    #pCents[ind_zero_pitch]= -5000   ###TODO remove this line, this is just to make it same as original version, but its not needed
    return pCents
    
    
def downsamplesPitchData(pCents,pHop,timeData, factor):
    
    if not np.round(factor)==factor:
        print 'Please provide a integral factor for downsampling, later we might support fractional'
        return 0,0,0
    
    pHop = pHop*factor
    ind_new = np.arange(0,timeData.size,factor)
    timeData = timeData[ind_new]
    pCents = pCents[ind_new]
    
    return pCents, pHop, timeData


