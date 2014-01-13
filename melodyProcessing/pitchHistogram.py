#! /usr/bin/python
import numpy as np
import os,sys, copy
import matplotlib.pyplot as plt
import math
import scipy.stats as sps
import scipy as sp


sys.path.append(os.path.join(os.path.dirname(__file__), '../basicFunctions/'))

import basicOperations as BOP
import DSPFunctions as DF

class PitchHistogram():
    
    hRange = -1
    hResolution = -1
    nBins = -1
    
    hist_Yval = []
    hist_Xval = []
    
    pitch = -1
    timeStamps = -1
    tonic = -1
    pCents = -1
    swarLocs = []
    
    def __init__(self, pitch, tonic, hRange = [-1200, 3600], hResolution = 5, timeStamps = -1):
        """
        Initialization of parameters
        pitch = pitch sequence of the audio
        tonic = tonic of the lead artist
        hRange = histogram range (x axis) in cents
        hresolution = resolution of the histogram (x axis) in cents. i.e. how many cents for 1 bin
        """
        # if parameters are provided in initialization just update them all
        if self.hRange == -1:
            self.hRange = hRange
        
        if self.hResolution == -1:
            self.hResolution = hResolution
            
        if self.nBins == -1:
            self.nBins = int((hRange[1]-hRange[0])/hResolution)
            
        if type(timeStamps) != int:
            self.timeStamps = timeStamps
        
        try:
            if len(pitch)>0:
                self.pitch = pitch
            else:
                print "Specify valid pitch data"
        except:
            print "Provide valid pitch data as an array of pitch sequence"
        
        try:
            if type(tonic) == int or type(tonic) == float:
                self.tonic = tonic
            else:
                print "Specify valid tonic value"
        except:
            print "Provide valid tonic value"
        
        
    def setPitch(self, pitch):
        self.pitch = pitch

    def setTonic(self, tonic):
        self.tonic = tonic
    def setTimeStamps(self,timeStamps):
        self.timeStamps = timeStamps
        
    def GetPitchHistogram(self):
        return self.hist_Yval, self.hist_Xval
    
    
    def SmoothPitchHistogram(self, Histogram=[-1], Variance=15):
        """This function performs low pass filtering (using a window of the shape normal distribution) of the histogram
        Input parameters:
        Variance = variance of the normal distribution used for smoothening of the histogram (default ==15) IN CENTS!!!
        """
        
        variance = float(Variance)/self.hResolution  #variance of the normal distribution, in bins
        wind_range = (np.array([-50, 50])/self.hResolution) # from index corresponding to -50 to 50 cents
        norm_win = sps.norm(0, variance).pdf(np.linspace(wind_range[0],wind_range[1],num=1+wind_range[1]-wind_range[0]))
        norm_win = norm_win/sum(norm_win)
        
        #convolving histogram withh  this window
        if (len(Histogram)==1):
            hist_smooth = sp.convolve(self.hist_Yval,norm_win, 'same')
            
        else:
            hist_smooth = sp.convolve(Histogram,norm_win, 'same')
            
        #normalizing both histograms
        hist_smooth = hist_smooth/max(hist_smooth)
        
        if (len(Histogram)==1):
            self.hist_Yval = hist_smooth
        else:
            return hist_smooth
    
    def PlotHistogram(self):
        """This function plots the computed histogram
        """    
        if ((len(self.hist_Yval)>0)&(len(self.hist_Xval)>0)):
            plt.plot(self.hist_Xval, self.hist_Yval)
            plt.show()
        
    def ComputePitchHistogram(self, pitch=-1, timeStamps=-1, tonic=-1, tRange=-1, Oct_fold =0):
        """This function computes pitch histogram 
        Input parameters:
        pitch = pitch sequence 
        tonic = tonic value of the lead artist
        tRange = time range within which pitch has to be considered for constructing pitch histogram
        timeStamps = time stamps needed if tRange is specified for histogram construction
        Oct_fold = (0 or 1); 0 for no octave folding, 1 for octave folding of the pitch histogram
        """
    
        ### reading values and throwing errors if any
        if type(pitch) != int:
            self.setPitch(pitch)        
        elif type(self.pitch)==int:
            print "Please provide a pitch file name, it was not provided during initialization"
            return -1

        if tonic!=-1:
            self.setTonic(tonic)
        elif self.tonic==-1:
            print "Please provide tonic information, it was not provided during initialization"
            return -1
    
        if type(timeStamps) != int:
            self.setTimeStamps(timeStamps)
        
        if (tRange!=-1) and (type(self.timeStamps)==int):
            print "For this option of using tRange for histogram computation, timeStamps information should also be provided"
            return -1
        
    
        ### before starting anything (and after updating all important parameters), lets convert pitch to cents
        self.pCents = BOP.PitchHz2Cents(self.pitch, self.tonic)
     
        ### Copying in local buffer to process pitch in this function and octave folding if specified
        
        #sil_loc_inds = np.where(self.pCents>-1200)[0] ###TODO remove this comment this is to reproduce the same error in original version
        #pCents_local = copy.deepcopy(self.pCents[sil_loc_inds])
        
        pCents_local = copy.deepcopy(self.pCents)   ###TODO after uncommenting above two lines comment this line
        
        if (Oct_fold==1):
            pCents_local = np.mod(pCents_local,1200)
        
        ### histogram computation        
        if tRange ==-1:
        
            histogram = np.histogram(pCents_local, bins = self.nBins, range = self.hRange)
        
        else:
            str_ind=find_nearest_element_ind(self.timeStamps,tRange[0])
            end_ind=find_nearest_element_ind(self.timeStamps,tRange[1])
        
            histogram = np.histogram(pCents_local[str_ind:end_ind], bins = self.nBins, range = self.hRange)
        
        ### assigning the obtained values to the class global variables
        hist_Yval = copy.deepcopy(histogram[0])
        hist_Yval = hist_Yval.astype(float)
    
        hist_Xval = copy.deepcopy(histogram[1][1:])
        hist_Xval = hist_Xval.astype(float)
    
        ### Normalization of the histogram
        hist_Yval = hist_Yval/max(hist_Yval)

        ### if Octave folding is performed, to avoid splitting of tonic note into two parts (think, why will it happen!!) we just copy paste small end part of histogram to negative values (very intuitive)
        if Oct_fold==1: # we need to cut the negative size of Sa (tonic) i.e. from <--1200 at the nearest valley to 1200. If we dont do it at valley then there can be a shart popping hump because of low pass filtering which will be detected as a valid swar. And if you apply a lot of mind you will find that this valley has to be detected from a smoothened version of histogram otherwise jittering can cause everything go wrong
            temp_smooth = self.SmoothPitchHistogram(Histogram=hist_Yval)  
            peak_ind , valley_ind = DF.PeakValleyPicking(temp_smooth)
            valley_location = hist_Xval[max(valley_ind)]
            ind2 = np.where((hist_Xval>valley_location)&(hist_Xval<=1200))
            ind1 = np.where((hist_Xval>-(1200-valley_location))&(hist_Xval<=0))
            hist_Yval[ind1]=hist_Yval[ind2]
            hist_Yval[ind2]=0
            
        #updating class variables
        self.hist_Yval = hist_Yval
        self.hist_Xval = hist_Xval        
                
        self.SmoothPitchHistogram()
        
    
    def ValidSwarLocEstimation(self, histogram=-1, peak_valley_thsld = 0.01, Oct_fold=1):

        ### if histogram is provided process that otherwise from the self instance
        no_hist_provided=0
        if (histogram==-1):
            no_hist_provided=1
            self.ComputePitchHistogram(Oct_fold=Oct_fold)
            if (len(self.hist_Yval)>0):
                histogram = self.hist_Yval
            else:            
                print "Please initialize class instance with pitch, tonic or setPitch and setTonic appropriately"
                return -1
    
        #finding peaks and valleys
        peak_ind , valley_ind = DF.PeakValleyPicking(histogram)

        # creating an array where even values are valleys and odds are peaks, this is used for easy access and thresholding of valid peaks
        peak_valley = np.append(peak_ind, valley_ind)

        peak_valley = np.sort(peak_valley)
        peak_ind = np.sort(peak_ind)
        valley_ind = np.sort(valley_ind)

        ind_str = np.where(peak_valley==peak_ind[0])[0]
        if(ind_str==0):
            peak_valley = np.append(0,peak_valley)
        else:
            ind_str = ind_str-1

        ind_end = np.where(peak_valley==peak_ind[-1])[0]
        if ind_end == len(peak_valley)-1:
            peak_valley = np.append(peak_valley,0)
            ind_end=ind_end+1
    
        else:
            ind_end=ind_end+1    
    
        peak_valley = peak_valley[ind_str:ind_end+1]


        #iterate over all peaks and based on threshold of minimum peak to valley height decide if its valid or not.
        peak_final=np.array([])
        for i, peak in enumerate(peak_ind):
            diff = min(abs(histogram[peak_valley[2*i]]-histogram[peak_valley[(2*i)+1]]), abs(histogram[peak_valley[(2*i)+1]]-histogram[peak_valley[(2*i)+2]]))
    
            if(diff>peak_valley_thsld):
                peak_final = np.append(peak_final,peak)
        


        peak_final = peak_final.astype(int)
    
        if (no_hist_provided==1):
            self.swarLocs = peak_final
        else:
            return peak_final
            
    def ExtendSwarOctaves(self):
        sl = copy.deepcopy(self.swarCents)
        self.swarCents = np.concatenate((sl-1200, sl, sl+1200, sl + 2400))
    
    def SwarLoc2Cents(self):
        
        if ((len(self.swarLocs)>0)&(len(self.hist_Xval)>0)):
            self.swarCents = self.hist_Xval[self.swarLocs]
            
    def PlotSwarOnHistogram(self):
        if((len(self.hist_Yval)>0)&(len(self.swarLocs)>0)):
            plt.hold('True')
            plt.plot(self.hist_Xval, self.hist_Yval)
            plt.plot(self.hist_Xval[self.swarLocs], self.hist_Yval[self.swarLocs],'ro')
            plt.show()
        else:
            print "Either histogram is not computed or swar locations are not obtained!! take care!! "  
        
