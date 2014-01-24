import numpy as np
import sys, os
import yaml
import copy

sys.path.append(os.path.join(os.path.dirname(__file__), '../melodyProcessing'))
sys.path.append(os.path.join(os.path.dirname(__file__), '../batchProcessing'))
import segmentation as MS
import basicOperations as BPO
import batchProcessing as BP


eps = np.finfo(np.float).eps


class motifDiscovertIO():
    
    params= {}
    
    def __init__(self, params):
        print 'Init motif generation step'
        self.params = params
        
    def generateSubsequenceDataset(self, pitchFile, tonicFile, outputCandidateFile):
        
        #reading pitch data
        tonic = np.loadtxt(open(tonicFile,"r"))
        pitchData,timeData,pHop = BPO.readPitchFile(pitchFile)
        pCents = BPO.PitchHz2Cents(pitchData, tonic)
        
        #preprocessing pitch data
        factor=3
        pCents, pHop, timeData = BPO.downsamplesPitchData(pCents,pHop,timeData, factor)
        
        
        if self.params.keys()[0]=='slidingWindow':
            segments = self.slidingWindowCandidates(pCents, pHop, self.params['slidingWindow']['windowLength'],pHop*3)
        dataset = [] 
        [dataset.append(pCents[segment[0]:segment[1]]) for segment in segments]
        np.savetxt(outputCandidateFile,np.array(dataset),fmt='%.2f')
        
        
    def slidingWindowCandidates(self, pCents, pHop, windowLength, windowHop):
        
        msObj = MS.melodySegmentation()
        bPhrases = msObj.ExtractBreathPhrases(pCents, pHop, 0.05)
        windowSamples = np.round(windowLength/pHop)
        hopSamples = np.round(windowHop/pHop)
        
        segments = []
        for bPhrase in bPhrases:
            if bPhrase[1]-bPhrase[0]>=windowSamples:
                
                for ii in np.arange(bPhrase[0], bPhrase[1]-windowSamples, hopSamples):
                    segments.append([ii,ii+windowSamples])
                    
        
        return segments
    
    def isflat(self, a):
        
        if abs(min(a)-max(a))<300 or min(a)<-1200:
            return True
        else:
            return False
            
    def ComputeLocalVariance(self, pCents, pHop):
        var_len = 100.0/1000.0    #in ms
        var_samples = int(round(var_len/pHop))
        
        pVar = np.zeros(pCents.shape[0])
        for i in range(0+var_samples,pCents.shape[0]-var_samples):
            pVar[i] = min(np.var(pCents[i:i+var_samples]),np.var(pCents[i-var_samples:i]))
            
        return pVar
    
    def nonFlatIndexes(self, pCents, pHop, binsPOctave):
        var_len = 100.0/1000.0    #in ms
        var_samples = int(round(var_len/pHop))
        threhsold = 1.5*1200.0/binsPOctave
        
        running_mean = np.zeros(pCents.shape[0])
        running_var = np.zeros(pCents.shape[0])
        running_mean[var_samples] = np.sum(pCents[0:2*var_samples+1])
        for ii in np.arange(var_samples+1,pCents.size-var_samples-1):
            running_mean[ii] = running_mean[ii-1]+ pCents[ii+var_samples]-pCents[ii-var_samples-1]
        
        running_mean= running_mean/(1+(var_samples*2))
        
        pCents_sq = np.square(pCents-running_mean)
        
        running_var[var_samples] = np.sum(pCents_sq[0:2*var_samples+1])
        for ii in np.arange(var_samples+1,pCents.size-var_samples-1):
            running_var[ii] = running_var[ii-1]+ pCents_sq[ii+var_samples]-pCents_sq[ii-var_samples-1]
        
        running_var = np.sqrt(running_var)
        
        flatRegion = np.zeros(pCents.shape[0])
        
        ind_Nonflat = np.where(running_var>threhsold)[0]
        
        return ind_Nonflat
        
    
    def generateSubsequenceDataset(self, root_dir, output_dir, pitchExt, tonicExt, downsampleFactor, windowLength, combineData = 0, writeBinary = 1, meanNormalize = 0, flatnessThreshold=0.8, binsPOctave=120, fixPointData=1):
        
        
        if combineData:
            timeInfo=np.array([])
            fileInfo={}
        
        #obtaining all the files in the directory
        filenames = BP.GetFileNamesInDir(root_dir,pitchExt)
        
        # iterating over each file
        for kk, filename in enumerate(filenames):
            
            #separate file name from extension
            fname, ext = os.path.splitext(filename)
            
            audiofile = fname.split('/')[-1]
            
            #if data is not to be combined, then output goes into individual directories (dirname = filename)
            if not combineData:
                out_dir = output_dir + '/' + audiofile
                #creating directory if doesn't exist
                if not os.path.isdir(out_dir):
                    os.makedirs(out_dir)
            
            #reading pitch and tonic data
            pitchData,timeStamps,pHop = BPO.readPitchFile(fname+pitchExt)
            tonic = np.loadtxt(open(fname+tonicExt,"r"))            
            
            #Convert the pitch values into cents, Note that we add a offset of 1 octave to make everything positive, assuming that pitch wont go below one octave to tonic
            pCents=np.round(binsPOctave*np.log2((eps+pitchData)/tonic)).astype(np.int) + binsPOctave 
            
            
            #downsampling pitch data
            factor=downsampleFactor
            pCents, pHop, timeStamps = BPO.downsamplesPitchData(pCents,pHop,timeStamps, factor)
            
            
            #removing silence regions from the pitch sequence
            ind_silence = np.where(pCents<0)[0] ###Please correct this silence condition once log eps is used
            pCents = np.delete(pCents,ind_silence)
            timeStamps = np.delete(timeStamps,ind_silence)
            
            #computing all the indices which have non flat region of pitch (binsPOctave is an input to decide threshold, thats it!!)
            nonFlatIndexes = self.nonFlatIndexes(pCents, pHop,binsPOctave)
            
            #making an array which will be 1 for nonflat and 0 for flat regions
            flatNonflat = np.zeros(pCents.shape[0])
            flatNonflat[nonFlatIndexes]=1
            
            #given the windowLength or motigLength compute how many samples do they correspond to
            windowSamples = int(np.round(windowLength/pHop))
            
            #to create a matlab buffer like thing, we do the following
            row = np.array([np.arange(windowSamples)])
            col = np.array([np.arange(pCents.size-windowSamples)])
            col = np.transpose(col)
            col2 = copy.deepcopy(col)
            col2[:]=1
            ind = row*col2 + col
            
            
            #so after creating matrix where each row is a subsequence, look which subsequence to reject based on flatness threshold
            mtx = flatNonflat[ind]            
            mean_array = np.mean(mtx,axis=1)            
            ind_Invalid = np.where(mean_array<flatnessThreshold)[0]
            
            ind = np.delete(ind,[ind_Invalid],axis=0)
            
            #finally obtain pitch matrix and timestamp array
            mtx = pCents[ind]
            timeStamps = timeStamps[ind[:,0]]
            
            ###TODO impement mean normalizatoin
            
            ### Adding a crucial check!!
            if timeStamps.shape[0] != mtx.shape[0]:
                print filename
            
            if combineData:#keep on appending the data to combine it
                
                if kk==0:
                    pitch = copy.deepcopy(mtx)
                    timeInfo = copy.deepcopy(timeStamps)
                else:
                    pitch = np.append(pitch, mtx,axis=0)
                    timeInfo = np.append(timeInfo, timeStamps)
                fileInfo[filename]= [timeInfo.size-timeStamps.size, timeInfo.size]
            
            
            else:#just dump for each file
                
                if writeBinary:
                    if fixPointData:
                        mtx.astype(np.uint32).tofile(out_dir+'/'+ audiofile +'.pitchSubDBbin')
                    else:
                        mtx.astype(np.float).tofile(out_dir+'/'+ audiofile +'.pitchSubDBbin')
                        
                    timeStamps.astype(np.float).tofile(out_dir+'/'+ audiofile +'.timeSubDBbin')
                else:
                    np.savetxt(out_dir+'/'+ audiofile +'.pitchSubDBtxt', mtx , fmt='%d')
                    np.savetxt(out_dir+'/'+ audiofile +'.timeSubDBtxt', timeStamps, fmt='%.3f')
                
            if combineData:
                #another crucial check at each step
                if pitch.shape[0] != timeInfo.shape[0]:
                    print filename
            
        if combineData:
            
            if writeBinary:
                if fixPointData:
                        pitch.astype(np.uint32).tofile(output_dir+'/'+'AggPitch.bin')
                else:
                        pitch.astype(np.float).tofile(output_dir+'/'+'AggPitch.bin')
                timeInfo.astype(np.float).tofile(output_dir+'/'+'AggTime.bin')
            else:
                np.savetxt(output_dir+'/'+'AggPitch.txt', pitch , fmt='%d')
                np.savetxt(output_dir+'/'+'AggTime.txt', timeInfo, fmt='%.3f')
            
            stream = file(output_dir+'/'+'fileInfo.yaml','w')
            yaml.dump(fileInfo, stream)
            
            print "Total number of time series : " + str(pitch.shape[0])
            print "Length of single Sub sequence : " + str(pitch.shape[1])
        
            

        
    
    def generateLinearDataset(self, root_dir, output_dir, pitchExt, tonicExt, downsampleFactor, min_nyas_dur=-1):
        
        pitch=np.array([])
        timeInfo=np.array([])
        fileInfo={}
        
        filenames = BP.GetFileNamesInDir(root_dir,pitchExt)
        
        for filename in filenames:
            fname, ext = os.path.splitext(filename)
            #reading pitch and tonic data
            pitchData,timeStamps,pHop = BPO.readPitchFile(fname+pitchExt)
            tonic = np.loadtxt(open(fname+tonicExt,"r"))
            pCents = BPO.PitchHz2Cents(pitchData, tonic)
            
            
            #some preprocessing
            
            #removing flat regions
            if (min_nyas_dur>0):
                msObj = MS.nyasSegmentation()
                msObj.ComputeNyasCandidates(pitchData, tonic.tolist(), pHop)
                msObj.FilterNyasCandidates(min_nyas_duration=min_nyas_dur)
            
                for swar in msObj.nyasInfo.keys():
                    for seg in msObj.nyasInfo[swar]:
                        pCents[seg[0]:seg[1]]=-5000
            
            
            #downsampling
            factor=downsampleFactor
            pCents, pHop, timeStamps = BPO.downsamplesPitchData(pCents,pHop,timeStamps, factor)
            
            
            
            #removing silence regions
            ind_silence = np.where(pCents<-4000)[0] ###Please correct this silence condition once log eps is used
            pCents = np.delete(pCents,ind_silence)
            timeStamps = np.delete(timeStamps,ind_silence)
            
            #accumulating
            pitch = np.append(pitch, pCents)
            timeInfo = np.append(timeInfo, timeStamps)
            fileInfo[filename]= [timeInfo.size-timeStamps.size, timeInfo.size]
            
        np.savetxt(output_dir+'/'+'AggPitch.txt', pitch, fmt='%.2f')
        np.savetxt(output_dir+'/'+'AggTime.txt', timeInfo, fmt='%.2f')
        stream = file(output_dir+'/'+'fileInfo.yaml','w')
        yaml.dump(fileInfo, stream)
        