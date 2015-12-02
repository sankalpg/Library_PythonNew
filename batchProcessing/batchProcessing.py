#!/usr/bin/env python


########################################################3
# Author: Sankalp Gulati (sankalp.gulati@upf.edu)
# Other Authors for individual functions are indicated inside the function!!
# Music Technology Group - Universitat Pompeu Fabra
# Created:     18/09/2012
# Last modified: --
########################################################3

try:
    import essentia.standard as ES
except:
    pass
import os, sys, copy, shutil
import numpy as np
#import eyed3
import scipy.interpolate as scpyinterp
from scipy.signal import medfilt
from mutagen import easyid3

try:
    from mutagen.mp3 import MP3
except:
    pass
sys.path.append(os.path.join(os.path.dirname(__file__), '../melodyProcessing/'))
import pitchHistogram as PH

def packTimePitch(time, pitch):
    
    if time.size != pitch.size:
        print "Please provide time and pitch arrays of the same length"
        return -1
    return np.vstack((time,pitch)).transpose()




def fineTuneTonicValue(tonicInFile, tonicOutFile, pitchFile):
  #reading pitch
  timePitch = np.loadtxt(pitchFile)
  tonic = np.loadtxt(tonicInFile)
  
  phObj = PH.PitchHistogram(timePitch[:,1], tonic.tolist(), hResolution=1)
  phObj.ComputePitchHistogram(Oct_fold=1)
  phObj.SmoothPitchHistogram()
  phObj.ValidSwarLocEstimation(Oct_fold=1)
  phObj.SwarLoc2Cents()
  phObj.ExtendSwarOctaves()
  swars = phObj.swarCents
  ind = np.argmin(abs(swars))
  if abs(swars[ind]) < 100:
    print "The tonic is off by %f cents"%swars[ind]
    offset = swars[ind]
    newTonic = tonic*np.power(2,offset/1200)
    print "Old Tonic %f, New tonic %f\n"%(tonic, np.array([newTonic]))
    np.savetxt(tonicOutFile, np.array([newTonic]))
  
def batchProcessTonicFineTune(root_dir, tonicExt = '.tonic', tonicExtOut = '.tonicFine', pitchExt = '.pitch'):
    
    audiofilenames = GetFileNamesInDir(root_dir,tonicExt)

    for audiofilename in audiofilenames:
      fname,ext = os.path.splitext(audiofilename)
      try:     
 	fineTuneTonicValue(fname+tonicExt, fname+tonicExtOut, fname+pitchExt)
      except:
	print fname

def fetchMBID(mp3File):
    try:
        mbid = easyid3.ID3(mp3File)['UFID:http://musicbrainz.org'].data
    except:
        print mp3File
        raise MBIDError('MBID not embedded')
    return mbid 

def computeDutationSongs(root_dir, FileExt2Proc='.mp3'): 
    
    audiofilenames = GetFileNamesInDir(root_dir, FileExt2Proc)

    totalLen = 0
    length = []
    mbid_dur = {}
    for audiofile in audiofilenames:
        if  FileExt2Proc=='.mp3':
            audio = MP3(audiofile)
            totalLen += audio.info.length
            length.append(audio.info.length)		
            mbid = fetchMBID(audiofile)
            if not mbid_dur.has_key(mbid):
                mbid_dur[mbid] = audio.info.length
        elif FileExt2Proc=='.wav':
            totalLen +=ES.MetadataReader(filename = audiofile)()[7]
            length.append(ES.MetadataReader(filename = audiofile)()[7]) 
    print "total files %d\n"%len(audiofilenames)
    print "Total length %d\n"%totalLen
    print "Max length %d\n"%np.max(length)
    print "Min length %d\n"%np.min(length)
    print "Mean length %d\n"%np.mean(length)
    print "median length %d\n"%np.median(length)
    
    return mbid_dur
                

# This function is to batch process 
def BatchProcess_PitchExtraction(RootDir, FileExt2Proc = ".mp3", HopSize = 128, FrameSize = 2048, BinResolution = 10, GuessUnvoiced=True, output = "Pitch", VoicingTolerance=0.2, MaxFrequency=20000, extension=".pit_justin.txt", PostProcess =0):
    
    
    audiofilenames = GetFileNamesInDir(RootDir, FileExt2Proc)
    
    for audiofilename in audiofilenames:
        
        outputfilename, ext = os.path.splitext(audiofilename) 
        
        outputfilename = outputfilename + extension
        
        if os.path.exists(outputfilename):
        
            print "file "  + outputfilename + " already exists"
            
        else:
            
            audio = ES.MonoLoader(filename = audiofilename)()
                
            pitch = ES.PredominantMelody(hopSize = HopSize, frameSize = FrameSize, binResolution = BinResolution, guessUnvoiced=GuessUnvoiced, voicingTolerance= VoicingTolerance, maxFrequency = MaxFrequency, minFrequency = 60)(audio)
            
            if output == "Pitch":
                pitch = pitch[0]
            else:
                pitch = pitch[1]        
            
            
            #generating time stamps (because its equally hopped)
            TStamps = np.array(range(0,len(pitch)))*np.float(HopSize)/44100.0
            
            #try to post process pitch to remove octave errors
            
            if(PostProcess==1):    
                pitch = PostProcessPitch(pitch)
            
            dump = np.array([TStamps, pitch]).transpose()
            
            np.savetxt(outputfilename, dump, delimiter = "\t")
            
def BatchConvertMp32Wav(root_dir):
        
        audiofilenames = GetFileNamesInDir(root_dir, '.mp3')
        
        for filename in audiofilenames:
                wavfilename = filename.replace('.mp3','.wav')
                cmd = "lame --decode "+ '"'+filename+'"' + " " + '"'+wavfilename+'"'
                os.system(cmd)
                
def BatchConvertWav2Mp3(root_dir):
        
        audiofilenames = GetFileNamesInDir(root_dir, '.wav')
        
        for filename in audiofilenames:
                mp3filename = filename.replace('.wav', '.mp3')
                cmd = "lame "+ '"'+filename+'"' + " " + '"'+mp3filename+'"'
                os.system(cmd)                
                
def PostProcessPitch(pitch):
    pitch_out = copy.deepcopy(pitch)
    for i in range(0,pitch.shape[0]-1):
        if((pitch[i]*pitch[i+1]!=0)):
            if((abs(((pitch[i+1]/pitch_out[i])-2))<.1)):
                pitch_out[i+1]=pitch[i+1]/2
            if((abs(((pitch[i+1]/pitch_out[i])-3))<.1)):
                    pitch_out[i+1]=pitch[i+1]/3
            if((abs(((pitch[i+1]/pitch_out[i])-4))<.1)):
                    pitch_out[i+1]=pitch[i+1]/4
            if((abs(((pitch[i+1]/pitch_out[i])-0.5))<.05)):
                    pitch_out[i+1]=pitch[i+1]*2                
    return pitch_out
    
    
def BatchProcess_TonicIdentification(RootDir, tonicExt = '.tonic', FileExt2Proc = ".wav"):
    
    audiofilenames = GetFileNamesInDir(RootDir, FileExt2Proc)
    
    for audiofilename in audiofilenames:
        
        path, fileN = os.path.split(audiofilename)
        
        tonic_file = open(path + '/'+ fileN.split('.')[0] + tonicExt, 'w')
        
        audio = ES.MonoLoader(filename = audiofilename)()
            
        tonic = ES.TonicIndianArtMusic()(audio)
        
        MBID = audiofilename.split("/")[-1].strip()
        print MBID
        tonic_file.write(str(tonic)+"\n")
        tonic_file.close()
    

def BatchProcess_PitchSalienceExtraction(RootDir, FileExt2Proc = ".mp3", ExeFile= -1):
  
    
    audiofilenames = GetFileNamesInDir(RootDir, FileExt2Proc)
    
    for audiofilename in audiofilenames:
        print "file being processed "+audiofilename
        #Converting mp3 to wav
        
        wavfilename, ext = os.path.splitext(audiofilename) 
        wavfilename = wavfilename + ".wav"

        ph_file_name, ext = os.path.splitext(audiofilename) 
        ph_file_name = ph_file_name +".mph.txt"
        
        if os.path.exists(ph_file_name):
            print "PH file already exist"
        else:
            cmd = "lame --decode "+ '"'+audiofilename+'"' + " " + '"'+wavfilename+'"'
            os.system(cmd)
            if os.path.exists(wavfilename):
                if ExeFile is not -1:
                    cmd =  ExeFile +" -m H -t I -i " + '"'+wavfilename +'"' + " -o "+ '"'+ph_file_name+'"'
                    print cmd
                    os.system(cmd)
                else:
                    print "Please specify the binary file for execution"
            
                os.remove(wavfilename)
            else:
                print "Hey file " + audiofilename + " didn't produce any wave file"
                 
def BatchProcess_Tonic_CCode(RootDir, FileExt2Proc = ".mp3", ExeFile= -1, outFileExt= ".tonic"):


    audiofilenames = GetFileNamesInDir(RootDir, FileExt2Proc)

    for audiofilename in audiofilenames:
        print "file being processed "+audiofilename


        filename, ext = os.path.splitext(audiofilename)
        tonicfilename = filename + outFileExt
        wavfilename =  audiofilename
        
        if os.path.isfile(tonicfilename):
            print "File already exist %s\n"%tonicfilename
            continue

        #Converting mp3 to wav if mp3 is to be processed\
        if FileExt2Proc == ".mp3":

            wavfilename = filename + ".wav"
            cmd = "lame --decode "+ '"'+audiofilename+'"' + " " + '"'+wavfilename+'"'
            os.system(cmd)

        if ExeFile is not -1:
            cmd =  ExeFile +" -m T -t V -i " + '"'+wavfilename +'"' + " -o "+ '"'+tonicfilename+'"'
            print cmd
            os.system(cmd)
        else:
            print "Please specify the binary file for execution"


        
def GetFileNamesInDir(dir_name, filter=".wav"):
    names = []
    for (path, dirs, files) in os.walk(dir_name):
        for f in files:
            #if filter in f.lower():
            if filter.split('.')[-1].lower() == f.split('.')[-1].lower():
                #print(path+"/"+f)
                #print(path)
                #ftxt.write(path + "/" + f + "\n")
                names.append(path + "/" + f)
                
    return names
      
      
def BatchPRocess_DeleteFiles(RootDir, FileExt2Proc = ".wav"):
    
    audiofilenames = GetFileNamesInDir(RootDir, FileExt2Proc)
    
    for audiofile in audiofilenames:
        print "File removed " + audiofile
        os.remove(audiofile)
    
    
def BatchPRocess_CheckDatabase(RootDir, FileExt2Proc = ".mp3"):
  #note that this function is heavily hardcoded
    
    audiofilenames = GetFileNamesInDir(RootDir, FileExt2Proc)
    
    for audiofile in audiofilenames:
        
        pitchfile, ext = os.path.splitext(audiofile)         
        pitchfile = pitchfile + ".pit.txt"        
        
        phfile, ext = os.path.splitext(audiofile) 
        phfile = phfile +".ph.txt"
        
        if not os.path.exists(phfile):
            print "phfile " + phfile+ "does not exist"
        if not os.path.exists(pitchfile):
            print "pitchfile " + pitchfile + "does not exist"
        


    


def PreprocessMelodicFeatures(RootDir, Ext2Proc = ".wav", FileExts2Proc = [(".tpe",2)], RefFeatExt = (".tpe",2), InterpExt = ".interp", SilVal = 0):
#This function preproces the melody features, by that we mean
#1) Converts all the features to the same sampling rate, determined by the feature with file extension RefFeatExt 
#2) Linearly Interpolate the unvoiced portions in pitch and other features. The reference for detecting unvoiced portions is feature with extension RefFeatExt
# It is important to take pitch as reference because in other features detecting unvoiced part reliably is almost impossible
# INPUTS
# RootDir - root directory which should be searched for reading features which are to be preprocessed
# Ext2Proc - all the files names with this extensions wwill be looked upon for preprocessing. This is just to know what filesnames to process (seed files)
# FileExts2Proc = list of tuples of file extensions (features) to process. The first element of the tuple is file extension and the second element is the index of column which represent this feature because ther emight be time stamps
# RefFeatExt = which feature should become our reference for selecting unvoiced regions (Pitch one ofcourse). So just provide the extension for these files and the index of column of pitch values as tuple
# InterpExt = extension which should be appended after the file name and before the original extension to mark the interpolated file
# SilVal = values in the feautres in RefFeatExt files which should be considered as values for unvoice regions.
# NOTE that big assumption here is that time stamps are always the first column of the filename
    
    #fetching all the audio file names which have extension RefFeatExt
    audiofilenames = GetFileNamesInDir(RootDir, Ext2Proc);
    
    # for each file performing resampling and interpolation
    for audiofile in audiofilenames:
        proc_filename, proc_fileext = os.path.splitext(audiofile);
        time_feature = np.loadtxt(open(proc_filename+RefFeatExt[0],"r"))
        
        ref_time = time_feature[:,0]
        ref_feat = time_feature[:,RefFeatExt[1]-1]
        
        for exts in FileExts2Proc:
            filename, fileext = os.path.splitext(audiofile);
            newfilename = filename + exts[0]
            local_time_feature = np.loadtxt(open(newfilename,"r"))
            local_time = local_time_feature[:,0]
            local_feature = local_time_feature[:,exts[1]-1]
            func= scpyinterp.interp1d(local_time, local_feature, fill_value = SilVal, bounds_error = False)
            #print local_time[-1],ref_time[-1]
            resample_local_feature = func(ref_time)
            sil_ind = np.where(ref_feat==SilVal)[0]
            resample_local_feature[sil_ind]=SilVal
            resample_interp_local_features = InterpolateSilence(resample_local_feature,SilVal)
            local_time_feature= time_feature[:,0:2]
            local_time_feature[:,1]=resample_interp_local_features
            np.savetxt(filename + InterpExt + exts[0], local_time_feature, delimiter = "\t", fmt='%1.6f',)
            

        
def interpolatePitchTracks(timeArrIn, pitchArrIn, timeArrOut, SilVal):
    
    pitchArrOut = np.interp(timeArrOut, timeArrIn, pitchArrIn)
    
    #now deadling with silence regions
    indSil = np.where(pitchArrIn<=SilVal)[0]
    
    interpFactor = float(timeArrIn[1]-timeArrIn[0])/(timeArrOut[1]-timeArrOut[0])
    
    indSilNew1 = np.floor(indSil*interpFactor).astype(np.int)
    indSilNew2 = np.ceil(indSil*interpFactor).astype(np.int)
    
    indSil = np.intersect1d(indSilNew1, indSilNew2)
    
    pitchArrOut[indSil] = SilVal
    
    
    return pitchArrOut
    
    
    
    

"""def renameMP3FilesUsingID3(root_dir, dst_dir, fileext='.mp3'):
    
    filenames = GetFileNamesInDir(root_dir, filter = fileext)
    
    for filename in filenames:
        mdata = eyed3.load(filename)
        artist = mdata.tag.artist
        album =  mdata.tag.album
        title = mdata.tag.title
        path = os.path.dirname(filename)
        outname = path + '/' + artist + '_' + album + '_' + title + fileext
        outname.replace("\s","_")
        shutil.move(filename, outname)
"""

def createValidityFiles(root_dir, searchExt, validFileList):
    """
    This function searches files with extention (searchExt) in the root_dir folder (recursively) and generate a validity file containing (1) for valid files and (0) for non valid files. 
    This kind of arrangement was needed for some codes to not consider file with validity 0.
    
    """
    
    #read the valid file list and generate a dictionary
    lines = open(validFileList, "r").readlines()
    validityInfo = {}
    for line in lines:
        linesplit = line.strip().split()
        validityInfo[linesplit[0]]={}
        validityInfo[linesplit[0]]['validity'] = linesplit[1]
        validityInfo[linesplit[0]]['done'] = 0
    
        
    
    filenames = GetFileNamesInDir(root_dir, searchExt)
    filesNotInList = []
    genCnt = 0
    for f in filenames:
        fn, fext = os.path.splitext(f)
        if validityInfo.has_key(fn):            
            fid=open(fn+'.validity', "w")
            fid.write(validityInfo[fn]['validity'])
            fid.close()
            validityInfo[fn]['done']=1
            genCnt = genCnt+1
        else:
            filesNotInList.append(f)
    filesNoMp3 = []     
    for key in validityInfo.keys():
        if validityInfo[key]['done']==0:
            filesNoMp3.append(key)
    
    print "total files generate %d" % (genCnt)
    print "total files not generate %d" % (len(filesNotInList))
    print "total mp3 files not found %d" % (len(filesNoMp3))
    

def testFileExist(fileList, Extension):
    
    #read the valid file list and generate a dictionary
    lines = open(fileList, "r").readlines()
    nonExist = [] 
    for line in lines:
        fname = line.strip()
        if not fname:
            continue
        fname = fname + Extension
        
        if not os.path.isfile(fname):
            nonExist.append(fname)
            
    print "total non existing files are %d"%(len(nonExist))
    print nonExist
    
def copyFileListForSearchingMotifs(root_dir, fileList, extOut, ext = '.wav'):
    
    filenames = GetFileNamesInDir(root_dir, ext)
    
    for f in filenames:
        filename, ext = os.path.splitext(f)
        filename = filename + extOut        
        shutil.copy(fileList, filename)

def generateFileList(root_dir,fileOut,ext = '.wav'):
    
    filenames = GetFileNamesInDir(root_dir, ext)
    
    fid = open(fileOut, "w")
    for f in filenames:
        filename, ext = os.path.splitext(f)
        fid.write("%s\n"%filename)
    
    fid.close()
        
def batchDeleteFiles(RootDir, delExts= [""]):
    
    for ext in delExts:
        filenames = GetFileNamesInDir(RootDir, ext)        
        for filename in filenames:
            print "File removed " + filename
            os.remove(filename)            
    
def batchRenameFiles(RootDir, renExts= [{}]):            
    """
    provide new and old extionsions like {'oldExt':'xx','newExt':'yy'}
    """
     
    for extExt in renExts:
        filenames = GetFileNamesInDir(RootDir, extExt['oldExt'])        
        for filename in filenames:
            fname ,ext = os.path.splitext(filename)
            os.rename(fname + extExt['oldExt'], fname + extExt['newExt'])
        
        
def validateMotifSearchDB(root_dir, fileout):
    
    filenames = GetFileNamesInDir(root_dir, '.mp3')
    
    ExtensionsRef = ['.flist', '.pitch', '.tonic', '.tablaSec']
    
    fid = open(fileout, "w")
    
    for filename in filenames:
        fname, ext = os.path.splitext(filename)
        
        for ext in ExtensionsRef:
            filecheck = fname + ext
            if not os.path.isfile(filecheck):
                fid.write("%s\t%s\n"%(filename, ext))
                #break
            
    fid.close()
    
def validateMotifSearchOutput(root_dir, fileout):
    
    filenames = GetFileNamesInDir(root_dir, '.2s25Motif_CONF1')
    
    ExtensionsRef = ['.2s25SearchMapp_CONF1', '.2s25MotifSearch_CONF1SqEuclidean', '.2s25MotifSearch_CONF1ShiftLinExp', '.2s25MotifSearch_CONF1ShiftCityBlock', '.2s25MotifSearch_CONF1CityBlock']
    
    fid = open(fileout, "w")
    
    for filename in filenames:
        fname, ext = os.path.splitext(filename)
        
        for ext in ExtensionsRef:
            filecheck = fname + ext
            if not os.path.isfile(filecheck):
                fid.write("%s\t%s\n"%(filename, ext))
                #break
            
    fid.close()    


def checkFileExistance(root_dir, outfile, fileList, Ext2Search = ['.mp3']):
    
    lines = open(fileList).readlines()
    
    fid = open(outfile, "w")
    strFMT = "FILE NAME\t"+"%s\t"*len(Ext2Search)+"\n"
    fid.write(strFMT%tuple(Ext2Search))
    totalCount = [0]*len(Ext2Search)
    for line in lines:
        line = line.strip()
        fid.write("%s\t"%line )
        for ii,ext in enumerate(Ext2Search):
            filename = line  + ext
            if os.path.isfile(filename) and os.path.getsize(filename) > 0:
                fid.write("%d\t"%1)
                totalCount[ii] = totalCount[ii]+1
            else:
                fid.write("%d\t"%0)
                
        fid.write("\n")
    strFMT = "Total count\t"+"%s\t"*len(totalCount)+"\n"
    fid.write(strFMT%tuple(totalCount))
    fid.close()
        
def countNumberOfRemovedPatterns(root_dir, Ext):
    
    filenames = GetFileNamesInDir(root_dir, Ext)
    
    totalPatterns = 0 
    validPatterns = 0
    
    for filename in filenames:
        out = np.loadtxt(filename)
        totalPatterns = totalPatterns + out.size
        validPatterns = validPatterns + out.size - np.sum(out)
        
        
    print "Total number of patterns are %d and total valid patterns are %d\n"%(totalPatterns, validPatterns)


    
def computeLoudness(audioFile, outputExt='.loudness', f0=-1, HopSize = 0.01, FrameSize = 0.04643990929, BinResolution = 10, GuessUnvoiced=True, VoicingTolerance=0.2, MaxFrequency=20000, interpolateLoudness=0, maxSilDurIntp=0.25, smoothLoudness=0):
  """
  This function computes loudness (represented by energy) of the predominant source assuming either you have provided pitch of the predominant melodic source or if f0=-1, it uses Essentia-Melodia to estimate pitch of the predominant melodic source and uses harmonic detection to compute energy (treated as loudness).
  Any sudden gap in the harmonic magnitudes (undetected harmonics) which span a continous time duration < maxSilDurIntp will be interpolated. You should set this value exactly the same you used for interpolating pitch sequence to accound for short intra pattern pauses.
  """
  #reading audio file
  fs = 44100.0#ES.AudioLoader(filename = audioFile)()[1]
  audio = ES.MonoLoader(filename = audioFile, sampleRate = fs)()
  
  
  #obtaining just the file name and splitting extionsion
  fname, ext = os.path.splitext(audioFile)
  
  frameNSamples = np.round(FrameSize*fs).astype(np.int)
  frameNSamples = frameNSamples + np.mod(frameNSamples,2)
  
  #checking the cases, possible types of input parameter f0 
  if type(f0)==int:
    #if its an integer (which essentially means the user has not provided any input and its -1), run the predominant melody estimation and obtain pitch estimate
    pitch = ES.PredominantMelody(hopSize = np.round(HopSize*fs).astype(np.int), frameSize = frameNSamples, binResolution = BinResolution, guessUnvoiced=GuessUnvoiced, voicingTolerance= VoicingTolerance, maxFrequency = MaxFrequency, minFrequency = 60)(audio)[0]
  if type(f0) == str:
    #if its a string that means a user has provided input file name of the pitch file stored in the format <time stamps><pitch value>
    pitch = np.loadtxt(f0)[:,1]
  if type(f0) == np.ndarray:
    # if its an ndarray, this means that the given sequence is the pitch sequence to be used for loudness computation
    pitch = f0
    
  #creating algorithm objects to be used for harmonic detection for each audio frame
  NFFT = (2**np.ceil(np.log2(frameNSamples)+1)).astype(np.int)
  WIN=ES.Windowing()
  SPECTRUM=ES.Spectrum()
  EQUALLOUD = ES.EqualLoudness()
  SPECPEAKS = ES.SpectralPeaks(sampleRate = fs, maxFrequency = 8000)
  HARMDET = ES.HarmonicDetection(nH=30, freqDevThsld=10, sampleRate = fs)#TODO: use the function that Dmitry implemented and not this one!!!
  
  audio_in = EQUALLOUD(audio)
  
  cnt = 0
  harmWghts = []
  for frame in ES.FrameGenerator(audio_in, frameSize = frameNSamples, hopSize = np.round(HopSize*fs).astype(np.int)):
    if cnt >= len(pitch):
        break
    spec = SPECTRUM(WIN(frame))
    peaks = SPECPEAKS(spec)
    wghtsLocal = HARMDET(peaks[0], peaks[1], pitch[cnt])
    harmWghts.append(wghtsLocal)
    cnt+=1

  if interpolateLoudness==1:
    #interpolating harmonic weights
    harmWghts = np.array(harmWghts)
    harmWghtsIntrp = np.zeros(harmWghts.shape)
    for ii in range(harmWghts.shape[1]):
       harmWghts_temp = InterpolateSilence(harmWghts[:,ii], 0, HopSize, maxSilDurIntp)
       harmWghtsIntrp[:,ii] = harmWghts_temp
  else:
    harmWghtsIntrp = harmWghts

  loudness=[]
  for wghtsLocal in harmWghtsIntrp:
    indValid = np.where(wghtsLocal>0)[0]
    loudness.append(np.sqrt(np.sum(np.power(wghtsLocal[indValid],2))))
  
    
  if interpolateLoudness ==1:
    loudness = InterpolateSilence(loudness, 0, HopSize, maxSilDurIntp)
  
  if smoothLoudness ==1:
    loudness = medfilt(loudness,np.round(50.0/(HopSize*1000)).astype(np.int))
  
  #generating time stamps (because its equally hopped)
  TStamps = np.array(range(0,len(loudness)))*np.float(HopSize)
  dump = np.array([TStamps, loudness]).transpose()
  np.savetxt(fname+outputExt, dump, delimiter = "\t")  
  
  
def batchProcessComputeLoudness(root_dir, searchExt = '.wav', outputExt='.loudness', pitchExt='.tpe', f0=-1, HopSize = 0.01, FrameSize = 0.04643990929, BinResolution = 10, GuessUnvoiced=True, VoicingTolerance=0.2, MaxFrequency=20000, interpolateLoudness=0, maxSilDurIntp=0.25, smoothLoudness=0):
  
    filenames = GetFileNamesInDir(root_dir, searchExt)
    
    for filename in filenames:
      fname, ext = os.path.splitext(filename)
      print "processing file: %s"%fname
      computeLoudness(filename, outputExt=outputExt, f0=fname+pitchExt, HopSize = HopSize, FrameSize = FrameSize, BinResolution = BinResolution, GuessUnvoiced=GuessUnvoiced, VoicingTolerance=VoicingTolerance, MaxFrequency=MaxFrequency, interpolateLoudness=interpolateLoudness, maxSilDurIntp=maxSilDurIntp, smoothLoudness=smoothLoudness)
      
    
