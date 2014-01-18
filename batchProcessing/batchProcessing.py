#!/usr/bin/env python


########################################################3
# Author: Sankalp Gulati (sankalp.gulati@upf.edu)
# Other Authors for individual functions are indicated inside the function!!
# Music Technology Group - Universitat Pompeu Fabra
# Created:     18/09/2012
# Last modified: --
########################################################3

import essentia.standard as ES
import os, sys, copy, shutil
import numpy as np
#import eyed3
import scipy.interpolate as scpyinterp

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
            if filter in f.lower():
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
        

def InterpolateSilence(array, silence_val):
    
    if type(array) ==list:
        array = np.array(array)    
    
    # interpolating zeros in middle of the time series
    array = array.astype('float')
    sil_ind = np.where(array==silence_val)[0]    
    last_sil_ind=sil_ind[0]
    sil_ind = np.append(sil_ind,0)
    length= len(array)
    for ii in range(0,len(sil_ind)-1):
        if sil_ind[ii] +1 != sil_ind[ii+1]:
            if last_sil_ind!=0 and sil_ind[ii]!=length-1:
                inter_data = np.linspace(array[last_sil_ind-1],array[sil_ind[ii]+1] , sil_ind[ii]-last_sil_ind+3)
                array[last_sil_ind-1:sil_ind[ii]+2] = inter_data
                last_sil_ind=sil_ind[ii+1]
            else:
                last_sil_ind=sil_ind[ii+1]
    
    # zeros at the beginning and at the end of the time series are replaced by the mean of the time series
    sil_ind = np.where(array==silence_val)[0]
    allind = np.array(range(0,length))
    ind_nonsilence = np.extract(np.invert(np.in1d(allind,sil_ind)),allind)
    nonsilvals = array[ind_nonsilence]
    #print np.sum(nonsilvals)
    array[sil_ind] = np.mean(nonsilvals)
    
    return array
        
def BatchProcessInterpPitchSilence(RootDir, FileExt2Proc = ".tpe", NewExt = "", PitchCol=2, SilVal=0):
    
    audiofilenames = GetFileNamesInDir(RootDir, FileExt2Proc)
    
    if len(NewExt)==0:
        NewExt = ".interp" + FileExt2Proc
    
    for audiofile in audiofilenames:
          time_pitch = np.loadtxt(open(audiofile,"r"))
          new_pitch = InterpolateSilence(time_pitch[:,PitchCol-1],SilVal)
          time_pitch[:,PitchCol-1]=new_pitch
          file,ext = os.path.splitext(audiofile)
          np.savetxt(file + NewExt, time_pitch, delimiter = "\t", fmt='%1.6f',)


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