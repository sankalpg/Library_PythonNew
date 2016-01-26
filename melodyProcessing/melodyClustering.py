import numpy as np
import scipy.signal as sig
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
import copy

min_pitch = 60	#Hz, humanly not possible below this 


def batchProc(root_dir, audioExt = '.mp3', pitchExt = '.pitchSilIntrpPP', transExt = '.transcription', outputExt = '.seg'):
    
    filenames = BP.GetFileNamesInDir(root_dir, '.mp3')
    
    for filename in filenames[:]:
        print "Processing file %s" %filename
        
        #======================
        ## This is done for all
        #======================
        
        fname, ext = os.path.splitext(filename)
        pitch_file = fname + pitchExt
        transcription_file = fname + transExt
        output_file = fname + outputExt
        
        # Dump segentation file
        #----------------------
        generate_segmentation_file(pitch_file, transcription_file, output_file)
        
        
        
        


def find_ind(time_stamps, time):
  ind = np.argmin(np.abs(time_stamps-time))
  return ind


def generate_segmentation_file(pitch_file, transcription_file, output_file):
  """
  
  """
  time_pitch = np.loadtxt(pitch_file)
  
  time_stamps = time_pitch[:,0]
  sec_labels = copy.deepcopy(time_pitch[:,0])*0
  sec_labels = sec_labels + -10000 	#-10000 is for non-silent regions
  
  sil_inds = np.where(time_pitch[:,1] < min_pitch)[0]
  sec_labels[sil_inds] = -100000	#-100000 is for silence regions
  
  trans = np.loadtxt(transcription_file)
  
  for t in trans:
    start_time = t[0]
    end_time = t[1]
    start_ind = find_ind(time_stamps, start_time)
    end_ind = find_ind(time_stamps, end_time)
    sec_labels[start_ind: end_ind +1 ] = t[2]	#  2 is for flat regions
    
  ind_sil = np.where(sec_labels == -100000)[0]
  ind_trans = np.where(sec_labels == -10000)[0]
  ind_flat = np.where(sec_labels >-1500)[0]
  
  seg_sil =  np.array(seg.groupIndices(ind_sil))
  seg_trans =  np.array(seg.groupIndices(ind_trans))
  seg_flat = np.array(seg.groupIndices(ind_flat))
  
  n_sil = seg_sil.shape[0]
  n_trans = seg_trans.shape[0]
  n_flat = seg_flat.shape[0]
  
  start_end_label = []#np.zeros((n_sil+n_trans+n_flat,4))
  
  for s in seg_sil:
    start_end_label.append([s[0], s[1], -100000, np.round(np.mean([s[0],s[1]]))])
  
  for s in seg_trans:
    start_end_label.append([s[0], s[1], -20000, np.round(np.mean([s[0],s[1]]))])
    
  for s in seg_flat:
    start_end_label.append([s[0], s[1], sec_labels[np.round(np.mean([s[0],s[1]]))], np.round(np.mean([s[0],s[1]]))])
  
  start_end_label = np.array(start_end_label)
  
  sort_inds = np.argsort(start_end_label[:,3])
  
  fid = open(output_file, 'w')
  for ii in sort_inds:
    fid.write("%f\t%f\t%s\n"%(time_stamps[start_end_label[ii,0]], time_stamps[start_end_label[ii,1]], start_end_label[ii,2]))
    
  
  
  
  
  
  
  
  
  
  