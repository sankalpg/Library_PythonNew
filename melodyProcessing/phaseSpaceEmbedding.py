import numpy as np
import os, sys
from scipy.signal import triang, filtfilt, gaussian
import matplotlib.pyplot as plt
import pickle
import scipy.ndimage as ndimage
sys.path.append(os.path.join(os.path.dirname(__file__), '../batchProcessing'))
sys.path.append(os.path.join(os.path.dirname(__file__), '../networkAnalysis'))

import batchProcessing as BP
import ragaRecognition as rr



def getDelayCoordinates(inp, m, t):
  """
  This function extract phase space embedding given an input signal (inp), number of delay 
  corrdinates (m), and delay per corrdinate (t) in samples.
  """
  if not type(t) == int:
    print "Provide delay per coordinate as an integer (samples)"
  
  tDelay = (m-1)*t
  nSamples = len(inp)
  outArray = np.zeros((nSamples-tDelay, m))
  for ii in range(m):
    outArray[:,ii] = inp[tDelay-(ii*t):nSamples-(ii*t)]
    
  return outArray

def getCentBinMapp(range_cents = (-1200, 3600), resolution = 10, oct_fold = 0):

  mapp = {}
  bins = []

  if oct_fold == 1:
    range_cents = (0, 1199)
  for r in range(range_cents[0], range_cents[1]+1):
    bin = int(-1*range_cents[0]/float(resolution)) + int(np.round(r/float(resolution)))
    mapp[r] = bin
    bins.append(bin)

  return mapp, np.unique(np.array(bins))

def get2DRagaRepresentation(delay_coord, tonic, resolution = 10, range_cents = (-1200, 3600), oct_fold = 0):
  """
  Given the delay coordinates in Hz and tonic, this function extracts a 2 dimensions mapping
  """
  mapp, bins = getCentBinMapp(range_cents, resolution, oct_fold)

  mtx = np.zeros((len(bins), len(bins)))
  if oct_fold == 1:
    for d in delay_coord:
      x_coord = int(1200*np.log2(d[0]/tonic))%1200
      y_coord = int(1200*np.log2(d[1]/tonic))%1200
      mtx[mapp[x_coord], mapp[y_coord]]+=1
  else:
    for d in delay_coord:
      x_coord = int(1200*np.log2(d[0]/tonic))
      y_coord = int(1200*np.log2(d[1]/tonic))
      mtx[mapp[x_coord], mapp[y_coord]]+=1

  return mtx


def getPhaseSpaceEmbPitchSignal(pitch_file, delay_tap, tonic, min_pitch = 60, resolution =5, range_cents = (-1200, 3600), out_ext = '.phasespace', oct_fold = 0):
  """
  This function computes phrase space embedding (delay coordinates) for pitch pitch_file
  """  

  time_pitch = np.loadtxt(pitch_file)
  hop_size = time_pitch[1,0] - time_pitch[0,0]

  delay_tap = int(np.round(delay_tap/hop_size))

  delay_coord = getDelayCoordinates(time_pitch[:,1], 2, delay_tap)

  ind_sil1 = np.where(delay_coord[:,0] < min_pitch)[0]
  ind_sil2 = np.where(delay_coord[:,1] < min_pitch)[0]

  ind_sil = np.unique(np.append(ind_sil1, ind_sil2))

  delay_coord = np.delete(delay_coord, ind_sil, axis=0)

  mtx = get2DRagaRepresentation(delay_coord, tonic, resolution = resolution, range_cents = range_cents, oct_fold = oct_fold)
  fname, ext = os.path.splitext(pitch_file)

  np.savetxt(fname + out_ext, mtx)


def BatchProcessGetPhaseSpaceEmbPitchSignal(root_dir, pitch_ext, out_ext, tonic_ext, delay_tap, oct_fold):


  filenames = BP.GetFileNamesInDir(root_dir, '.mp3')

  for filename in filenames:
    fname, ext = os.path.splitext(filename)
    tonic = np.loadtxt(fname + tonic_ext)
    getPhaseSpaceEmbPitchSignal(fname + pitch_ext, delay_tap, tonic, min_pitch = 60, resolution =10, range_cents = (0, 1200), out_ext = out_ext, oct_fold = oct_fold)
  

def phase_space_dist(mtx1, mtx2):
  return np.sqrt(np.sum(np.power(mtx1-mtx2, 2)))

  
def ragaRecognitionPhaseSpaceKNN(file_list, database = '', user= '', phase_ext = '.phasespace', smooth_gauss_sigma = -1, compression = -1, normalize = -1, dist_metric = 'Euclidean', KNN=1):
  """
  This function performs raga recognition using phase space embeddings and KNN classifier.

  There are options for different types of distane measures and pre-processing parameters of the embeddings.
  """

  #available distance measure
  dist_metrics = {'Euclidean': phase_space_dist}

  if not dist_metric in dist_metrics.keys():
    print "Please choose a valid distance metric that is supported by this function"
    return False

  #fetching mbid and raga list for collection of files in file_list (ragas are fetched from the database)
  raga_mbid = rr.get_mbids_raagaIds_for_collection(file_list, database = database, user= user)
  mbid2raga = {}  #creading an mbid to raga id mapping
  for r in raga_mbid:
    mbid2raga[r[1]] = r[0]

  #reading files listed in file_list
  lines = open(file_list, 'r').readlines()

  phase_space_embeds = []  #matrix to store 
  ind2mbid = {} #index (in file_list) to mbid mapping
  for ii, line in enumerate(lines):
    phase_file = line.strip() + phase_ext
    mtx = np.loadtxt(phase_file)

    #pre-processing phase space embedding based reprsentation
    
    #smoothening operation?
    if smooth_gauss_sigma >0:
      mtx = ndimage.gaussian_filter(mtx, smooth_gauss_sigma)
    
    #dymanic compression
    if compression > 0:
      mtx = np.power(mtx, compression)
    
    #Normalization (max val = 1)
    if normalize > 0:
      mtx = mtx/np.max(mtx)

    #appending mtx in the order of filelists
    phase_space_embeds.append(mtx)
    mbid = rr.get_mbid_from_mp3(line.strip()+'.mp3')
    ind2mbid[ii] = mbid

  dist = np.finfo(np.float).max*np.ones((levalFeaturesen(phase_space_embeds),len(phase_space_embeds)))
  for ii in range(len(phase_space_embeds)):
    for jj in range(len(phase_space_embeds)):
      if ii == jj:
        continue
      dist[ii,jj] = dist_metrics[dist_metric](phase_space_embeds[ii], phase_space_embeds[jj])

  nearest_dists = []
  cnt =0
  for ii in range(len(phase_space_embeds)):
    nearest_dists.append(np.argsort(dist[ii,:])[0])
    if mbid2raga[ind2mbid[ii]] == mbid2raga[ind2mbid[nearest_dists[-1]]]:
      cnt+=1
  print cnt

  np.savetxt('dist.txt', nearest_dists)

