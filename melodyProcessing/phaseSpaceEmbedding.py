import numpy as np
import os, sys
from scipy.signal import triang, filtfilt, gaussian
import matplotlib.pyplot as plt
import pickle
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
  return np.sqrt(np.sum(np.power(np.power(mtx1,0.1)-np.power(mtx2,0.1),2)))

  
def ragaRecognition(file_list, database = '', user= '', phase_ext = '.phasespace'):

  raga_mbid = rr.get_mbids_raagaIds_for_collection(file_list, database = database, user= user)
  mbid2raga = {}
  for r in raga_mbid:
    mbid2raga[r[1]] = r[0]


  lines = open(file_list, 'r').readlines()

  mtx = []
  ind2mbid = {}
  for ii, line in enumerate(lines):
    phase_file = line.strip() + phase_ext
    mtx.append(np.loadtxt(phase_file))
    mbid = rr.get_mbid_from_mp3(line.strip()+'.mp3')
    ind2mbid[ii] = mbid


  dist = 100000000000*np.ones((len(mtx),len(mtx)))
  for ii in range(len(mtx)):
    for jj in range(len(mtx)):
      if ii == jj:
        continue
      dist[ii,jj] = phase_space_dist(mtx[ii], mtx[jj])

  nearest_dists = []
  cnt =0
  for ii in range(len(mtx)):
    nearest_dists.append(np.argsort(dist[ii,:])[0])
    if mbid2raga[ind2mbid[ii]] == mbid2raga[ind2mbid[nearest_dists[-1]]]:
      cnt+=1
  print cnt

  np.savetxt('dist.txt', nearest_dists)

