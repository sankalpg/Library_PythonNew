#!/usr/bin/python
# -*- coding: utf-8 -*-

import psycopg2 as psy
import sys, os
from mutagen import easyid3
import numpy as np
import pickle

sys.path.append(os.path.join(os.path.dirname(__file__), '../../batchProcessing'))

import batchProcessing as BP
import matplotlib.pyplot as plt

try:
    from mutagen.mp3 import MP3
except:
    pass

myUser = 'sankalp'
myDatabase = 'motifDB_CONF1'

root_path = '/media/Data/Datasets/MotifDiscovery_Dataset/CompMusic/'


def getAudioDurationForDataSet():
    
    cmd1 = "select filename from file where hasseed=1"
    
    try:
        con = psy.connect(database=myDatabase, user=myUser) 
        cur = con.cursor()
        cur.execute(cmd1)
        audiofiles = cur.fetchall()
        audiofiles = [x[0] for x in audiofiles]
        
    except psy.DatabaseError, e:
        print 'Error %s' % e
        if con:
            con.rollback()
            con.close()
        sys.exit(1)
    
    if con:
        con.close()
        
    total_duration = 0
    for audiofile in audiofiles:
        total_duration += computeDutationSong(root_path+audiofile)
        
    print "total duration is %d"%total_duration
    print "total number of files is %d"%len(audiofiles)
    
    return total_duration

def computeDutationSong(audiofile): 
    
    filename, ext = os.path.splitext(audiofile)

    if  ext=='.mp3':
        audio = MP3(audiofile)
        duration = audio.info.length
    elif ext=='.wav':
        duration = ES.MetadataReader(filename = audiofile)()[7]
            
    return duration

def getSeedPatternDistancesTXT(root_dir, seedExt = '.2s25Motif_CONF1', distFile = 'seedDistances'):
    
    filenames = BP.GetFileNamesInDir(root_dir,seedExt)
    distArray = []
    for ii, filename in enumerate(filenames):
        #print "processing %d of %d files"%(ii+1, len(filenames))
        seedData = np.loadtxt(filename)
        indValid = np.where(seedData[:,4]<99999999999999999999)[0]
        distArray.extend(seedData[indValid,4].tolist())
        
    #np.save(distFile, distArray)
    
    return distArray

def computeSeedPatternDistHistogramTXT(root_dir, seedExt = '.2s25Motif_CONF1', nBins=100, plotHist=0):
    
    dist = getSeedPatternDistances(root_dir, seedExt)
    
    dist=np.log10(dist)
    
    min_val = np.min(dist)
    max_val = np.max(dist)
    
    #bins = np.arange(0,max_val, 10000)
    bins = np.linspace(min_val, max_val, num=nBins+1)
    
    hist = np.histogram(dist,bins=bins)
    
    if plotHist:
        fig = plt.figure()
        fsize=14
        plt.plot(hist[1][:-1], hist[0])
        plt.ylabel("Frequency", fontsize=fsize)
        plt.xlabel("Log distance", fontsize=fsize)
        fig.savefig('seedDistanceDistribution.pdf')
    
    return hist[0], hist[1]

def getSeedPatternDistancesDB():
    
    cmd1 = "select distance from match where version <0"
    
    distArray = []
    
    try:
        con = psy.connect(database=myDatabase, user=myUser) 
        cur = con.cursor()
        cur.execute(cmd1)
        distArray = cur.fetchall()
       
    except psy.DatabaseError, e:
        print 'Error %s' % e
        if con:
            con.rollback()
            con.close()
        sys.exit(1)
    
    if con:
        con.close()
    
    return distArray    
    
def computeSeedPatternDistHistogramDB(nBins=100, plotHist=0):
    
    dist = getSeedPatternDistancesDB()
    
    dist=np.log10(dist)
    
    min_val = np.min(dist)
    max_val = np.max(dist)
    
    #bins = np.arange(0,max_val, 10000)
    bins = np.linspace(min_val, max_val, num=nBins+1)
    
    hist = np.histogram(dist,bins=bins)
    
    if plotHist:
        fig = plt.figure()
        fsize=14
        plt.plot(hist[1][:-1], hist[0])
        plt.ylabel("Frequency", fontsize=fsize)
        plt.xlabel("Log distance", fontsize=fsize)
        fig.savefig('seedDistanceDistribution.pdf')
    
    return hist[0], hist[1]




def createEvaluationSubset(nBins=10, totalPatterns=200, nSearchItems=10, nVersions=4):
    """
    this method will sample the seed pattern space and will select a subset such that they are equally disctibuted over stratas of 
    distances, where each strata is basically equi-spaced bin from min to max in log distance domain.
    """    

    cmd1 = "select source_id, distance from match where version =-1 order by distance"

    try:
        con = psy.connect(database=myDatabase, user=myUser) 
        cur = con.cursor()
        cur.execute(cmd1)
        seedData = cur.fetchall()
        
    except psy.DatabaseError, e:
        print 'Error %s' % e
        if con:
            con.rollback()
            con.close()
        sys.exit(1)

    if con:
        con.close()

    seedData = np.array(seedData)
    distArray = seedData[:,1]
    
    distArrayLOG=np.log10(distArray)

    min_val = np.min(distArrayLOG)
    max_val = np.max(distArrayLOG)

    bins = np.linspace(min_val, max_val, num=nBins+1)
    
    #just store indexes for each strata
    indStratas=[]
    for ii in range(len(bins)-1):        
        indMore = np.where(distArrayLOG>=bins[ii])[0]
        indLess = np.where(distArrayLOG<bins[ii+1])[0]
        inds = np.intersect1d(indMore, indLess)
        indStratas.append(inds)
    
    #performing a greedy sampling to cover all bins, selecting one from every bin unless total cnt is reached
    indexSample = []
    category = []
    cnt = totalPatterns
    binIndex = 0
    perBinSamples = np.zeros(nBins)
    while(cnt>0):
        binIndex = np.mod(binIndex,nBins)
        leftSize = indStratas[binIndex].size
        
        if leftSize>0:
            indRand = np.random.randint(leftSize)
            indexSample.extend([indStratas[binIndex][indRand]])
            category.extend([binIndex])
            indStratas[binIndex] = np.delete(indStratas[binIndex],indRand)
            perBinSamples[binIndex]+=1
            cnt-=1
        binIndex+=1
            
    seedSubset = seedData[indexSample,0]
    distSubset = seedData[indexSample,1]
    
    #ok so after subsampling now lets fetch the searched patterns
    evalSubSet = -1*np.ones((1+(nVersions*nSearchItems), totalPatterns))
    
    cmd1 = "select target_id from match where source_id=%d and version= %d order by distance"

    try:
        con = psy.connect(database=myDatabase, user=myUser) 
        cur = con.cursor()
        
        for ii,seed in enumerate(seedSubset):
            evalSubSet[0,ii]=seedSubset[ii]
            for jj in range(nVersions):
                cur.execute(cmd1%(seed,jj))   
                searchData = cur.fetchall()
                searchData = [x[0] for x in searchData]
                searchData = np.array(searchData[:nSearchItems])
                evalSubSet[1+(jj*nSearchItems) : 1 + ((jj+1)*nSearchItems),ii] = searchData
 
        
    except psy.DatabaseError, e:
        print 'Error %s' % e
        if con:
            con.rollback()
            con.close()
        sys.exit(1)
        
    if con:
       con.close()
    
    data = {'seedSubset':seedSubset, 'distSubset':distSubset, 'category':category, 'bins':bins,  'perBinSamples':perBinSamples}
    fid = open("evaluationFullData.pkl", "w")
    pickle.dump(data, fid)
    fid.close()
    
    np.savetxt("patternInfo.txt", evalSubSet.astype(np.int32), fmt="%ld")
    np.savetxt("evalationInfo.txt", 0*evalSubSet.astype(np.int32), fmt="%ld")
    
    return 1

    
    

