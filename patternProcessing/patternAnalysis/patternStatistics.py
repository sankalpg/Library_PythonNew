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




def createISMIR2014EvaluationSubset(nBins=10, totalPatterns=200, nSearchItems=10, nVersions=4, splitLogic = 1):
    """
    this method will sample the seed pattern space and will select a subset such that they are equally disctibuted over stratas of 
    distances, where each strata is basically equi-spaced bin from min to max in log distance domain.
    
    In addtion for each seed pattern top N searched patterns are selected for different versions of the rank refinement.
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
    
    #constructing bins with different different logic, very crucial step.
    
    ### Logic 1: Equal width bins from min_val to max_val
    if splitLogic ==1:
        bins = np.linspace(min_val, max_val, num=nBins+1)
    
    ### Logic 2: for the case of three bins taking  min_val < bin1 < mean-2*std < bin2 < mean+2*std < bin3 < max_val
    if splitLogic ==2:
        if nBins ==3:
            mean = np.mean(distArrayLOG)
            std = np.std(distArrayLOG)
            threshold1 = mean - 2*std
            threshold2 = mean + 2*std
            bins = [min_val, threshold1, threshold2, max_val]
    
    
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
    evalSubSet = -1*np.ones((2+(nVersions*nSearchItems), totalPatterns))
    
    cmd1 = "select target_id from match where source_id=%d and version= %d order by distance"
    cmd2 = "select pair_id from pattern where id=%d"

    try:
        con = psy.connect(database=myDatabase, user=myUser) 
        cur = con.cursor()
        
        for ii,seed in enumerate(seedSubset):
            evalSubSet[1,ii]=seed
            cur.execute(cmd2%seed)
            seedPair = int(cur.fetchone()[0])
            evalSubSet[0,ii]=seedPair
            
            for jj in range(nVersions):
                cur.execute(cmd1%(seed,jj))   
                searchData = cur.fetchall()
                searchData = [x[0] for x in searchData]
                searchData = np.array(searchData[:nSearchItems])
                evalSubSet[2+(jj*nSearchItems) : 2 + ((jj+1)*nSearchItems),ii] = searchData
 
        
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


def generateHeatMapPlotForSeedVsSearchEvalSubSetISMIR(patternInfoFile, versionSet= [0]):
    
    #reading the pattern ids from info filename
    patternData = np.loadtxt(patternInfoFile)
    
    seedPatterns = patternData[1,:]
    seedPatternPairs = patternData[0,:]
    
    searchPatternsAll = patternData[2:,:]
    
    cmd1 = "select distance in match where source_id=%d and target_id=%d"
    
    try:
        con = psy.connect(database=myDatabase, user=myUser) 
        cur = con.cursor()
        distVals = np.array([[0,0]])
        for ii, seed in enumerate(seedPatterns):
            for version in versionSet:
                searchPatterns = searchPatternsAll[ii,version*10:(version+1)*10]
                
                for searches in searchPatterns:
                    cur.execute(cmd1%(seed, searches))
                    distVals.append(np.array([[]]))
                
                
                    
            
        
        
        
        seedData = cur.fetchall()
        
    except psy.DatabaseError, e:
        print 'Error %s' % e
        if con:
            con.rollback()
            con.close()
        sys.exit(1)

    if con:
        con.close()
    
    
        
    

    
    

