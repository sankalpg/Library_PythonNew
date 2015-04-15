#!/usr/bin/python
# -*- coding: utf-8 -*-

import psycopg2 as psy
import numpy as np
import os, sys
from mutagen import easyid3
import json

sys.path.append(os.path.join(os.path.dirname(__file__), '../batchProcessing'))

import batchProcessing as BP

myUser = 'sankalp'
#myDatabase = 'motifCarnatic_CONF2'


#root_path = '/media/Data/Datasets/MotifDiscovery_Dataset/CompMusic/'

    
        
def fetchMBID(mp3File):
    try:
        mbid = easyid3.ID3(mp3File)['UFID:http://musicbrainz.org'].data
    except:
        print mp3File
        raise MBIDError('MBID not embedded')
    return mbid

        
def getPatternsPerFile(root_path, myDatabase = '', outExt = '.allPatts'):
    
    cmd1 = "select id, filename, mbid from file where hasseed=1"
    cmd2 = "select id, start_time, end_time from pattern where file_id = %d"

    
    try:
        con = psy.connect(database=myDatabase, user=myUser) 
        cur = con.cursor()
        print "Successfully connected to the server"
        cur.execute(cmd1)
        output = cur.fetchall()
        fileIds = [x[0] for x in output]
        audiofiles = [x[1] for x in output]
        mbids = [x[2] for x in output]
        
        for ii, fileId in enumerate(fileIds):
            filename, ext = os.path.splitext(audiofiles[ii])
            print "processing %d out of %d files"%(ii+1, len(fileIds))
            
            cur.execute(cmd2%(fileId))
            allPatterns = cur.fetchall()
            allPatterns = np.array(allPatterns)
            np.savetxt(os.path.join(root_path, filename+outExt), allPatterns, fmt = ['%ld', '%f', '%f'], delimiter = '\t')
    
    except psy.DatabaseError, e:
        print 'Error %s' % e
        if con:
            con.rollback()
            con.close()
        sys.exit(1)
    
    if con:
        con.close()
        
        
def batchProcessPatternDataDump(root_dir, exePath):
    """
    this is a batch function for dumping data for the extracted pattern for each file. it uses a C code which does exactly the same preprocessing as was done at the time of pattern extraction.
    """
    
    filenames = BP.GetFileNamesInDir(root_dir, '.allPatternsInfoConf2')
    
    for ii, filename in enumerate(filenames):
        print "Processing %d out of %d file\n"%(ii+1, len(filenames))
        print "%s\n"%(filename)
        fname,ext = os.path.splitext(filename)
        cmd = exePath+'/'+'DumpTimeSeries_O3 ' + "\""+ fname + "\""+ " '.pitch' '.tonic' '.allPatternsInfoConf2' '.patternDataConf1' 2.0 5 5 1"
        os.system(cmd)    
            
        
def combinePatternInfoAndDataFiles(root_dir, infoExt, dataExt, outInfoFile, outDataFile):
    
    filenames = BP.GetFileNamesInDir(root_dir, infoExt)
        
    for ii, filename in enumerate(filenames):
        print "Processing %d out of %d file\n"%(ii+1, len(filenames))
        print "%s\n"%(filename)
        
        fid1 = open(outInfoFile, "a")
        fid2 = open(outDataFile, "ab")
        
        fname, ext  = os.path.splitext(filename)
        
        Info = np.loadtxt(fname + infoExt)
        np.savetxt(fid1, Info, fmt = ['%ld', '%f', '%f'], delimiter = '\t')
        
        data = open(fname + dataExt, 'rb').read()
        
        fid2.write(data)
        
        fid1.close()
        fid2.close()
        
        
def getMBID_RaagaID_For_PatternID(outputFile, myDatabase= ''):
    """
    This function dumps all the mbid and raagaid associated with a pattern id. 
    The idea is to use this info while community ranking. Instead of fetching it from the datase,
    I want to load it offfline. Database query is taking time.
    """
    info = {}
    con = None
    try:
        con = psy.connect(database=myDatabase, user=myUser) 
        cur = con.cursor()
        query = "select raagaid, mbid from file where id = (select file_id from pattern where id = %s)"
        cur.execute("select count(*) from pattern")
        nPatterns = cur.fetchone()[0]
        cnt=0
        print "Total number of patterns are %d\n"%nPatterns
        for ii in range(1,nPatterns+1):
            cur.execute(query%ii)
            ragaId, mbid = cur.fetchone()
            info[ii]={}
            info['ragaId'] = ragaId
            info['mbid'] = mbid
            cnt+=1
            if ii%100==0:
                print cnt
            
    except psy.DatabaseError, e:
        print 'Error %s' % e
        if con:
            con.rollback()
            con.close()
        sys.exit(1)
    
    if con:
        con.close()  
        
    json.dump(info, open(outputFile,'w'))
        
        
        
        
        
        
        
        
    