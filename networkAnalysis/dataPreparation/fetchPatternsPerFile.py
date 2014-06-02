#!/usr/bin/python
# -*- coding: utf-8 -*-

import psycopg2 as psy
import numpy as np
import os, sys
from mutagen import easyid3

sys.path.append(os.path.join(os.path.dirname(__file__), '../../batchProcessing'))

import batchProcessing as BP

myUser = 'sankalp'
myDatabase = 'motifCarnatic_CONF2'


root_path = '/media/Data/Datasets/MotifDiscovery_Dataset/CompMusic/'

    
        
def fetchMBID(mp3File):
    try:
        mbid = easyid3.ID3(mp3File)['UFID:http://musicbrainz.org'].data
    except:
        print mp3File
        raise MBIDError('MBID not embedded')
    return mbid

        
def getPatternsPerFile():
    
    extension = '.allPatternsInfoConf2'
    
    cmd1 = "select id, filename, mbid from file where hasseed=1"
    cmd2 = "select id, start_time, end_time from pattern where file_id = %d"

    
    try:
        con = psy.connect(database=myDatabase, user=myUser) 
        cur = con.cursor()
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
            np.savetxt(os.path.join(root_path, filename+extension), allPatterns, fmt = ['%ld', '%f', '%f'], delimiter = '\t')
    
    except psy.DatabaseError, e:
        print 'Error %s' % e
        if con:
            con.rollback()
            con.close()
        sys.exit(1)
    
    if con:
        con.close()
        
        
def batchProcessPatternDataDump(root_dir, exePath):
    
    filenames = GetFileNamesInDir(root_dir, '.allPatternsInfoConf2')
    
    for ii, filename in enumerate(filenames):
        print "Processing %d out of %d file\n"%(ii+1, len(filenames))
        print "%s\n"%(filename)
        fname,ext = os.path.splitext(filename)
        cmd = exePath+'/'+'DumpTimeSeries_O3 ' + "\""+ fname + "\""+ " '.pitch' '.tonic' '.allPatternsInfoConf2' '.patternDataConf1' 2.0 5 5 1"
        os.system(cmd)    
            
        
        
        
        
        
        
        
        
        
        
        
    