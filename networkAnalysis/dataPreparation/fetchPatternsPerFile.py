#!/usr/bin/python
# -*- coding: utf-8 -*-

import psycopg2 as psy
import numpy as np
import os, sys
from mutagen import easyid3

sys.path.append(os.path.join(os.path.dirname(__file__), '../../batchProcessing'))

import batchProcessing as BP

myUser = 'sankalp'
myDatabase = 'motifCarnatic_CONF1'


root_path = '/media/Data/Datasets/MotifDiscovery_Dataset/CompMusic/'

    
        
def fetchMBID(mp3File):
    try:
        mbid = easyid3.ID3(mp3File)['UFID:http://musicbrainz.org'].data
    except:
        print mp3File
        raise MBIDError('MBID not embedded')
    return mbid

        
def getPatternsPerFile():
    
    versionIndexSelected = 1
    nTopSearch = 20
    extension = '.seedPlus20Search'
    
    cmd1 = "select id, filename, mbid from file where hasseed=1"
    cmd2 = "select id, start_time, end_time from pattern where file_id = %d and isseed=1"
    cmd3 = "select target_id from match where source_id = %d and version = %d order by distance limit %d "
    cmd4 = "select id, start_time, end_time from pattern where id=%d"
    
    try:
        con = psy.connect(database=myDatabase, user=myUser) 
        cur = con.cursor()
        cur.execute(cmd1)
        output = cur.fetchall()
        fileIds = [x[0] for x in output]
        audiofiles = [x[1] for x in output]
        mbids = [x[2] for x in output]
        
        for ii, fileId in enumerate(fileIds):
            filePatterns=[]
            filename, ext = os.path.splitext(audiofiles[ii])
            
            print "processing %d out of %d files"%(ii+1, len(fileIds))
            
            cur.execute(cmd2%(fileId))
            seedPatterns = cur.fetchall()
            
            for seedPattern in seedPatterns:
                filePatterns.append(seedPattern)
                cur.execute(cmd3%(seedPattern[0], versionIndexSelected, nTopSearch))
                searchPatterns = cur.fetchall()

                for searchPattern in searchPatterns:
                    cur.execute(cmd4%(searchPattern[0]))
                    filePatterns.append(cur.fetchone())

            filePatterns = np.array(filePatterns)
            np.savetxt(os.path.join(root_path, filename+extension), filePatterns, fmt = ['%ld', '%f', '%f'], delimiter = '\t')
        
    
    except psy.DatabaseError, e:
        print 'Error %s' % e
        if con:
            con.rollback()
            con.close()
        sys.exit(1)
    
    if con:
        con.close()
        
        
    
            
        
        
        
        
        
        
        
        
        
        
        
    