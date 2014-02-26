#!/usr/bin/python
# -*- coding: utf-8 -*-

import psycopg2 as psy
import sys, os
from mutagen import easyid3
import numpy as np

sys.path.append(os.path.join(os.path.dirname(__file__), '../../batchProcessing'))


import batchProcessing as BP


myUser = 'sankalp'
myDatabase = 'motif_local'

serverPrefix = "/homedtic/sgulati/motifDiscovery/dataset/carnatic/compMusic/"
localPrefix = "/media/Data/Datasets/MotifDiscovery_Dataset/CompMusic/"


def resetAllTables():
    
    cmd1 = "DELETE FROM match"
    cmd2 = "DELETE FROM pattern"
    cmd3 = "DELETE FROM file"
    cmd4 = "ALTER SEQUENCE match_id_seq RESTART WITH 1"
    cmd5 = "ALTER SEQUENCE pattern_id_seq RESTART WITH 1"
    cmd6 = "ALTER SEQUENCE file_id_seq RESTART WITH 1"
    
    try:
        con = psy.connect(database=myDatabase, user=myUser) 
        cur = con.cursor()
        cur.execute(cmd1)
        cur.execute(cmd2)
        cur.execute(cmd3)
        cur.execute(cmd4)
        cur.execute(cmd5)
        cur.execute(cmd6)
        con.commit()
        print "Successfully updated file table in %s database"%(myDatabase)

    except psycopg2.DatabaseError, e:
        print 'Error %s' % e
        if con:
            con.rollback()
            con.close()
        sys.exit(1)
    
    if con:
        con.close()
        


def createFileTable(root_dir, filterExt = '.mp3'):
    
    con = None
    audioFiles = BP.GetFileNamesInDir(root_dir, filterExt)
    dataDump = []
    query = "INSERT INTO file (filename, mbid) VALUES (%s, %s)"
    
    
    for audiofile in audioFiles:
        
        try:
            mbid = fetchMBID(audiofile)
        except:
            print "MBID not embedded in file %s"%audiofile
        else:
            #removing prefix
            if audiofile.count(serverPrefix):
                audiofile_WOPre = audiofile.split(serverPrefix)[1]
            elif audiofile.count(localPrefix):
                audiofile_WOPre = audiofile.split(localPrefix)[1]
            else:
                print "please provide files with known prefixes (paths)"
                sys.exit(1)
                
            dataDump.append((audiofile_WOPre, mbid))
        
    print "Dump array successfully created"
    try:
        con = psy.connect(database=myDatabase, user=myUser) 
        cur = con.cursor()
        cur.executemany(query, dataDump)
        con.commit()
        print "Successfully updated file table in %s database"%(myDatabase)

    except psycopg2.DatabaseError, e:
        print 'Error %s' % e
        if con:
            con.rollback()
            con.close()
        sys.exit(1)
    
    if con:
        con.close()


def createPatternMatchTable(root_dir, motifDiscExt, motifSearchExt, motifSearchMappExt, version, resetPattern= 0 , resetMatch = 0):
    
    fileNames = BP.GetFileNamesInDir(root_dir, motifDiscExt)
    
    
    #commands for doing different tasks
    cmd1 = "INSERT INTO pattern (file_id, start_time, end_time, version) VALUES (%ld, %f, %f, %d)"#storing seed motifs
    cmd2 = "SELECT id FROM file WHERE mbid = '%s'"#searching file id
    cmd3 = "INSERT INTO pattern (id, pair_id) VALUES (%ld, %ld)"#storing seed motifs
    cmd4 = "SELECT currval('pattern_id_seq')"
    cmd5 = "UPDATE pattern SET pair_id = %ld WHERE id = %ld"
    cmd6 = "INSERT INTO match (source_id, target_id, distance) VALUES (%ld, %ld, %f)"
    cmd7 = "SELECT id FROM file WHERE filename = $$%s$$"#searching file id
    
    
    con = None
    try:
        con = psy.connect(database=myDatabase, user=myUser) 
        cur = con.cursor()
        print "Successfully connected to %s database"%(myDatabase)
        
    except psycopg2.DatabaseError, e:
        print 'Error %s' % e
        if con:
            con.close()
        sys.exit(1)

    try:
        for filename in fileNames:
            
            fname, ext = os.path.splitext(filename)
            motifSeedFile = fname + motifDiscExt
            motifSearchFile = fname + motifSearchExt
            motifMappFile = fname + motifSearchMappExt
            audiofile = fname + '.mp3'
            
            
            #reading discovered motifs file
            seedMotifData = np.loadtxt(motifSeedFile)
            
            #reading searched motifs file
            searchMotifData = np.loadtxt(motifSearchFile)
            searchMappData = open(motifMappFile,'r').readlines()
            
            #storing seed patterns
            try:
                mbid_seed_file = fetchMBID(audiofile)
            except:
                print "MBID not embedded in file %s"%audiofile
                if con:
                    con.rollback()
                    con.close()
                sys.exit(1)
            
            #storing seed + searched patterns in database
            #array to store the pattern ids
            seedPatternIds = []
            for ii in range(0, seedMotifData.shape[0]):
                
                if seedMotifData[ii][4]>9999999999:
                    break
                
                cur.execute(cmd2%(mbid_seed_file))
                file_id = cur.fetchone()[0]
                
                #entering in table pattern the first instance of seed pair
                cur.execute(cmd1%(file_id, seedMotifData[ii][0], seedMotifData[ii][1], version))
                cur.execute(cmd4)
                pattern_id1 = cur.fetchone()[0]
                seedPatternIds.append(pattern_id1)
                
                #entering in table pattern the second instance of seed pair
                cur.execute(cmd1%(file_id, seedMotifData[ii][2], seedMotifData[ii][3], version))
                cur.execute(cmd4)
                pattern_id2 = cur.fetchone()[0]
                seedPatternIds.append(pattern_id2)
                
                
                #cross referencing these two ids
                cur.execute(cmd5%(pattern_id2,pattern_id1))
                cur.execute(cmd5%(pattern_id1,pattern_id2))
                
                #also inserting match info about the seed pair
                cur.execute(cmd6%(pattern_id1, pattern_id2, seedMotifData[ii][4]))
                cur.execute(cmd6%(pattern_id2, pattern_id1, seedMotifData[ii][4]))
            
            con.commit()
            
            #start inserting information about the searched patterns.
            for line in searchMappData:
                fileSearched, start, end = line.split('\t')
                fileSearched = fileSearched.strip() + '.mp3'
                start = int(start.strip())-1
                end = int(end.strip())
                
                #removing prefix
                if fileSearched.count(serverPrefix):
                    fileSearched_WOPre = fileSearched.split(serverPrefix)[1]
                elif audiofile.count(localPrefix):
                    fileSearched_WOPre = fileSearched.split(localPrefix)[1]
                else:
                    print "please provide files with known prefixes (paths)"
                    if con:
                        con.rollback()
                        con.close()
                    sys.exit(1)
                
                cur.execute(cmd7%(fileSearched_WOPre))
                file_id_Searched = cur.fetchone()[0]
                
                for ii in range(len(seedPatternIds)):
                    colInd = 5*ii                    
                    for jj in range(start, end):                    
                        if searchMotifData[jj][colInd+4] != -1:
                            cur.execute(cmd1%(file_id_Searched, searchMotifData[jj][colInd+2], searchMotifData[jj][colInd+3], version))
                            cur.execute(cmd4)
                            pattern_id3 = cur.fetchone()[0]
                            cur.execute(cmd6%(seedPatternIds[ii], pattern_id3, searchMotifData[jj][colInd+4]))
                    
            con.commit()
    except psy.DatabaseError, e:
        print 'Error %s' % e
        if con:
            con.rollback()
            con.close()
        sys.exit(1)
                
        
    if con:
        con.close()

def fetchMBID(mp3File):
    try:
        mbid = easyid3.ID3(mp3File)['UFID:http://musicbrainz.org'].data
    except:
        print mp3File
        raise MBIDError('MBID not embedded')
    return mbid


def MBIDError(a):
    return -1
        

"""
CREATE TABLE file (
    id bigserial NOT NULL primary key,
    filename character varying(1000),
    mbid uuid not null
);

CREATE TABLE pattern (
   id bigserial not null primary key,
   file_id bigint not null references file(id),
   start_time double precision not null,
   end_time double precision not null,
   pair_id bigint references pattern(id),
   version int not null
);

CREATE TABLE match (
    id bigserial not null primary key,
    source_id bigint not null references pattern(id),
    target_id bigint not null references pattern(id),
    distance double precision NOT NULL
);

sudo -u postgres dropdb motif_local
sudo -u postgres createdb motif_local -O sankalp


"""