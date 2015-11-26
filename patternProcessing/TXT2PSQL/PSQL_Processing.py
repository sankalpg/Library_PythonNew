#!/usr/bin/python
# -*- coding: utf-8 -*-

import psycopg2 as psy
import sys, os
from mutagen import easyid3
import numpy as np
import time 
from compmusic import dunya as dn
from compmusic.dunya import hindustani as hn
from compmusic.dunya import carnatic as ca

sys.path.append(os.path.join(os.path.dirname(__file__), '../../batchProcessing'))


import batchProcessing as BP


myUser = 'sankalp'
#myDatabase = 'motifHindustani_CONF1'

serverPrefix = "/homedtic/sgulati/motifDiscovery/dataset/PatternProcessing_DB/unsupervisedDBs/carnaticDB/Carnatic40RagaICASSP2016/audio/"
localPrefix = "/media/Data/Datasets/PatternProcessing_DB/unsupervisedDBs/carnaticDB/Carnatic40RagaICASSP2016/audio/"


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

    except psy.DatabaseError, e:
        print 'Error %s' % e
        if con:
            con.rollback()
            con.close()
        sys.exit(1)
    
    if con:
        con.close()
        
def stripPrefix(audiofile):
    
    if audiofile.count(serverPrefix):
        audiofile_WOPre = audiofile.split(serverPrefix)[1]
    elif audiofile.count(localPrefix):
        audiofile_WOPre = audiofile.split(localPrefix)[1]
    else:
        print "please provide files with known prefixes (paths)"
        audiofile_WOPre = audiofile;
        
    return audiofile_WOPre


 #raagaBins = {'bf4662d4-25c3-4cad-9249-4de2dc513c06':0,
 #'a9413dff-91d1-4e29-ad92-c04019dce5b8':1,
 #'f0866e71-33b2-47db-ab75-0808c41f2401':2,
 #'123b09bd-9901-4e64-a65a-10b02c9e0597':3,
 #'d5285bf4-c3c5-454e-a659-fec30075990b':4,
 #'09c179f3-8b19-4792-a852-e9fa0090e409':5,
 #'700e1d92-7094-4c21-8a3b-d7b60c550edf':6,
 #'77bc3583-d798-4077-85a3-bd08bc177881':7,
 #'aa5f376f-06cd-4a69-9cc9-b1077104e0b0':8,
 #'85ccf631-4cdf-4f6c-a841-0edfcf4255d1':9}

def get_db_raga_dunya_raga(outputfile, myDatabase, collection = 'carnatic'):
	"""
	Due to changes in raga labels there are cases where the raga labels in the dunya are changed.This function just dumps a file 
	with mbid, raga in the curent db and in the dunya (from musicbrainz)
	"""
	dn.set_token("60312f59428916bb854adaa208f55eb35c3f2f07")
	if collection == 'hindustani':
		tradition = hn
	elif collection == 'carnatic':
		tradition = ca
		
	getMBIDs = "select mbid from file"
	cmd1 = "select raagaId from file where mbid='%s'"
		
	try:
		con = psy.connect(database=myDatabase, user=myUser) 
		cur = con.cursor()
		cur.execute(getMBIDs)
		mbids = cur.fetchall()
		
		for ii,mbid in enumerate(mbids):
			mbid = mbid[0]
			try:
				recData = tradition.get_recording(mbid)
				ragaId_dunya = recData['raaga'][0]['uuid']
				cur.execute(cmd1%(mbid))
				raaga_id = cur.fetchone()[0]
				print mbid, ragaId_dunya, raaga_id
			except:
				pass
	except psy.DatabaseError, e:
		print 'Error %s' % e	
		if con:
			con.rollback()
			con.close()
			sys.exit(1)

	if con:
		con.close()
	
def update_tonic_in_db(root_dir, tonicExt = '.tonicFine', myDatabase= '', myUser=''):	
    
    audioFiles = BP.GetFileNamesInDir(root_dir, '.mp3')
    cmd = u"update file set tonic = %f where mbid='%s'"
    try:
        con = psy.connect(database=myDatabase, user=myUser) 
        cur = con.cursor()        
        for audiofile in audioFiles:
            print audiofile
            mbid = fetchMBID(audiofile)
            fname, ext = os.path.splitext(audiofile)
            tonic = np.loadtxt(fname+tonicExt)
            
            cur.execute(cmd%(tonic, mbid))
    
        con.commit()
    except:
        pass
    
        

def fillFileTableWithRagaIds(myDatabase, collection = 'carnatic'):
    
    dn.set_token("60312f59428916bb854adaa208f55eb35c3f2f07")
    if collection == 'hindustani':
        tradition = hn
    elif collection == 'carnatic':
        tradition = ca
        
    getMBIDs = "select mbid from file"
    cmd1 = "update file set raagaId ='%s' where mbid='%s'"
        
    try:
        con = psy.connect(database=myDatabase, user=myUser) 
        cur = con.cursor()
        cur.execute(getMBIDs)
        mbids = cur.fetchall()
        
        for ii,mbid in enumerate(mbids):
            mbid = mbid[0]
            print ii
            try:
                recData = tradition.get_recording(mbid)
            except:
                print "Api didn't respond to this mbid = %s"%mbid
            
            if len(recData['raaga'])==1:
                ragaId = recData['raaga'][0]['uuid']
                cur.execute(cmd1%(ragaId, mbid))
            else:
                print "Can't fetch correct raga for mbid = %s"%mbid
        con.commit()
        print "Successfully updated file table in %s database"%(myDatabase)

    except psy.DatabaseError, e:
        print 'Error %s' % e    
        if con:
            con.rollback()
            con.close()
            sys.exit(1)

    if con:
        con.close()
        
        
def fillFileTableWithRagaIds_withoutAPI(myDatabase, mbid_ragaid_file, collection = 'carnatic'):
    
    mbid_raga = {}
    lines  = open(mbid_ragaid_file,'r').readlines()
    for line in lines:
        sline = line.split('\t')
        mbid = sline[0].strip()
        ragaid = sline[1].strip()
        mbid_raga[mbid] = ragaid
    dn.set_token("60312f59428916bb854adaa208f55eb35c3f2f07")
    if collection == 'hindustani':
        tradition = hn
    elif collection == 'carnatic':
        tradition = ca
        
    getMBIDs = "select mbid from file"
    cmd1 = "update file set raagaId ='%s' where mbid='%s'"
        
    try:
        con = psy.connect(database=myDatabase, user=myUser) 
        cur = con.cursor()
        cur.execute(getMBIDs)
        mbids = cur.fetchall()
        
        for ii,mbid in enumerate(mbids):
            mbid = mbid[0]
            if mbid_raga.has_key(mbid):
                cur.execute(cmd1%(mbid_raga[mbid], mbid))
            else:
                print "problem with mbid = %s"%(mbid)
            
        con.commit()
        print "Successfully updated file table in %s database"%(myDatabase)

    except psy.DatabaseError, e:
        print 'Error %s' % e    
        if con:
            con.rollback()
            con.close()
            sys.exit(1)

    if con:
        con.close()

def createFileTable(root_dir, myDatabase, filterExt = '.mp3'):
    
    con = None
    audioFiles = BP.GetFileNamesInDir(root_dir, filterExt)
    dataDump = []
    query = "INSERT INTO file (filename, mbid) VALUES (%s, %s)"
    cmd1 = "create index on file(filename)"
    cmd2 = "create index on file(mbid)"
    cmd3 = "create index on file(id)"
    
    
    for audiofile in audioFiles:
        filename, ext = os.path.splitext(audiofile)
        audiofile = filename + '.mp3'
        try:
            mbid = fetchMBID(audiofile)
        except:
            print "MBID not embedded in file %s"%audiofile
        else:
            audiofile_WOPre = stripPrefix(audiofile)
            dataDump.append((audiofile_WOPre, mbid))
        
    print "Dump array successfully created"
    try:
        con = psy.connect(database=myDatabase, user=myUser) 
        cur = con.cursor()
        cur.executemany(query, dataDump)
        cur.execute(cmd1)
        cur.execute(cmd2)
        cur.execute(cmd3)
        con.commit()
        print "Successfully updated file table in %s database"%(myDatabase)

    except psy.DatabaseError, e:
        print 'Error %s' % e
        if con:
            con.rollback()
            con.close()
        sys.exit(1)
    
    if con:
        con.close()

def addHasSeedValeInFileTable(myDatabase = ''):
    
    cmd1 = "drop index if exists file_id_idx"
    cmd2 = "drop index if exists file_filename_idx"
    cmd3 = "drop index if exists file_mbid_idx"
    
    cmd4 = "select id from file"
    cmd5 = "select * from pattern where file_id=%d and isseed=1"
    cmd6 = "update file  set hasseed =%d where id=%d"
    
    cmd7 = "create index on file(id)"
    cmd8 = "create index on file(filename)"
    cmd9 = "create index on file(mbid)"
    cmd10 = "create index on file(hasseed)"
    
    try:
        con = psy.connect(database=myDatabase, user=myUser) 
        cur = con.cursor()
        print "Successfully connected to %s database"%(myDatabase)
        
        cur.execute(cmd4)
        fileIds = cur.fetchall()
        cur.execute(cmd1)
        cur.execute(cmd2)
        cur.execute(cmd3)
        con.commit()
        
        for fileId in fileIds:
            fileId = fileId[0]
            cur.execute(cmd5%fileId)
            patterns = cur.fetchall()
            if len(patterns)>0:
                cur.execute(cmd6%(1,fileId))
            else:
                cur.execute(cmd6%(0,fileId))
        
        cur.execute(cmd7)
        cur.execute(cmd8)
        cur.execute(cmd9)
        cur.execute(cmd10)
        con.commit()
        
    except psy.DatabaseError, e:
        print 'Error %s' % e
        if con:
            con.close()
        sys.exit(1)
        
    if con:
        con.close()
        

def createPatternMatchTable(root_dir, logFile, myDatabase=''):
    
    t1 = time.time()
    #important stuff
    motifDiscExt = '.disPatt_2s'
    #motifSearchExt = ['.srhPatts1SqEuclidean', '.srhPatts1CityBlock', '.srhPatts1ShiftCityBlock', '.srhPatts1ShiftLinExp']
    motifSearchExt = ['.srhPatt_2sSqEuclidean']
    mainDistanceIndex = 0           # this is basically index of extions stored in motifSearchExt which is the main distance files
    motifSearchMappExt = '.srhMap_2s'
    
    #lets get all the files in the root_dir
    fileNames = BP.GetFileNamesInDir(root_dir, motifDiscExt)
    
    
    #commands for doing different tasks with psql
    
    cmd1 = "INSERT INTO pattern (file_id, start_time, end_time, isseed) VALUES (%ld, %f, %f, %d) RETURNING id"#storing seed motifs
    cmd6 = "INSERT INTO match (source_id, target_id, distance, version) VALUES (%ld, %ld, %f, %d)"
    
    #cmd3 = "INSERT INTO pattern (id, pair_id) VALUES (%ld, %ld)"#storing seed motifs
    #cmd4 = "SELECT currval('pattern_id_seq')"
    cmd5 = "UPDATE pattern SET pair_id = %ld WHERE id = %ld"
    
    
    cmd2 = "SELECT id FROM file WHERE mbid = '%s'"#searching file id
    cmd7 = "SELECT id FROM file WHERE filename = $$%s$$"#searching file id
    
    cmd8 = "create index on pattern(file_id)"
    cmd9 = "create index on pattern(pair_id)"
    cmd10 = "create index on pattern(isseed)"
    cmd11 = "create index on match(source_id)"
    cmd12 = "create index on match(target_id)"
    cmd13 = "create index on match(version)"
    
    #rety the connection
    con = None
    try:
        con = psy.connect(database=myDatabase, user=myUser) 
        cur = con.cursor()
        print "Successfully connected to %s database"%(myDatabase)
        
    except psy.DatabaseError, e:
        print 'Error %s' % e
        if con:
            con.close()
        sys.exit(1)

    try:
        for ii,filename in enumerate(fileNames):
            
            ####################### First inserting discovered seed patterns in the table ###################
            fname, ext = os.path.splitext(filename)
            
            print "Processing file %d out of %d\n"%(ii+1,len(fileNames))
            print "Processing %s\n"%fname
            
            motifSeedFile = fname + motifDiscExt
            audiofile = fname + '.mp3'
            
            audiofile = stripPrefix(audiofile)
            cur.execute(cmd7%(audiofile))
            file_id = cur.fetchone()[0]
            
            #reading discovered motifs file
            try:
                seedMotifData = np.loadtxt(motifSeedFile)
            except:
                print ("seed file couldnt be opened%s\n"%motifSeedFile)
                fid = open(logFile,'a')
                fid.write('###########BUG3###########\n')
                fid.write("seed file couldnt be opened%s\n"%motifSeedFile)
                fid.close()
                continue
            
            #array to store the pattern ids since we need this while storing searched patterns
            seedPatternIds = []
            for ii in range(0, seedMotifData.shape[0]): 
                
                if seedMotifData[ii][4]>9999999999999999:
                    break
                
                #entering in table pattern the first instance of seed pair
                cur.execute(cmd1%(file_id, seedMotifData[ii][0], seedMotifData[ii][1], 1))
                pattern_id1 = cur.fetchone()[0]
                seedPatternIds.append(pattern_id1)
                
                #entering in table pattern the second instance of seed pair
                cur.execute(cmd1%(file_id, seedMotifData[ii][2], seedMotifData[ii][3], 1))
                pattern_id2 = cur.fetchone()[0]
                seedPatternIds.append(pattern_id2)
                
                #cross referencing these two ids
                cur.execute(cmd5%(pattern_id2,pattern_id1))
                cur.execute(cmd5%(pattern_id1,pattern_id2))
                
                cur.execute(cmd6%(pattern_id1, pattern_id2, seedMotifData[ii][4], -1))
                cur.execute(cmd6%(pattern_id2, pattern_id1, seedMotifData[ii][4], -2))
                
            con.commit()
            
            ####################### Now inserting searched patterns in the table ###################
            ### NOTE that since searched patterns for a song remains same across different versions (only distance changes not the position) the differnce should only be in match table and not in pattenr table
            
            #managing the mapp file to get file ids and their corresponding names in the mapp dump
            motifMappFile = fname + motifSearchMappExt
            try:
                searchMappData = open(motifMappFile,'r').readlines()
            except:
                print ("Mapp file couldnt be opened%s\n"%motifMappFile)
                fid = open(logFile,'a')
                fid.write('###########BUG3###########\n')
                fid.write("Mapp file couldnt be opened%s\n"%motifMappFile)
                fid.close()
                continue
            
            searchedFileArray = np.array([])
            searchedFileIDArray = np.array([]).astype(np.int)
            
            for line in searchMappData:
                fid,fileSearched = line.split('\t')
                
                fileSearched = fileSearched.strip() + '.mp3'
                fileSearched = stripPrefix(fileSearched)
                searchedFileArray = np.append(searchedFileArray,fileSearched)
                # lets search in the database what is the id for this file
                cur.execute(cmd7%(fileSearched))
                file_id = cur.fetchone()[0]
                searchedFileIDArray = np.append(searchedFileIDArray, file_id)
            
            ### Insert first data for one file so that we know what patterns have what id in the "pattern" table and then for rest of the version we can directly get these ids for storing in match table
            
            #now inserting searched data
            SearchExt = motifSearchExt[mainDistanceIndex]
            version=1
            motifSearchFile = fname + SearchExt
            #reading searched motifs file
            try:
                searchMotifData = np.loadtxt(motifSearchFile)
            except:
                print ("Search file couldnt be opened%s\n"%motifSearchFile)
                fid = open(logFile,'a')
                fid.write('###########BUG3###########\n')
                fid.write("Search file couldnt be opened%s\n"%motifSearchFile)
                fid.close()
                continue
            
            patternIdMTX = np.zeros((searchMotifData.shape[0], len(seedPatternIds))).astype(np.int64)
            
            #start inserting information about the searched patterns.
            for ii in range(len(seedPatternIds)):
                colOffset = 6*ii
                for rowInd in np.arange(searchMotifData.shape[0]):
                    cur.execute(cmd1%(searchedFileIDArray[searchMotifData[rowInd][colOffset+5]], searchMotifData[rowInd][colOffset+2], searchMotifData[rowInd][colOffset+3], 0))
                    patternID = cur.fetchone()[0]
                    cur.execute(cmd6%(seedPatternIds[ii], patternID, searchMotifData[rowInd][colOffset+4], version))
                    patternIdMTX[rowInd,ii] = patternID
                    
            
            #now inserting searched data
            for version,SearchExt in enumerate(motifSearchExt):
                
                if SearchExt == motifSearchExt[mainDistanceIndex]:
                    continue
                motifSearchFile = fname + SearchExt
                #reading searched motifs file
                try:
                    searchMotifDataNEW = np.loadtxt(motifSearchFile)
                except:
                    print ("Search file couldnt be opened%s\n"%searchMotifDataNEW)
                    fid = open(logFile,'a')
                    fid.write('###########BUG3###########\n')
                    fid.write("Search file couldnt be opened%s\n"%searchMotifDataNEW)
                    fid.close()
                    continue
                
                #start inserting information about the searched patterns.
                for ii in range(len(seedPatternIds)):
                    colOffset = 6*ii
                    for rowInd in np.arange(searchMotifDataNEW.shape[0]):
                        
                        indSameFiles = np.where(searchMotifData[:,colOffset+5]==searchMotifDataNEW[rowInd,colOffset+5])[0]
                        indSameStart = np.where(searchMotifData[indSameFiles,colOffset+2]==searchMotifDataNEW[rowInd,colOffset+2])[0]
                        indMatch = indSameFiles[indSameStart]
                        if indMatch.size>1:
                            fid = open(logFile,'a')
                            fid.write('###########BUG1###########\n')
                            fid.write('Seed file name %s\n'%motifSeedFile)
                            fid.write('Search file name %s\n'%motifSearchFile)
                            fid.write('Seed pattern Index %d\n'%ii)
                            fid.write('Search pattern index (row)%d\n'%rowInd)
                            fid.close()
                            
                            indSameSeedEnd = np.where(searchMotifData[indSameFiles,colOffset+1]==searchMotifDataNEW[rowInd,colOffset+1])[0]
                            indSameStartEnd = np.intersect1d(indSameStart, indSameSeedEnd)
                            indMatch = indSameFiles[indSameStartEnd]
                            if indMatch.size>1:
                                fid = open(logFile,'a')
                                fid.write('###########BUG2###########\n')
                                fid.write('Seed file name %s\n'%motifSeedFile)
                                fid.write('Search file name %s\n'%motifSearchFile)
                                fid.write('Seed pattern Index %d\n'%ii)
                                fid.write('Search pattern index (row)%d\n'%rowInd)
                                fid.close()
                                indMatch=indMatch[0]
                                print "THERE IS A BIG PROBLEM HERE, YOUR ASSUMPTION IS WRONG SIR"
                        cur.execute(cmd6%(seedPatternIds[ii], patternIdMTX[indMatch,ii], searchMotifDataNEW[rowInd][colOffset+4], version))
            con.commit()
            
        cur.execute(cmd8)
        cur.execute(cmd9)
        cur.execute(cmd10)
        cur.execute(cmd11)
        cur.execute(cmd12)
        cur.execute(cmd13)
        con.commit()
            
    except psy.DatabaseError, e:
        print 'Error %s' % e
        if con:
            con.rollback()
            con.close()
        sys.exit(1)
                
        
    if con:
        con.close()
    
    print "Everything went fine, enjoy!! Time taken for inserting all results into DB %f\n"%(time.time()-t1)


def upload_network_data(fileListFile, pattDistExt = '', myDatabase = '', myUser = ''):
    """
    This function uploads the distances obtained from all the non-overlapping phrases to all the 
    other such phrases. This is the final data obtained for building the network.
    
    """
    cmd = "INSERT INTO network (source_id, target_id, distance) VALUES (%s, %s, %s)"
    
    try:
        con = psy.connect(database=myDatabase, user=myUser) 
        cur = con.cursor()
        print "Successfully connected to %s database"%(myDatabase)
        
        filenames = open(fileListFile).readlines()
        for ii, filename in enumerate(filenames):
            print "uploading phrases for file number %d, name: %s"%(ii, filename)
            fname = filename.strip() + pattDistExt
            data = np.loadtxt(fname)
            dump = []
            #for ii in range(data.shape[0]):
                #dump.append((int(data[ii][0]), int(data[ii][1]), float(data[ii][2])))
                
            cur.executemany(cmd, data)
            con.commit()
        
    except psy.DatabaseError, e:
        print 'Error %s' % e
        if con:
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
    mbid uuid not null,
    hasseed int,
    raagaid uuid,
    tonic double precision
);

CREATE TABLE pattern (
   id bigserial not null primary key,
   file_id bigint not null references file(id),
   start_time double precision not null,
   end_time double precision not null,
   pair_id bigint references pattern(id),
   isSeed int not null
);

CREATE TABLE match (
    id bigserial not null primary key,
    source_id bigint not null references pattern(id),
    target_id bigint not null references pattern(id),
    distance double precision NOT NULL,
    version int not null
);

CREATE TABLE network (
    source_id bigint not null references pattern(id),
    target_id bigint not null references pattern(id),
    distance double precision NOT NULL
);


sudo -u postgres dropdb motif_local
sudo -u postgres createdb motif_local -O sankalp
create index on file(filename);

select last_value from match_id_seq; ALTERNATIVE OF COUNT!!!

alter table file add hasseed int;





Useful tips
1) Use truncate instead of delete
2) Use escape function to handle escape characters
3) Use return id statement instead of querrying for currval
4) index columns once the dataset is frozen. Beware to remove indexing when you want to add new entries
5) Also we can use tables without references and we can update references after the entries are done
6) ALWAYS remember to remove indexes before adding new entries



"""
