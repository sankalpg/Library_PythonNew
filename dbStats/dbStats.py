from __future__ import unicode_literals
import codecs
import json, os, sys
import numpy as np
import compmusic
import json
from compmusic import dunya as dn
from compmusic.dunya import hindustani as hi
from compmusic.dunya import carnatic as ca
from compmusic.dunya import docserver as ds
from compmusic import musicbrainz
import codecs

sys.path.append('/home/sankalp/Work/CompanyWork/camut/saragautils')
import metadata as meta

import psycopg2 as psy	


def getDatasetStats(mbids, output_file, music_tradition = ''):
	"""
	This function obtains a set of statistics/numbers on the dataset. This is needed to report in papers/thesis. Basically number of unique artists, duration, releases, recordings, raga, tala etc etc

	mbids: list of mbids
	music_tradition: 'hindustani' or 'carnatic'
	output_file: json file path to write. This file will contain the set of extracted number fields
	"""
	entities = []
	failure = 0
	if music_tradition == 'hindustani':
		tradition = hi
		concert = 'release'  # in hindustani album level items are referred by 'release'
		work = 'works'
		raga = 'raags'
		tala = 'taals'
		form = 'forms'
		laya = 'layas'
		lead_artists = 'album_artists'
		artists = 'artists'
		entities = [concert, work, raga, tala, form, laya, artists, lead_artists, 'length']
		object_to_fetch = {concert: 'mbid', work: 'mbid', raga: 'uuid', tala: 'uuid', form: 'name', laya: 'uuid', lead_artists: 'mbid', artists: 'mbid'}

	elif music_tradition == 'carnatic':
		tradition = ca
		concert = 'concert'  # in carnatic album level items are referred by 'concerts'
		work = 'work'
		raga = 'raaga'
		tala = 'taala'
		form = 'form'
		laya = 'laya'        
		lead_artists = 'album_artists'
		artists = 'artists'
		entities = [concert, work, raga, tala, form, artists, lead_artists, 'length']
		object_to_fetch = {concert: 'mbid', work: 'mbid', raga: 'uuid', tala: 'uuid', form: 'name', lead_artists: 'mbid', artists: 'mbid'}
	else:
	    print "Please specify a valid music tradition"


	stats = {}
	for e in entities:
		stats[e] = []


	for mbid in mbids:
		try:
			rec_info = tradition.get_recording(mbid)
		except:
			failure+=1
			print "Failed to fetch info for file %s"%mbid
		for e in entities:
			if rec_info.has_key(e):
				#special case for parsing artist field
				if e == 'artists':
					rec_info[e] = [a['artist'] for a in rec_info[e]]
				if isinstance(rec_info[e], int):
					stats[e].append(rec_info[e])
				elif isinstance(rec_info[e], list):
					temp = []
					for item in rec_info[e]:
						temp.append(item[object_to_fetch[e]])
					if len(temp)>0:
						stats[e].append(temp)
	
	output_stats = {}
	for e in entities:
		if e == 'length':
			output_stats[e] = {'total_length': np.sum(stats[e]), 'total_recs': len(stats[e])}
		else:
			output_stats[e] = {'total_unique': len(np.unique(sum(stats[e], []))), 'unique_elems': np.unique(sum(stats[e], [])).tolist(), 'total_rels': len(sum(stats[e], [])), 'total_recs': len(stats[e])}
	json.dump(output_stats, open(output_file, 'w'))



def getStatsDunyaCorpus(collectionId, output_file, music_tradition = ''):
	"""
	This function will fetch all the stats/numbers for a given collection/corpus/set using function 'getDatasetStats'
	collectionId: list of mbids of the collections
	music_tradition: 'hindustani' or 'carnatic'
	output_file: json output file to dump the stats
	"""

	if music_tradition == 'hindustani':
		tradition = hi
	elif music_tradition == 'carnatic':
		tradition = ca        
	else:
		print "Please specify a valid music tradition"

	dn.set_token("60312f59428916bb854adaa208f55eb35c3f2f07")
	tradition.set_collections(collectionId)
	recs = tradition.get_recordings()
	mbids = [r['mbid'] for r in recs]
	getDatasetStats(mbids, output_file, music_tradition)



def getMBIDsInCollections(collection_id, music_tradition, output_file):

	if music_tradition == 'hindustani':
		tradition = hi
	elif music_tradition == 'carnatic':
		tradition = ca        
	else:
		print "Please specify a valid music tradition"

	dn.set_token("60312f59428916bb854adaa208f55eb35c3f2f07")
	tradition.set_collections([collection_id])
	recs = tradition.get_recordings()
	mbids = [r['mbid'] for r in recs]
	json.dump(mbids, open(output_file, 'w'))
	return len(mbids)


def getMBIDsInCollectionsBatchRun():
	"""
	Enough of writing generic functions, 
	"""
	## Dumping data for carnatic music corpus
	collection_id = 'f96e7215-b2bd-4962-b8c9-2b40c17a1ec6'
	music_tradition = 'carnatic'
	output_file = '/home/sankalp/Work/Work_PhD/library_pythonnew/dbStats/DunyaCorpusCarnatic/CarnaticCorpus_MBIDS.json'
	n_recs = getMBIDsInCollections(collection_id, music_tradition, output_file)
	print "For %s collection we got %d number of recordings"%("Carnatic corpus", n_recs)


	## Dumping data for carnatic bootleg music corpus
	collection_id = 'f8bf7d1e-70d2-44f6-a3cb-5a6ded00be1f'
	music_tradition = 'carnatic'
	output_file = '/home/sankalp/Work/Work_PhD/library_pythonnew/dbStats/DunyaCorpusBootledCarnatic/CarnaticCorpusBootleg_MBIDS.json'
	n_recs = getMBIDsInCollections(collection_id, music_tradition, output_file)
	print "For %s collection we got %d number of recordings"%("Carnatic bootleg corpus", n_recs)

	## Dumping data for carnatic music corpus CC
	collection_id = 'a163c8f2-b75f-4655-86be-1504ea2944c2'
	music_tradition = 'carnatic'
	output_file = '/home/sankalp/Work/Work_PhD/library_pythonnew/dbStats/DunyaCorpusCC/Carnatic/CarnaticCorpusCC_MBIDS.json'	
	n_recs = getMBIDsInCollections(collection_id, music_tradition, output_file)
	print "For %s collection we got %d number of recordings"%("Carnatic CC corpus", n_recs)


	## Dumping data for hindustani music corpus
	collection_id = '213347a9-e786-4297-8551-d61788c85c80'
	music_tradition = 'hindustani'
	output_file = '/home/sankalp/Work/Work_PhD/library_pythonnew/dbStats/DunyaCorpusHindustani/HindustaniCorpus_MBIDS.json'
	n_recs = getMBIDsInCollections(collection_id, music_tradition, output_file)
	print "For %s collection we got %d number of recordings"%("Hindustani corpus", n_recs)

	## Dumping data for hindustani music corpus CC
	collection_id = '6adc54c6-6605-4e57-8230-b85f1de5be2b'
	music_tradition = 'hindustani'
	output_file = '/home/sankalp/Work/Work_PhD/library_pythonnew/dbStats/DunyaCorpusCC/Hindustani/HindustaniCorpusCC_MBIDS.json'
	n_recs = getMBIDsInCollections(collection_id, music_tradition, output_file)
	print "For %s collection we got %d number of recordings"%("Hindustani CC corpus", n_recs)



def getMBIDSInDatasets():

	#gathering mbids of carnatic raga recognition datasrt
	database = 'Raga_Rec_Carnatic_40Raga_Config0'
	user = 'sankalp'
	cmd = 'select mbid from file'
	output_file = '/home/sankalp/Work/Work_PhD/library_pythonnew/dbStats/RagaRecognitionDS/carnatic/CarnaticRagaRec40_MBIDS.json'
	try:
		con = psy.connect(database=database, user=user) 
		cur = con.cursor()
		print "Successfully connected to the server"
		cur.execute(cmd)
		results = cur.fetchall()

	except psy.DatabaseError, e:
		print 'Error %s' % e
		if con:
		    con.rollback()
		    con.close()
		sys.exit(1)

	if con:
		con.close()

	mbids = [res[0] for res in results]
	json.dump(mbids, open(output_file, 'w'))

	database = 'Raga_Rec_Hindustani_30Raga_Config0'
	output_file = '/home/sankalp/Work/Work_PhD/library_pythonnew/dbStats/RagaRecognitionDS/hindustani/HindustaniRagaRec30_MBIDS.json'
	try:
		con = psy.connect(database=database, user=user) 
		cur = con.cursor()
		print "Successfully connected to the server"
		cur.execute(cmd)
		results = cur.fetchall()

	except psy.DatabaseError, e:
		print 'Error %s' % e
		if con:
		    con.rollback()
		    con.close()
		sys.exit(1)

	if con:
		con.close()

	mbids = [res[0] for res in results]
	json.dump(mbids, open(output_file, 'w'))


def getStatsAllCollections():

	# for corpus first 

	mbids_file = '/home/sankalp/Work/Work_PhD/library_pythonnew/dbStats/DunyaCorpusCarnatic/CarnaticCorpus_MBIDS.json'
	music_tradition = 'carnatic'
	output_file = '/home/sankalp/Work/Work_PhD/library_pythonnew/dbStats/DunyaCorpusCarnatic/CarnaticCorpus_Stats.json'
	mbids = json.load(open(mbids_file, 'r'))
	getDatasetStats(mbids, output_file, music_tradition)

	mbids_file = '/home/sankalp/Work/Work_PhD/library_pythonnew/dbStats/DunyaCorpusBootledCarnatic/CarnaticCorpusBootleg_MBIDS.json'
	music_tradition = 'carnatic'
	output_file = '/home/sankalp/Work/Work_PhD/library_pythonnew/dbStats/DunyaCorpusBootledCarnatic/CarnaticCorpusBootleg_Stats.json'
	mbids = json.load(open(mbids_file, 'r'))
	getDatasetStats(mbids, output_file, music_tradition)


	mbids_file = '/home/sankalp/Work/Work_PhD/library_pythonnew/dbStats/DunyaCorpusCC/Carnatic/CarnaticCorpusCC_MBIDS.json'
	music_tradition = 'carnatic'
	output_file = '/home/sankalp/Work/Work_PhD/library_pythonnew/dbStats/DunyaCorpusCC/Carnatic/CarnaticCorpusCC_Stats.json'
	mbids = json.load(open(mbids_file, 'r'))
	getDatasetStats(mbids, output_file, music_tradition)


	mbids_file = '/home/sankalp/Work/Work_PhD/library_pythonnew/dbStats/DunyaCorpusHindustani/HindustaniCorpus_MBIDS.json'
	music_tradition = 'hindustani'
	output_file = '/home/sankalp/Work/Work_PhD/library_pythonnew/dbStats/DunyaCorpusHindustani/HindustaniCorpus_Stats.json'
	mbids = json.load(open(mbids_file, 'r'))
	getDatasetStats(mbids, output_file, music_tradition)

	mbids_file = '/home/sankalp/Work/Work_PhD/library_pythonnew/dbStats/DunyaCorpusCC/Hindustani/HindustaniCorpusCC_MBIDS.json'
	music_tradition = 'hindustani'
	output_file = '/home/sankalp/Work/Work_PhD/library_pythonnew/dbStats/DunyaCorpusCC/Hindustani/HindustaniCorpusCC_Stats.json'
	mbids = json.load(open(mbids_file, 'r'))
	getDatasetStats(mbids, output_file, music_tradition)


	# for datasets now
	mbids_file = '/home/sankalp/Work/Work_PhD/library_pythonnew/dbStats/RagaRecognitionDS/carnatic/CarnaticRagaRec40_MBIDS.json'
	music_tradition = 'carnatic'
	output_file = '/home/sankalp/Work/Work_PhD/library_pythonnew/dbStats/RagaRecognitionDS/carnatic/CarnaticRagaRec40_Stats.json'
	mbids = json.load(open(mbids_file, 'r'))
	getDatasetStats(mbids, output_file, music_tradition)

	mbids_file = '/home/sankalp/Work/Work_PhD/library_pythonnew/dbStats/RagaRecognitionDS/hindustani/HindustaniRagaRec30_MBIDS.json'
	music_tradition = 'hindustani'
	output_file = '/home/sankalp/Work/Work_PhD/library_pythonnew/dbStats/RagaRecognitionDS/hindustani/HindustaniRagaRec30_Stats.json'
	mbids = json.load(open(mbids_file, 'r'))
	getDatasetStats(mbids, output_file, music_tradition)	




def getStatsDunyaCarnaticCorpusCC():
	collectionId = ['a163c8f2-b75f-4655-86be-1504ea2944c2']
	music_tradition = 'carnatic'
	output_file = '/home/sankalp/Work/Work_PhD/library_pythonnew/dbStats/DunyaCorpusCC/Carnatic/stats.json'
	getStatsDunyaCorpus(collectionId, output_file, music_tradition)


def getStatsDunyaHindustaniCorpusCC():
	collectionId = ['6adc54c6-6605-4e57-8230-b85f1de5be2b']
	music_tradition = 'hindustani'
	output_file = '/home/sankalp/Work/Work_PhD/library_pythonnew/dbStats/DunyaCorpusCC/Hindustani/stats.json'
	getStatsDunyaCorpus(collectionId, output_file, music_tradition)

def generatePretyReport(stats_file, out_file):

	fid = open(out_file, 'w')

	data = json.load(open(stats_file, 'r'))
	for key1 in data.keys():
		fid.write('------------ %s ------------\n'%str(key1))
		if key1 == 'length':
			for key2 in data[key1].keys():
				fid.write('%s\t%f\n'%(str(key2), float(data[key1][key2])/(1000.0*3600.0)))
		else:
			for key2 in data[key1].keys():
				if key2 == 'unique_elems':
					fid.write('%s\t%d\n'%(str(key2), len(data[key1][key2])))
				else:
					fid.write('%s\t%d\n'%(str(key2), data[key1][key2]))
		fid.write('\n')
	fid.close()

def generatePrettyReportAllStats():

	input_file = '/home/sankalp/Work/Work_PhD/library_pythonnew/dbStats/DunyaCorpusCarnatic/CarnaticCorpus_Stats.json'
	output_file = '/home/sankalp/Work/Work_PhD/library_pythonnew/dbStats/DunyaCorpusCarnatic/CarnaticCorpus_StatsPretty.txt'
	generatePretyReport(input_file, output_file)
	
	input_file = '/home/sankalp/Work/Work_PhD/library_pythonnew/dbStats/DunyaCorpusBootledCarnatic/CarnaticCorpusBootleg_Stats.json'
	output_file = '/home/sankalp/Work/Work_PhD/library_pythonnew/dbStats/DunyaCorpusBootledCarnatic/CarnaticCorpusBootleg_StatsPretty.txt'
	generatePretyReport(input_file, output_file)

	input_file = '/home/sankalp/Work/Work_PhD/library_pythonnew/dbStats/DunyaCorpusCC/Carnatic/CarnaticCorpusCC_Stats.json'
	output_file = '/home/sankalp/Work/Work_PhD/library_pythonnew/dbStats/DunyaCorpusCC/Carnatic/CarnaticCorpusCC_StatsPretty.txt'
	generatePretyReport(input_file, output_file)

	input_file = '/home/sankalp/Work/Work_PhD/library_pythonnew/dbStats/DunyaCorpusHindustani/HindustaniCorpus_Stats.json'
	output_file = '/home/sankalp/Work/Work_PhD/library_pythonnew/dbStats/DunyaCorpusHindustani/HindustaniCorpus_StatsPretty.txt'
	generatePretyReport(input_file, output_file)

	input_file = '/home/sankalp/Work/Work_PhD/library_pythonnew/dbStats/DunyaCorpusCC/Hindustani/HindustaniCorpusCC_Stats.json'
	output_file = '/home/sankalp/Work/Work_PhD/library_pythonnew/dbStats/DunyaCorpusCC/Hindustani/HindustaniCorpusCC_StatsPretty.txt'
	generatePretyReport(input_file, output_file)

	input_file = '/home/sankalp/Work/Work_PhD/library_pythonnew/dbStats/RagaRecognitionDS/carnatic/CarnaticRagaRec40_Stats.json'
	output_file = '/home/sankalp/Work/Work_PhD/library_pythonnew/dbStats/RagaRecognitionDS/carnatic/CarnaticRagaRec40_StatsPretty.txt'
	generatePretyReport(input_file, output_file)

	input_file = '/home/sankalp/Work/Work_PhD/library_pythonnew/dbStats/RagaRecognitionDS/hindustani/HindustaniRagaRec30_Stats.json'
	output_file = '/home/sankalp/Work/Work_PhD/library_pythonnew/dbStats/RagaRecognitionDS/hindustani/HindustaniRagaRec30_StatsPretty.txt'
	generatePretyReport(input_file, output_file)



def getAnotsLenPerFile(filename):

	lines = open(filename, 'r')
	anots = 0
	for line in lines:
		if len(line.strip())>0:
			anots+=1
	return anots

def dumpSaragaAnotStats(base_dir, mbids_file, output_file, tradition, token='60312f59428916bb854adaa208f55eb35c3f2f07', exts = []):
	
	mbids = json.load(open(mbids_file, 'r'))

	stats = {}
	for k in exts:
		stats[k] = {}
		stats[k]['files'] = 0
		stats[k]['lines'] = 0

	for mbid in mbids:
		concert_dir, con_name, recording_dir, rec_name = meta.get_recording_path(mbid, tradition, token)
		for ext in exts:
			filename = os.path.join(base_dir, recording_dir, rec_name + ext)
			if os.path.isfile(filename):
				stats[ext]['files'] +=1
				if ext == '.pitch' or ext == '.mpitch':
					pass
				else:
					stats[ext]['lines'] += getAnotsLenPerFile(filename)
	
	json.dump(stats, open(output_file, 'w'))


def getSaragaAnnotationStats():

	input_file = '/home/sankalp/Work/Work_PhD/library_pythonnew/dbStats/DunyaCorpusCC/Carnatic/CarnaticCorpusCC_MBIDS.json'
	output_file = '/home/sankalp/Work/Work_PhD/library_pythonnew/dbStats/SaragaDBStats/Saraga_Carnatic_Anots.txt'
	base_dir = '/media/Data/Dropbox/CAMUT/SaragaV1/Carnatic/Concerts'
	tradition = 'carnatic'
	exts = ['.sama', '.sections_p', '.mphrases', '.pitch', '.tonic', '.mpitch', '.bpm'] 
	dumpSaragaAnotStats(base_dir, input_file, output_file, tradition, token='60312f59428916bb854adaa208f55eb35c3f2f07', exts = exts)



	input_file = '/home/sankalp/Work/Work_PhD/library_pythonnew/dbStats/DunyaCorpusCC/Hindustani/HindustaniCorpusCC_MBIDS.json'
	output_file = '/home/sankalp/Work/Work_PhD/library_pythonnew/dbStats/SaragaDBStats/Saraga_Hindustani_Anots.txt'
	base_dir = '/media/Data/Dropbox/CAMUT/SaragaV1/Hindustani/Concerts'
	tradition = 'hindustani'
	exts = ['.sama', '.sections_p', '.mphrases', '.pitch', '.tonic', '.mpitch', '.bpm'] 
	dumpSaragaAnotStats(base_dir, input_file, output_file, tradition, token='60312f59428916bb854adaa208f55eb35c3f2f07', exts = exts)


def getRagaWiseStats(raga_mbid_file, output_dir, music_tradition):

	#Hindustani: dbs.getRagaWiseStats('RagaRecognitionDS/hindustani/PerRagaStats/Hindustani30Raga300_FILE_MBID_RAGA.txt', 'RagaRecognitionDS/hindustani/PerRagaStats/', 'hindustani')
	#Carnatic: dbs.getRagaWiseStats('RagaRecognitionDS/carnatic/PerRagaStats/Carnatic40Raga480_FILE_MBID_RAGA.txt', 'RagaRecognitionDS/carnatic/PerRagaStats/', 'carnatic')

	lines = open(raga_mbid_file, 'r').readlines()
	mbid_ragaid = []

	dn.set_token("60312f59428916bb854adaa208f55eb35c3f2f07")

	for line in lines:
		lsplit = line.split('\t')
		mbid_ragaid.append([lsplit[1].strip(), lsplit[2].strip()])

	mbid_ragaid = np.array(mbid_ragaid)
	mbids = mbid_ragaid[:,0]
	ragaids = mbid_ragaid[:,1]

	uragaids = np.unique(ragaids)

	for ragaid in uragaids:
		ind_raga = np.where(ragaids == ragaid)[0]
		mbid_selected = mbids[ind_raga]
		output_file = os.path.join(output_dir, ragaid + '.json')
		getDatasetStats(mbid_selected, output_file , music_tradition = music_tradition)


def generatePerRagaPrettyReport(root_dir, raga_mbid_file, raga_map, output_file):

	#Hindustani: dbs.generatePerRagaPrettyReport('RagaRecognitionDS/hindustani/PerRagaStats/', 'RagaRecognitionDS/hindustani/PerRagaStats/Hindustani30Raga300_FILE_MBID_RAGA.txt', 'raga_name_mapping.json', 'RagaRecognitionDS/hindustani/PerRagaStats/Hindustani30Raga300_PerRagaStats.txt')
	#Carnatic: dbs.generatePerRagaPrettyReport('RagaRecognitionDS/carnatic/PerRagaStats/', 'RagaRecognitionDS/carnatic/PerRagaStats/Carnatic40Raga480_FILE_MBID_RAGA.txt', 'raga_name_mapping.json', 'RagaRecognitionDS/carnatic/PerRagaStats/Carnatic40Raga480_PerRagaStats.txt')

	lines = open(raga_mbid_file, 'r').readlines()
	mbid_ragaid = []

	mapping = json.load(open(raga_map, 'r'))

	for line in lines:
		lsplit = line.split('\t')
		mbid_ragaid.append([lsplit[1].strip(), lsplit[2].strip()])

	mbid_ragaid = np.array(mbid_ragaid)
	mbids = mbid_ragaid[:,0]
	ragaids = mbid_ragaid[:,1]

	uragaids = np.unique(ragaids)

	fid = codecs.open(output_file, 'w', encoding = 'utf8')
	for ragaid in uragaids:
		ind_raga = np.where(ragaids == ragaid)[0]
		mbid_selected = mbids[ind_raga]
		stat_file = os.path.join(root_dir, ragaid + '.json')
		stats = json.load(open(stat_file,'r'))
		try:
			fid.write("%s\t%s\t%0.2f\t%d\t%d\n"%(mapping[ragaid], ragaid, stats['length']['total_length']/(1000.0*3600.0),  stats['album_artists']['total_unique'], stats['release']['total_unique']))
		except:
			fid.write("%s\t%s\t%0.2f\t%d\t%d\n"%(mapping[ragaid], ragaid, stats['length']['total_length']/(1000.0*3600.0),  stats['album_artists']['total_unique'], stats['concert']['total_unique']))

	fid.close()