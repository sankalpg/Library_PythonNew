import json, os, sys
import numpy as np
import compmusic
import json
from compmusic import dunya as dn
from compmusic.dunya import hindustani as hn
from compmusic.dunya import carnatic as ca
from compmusic.dunya import docserver as ds
from compmusic import musicbrainz
import fixpath

ISMIR2015_10RAGASCarnatic = ['55', '10', '27', '8', '9', '137', '159', '20', '13', '210']

dn.set_token("60312f59428916bb854adaa208f55eb35c3f2f07")

####
#TODO This is a temporary code becasue audio download was not functioning properly
lines = open('CarnaticInfo/carnaticMBIDLocationMapping_3_March_2015.txt','r').readlines()
location = {}
for line in lines:
	splitLine = line.split()
	location[splitLine[1].strip()] = splitLine[0].strip()
	
####

def getRagaNameMapping(outFile, ragas = ISMIR2015_10RAGASCarnatic, collection = 'carnatic'):
    
    if collection == 'hindustani':
        tradition = hn
    elif collection == 'carnatic':
        tradition = ca
    fid = open(outFile, 'w')
    for raga in ragas:
        ragaName = tradition.get_raaga(raga)['name']
        fid.write("%s\t%s\n"%(raga,unicode(ragaName).encode('utf8')))
    fid.close()
        

def getRagasWithRecordings(outputFile, collection = 'hindustani', with_bootlegs=False):
	ragas= {}
	if collection == 'hindustani':
		tradition = hn
	elif collection == 'carnatic':
		tradition = ca
	#fetch all the recordings in the collection
	recordings = tradition.get_recordings(with_bootlegs=with_bootlegs)
	for ii,rec in enumerate(recordings):
		if ii%100.0==0:
			print ii, rec['mbid']
		try:
			recData=None
			recData = tradition.get_recording(rec['mbid'])
		except:
			print "ERROR Encountered for file %s"%rec['mbid']
			continue
		#print recData
		if type(recData['raaga']) == dict:
			ragaId = recData['raaga']['id']
			if not ragas.has_key(ragaId):
				ragas[ragaId]={}
				ragas[ragaId]['name'] = tradition.get_raaga(ragaId)['name']
				ragas[ragaId]['recs'] = []
			recData['concert_artist'] = tradition.get_concert(recData['concert']['mbid'])['concert_artists'][0]
			ragas[ragaId]['recs'].append(recData)
	json.dump(ragas,open(outputFile,'w'))
		
def countAndSortRagas(ragafile):
	ragas = json.load(open(ragafile,'r'))
	count = []
	ragaIds = []
	for r in ragas.keys():
		count.append(len(ragas[r]['recs']))
		ragaIds.append(r)
		#print r, len(ragas[r])
	
	sortInd = np.argsort(count)
	for ind in sortInd:
		print ragaIds[ind], count[ind]
		print ca.get_raaga(ragaIds[ind])['name']
		

def getAristMBIDPerRec(recObj):
	return recObj['concert_artist']['mbid']
def getDefaultIDPerRec(recObj):
	return '0000-0000-0000-0000'

def getDutationForAllRecordings(outputDurationFile, collection='carnatic', with_bootlegs = False):
	
	musicbrainz.mb.auth('compmusic','compmusic321')
	
	durationData = {}
	
	if collection == 'hindustani':
		tradition = hn
	elif collection == 'carnatic':
		tradition = ca
		
	recordings = tradition.get_recordings(with_bootlegs=with_bootlegs)
	
	for ii,rec in enumerate(recordings):
		if ii%100.0==0:
			print ii, rec['mbid']
		try:
			duration = int(musicbrainz.mb.get_recording_by_id(rec['mbid'])['recording']['length'][:-3])
			durationData[rec['mbid']]=duration
		except:
			print rec['mbid']
	
	json.dump(durationData, open(outputDurationFile,'w'))
	
	
	
	

def generateRagadataset(allRagaRecordingfile, durationFile, selectedRecordingOutputFile, outputTXTFILE, raagas=[], samplingMethod= 'maxArtist', filterBy = 'year', sortBy = 'duration', collection ='carnatic'):
	"""
	This function selects recordings of each raga according to specified criterions to build a raga dataset. Since selection of raagas for the dataset also involves consultation with musician, in addition to just looking
	at number of recordings per raga, the selected ragas are given as input after consultation with a musician.
	
	samplingMethod is basically how we would like to sample the recordings so that we have same number of selected reocrdings in each raga. The final number of recordings per raga is basically the minimum number of 
	recordings in the provided list of ragas.
	
	TODO: Additionally one can filter the recordings based on critrion. For example he can specify year threshold and can filter recordings < that year etc.
	
	"""
	
	if collection == 'hindustani':
		tradition = hn
	elif collection == 'carnatic':
		tradition = ca
			
	musicbrainz.mb.auth('compmusic','compmusic321')
	
	ragaData = json.load(open(allRagaRecordingfile,'r'))
	
	durationData = json.load(open(durationFile,'r'))
	
	#Computing what is the minimum number of recordings per raga
	minRecording = 1000000000000000
	for raga in raagas:
		nRecs = len(ragaData[raga]['recs'])
		print "Processing raga %s which has %d number of recordings"%(ragaData[raga]['name'], nRecs)
		if minRecording > nRecs:
			minRecording = nRecs
	
	print "Number of recordings per raga to be selected for the dataset is %d"%minRecording
	
	# applying sampling methods
	if samplingMethod == 'maxArtist':
		samplingMethod = getAristMBIDPerRec
	else:
		samplingMethod = getDefaultIDPerRec
		
	data={}
	dataSelected={}
	for raga in raagas:
		defaultDuration = 0
		if not data.has_key(raga):
			data[raga]={}
		if not dataSelected.has_key(raga):
			dataSelected[raga]={}
		
		for rec in ragaData[raga]['recs']:
			
			#Here comes the filtering stage
			
			#Applying sampling method
			samplingKey = samplingMethod(rec)
			try:
				if sortBy =='duration':
					#duration = int(musicbrainz.mb.get_recording_by_id(rec['mbid'])['recording']['length'][:-3])
					duration = durationData[rec['mbid']]
				else:
					duration = defaultDuration
					defaultDuration+=0.1
			except:
				duration = defaultDuration
				defaultDuration+=0.1
			if not data[raga].has_key(samplingKey):
				data[raga][samplingKey] = []
			
			data[raga][samplingKey].append({duration:rec['mbid']})
		
		#Here is the sorting stage.
		if sortBy == 'duration':
			for k in data[raga].keys():
				data[raga][k].sort(reverse=True)
		
		selectedRecs =0
		iterationCount =0
		run = True
		while (run):
			for k in data[raga].keys():
				#print k, selectedRecs, iterationCount
				if len(data[raga][k]) > iterationCount:
					key = data[raga][k][iterationCount].keys()[0]
					val = data[raga][k][iterationCount][key]
					if not dataSelected[raga].has_key(val):
						dataSelected[raga][val]={}
					selectedRecs+=1
					if selectedRecs >= minRecording:
						run = False
						break
			iterationCount+=1

	#writing selected material in a human readable file to be verified by a musician
	fid = open(outputTXTFILE,'w')
	for raga in raagas:
		fid.write("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n")
		fid.write("Files selected for raga %s\n"%unicode(ragaData[raga]['name']).encode('utf8'))
		fid.write("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n")
		
		for mbid in dataSelected[raga].keys():
			recData = tradition.get_recording(mbid)
			recData['concert_artist'] = tradition.get_concert(recData['concert']['mbid'])['concert_artists'][0]
			fid.write("%s\t%s\t%s\t%s\n"%(recData['concert_artist']['name'],recData['concert']['title'], recData['title'], recData['work']['title']))
			dataSelected[raga][mbid]={'artist':recData['concert_artist']['name'], 'concert':recData['concert']['title'], 'title':recData['title']}
		fid.write("\n")
		fid.write("\n")
	fid.close()
	
	#writing selected mbids to a file
	json.dump(dataSelected, open(selectedRecordingOutputFile, 'w'))	
		


def getDataPerMBID(MBID, path, features=[{'name':'mp3', 'extension':'.mp3'}], collection = 'carnatic'):
	
	if collection == 'hindustani':
		tradition = hn
	elif collection == 'carnatic':
		tradition = ca
		
	for ftr in features:
		if ftr['name']=='mp3':
			src = location[MBID]
			cmd = "scp -r kora:'%s' '%s'"%(src, path)
			os.system(cmd)
			#tradition.download_mp3(MBID, path+ftr['extension'])
			#tradition.download_mp3(MBID, path)


def downloadRagaDataset(selectedRagaDBFile, outputDir, logFile):
	
	data = json.load(open(selectedRagaDBFile,'r'))
	fid = open(logFile,'w')
	for raga in data.keys():
		for mbid in data[raga].keys():
			print "Processing", mbid, data[raga][mbid]
			#outPath = fixpath.fixpath(os.path.join(outputDir,raga, data[raga][mbid]['artist'],data[raga][mbid]['concert'], data[raga][mbid]['title']))
			outPath = fixpath.fixpath(os.path.join(outputDir,raga, data[raga][mbid]['artist'],data[raga][mbid]['concert']))
			if not os.path.exists(outPath):
				os.makedirs(outPath)
			try:
				#getDataPerMBID(mbid, outPath)
				pass
			except:
				print mbid
			
			has_src = 0
			try:
				src = location[mbid]
				has_src = 1
			except:
				print "CAnnot find file for this MBID: %s\n"%mbid
				
			if has_src ==1 and os.path.isfile(outPath+'/'+src.split('/')[-1]):
				fid.write("%s\t%s\t%d\n"%(mbid,outPath+'/'+src.split('/')[-1], 1))
			else:
				fid.write("%s\t%s\t%d\n"%(mbid,outPath+'/'+	src.split('/')[-1], 0))
			
				
	fid.close()
	
		
		
			

"""
Main functinos needed
General Purpose
getDataPerMBID(MBID, path, features=[{'name':'pitch', 'extension':'.pitch'}])


dataset creation generic codes

getAllRagasWithRecordings(collectionId)


DatabaseSpecific
generateTopNRagaDataset(outDir)


"""
		
		
		

