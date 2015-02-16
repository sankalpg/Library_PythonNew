# -*- coding: UTF-8 -*-
import compmusic
import sys
import os
import pickle
import copy 



def getCollectionInfoAll(collID='213347a9-e786-4297-8551-d61788c85c80'):
    """
    This method should return all the artists + ragas + talas + works in the releases (in all the tracks)
    """
    compmusic.mb.auth('compmusic','compmusic321')
    
    #get all the releases in a collection
    releaseIds = compmusic.musicbrainz.get_releases_in_collection(collID)
    print "got list of all release ids"
    
   # set local host (will be much faster)
    compmusic.mb.set_hostname('musicbrainz.s.upf.edu')
    
    artistReleaseName = []
    artistReleaseId = []
    totalTracks = []
    artistsIds=[]
    artistsNames=[]
    workIds=[]
    workNames=[]
    
    workCompList = []
    ragaCompList = []
    talaCompList = []
    accArtistCompList = []
    
    totalLength = []
    
    for ii, releaseId in enumerate(releaseIds):
        
        hasRelAccomp=0
        
        print "processing release id %s %d out of %d"%(releaseId, ii+1, len(releaseIds))
        
        #get the release
        rel = compmusic.mb.get_release_by_id(releaseId, includes=["artists", "recordings", "artist-rels"])
        rel = rel["release"]
        
        releaseArtists=[]
        
        totalTracks.append(rel['medium-list'][0]['track-count'])
        #artist credited directly with release
        creditList = rel.get("artist-credit",[])
        for credit in creditList:
            if isinstance(credit, dict):
                artistReleaseName.append(credit['artist']['name'])
                artistReleaseId.append(credit['artist']['id'])
                #neede separately for coverage analysis
                releaseArtists.append(credit['artist']['id'])

    
        #get list of recordind ids in this release
        recordings = []
        for medium in rel["medium-list"]:
            for track in medium["track-list"]:
                recordings.append(track["recording"]["id"])
                totalLength.append(int(track["recording"]["length"]))
                
        
        
        releaseRelList = rel.get("artist-relation-list", [])
        releaseRelArtists=[]
        if len(releaseRelList)>0:
            performanceInfo = _get_artist_performances(releaseRelList)
            for perf in performanceInfo:
                    artistsIds.append(perf[1])
                    artistsNames.append(perf[0])
                    releaseRelArtists.append(perf[1])
                    temp1 = copy.deepcopy(releaseArtists)
                    temp1.extend(releaseRelArtists)
                    if not len(list(set(temp1))) == len(list(set(releaseArtists))):
                        hasRelAccomp=1

    
        #iterate over recordinds ids to fetch all the data
        for recid in recordings:
            hasAccomp=0
            hasRag=0
            hasTal=0
            hasWork=0
            
            recordingArtists=[]
            mbrec = compmusic.mb.get_recording_by_id(recid, includes=["tags", "work-rels", "artist-rels"])
            mbrec = mbrec["recording"]
            
            relList = mbrec.get("artist-relation-list",[])
            if len(relList)>0:
                performanceInfo = _get_artist_performances(relList)
                for perf in performanceInfo:
                    artistsIds.append(perf[1])
                    artistsNames.append(perf[0])
                    recordingArtists.append(perf[1])
                
                temp1 = copy.deepcopy(releaseArtists)
                temp1.extend(recordingArtists)
                if not len(list(set(temp1))) == len(list(set(releaseArtists))):
                    hasAccomp=1
                else:
                    print recid
            
            workList =  mbrec.get("work-relation-list", [])
            
            if(len(workList))>0:
                hasWork=1
                for workinfo in workList:
                    if workinfo["type"]=='performance':
                        workIds.append(workinfo["work"]["id"])
                        workNames.append(workinfo["work"]["title"])
                        
            tagList = mbrec.get("tag-list", [])
            for t in tagList:
                if t["name"].lower().count('rag')>0 or t["name"].lower().count('raag')>0:
                    hasRag=1
                if t["name"].lower().count('tal')>0 or t["name"].lower().count('taal')>0:
                    hasTal=1
                    
            
            #completeness 
            workCompList.append(hasWork)
            ragaCompList.append(hasRag)
            talaCompList.append(hasTal)
            accArtistCompList.append(hasAccomp or hasRelAccomp)            
        
    print "total artist %d"%len(artistsIds)
    print "total works %d"%len(workIds)
    
    info = {'artistRelNames': artistReleaseName, 'artistRelIds': artistReleaseId,'artistIds': artistsIds,'artistNames': artistsNames,'workIds': workIds,'workNames':workNames, 'nTracks': totalTracks, 'workCompList': workCompList, 'ragaCompList':ragaCompList, 'talaCompList':talaCompList, 'accArtistCompList':accArtistCompList, 'totalLength':totalLength}
    
    fid = open('HindustaniInfo.pkl', "w")
    pickle.dump(info, fid)
    fid.close()
    
    
        

def _get_artist_performances(artistrelationlist):
        performances = []
        for perf in artistrelationlist:
            if perf["type"] in ["vocal", "instrument"]:
                artistid = perf["target"]
                artistName = perf["artist"]["name"]
                attrs = perf.get("attribute-list", [])
                is_lead = False
                for a in attrs:
                    if "lead" in a:
                        is_lead = True
                if perf["type"] == "instrument":
                    inst = perf["attribute-list"][0]
                else:
                    inst = "vocal"
                performances.append((artistName, artistid, inst, is_lead))
        return performances
    
def getStatsHindustani(pickleDump):
    fid = open(pickleDump,'r')
    data = pickle.load(fid)
    fid.close()
    
    works = data['workIds']
    artist = data['artistIds']
    
    accArtistComp = data['accArtistCompList']
    accRagaComp = data['ragaCompList']
    accTalaComp = data['talaCompList']
    accWorkComp = data['workCompList']
    
    print "Total number of tracks \t\t: %d"%len(accTalaComp)
    print "Total number of works \t\t: %d"%len(list(set(works)))
    print "Total number of artists \t\t: %d"%len(list(set(artist)))
    
    print "Number of tracks with raga tag \t\t: %d"%sum(accRagaComp)
    print "Number of tracks with tala tag \t\t: %d"%sum(accTalaComp)
    print "Number of tracks with work relation \t\t: %d"%sum(accWorkComp)
    print "Number of tracks with accompanist relation \t\t: %d"%sum(accArtistComp)
    
 
