import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import networkx as nx
import snap
import community
import psycopg2 as psy
import json


myUser = 'sankalp'


def qsort1(list):
    """Quicksort using list comprehensions"""
    if list == []: 
        return []
    else:
        pivot = list[0]
        lesser = qsort1([x for x in list[1:] if x['ragaRank'] < pivot['ragaRank']])
        greater = qsort1([x for x in list[1:] if x['ragaRank'] >= pivot['ragaRank']])
        return lesser + [pivot] + greater
    
def peakiness(d):
    
    am = float(np.mean(d))
    gm = float(np.power(np.product(np.array(d).astype(np.float)), 1.0/len(d)))
    if gm==0:
        print gm
    return am/gm

def cmpComm(a,b):
    
    if a['rank'] > b['rank']:
        return 1
    elif a['rank'] < b['rank']:
        return -1
    else:
        return 0
    
def cmpBC_Nodes(a,b):
    if a['bc'] > b['bc']:
        return 1
    elif a['bc'] < b['bc']:
        return -1
    else:
        if a['deg'] > b['deg']:
            return 1
        elif a['deg'] < b['deg']:
            return -1
        else:
            return 0

def rankCommunities(communityInfoFile, outputFile, myDatabase= ''):
    """
    This function assigns ranks to detected communities. In order to obtain the best set of parameters 
    for the network analysis we plan to rank communities and empirically come up with a way to obtain 
    best set of prameters based on predefined criterions.
    
    """
    
    communityData = np.loadtxt(communityInfoFile)
    
    community = {}
    #creating dictionary of communities
    for row in communityData:
        if not community.has_key(row[1]):
            community[row[1]] = []
        community[row[1]].append({'nId': row[0], 'ragaId':'', 'mbid':''})
    
    print "Total number of communities are: %d\n"%(len(community.keys()))

    #obtaining raga and file id for every node in every community
    try:
        con = psy.connect(database=myDatabase, user=myUser) 
        cur = con.cursor()
        query = "select raagaid, mbid from file where id = (select file_id from pattern where id = %s)"
        cnt=0
        for comId in community.keys():
            for nodeInfo in community[comId]:
                cur.execute(query%int(nodeInfo['nId']))
                ragaId, mbid = cur.fetchone()
                nodeInfo['ragaId'] = ragaId
                nodeInfo['mbid'] = mbid
                cnt+=1
                #if cnt%100==0:
                    #print cnt
       
    except psy.DatabaseError, e:
        print 'Error %s' % e
        if con:
            con.rollback()
            con.close()
        sys.exit(1)
    
    if con:
        con.close()    
    
    communityRank = []
    #ranking communities based on their number of nodes, number of different ragaIds and number of different mbids
    # We rank a community high-> 1) Least number of unique ragas associated with it, 2) Max number of nodes 3) Max number of unique mbids
    for ii, comId in enumerate(community.keys()):
        
        ragaIds = [r['ragaId']  for r in community[comId]]
        nIds = [r['nId']  for r in community[comId]]
        mbids = [r['mbid']  for r in community[comId]]
        
        uniqueRagaIds = np.unique(ragaIds)
        ragaIdHist = []
        for uniqueRagaId in uniqueRagaIds:
            ragaIdHist.append(ragaIds.count(uniqueRagaId))
        
        sortInd = np.argsort(ragaIdHist)
        
        uniqueMBIDs = np.unique(mbids)
        mbidHist = []
        for mbid in uniqueMBIDs:
            mbidHist.append(mbids.count(mbid))
        
        communityRank.append({'rank':0, 'comId':-1 , 'nNodes':-1, 'nFiles':-1, 'nRagas':-1})
        ### TODO Criterion1: communityRank[-1]['rank'] = len(uniqueMBIDs)*len(np.unique(nIds))*(float(ragaIdHist[sortInd[-1]])/float(np.sum(ragaIdHist)))/peakiness(mbidHist)
        ### TODO Ctriterion2:
        mbidHist = np.array(mbidHist).astype(np.float)
        communityRank[-1]['rank'] = np.round(len(np.unique(nIds))*((float(ragaIdHist[sortInd[-1]])/float(np.sum(ragaIdHist)))**4)*fileCentroid(mbidHist),2)
        communityRank[-1]['comId'] = comId
        communityRank[-1]['nNodes'] = len(np.unique(nIds))
        communityRank[-1]['nFiles'] = len(uniqueMBIDs)
        communityRank[-1]['nRagas'] = len(uniqueRagaIds)
        
        #hard thresholds
        if len(uniqueRagaIds) >=2:
            if (float(ragaIdHist[sortInd[-1]])/float(np.sum(ragaIdHist))) < 0.8:
                 communityRank[-1]['rank'] = -100
                 
        """
        communityRank.append({'ragaRank':0, 'nodesRank':0, 'filesRank':0, 'comId':-1 })
        ragaIds = [r['ragaId']  for r in community[comId]]
        nIds = [r['nId']  for r in community[comId]]
        mbids = [r['mbid']  for r in community[comId]]
        
        uniqueRagaIds = np.unique(ragaIds)
        ragaIdHist = []
        for uniqueRagaId in uniqueRagaIds:
            ragaIdHist.append(ragaIds.count(uniqueRagaId))
        sortInd = np.argsort(ragaIdHist)
        ragaRank = float(ragaIdHist[sortInd[-1]])/float(np.sum(ragaIdHist))
        communityRank[-1]['ragaRank'] = ragaRank 
        communityRank[-1]['nodesRank'] = len(np.unique(nIds))
        communityRank[-1]['filesRank'] = len(np.unique(mbids))
        communityRank[-1]['comId'] = comId
        """
    
    communityRank.sort(cmp=cmpComm, reverse=True)
    
    json.dump({'comData':community, 'comRank': communityRank}, open(outputFile,'w'))

def fileCentroid(d):
    d = np.array(d)
    d = d[np.argsort(d)[::-1]]
    
    inds  = range(1,len(d)+1)
    centroid = float(np.sum(d*inds))/float(np.sum(d))
    
    return centroid

def reRankCommunity(inpComm, outComm):
    
    com1 = json.load(open(inpComm))
    community = com1['comData']
    
    communityRank = []
    #ranking communities based on their number of nodes, number of different ragaIds and number of different mbids
    # We rank a community high-> 1) Least number of unique ragas associated with it, 2) Max number of nodes 3) Max number of unique mbids
    for ii, comId in enumerate(community.keys()):
        
        ragaIds = [r['ragaId']  for r in community[comId]]
        nIds = [r['nId']  for r in community[comId]]
        mbids = [r['mbid']  for r in community[comId]]
        
        uniqueRagaIds = np.unique(ragaIds)
        ragaIdHist = []
        for uniqueRagaId in uniqueRagaIds:
            ragaIdHist.append(ragaIds.count(uniqueRagaId))
        
        sortInd = np.argsort(ragaIdHist)
        
        uniqueMBIDs = np.unique(mbids)
        mbidHist = []
        for mbid in uniqueMBIDs:
            mbidHist.append(mbids.count(mbid))
        
        communityRank.append({'rank':0, 'comId':-1 , 'nNodes':-1, 'nFiles':-1, 'nRagas':-1})
        
        communityRank[-1]['rank'] = np.round(len(np.unique(nIds))*((float(ragaIdHist[sortInd[-1]])/float(np.sum(ragaIdHist)))**2)*fileCentroid(mbidHist),2)
        communityRank[-1]['comId'] = comId
        communityRank[-1]['nNodes'] = len(np.unique(nIds))
        communityRank[-1]['nFiles'] = len(uniqueMBIDs)
        communityRank[-1]['nRagas'] = len(uniqueRagaIds)
        
        #hard thresholds
        if len(uniqueRagaIds) >=2:
            if (float(ragaIdHist[sortInd[-1]])/float(np.sum(ragaIdHist))) <.75:
                 communityRank[-1]['rank'] = -100
            
    
    communityRank.sort(cmp=cmpComm, reverse=True)
    
    json.dump({'comData':community, 'comRank': communityRank}, open(outComm,'w'))
    
def convertFormat(sec):
        
    hours = int(np.floor(sec/3600))
    minutes = int(np.floor((sec - (hours*3600))/60))
    
    seconds = sec - ( hours*3600 + minutes*60)
    
    return str(hours) + ':' + str(minutes) + ':' + str(seconds)  


def extractRagaPhrases(audio_root, root_name, N, comRankExt = '', netExt= '', myDatabase= '', betcenExt = '', degreeExt = '', audioGutter = 0):
    """
    This function parses community rank file and select best N community from each raga.
    The phrases corresponding to these selected communities are dumped in mp3 format.
    
    NOTE that this function expects community data in a sortedfashion
    """
    
    cmd1 = "select file_id, start_time, end_time from pattern where id = %s"
    cmd2 = "select filename from file where id = %s"
    
    #G = nx.read_pajek(root_name + netExt)
    #bc = nx.betweenness_centrality(G)
    bc = json.load(open(root_name+betcenExt))
    degree = json.load(open(root_name+degreeExt))
    
    phrases = {}
    
    data = json.load(open(root_name+comRankExt))
    comData = data['comData']
    comRank = data['comRank']
    comRank.sort(cmp=cmpComm, reverse=True)
    
    for c in comRank:
        community = comData[str(c['comId'])]
        ragaIds = [r['ragaId']  for r in community]
        
        uniqueRagaIds = np.unique(ragaIds)
        ragaIdHist = []
        for uniqueRagaId in uniqueRagaIds:
            ragaIdHist.append(ragaIds.count(uniqueRagaId))
        
        sortInd = np.argsort(ragaIdHist)
        raagaClass = uniqueRagaIds[sortInd[-1]]
        
        if not phrases.has_key(raagaClass):
            phrases[raagaClass]= []
        
        if len( phrases[raagaClass])<= N:
            phrases[raagaClass].append({'rank':c['rank'], 'comId':c['comId']})
        
    if not os.path.isdir(root_name):
        os.makedirs(root_name)
        
    try:
        con = psy.connect(database=myDatabase, user='sankalp') 
        cur = con.cursor()
        
        for ii, ragaId in enumerate(phrases.keys()):
            
            ragaPath = os.path.join(root_name, ragaId)
            if not os.path.isdir(ragaId):
                os.makedirs(ragaId)
            
            for jj, phraseInfo in enumerate(phrases[ragaId]):
                
                phrasePath = os.path.join(ragaPath, 'phrase_%d_rank%f'%((jj+1), phrases[ragaId][jj]['rank']))
                if not os.path.isdir(phrasePath):
                    os.makedirs(phrasePath)
                
                comId = str(phrases[ragaId][jj]['comId'])
                
                nIds = [r['nId']  for r in comData[comId]]
                nodeProp = [{'bc': bc[str(int(jj))], 'deg':degree[str(int(jj))], 'index':ii} for ii,jj in enumerate(nIds)]
                
                nodeProp.sort(cmp=cmpBC_Nodes, reverse = True)
                
                #ind_sort = np.argsort(bc_vals)
                
                #ind_sort = range(len(nIds))
                ind_sort = [ii['index'] for ii in nodeProp]
                
                for kk, nodeInd in enumerate(ind_sort):
                    cur.execute(cmd1%comData[comId][nodeInd]['nId'])
                    fileID, start_time, end_time = cur.fetchone()
                    cur.execute(cmd2%fileID)
                    filePath = cur.fetchone()[0]
                    audioFile = os.path.join(audio_root, filePath)                    
                    outfile = os.path.join(phrasePath, '%d_%s_%s.mp3'%(kk,comData[comId][nodeInd]['nId'], comData[comId][nodeInd]['ragaId']))
        
                    cmd = "sox \"%s\" \"%s\" trim %s =%s"%(audioFile, outfile, convertFormat(start_time-audioGutter), convertFormat(end_time+audioGutter))
                    os.system(cmd)

    except psy.DatabaseError, e:
        print 'Error %s' % e    
        if con:
            con.rollback()
            con.close()
            sys.exit(1)

    if con:
        con.close()
    

            
            
        
        
    
    
    
    
    
    