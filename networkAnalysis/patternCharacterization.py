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
    
    

def rankCommunitiesRagaUniquePhrases(communityInfoFile, outputFile, myDatabase= ''):
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
        for comId in community.keys():
            for nodeInfo in community[comId]:
                cur.execute(query%nodeInfo['nId'])
                ragaId, mbid = cur.fetchone()
                nodeInfo['ragaId'] = ragaId
                nodeInfo['mbid'] = mbid
       
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
    
    communityRank = qsort1(communityRank)
    
    json.dump({'comData':community, 'comRank': communityRank}, open(outputFile,'w'))
    
            
    
    
