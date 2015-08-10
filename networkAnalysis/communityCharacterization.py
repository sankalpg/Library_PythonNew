import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import networkx as nx
import snap
import community
import psycopg2 as psy
import json
from collections import Counter


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


def fetch_phrase_attributes(community_dict, database = '', user= ''):
    """
    This function fetches mbid raagaId and compId for all the phrases in the file community_dict.
    These attributes are fetched from the database.
    community_dict: is a dict with format {'community_id':{nId:<>,'ragaId':<>, 'mbid':<>, 'compId':<>}}
    """
    #obtaining raga and file id for every node in every community_dict
    con = psy.connect(database=database, user=user) 
    try:
        cur = con.cursor()
        query = "select raagaid, mbid from file where id = (select file_id from pattern where id = %s)"
        for comId in community_dict.keys():
            for nodeInfo in community_dict[comId]:
                cur.execute(query%int(nodeInfo['nId']))
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
    return community_dict

def get_histogram_sorted(data_list):
    """
    given a list of data points this function creates a sorted histogram
    """
    
    cntr = Counter(data_list)
    sort_data = cntr.most_common()
    
    hist = [s[1] for s in sort_data]
    items = [s[0] for s in sort_data]
    
    return hist, items
    

def rank_community_raga_phrases(comm_data):
    """
    This function ranks a community given the data (comm_dat) for that community. comm_data comprises
    phrases in the community along with their attributes.
    
    The rank given by this function is such that higher the rank more likely is that the community 
    represent a cluster of raga characteristic phrases.    
    """
    
    ragaIds = [r['ragaId']  for r in comm_data]
    nIds = [r['nId']  for r in comm_data]
    mbids = [r['mbid']  for r in comm_data]
    
    raga_hist, raga_val = get_histogram_sorted(ragaIds)
    mbid_hist, mbid_val = get_histogram_sorted(mbids)
    
    single_raga_fraction = float(raga_hist[0])/float(np.sum(raga_hist))
    
    raga_phrase_rank = np.round(len(np.unique(nIds))*(single_raga_fraction**4)*fileCentroid(mbid_hist),2)
    
    #NOTE this is a hard threshold to avoid big hubs with lots of nodes and junky info
    if len(raga_val) >=2:
            if (single_raga_fraction) < 0.8:
                 raga_phrase_rank = -100
    
    return raga_phrase_rank, len(np.unique(nIds)), len(mbid_val), len(raga_val)
    
    

def rankCommunities(communityInfoFile, outputFile, myDatabase= '', myUser = ''):
    """
    This function parses the community file and assigns a rank to each community based on a criterion.
    
    The criterion for ISMIR2015-initial submission paper was: "np.round(len(np.unique(nIds))*((float(ragaIdHist[sortInd[-1]])/float(np.sum(ragaIdHist)))**4)*fileCentroid(mbidHist),2"
    
    This is set empirically looking at the data and distribution of various attributes within the community.
    """
    
    #reading community data
    community = json.load(open(communityInfoFile, 'r'))
    #print "Total number of communities are: %d\n"%(len(community.keys()))
    
    #fetching attributes from database needed for assigning ranks to communities
    community = fetch_phrase_attributes(community, myDatabase, myUser)    

    communityRank = []
    #ranking communities based on their number of nodes, number of different ragaIds and number of different mbids
    # We rank a community high-> 1) Least number of unique ragas associated with it, 2) Max number of nodes 3) Max number of unique mbids
    for ii, comId in enumerate(community.keys()):
        
        raga_rank, n_nodes, n_files, n_ragas = rank_community_raga_phrases(community[comId])
        
        communityRank.append({'rank':raga_rank, 'comId':comId , 'nNodes':n_nodes, 'nFiles':n_files, 'nRagas':n_ragas})
    
    communityRank.sort(cmp=cmpComm, reverse=True)
    
    json.dump({'comData':community, 'comRank': communityRank}, open(outputFile,'w'))
    
    


def fileCentroid(d):
    d = np.array(d)
    d = d[np.argsort(d)[::-1]]
    
    inds  = range(1,len(d)+1)
    centroid = float(np.sum(d*inds))/float(np.sum(d))
    
    return centroid


    
def convertFormat(sec):
        
    hours = int(np.floor(sec/3600))
    minutes = int(np.floor((sec - (hours*3600))/60))
    
    seconds = sec - ( hours*3600 + minutes*60)
    
    return str(hours) + ':' + str(minutes) + ':' + str(seconds)  


def select_topN_community_per_raaga(comm_rank_file, N):
    """
    Given a community rank file, this function selectes the top N community according to the rank.
    """
    
    data = json.load(open(comm_rank_file))
    comData = data['comData']
    comRank = data['comRank']
    
    top_comms = {}
    #iterating over all the communities and selecting top N per raga
    for c in comRank:
        comm = comData[str(c['comId'])]
        ragaIds = [r['ragaId']  for r in comm]
        #assigning a raga to the community based on the majority voting of each node's raga
        raga_hist, raga_names = get_histogram_sorted(ragaIds)
        comm_raga = raga_names[0]
        
        if not top_comms.has_key(comm_raga):
            top_comms[comm_raga]= []
        
        if len(top_comms[comm_raga])< N:
            top_comms[comm_raga].append({'rank':c['rank'], 'comId':c['comId']})

    return top_comms
        
    


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

    data = json.load(open(root_name+comRankExt))
    comData = data['comData']
    comRank = data['comRank']

    phrases = select_topN_community_per_raaga(root_name+comRankExt, N)
        
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
    

            
            
        
if  __name__ == "__main__":       
    #pass
# Refactor this file. The code is written in a quite messy way. Split functions into smaller functions.    

#Unit test stuff    
#performing community detection   
    
    #rankCommunities('unitTests/Weight_-1___DTshld_10___PostFilt___C.community', 'unitTests/Weight_-1___DTshld_10___PostFilt___C.communityRank', myDatabase= 'ISMIR2015_10RAGA_TONICNORM', myUser = myUser)
    
    extractRagaPhrases('/media/Data/Datasets/PatternProcessing_DB/unsupervisedDBs/carnaticDB/Carnatic10RagasISMIR2015DB/audio', 'unitTests/Weight_-1___DTshld_10___PostFilt___C', 2, comRankExt = '.communityRank', netExt= '', myDatabase= 'ISMIR2015_10RAGA_TONICNORM', betcenExt = '.bwcen', degreeExt = '.degree', audioGutter = 0)
    
    
    