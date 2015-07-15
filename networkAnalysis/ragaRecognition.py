import sys
import os
import numpy as np
import snap
import networkx as nx
import matplotlib.pyplot as plt
import constructNetwork as cons_net
import psycopg2 as psy
from sklearn import cross_validation as cross_val
import time
import networkProcessing as net_pro

######## TASKS TO BE DONE IN V1 OF RAGA RECOGNITION

#Build a network with a threshold
#Split the dataset iteratively into training and testing set
#For training do NON_OVERLAPPING community detection
#Characterize detected communities and identify pakad communities per raga
#For testing set for each file obtain N nearest communities
#Classify testing set based on KNN kind of a classifier
#Evaluation of the raga recognition task
#Randomization and multiple iterations of the entire experiment


myUser = 'sankalp'
myDatabase = 'ISMIR2015_10RAGA_TONICNORM'


def get_mbids_raagaIds_for_collection(database = '', user= ''):
    
    cmd = "select raagaid, mbid from file"
    
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
    
    raga_mbid= []
    for r in results:
        raga_mbid.append([r[0], r[1]])
    
    return raga_mbid

def generate_raga_mapping(raga_ids):
    
    raga_map = {}
    map_raga = {}
    for ii, r in enumerate(np.unique(raga_ids)):
        raga_map[r] = ii
        map_raga[ii] = r
        
    return raga_map, map_raga
                

def get_phrase_ids_for_files(mbids, database = '', user= ''):
    
    cmd = "select id from pattern where file_id in (select id from file where mbid in %s)"
    
    try:
        con = psy.connect(database=database, user=user) 
        cur = con.cursor()
        print "Successfully connected to the server"
        cur.execute(cmd%str(tuple(mbids)))
        results = cur.fetchall()
        
    except psy.DatabaseError, e:
        print 'Error %s' % e
        if con:
            con.rollback()
            con.close()
        sys.exit(1)
    
    if con:
        con.close()
    
    phrases = []
    for r in results:
        phrases.append(int(r[0]))
    
    return phrases

def remove_nodes_graph(g, node_ids):
    for n in node_ids:
        if g.has_node(str(n)):
            g.remove_node(str(n))
    
    return g

def raga_recognition_V1(fileListFile, thresholdBin, pattDistExt, n_fold = 16):
    """
    This is a wrapper function which performs raga recognition (V1).
    """
    
    
    #constructing the network
    t1 = time.time()
    #full_net = cons_net.constructNetworkSNAP(fileListFile, 'graph_temp', thresholdBin, pattDistExt)
    t2 = time.time()
    print "time taken = %f"%(t2-t1)
    
    ##########Loop for N_Fold cross validataion##############
    raga_mbid = get_mbids_raagaIds_for_collection(myDatabase, myUser)
    raga_list = [r[0] for r in raga_mbid]
    raga_map, map_raga = generate_raga_mapping(raga_list)
    label_list = [raga_map[r] for r in raga_list]
    
    #initializing crossfold object
    cval = cross_val.StratifiedKFold(label_list, n_folds=n_fold)
    
    #splitting folds.
    for train_ind, test_ind in cval:
        mbid_list = [raga_mbid[i][1] for i in test_ind]
        phrases_remove = get_phrase_ids_for_files(mbid_list, myDatabase, myUser)
        
        #reading the original graph from the file
        g = nx.read_pajek('graph_temp')
        #removing the edges and nodes which corresponding to the testing data
        g = remove_nodes_graph(g, phrases_remove)
        
        nx.write_pajek(g, 'graph_temp_training')
        
        net_pro.detectCommunitiesInNewtworkNX('graph_temp_training', 'comm')
    
    
    