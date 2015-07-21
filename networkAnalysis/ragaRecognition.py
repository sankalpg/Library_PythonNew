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
import communityCharacterization as comm_char
import json

INF = np.finfo(np.float).max

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
    """
    Function to map raga_ids to integers
    """
    raga_map = {}
    map_raga = {}
    for ii, r in enumerate(np.unique(raga_ids)):
        raga_map[r] = ii
        map_raga[ii] = r
        
    return raga_map, map_raga
                

def get_phrase_ids_for_files(mbids, database = '', user= ''):
    
    cmd = "select id from pattern where file_id in (select id from file where mbid in %s)"
    cmd1 = "select id from pattern where file_id in (select id from file where mbid = '%s')"
    
    try:
        con = psy.connect(database=database, user=user) 
        cur = con.cursor()
        print "Successfully connected to the server"
        if len(mbids) ==1:
            cur.execute(cmd1%str(mbids[0]))
        else:
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

def find_1NN_comms_per_phrase(G, phrases, comm_data, raga_comms, dthresh=-1, comb_criterion = 'mean', myDatabase = '', myUser = ''):
    """
    This function finds the nearest community in "raga_comms" to the phrases inputted. 
    The community details are in the comm_datae. The distance of a phrase to another phrase 
    of a community is determined based on the input network G. Note that since G already has a threshold
    applied for some phrases there might not even be a connection with any community. There are different
    ways to combine the distances of a phrase to phrases of a community in order to obtain a single score
    for each community. 
    
    Input:
        G: input graph based on which the distances are to be obtained
        phrase: phrases for which the nearest communities are to be obtained
        comm_data: dictionary containing details of the community
        raga_comms: communities corresponding to various ragas
        comb_criterion: criterion to be used to combine the distances of all the phrases of a community to a single number
        dthresh : is the distance threshold applied while computing distance of a phrase from the community.
    Output:
        [(dist, comm_id, raga_id), (), ()]
    """
    cmd = "select distance from network where source_id = %s and target_id = %s"
    
    try:
        con = psy.connect(database=myDatabase, user=myUser) 
        cur = con.cursor()
        print "Successfully connected to the server"
        
        nn_comms = []
        for ii, phrase in enumerate(phrases):
            print "processing %d of total %d phrases"%(ii+1, len(phrases))
            dist_overall = []
            for raga in raga_comms.keys():
                dist_across_comms = []
                for raga_comm in raga_comms[raga]:
                    comId = raga_comm['comId']
                    comm_phrases = [r['nId']  for r in comm_data[comId]]
                    dist_within_comm = []
                    for ii, p in enumerate(comm_phrases):
                        cur.execute(cmd%(phrase,p))
                        dist = cur.fetchone()
                        if dist == None:
                            dist_within_comm.append(cons_net.ThshldArray[-1])
                        else:
                            dist_within_comm.append(dist[0])
                    if len(dist_within_comm)==0:
                        dist_comm = cons_net.ThshldArray[-1]
                    else:
                        if comb_criterion == 'min':
                            dist_comm = np.min(dist_within_comm)
                        elif comb_criterion == 'mean':
                            dist_comm = np.mean(dist_within_comm)
                        
                    
                    dist_across_comms.append(dist_comm)
                
                dist_raga = np.min(dist_across_comms)
                dist_overall.append(dist_raga)
            
            ind_min = np.argmin(dist_overall)
            nn_comms.append((dist_overall[ind_min], -1, raga_comms.keys()[ind_min]))
    
    except psy.DatabaseError, e:
        print 'Error %s' % e
        if con:
            con.rollback()
            con.close()
        sys.exit(1)
    
    if con:
        con.close()
            
        
    return nn_comms
        
        
def classify_file(nn_comms):
    
    dists = [r[0] for r in nn_comms]
    ind_min = np.argmin(dists)
    
    return nn_comms[ind_min][2]
    
    
def compute_raga_phrase_distance_distribution(nn_comms, unique_raga_ids):
    """
    This function computed ragas #nodes and distance distribution for a given recording.
    nn_comms contains distances and ragaids for all the phrases in a file.
    """
    
    for rid in unique_raga_ids:
        ind_rids = [i]
        return True
    
    
    
    
    
    

def raga_recognition_V1(fileListFile, thresholdBin, pattDistExt, n_fold = 16, top_N_com = 10, force_build_network=1):
    """
    This is a wrapper function which performs raga recognition (V1).
    """
    
    
    #constructing the network
    t1 = time.time()
    wghtd_graph_filename = 'graph_temp'+'_'+str(thresholdBin)
    
    if force_build_network or not os.path.isfile(wghtd_graph_filename):
        cons_net.constructNetwork_Weighted_NetworkX(fileListFile, wghtd_graph_filename , thresholdBin, pattDistExt, 0 , -1)
    
    full_net = nx.read_pajek(wghtd_graph_filename)
    
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
    fold_cnt = 1
    predicted_raga = ['' for r in range(len(raga_mbid))]
    
    for train_ind, test_ind in cval:
        mbid_list_test = [raga_mbid[i][1] for i in test_ind]
        raga_list_test = [raga_mbid[i][0] for i in test_ind]
        phrases_remove = get_phrase_ids_for_files(mbid_list_test, myDatabase, myUser)
        
        #reading the original graph from the file
        g = nx.read_pajek(wghtd_graph_filename)
        #removing the edges and nodes which corresponding to the testing data
        g = remove_nodes_graph(g, phrases_remove)
        
        training_graph_filename = 'graph_training'
        nx.write_pajek(g, training_graph_filename)
        
        comm_filename = 'comm'+'_'+str(fold_cnt)+'.community'
        net_pro.detectCommunitiesInNewtworkNX(training_graph_filename, comm_filename)
        
        comm_rank_filename  = 'comm'+'_'+str(fold_cnt)+'.communityRank'
        comm_char.rankCommunities(comm_filename, comm_rank_filename, myDatabase = myDatabase, myUser = myUser)
        
        raga_comm  = comm_char.select_topN_community_per_raaga(comm_rank_filename, top_N_com)
        
        comm_data = json.load(open(comm_filename,'r'))
        comm_char.fetch_phrase_attributes(comm_data, database = myDatabase, user= myUser)
    
        for ii, mbid in enumerate(mbid_list_test):
            phrases_recording = get_phrase_ids_for_files([mbid], myDatabase, myUser)
            
            nn_comms = find_1NN_comms_per_phrase(full_net, phrases_recording, comm_data, raga_comm, myDatabase = myDatabase, myUser = myUser)
            
            raga_id_classify = classify_file(nn_comms)
            
            predicted_raga[test_ind[ii]] = raga_id_classify
            print raga_list_test[ii], raga_id_classify
            
    cnt = 0
    for i in range(len(predicted_raga)):
        if raga_list[i] == predicted_raga[i]:
            cnt+=1
    
    print "You got %d number of ragas right for a total of %d number of recordings"%(cnt, len(predicted_raga))
       
    
    return
        
        
        
        
        
        
        
        
    
    
    