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
from collections import Counter
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.feature_extraction.text import TfidfTransformer
from sklearn.naive_bayes import MultinomialNB
from sklearn.linear_model import SGDClassifier

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

def find_1NN_comms_per_phrase_old(G, phrases, comm_data, raga_comms, comb_criterion):
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
    Output:
        [(dist, comm_id, raga_id), (), ()]
    """
    
    nn_comms = []
    for ii, phrase in enumerate(phrases):
        if not G.has_node(str(phrase)):
            nn_comms.append((-1*INF, -1, ''))
            continue
        dist_overall = []
        for raga in raga_comms.keys():
            dist_raga = []
            for raga_comm in raga_comms[raga]:
                comId = raga_comm['comId']
                comm_phrases = [r['nId']  for r in comm_data[comId]]
                dist_comm = []
                for ii, p in enumerate(comm_phrases):
                    if G.has_edge(str(phrase),str(p)):
                        print "Aaya", comm_data[comId][ii]['ragaId']
                        dist_comm.append(G.get_edge_data(str(phrase),str(p))[0]['weight'])
                if len(dist_comm)==0:
                    dist_comm = -1*INF
                else:
                    dist_comm = np.max(dist_comm)
                dist_raga.append(dist_comm)
            
            dist_raga = np.max(dist_raga)
            dist_overall.append(dist_raga)
        
        ind_max = np.argmax(dist_overall)
        nn_comms.append((dist_overall[ind_max], -1, raga_comms.keys()[ind_max]))
        
        
    return nn_comms


def find_1NN_comms_per_phrase(G, phrases, comm_data, raga_comms, comb_criterion = 'min', myDatabase = '', myUser = ''):
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
    threshold = cons_net.ThshldArray[-1]
    
    if comb_criterion == 'min':
        combine_distance = np.min
    elif comb_criterion == 'mean':
        combine_distance = np.mean
    
    #since reading nids for each community for every phrase in a file is too much, we make a datastructure
    # in which we have comids and nodes in them
    comm_infos = []
    for raga in raga_comms.keys():
        for raga_comm in raga_comms[raga]:
            comm_infos.append({'comId':raga_comm['comId']})
            comm_infos[-1]['nodes'] = [r['nId']  for r in comm_data[raga_comm['comId']]]
            comm_infos[-1]['ragaId'] = raga

    cmd1 = "select target_id, distance from network where source_id = %s"
    
    try:
        con = psy.connect(database=myDatabase, user=myUser) 
        cur = con.cursor()
        print "Successfully connected to the server"
        
        nn_comms = []
        for ii, phrase in enumerate(phrases):
            cur.execute(cmd1%(phrase))
            dist_dict = dict(cur.fetchall())
            if len(dist_dict.keys()) ==0:
                continue
            #print "processing %d of total %d phrases"%(ii+1, len(phrases))
            dist_across_comms = []
            for comm_info in comm_infos:
                dist_within_comm = []
                for p in comm_info['nodes']:
                    if not dist_dict.has_key(p):
                        dist_within_comm.append(threshold)
                    else:
                        dist_within_comm.append(dist_dict[p])
                dist_across_comms.append(combine_distance(dist_within_comm))
            ind_min = np.argmin(dist_across_comms)
            #print ind_min, comm_infos[ind_min]['ragaId'], dist_across_comms[ind_min]
            nn_comms.append((dist_across_comms[ind_min], comm_infos[ind_min]['comId'], comm_infos[ind_min]['ragaId']))
    
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

def classify_file_old(nn_comms):
    
    dists = [r[0] for r in nn_comms]
    ind_min = np.argmax(dists)
    
    return nn_comms[ind_min][2]
    
    
def compute_raga_phrase_distance_distribution(nn_comms, unique_raga_ids, gt_raga):
    """
    This function computes ragas #nodes and distance distribution for a given recording.
    nn_comms contains distances and ragaids of all the nearest neighbors of all the phrases in a file.
    """
    print "total number of phrases = %d"%len(nn_comms)
    dist_MTX = np.zeros((len(unique_raga_ids), len(cons_net.ThshldArray)))
    
    for ii, rid in enumerate(unique_raga_ids):
        ind_rids = [i for i,r in enumerate(nn_comms) if r[2] == rid]
        dists = [nn_comms[d][0] for d in ind_rids]
        disttribution = np.histogram(np.array(dists), bins = cons_net.ThshldArray)[0]
        dist_MTX[ii,1:] = disttribution
        if rid == gt_raga:
            dist_MTX[ii,0] = 1
            
    return dist_MTX
    
def dump_raga_phrase_distance_disribution(dist_MTX, out_dir, raga_name, mbid):
    """
    This function dumps the pdf/png corresponding to the distribution (raga, phrase, distance) for one recording. 
    It places the image in an appropriate folder inside out_dir, named after raga_name
    """
    dir_name = os.path.join(out_dir, raga_name)
    
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)
        
    filename = os.path.join(dir_name, mbid+'.pdf')
    
    #saving the plots
    fig = plt.figure()
    ax = fig.add_subplot(111)
    fsize = 22
    fsize2 = 16
    font="Times New Roman"
    
    
    plt.imshow(dist_MTX, interpolation="nearest")
    plt.xlabel("Distance", fontsize = fsize, fontname=font)
    plt.ylabel("Ragas", fontsize = fsize, fontname=font, labelpad=fsize2)
    fig.savefig(filename, bbox_inches='tight')
    

def raga_recognition_V1(fileListFile, thresholdBin, pattDistExt, n_fold = 16, top_N_com = 10, force_build_network=1, out_dir_dump = '/home/sankalp/Work/Work_PhD/experiments/MotifDiscovery/networkAnalysis/ragaRecognition/Raga_Recognition_debug/ISMIR2015_10RAGA_TONICNORM'):
    """
    This is a wrapper function which performs raga recognition (V1). The chain for raga recognitino being
    Ncrossfold split->training set->clustering->topNragaphrase->testing set->nearest raga communities->classification
    
    Essentially the logic is that one recording can be classified looking at the phrases in that file
    and the nearest neighbors of those phrases in terms of raga communities. A piece in a given raga
    is expected to have nearest communities as the ones which were selected for that specific raga. This is a very basic approach
    because how do we combine multiple phrases is a big question. Also it doesn't take into account occurrences of phrases in a raga. 
    In addition if we have several communities for each raga, we would like to exploit all of them for raga id. Also this method needs 
    a algorithm to select top N raga communities and currently its based on a heuristic formulation. This optimal set of communitie should
    be automatically obtained looking at the communities that best describe a raga. 
    
    NOTE: this task of automatically selecting relevant communities and then classifying a recording is very similar to document classification.
    There is a full theory to select meaningful words (phrase) relevant for a concept (ragas). raga_recognition_V2 function intend to do that.
    """
    
    append_dir = 'DUMP_DIST_Thsldbin_%d_pattDistExt_%s_NComms_%d'%(thresholdBin, pattDistExt, top_N_com)
    
    
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
            
            dist_MTX = compute_raga_phrase_distance_distribution(nn_comms, np.unique(raga_list), raga_list_test[ii])
            
            dump_raga_phrase_distance_disribution(dist_MTX, os.path.join(out_dir_dump, append_dir), raga_list_test[ii], mbid)
            
            
    cnt = 0
    for i in range(len(predicted_raga)):
        if raga_list[i] == predicted_raga[i]:
            cnt+=1
    
    print "You got %d number of ragas right for a total of %d number of recordings"%(cnt, len(predicted_raga))
       
    
    return
        
def get_per_recording_data(comm_data):
    """
    This function simply transforms the incoming community data which is in the form:
    {'comm_index': [u'compId': u'', u'mbid': u'', u'nId': <int>, u'ragaId': u''}]}
    
    this is transformed to a data type
    {'mbid': [(comm_id, nId)]}
    """
    
    per_rec_data = {}
    for comm_id in comm_data:
        for node_data in comm_data[comm_id]:
            if not per_rec_data.has_key(node_data['mbid']):
                per_rec_data[node_data['mbid']]=[]
            per_rec_data[node_data['mbid']].append((comm_id, node_data['nId']))
            
    return per_rec_data
        
def raga_recognition_V2(fileListFile, thresholdBin, pattDistExt, n_fold = 16, force_build_network=0):
    """
    Raga recognition system using document classification and topic modelling techniques.
    In this approach we treat phrases of a recording as words (basically cluster id). 
    """
    
    #constructing the network
    t1 = time.time()
    wghtd_graph_filename = 'graph_temp'+'_'+str(thresholdBin)
    
    if force_build_network or not os.path.isfile(wghtd_graph_filename):
        cons_net.constructNetwork_Weighted_NetworkX(fileListFile, wghtd_graph_filename , thresholdBin, pattDistExt, 0 , -1)
    
    full_net = nx.read_pajek(wghtd_graph_filename)
    
    #performing community detection on the built network
    comm_filename = 'comm_thresholdBin'+'_'+str(thresholdBin)+ pattDistExt +'.community'
    net_pro.detectCommunitiesInNewtworkNX(wghtd_graph_filename, comm_filename)
    
    #fetching relevant data for the community
    comm_data = json.load(open(comm_filename,'r'))
    comm_char.fetch_phrase_attributes(comm_data, database = myDatabase, user= myUser)
    
    #getting per document (recording) words (community index)
    per_rec_data = get_per_recording_data(comm_data)

    t2 = time.time()
    print "time taken = %f"%(t2-t1)
    
    ##########Loop for N_Fold cross validataion##############
    raga_mbid = get_mbids_raagaIds_for_collection(myDatabase, myUser)
    raga_list = [r[0] for r in raga_mbid]
    mbid_list = [r[1] for r in raga_mbid]
    raga_map, map_raga = generate_raga_mapping(raga_list)
    label_list = np.array([raga_map[r] for r in raga_list])
    
    #initializing crossfold object
    cval = cross_val.StratifiedKFold(label_list, n_folds=n_fold)
    
    #splitting folds.
    fold_cnt = 1
    predicted_raga = ['' for r in range(len(raga_mbid))]
    
    count_vect = CountVectorizer()
    tfidf_transformer = TfidfTransformer()
    for train_inds, test_inds in cval:
        docs_train = []
        for train_ind in train_inds:
            #if per_rec_data.has_key(mbid_list[test_ind]):
                #print mbid_list[test_ind], len(per_rec_data[mbid_list[test_ind]])
                #per_rec_words = [p[0] for p in per_rec_data[mbid_list[test_ind]]]
                #print per_rec_words
            #else:
                #print mbid_list[test_ind], 0
            if per_rec_data.has_key(mbid_list[train_ind]):
                per_rec_words = ' '.join([p[0] for p in per_rec_data[mbid_list[train_ind]]])
            else:
                per_rec_words = ''
            docs_train.append(per_rec_words)
        
        X_train_counts = count_vect.fit_transform(docs_train)
        X_train_tfidf = tfidf_transformer.fit_transform(X_train_counts)
        #clf = MultinomialNB().fit(X_train_tfidf, label_list[train_inds])
        clf = SGDClassifier(loss='hinge', penalty='l2', alpha=1e-3, n_iter=5, random_state=42).fit(X_train_tfidf, label_list[train_inds])
        
        docs_test = []
        for test_ind in test_inds:
            #if per_rec_data.has_key(mbid_list[test_ind]):
                #print mbid_list[test_ind], len(per_rec_data[mbid_list[test_ind]])
                #per_rec_words = [p[0] for p in per_rec_data[mbid_list[test_ind]]]
                #print per_rec_words
            #else:
                #print mbid_list[test_ind], 0
            if per_rec_data.has_key(mbid_list[test_ind]):
                per_rec_words = ' '.join([p[0] for p in per_rec_data[mbid_list[test_ind]]])
            else:
                per_rec_words = ''
            docs_test.append(per_rec_words)
        
        X_test_counts = count_vect.transform(docs_test)
        X_test_tfidf = tfidf_transformer.transform(X_test_counts)
        
        predicted = clf.predict(X_test_tfidf)
        for ii, pred_val in enumerate(predicted):
            predicted_raga[test_inds[ii]] = map_raga[pred_val]
        
            
    cnt = 0
    for i in range(len(predicted_raga)):
        if raga_list[i] == predicted_raga[i]:
            cnt+=1
    
    print "You got %d number of ragas right for a total of %d number of recordings"%(cnt, len(predicted_raga))
       
    
    return
    
    
        
        
        
        
    
    
    