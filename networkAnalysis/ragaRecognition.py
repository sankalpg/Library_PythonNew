import sys
import pickle
import inspect
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
from sklearn.decomposition import TruncatedSVD
import uuid
from sklearn.metrics import confusion_matrix
sys.path.append(os.path.join(os.path.dirname(__file__), '../machineLearning'))
import mlWrapper as ml
from sklearn import preprocessing
from scipy.sparse import csc
import codecs
from mutagen.mp3 import MP3
from mutagen import easyid3


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


#myUser = 'sankalp'
#myDatabase = 'ISMIR2015_10RAGA_TONICNORM'

def get_mbid_from_mp3(mp3_file):
    """
    fetch MBID form an mp3 file
    """
    try:
        mbid = easyid3.ID3(mp3_file)['UFID:http://musicbrainz.org'].data
    except:
        print "problem reading mbid for file %s\n" % mp3_file
        raise
    return mbid


def get_mbids_raagaIds_for_collection(file_list, database = '', user= ''):
    
    cmd = "select raagaid, mbid from file"
    #first lets get all the raga mbid mapp in a single shot
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
    
    #now lets filter this list by the files mentioned in the file_list
    lines = codecs.open(file_list, 'r', encoding = 'utf-8')
    mbids = []
    for line in lines:
        if line.strip() == '':
            continue
        try:
            mbids.append(get_mbid_from_mp3(line.strip()+'.mp3'))
        except:
            print "File %s coundn't fetch any mbid"%line.strip()
    
    mbid2raga = {}
    for r in results:
        mbid2raga[r[1]] = r[0]

    raga_mbid= []
    for mbid in mbids:
        if mbid2raga.has_key(mbid):
            raga_mbid.append([mbid2raga[mbid], mbid])

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
        
        
        
def raga_recognition_V2(out_dir, scratch_dir, fileListFile, thresholdBin, pattDistExt, network_wght_type = -1, force_build_network=0, feature_type = 'tf-idf', pre_processing = -1, norm_tfidf = None, smooth_idf = False, classifier = ('nbMulti', "default"), n_expts = 10, var1 = True, var2 = True, myDatabase = '', myUser = '', type_eval = ("kStratFoldCrossVal", 10), balance_classes =1):
    """
    Raga recognition system using document classification and topic modelling techniques.
    In this approach we treat phrases of a recording as words (basically cluster/community id). 
    phrases->communities(Words)->word frequencies per file->tf-idf kind of features->classification
    Input:
        fileListFile: file which lists all the files to be considered for this anlyasis (there is relevant data extracted and stored in a structed manner from these files)
        thresholdBin: distance threshold (in bins) which is applied the network
        pattDistExt: extension of the file that stores pattern distances
        n_fold: number of cross fold validations
        force_build_network: if 0 the network is not rebuild if it exists on the disk, if 1 its built again no matter what
        feature_type: the type of feature to be used for the classification. Options are
                      'tf': term frequency
                      'tp': term presence (binary value to indicate if the term is present or not)
                      'tf-idf': the typicall term frequency * inverse document frequency
                      'tf-idf_pp1': this is normal tf-idf but with a preprocessing to explicitely remove crappy phrases (gamakas). This is like removing stop words from word count computation.
                                    
                      'tf-idf_pp2': Along with the gamaka phrases (communities), we also remove (as stop word) the communities which are consists of only one mbid. 
        
        classifier: the classifier to be used for the classification task. Options are:
                    'NB': Naive naive bayes
                    'SVM': support vector machines
                    'SGD': svms with SGD training. Somewhere it was recommended for text classification.
                    
        pre_processing: -1 for no preprocessing
                         1 for removing gamaka communities from the analysis (treating them as stop words)
                         NOTE: It feels like this should be taken care of by the IDF computation, but for a small corpus if there is a lot of frequency, the weight is high no matter what. Just to try it out, brain worms!!
                         2: for removing communities which have only one mbid in them. 
                         3: for removing communities for option 1 and 2
        norm_feature: Normalize the final feature vector or not (NOTE: when its on the result seems to be affected a lot by the presence of the gamaka communities)
        network_wght_type: the schema used for weighting edges of the network. Either 0, 1 or -1 for unity weight
                         
    
    """
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    
    
    # path to store all the temporary files
    #scratch_dir = '/home/sankalp/Work/Work_PhD/library_pythonnew/networkAnalysis/scratch_raga_recognition'
    fileListFile_basename  = os.path.basename(fileListFile)
    root_filename = os.path.join(scratch_dir, 'network'+'_'+ str(fileListFile_basename.replace('.','_')) +'_'+ myDatabase+'_'+str(thresholdBin)+'_'+pattDistExt.replace('.',''))
    #root_filename = os.path.join(scratch_dir, 'network'+'_'+ myDatabase+'_'+str(thresholdBin)+'_'+pattDistExt.replace('.',''))
    
    #constructing the network
    t1 = time.time()
    wghtd_graph_filename = root_filename + '.net'
    #building network only when its not already present on the disk
    if force_build_network or not os.path.isfile(wghtd_graph_filename):
        cons_net.constructNetwork_Weighted_NetworkX(fileListFile, wghtd_graph_filename, thresholdBin, pattDistExt, network_wght_type, -1)   # we do not do any significance filtering
    
    #reading the network on the disk (either build in the current call or from previous ones)
    full_net = nx.read_pajek(wghtd_graph_filename)
    
    #performing community detection on the built network
    comm_filename = root_filename +'.community'
    net_pro.detectCommunitiesInNewtworkNX(wghtd_graph_filename, comm_filename)
    
    #fetching relevant data for the community (raga names and file names to analysis)
    comm_data = json.load(open(comm_filename,'r'))
    comm_char.fetch_phrase_attributes(comm_data, database = myDatabase, user= myUser)
    
    #since all the text mining tools consider 0-10 integers as stop words, for making system robust to any unexpected hiccups
    # all the the comm_ids are mapped to uuids
    com_id_2_uuid = {}
    for com_id in comm_data.keys():
		com_id_2_uuid[com_id] = uuid.uuid1().hex
    
    #getting per document (recording) words (community index, phrase instanes)
    per_rec_data = get_per_recording_data(comm_data)

    t2 = time.time()
    print "time taken = %f"%(t2-t1)    
    
    ##########Loop for N_Fold cross validataion##############
    raga_mbid = get_mbids_raagaIds_for_collection(fileListFile, myDatabase, myUser)
    raga_list = [r[0] for r in raga_mbid]
    mbid_list = [r[1] for r in raga_mbid]
    raga_map, map_raga = generate_raga_mapping(raga_list)
    label_list = np.array([raga_map[r] for r in raga_list])
    
    #if there has to be a pre-processing done to remove specific communities, estimating them to put them as stop words
    stop_words = []
    if pre_processing == 1:
        stop_words.extend(comm_char.find_gamaka_communities(comm_data, max_mbids_per_comm =label_list.size/np.unique(label_list).size)[0])
    if pre_processing == 2:
        stop_words.extend(comm_char.get_comm_1MBID(comm_data))
    if pre_processing == 3:
        stop_words.extend(comm_char.find_gamaka_communities(comm_data, max_mbids_per_comm =label_list.size/np.unique(label_list).size)[0])
        stop_words.extend(comm_char.get_comm_1MBID(comm_data))
    
    stop_words = [com_id_2_uuid[s] for s in stop_words]
    accuracy_var1=-1
    if var1:
        ########################### Performing cross fold train testing Variant 1 ###################################
        #In this variant for each fold we generate a training tf-idf vector. Meaning, vocabulary for each fold is solely determined by the training examples#    
        
        #initializing crossfold object
        if type_eval[0] == 'kStratFoldCrossVal':
            cval = cross_val.StratifiedKFold(label_list, n_folds=type_eval[1], shuffle = True, random_state=np.random.randint(100))
        elif  type_eval[0] == 'kFoldCrossVal':
            cval = cross_val.KFold(len(label_list), n_folds=type_eval[1], shuffle = True, random_state=np.random.randint(100))
        elif type_eval[0] == 'LeaveOneOut':
            cval = cross_val.LeaveOneOut(len(label_list))

        mlObj_var1  = ml.experimenter()
        sca=preprocessing.StandardScaler()
        
        # arrays for storing predicted labels and their names
        label_list_pred = -1*np.ones(label_list.shape)    
        predicted_raga = ['' for r in range(len(raga_mbid))]    #placeholder for storing predicted ragas
        
        #starting crossfold validation loop (NOTE: in this variant we only perform a single experiment)
        for mm, (train_inds, test_inds) in enumerate(cval):
            print "Processing fold %d\n"%mm
            
            #initializers needed for analysis of words (community indexes), WE DO IT IN EVERY FOLD TO MAK SURE EVERYTHING FROM THE PREV FOLD IS REMOVED
            count_vect = CountVectorizer(stop_words = stop_words)
            tfidf_transformer = TfidfTransformer(norm=norm_tfidf, smooth_idf=smooth_idf)
            
            docs_train = []                 #storing documents (phrases per recording)
            #preparing tf-idf matrix for the training data
            for train_ind in train_inds:
                if per_rec_data.has_key(mbid_list[train_ind]):  #not every file has phrases found!! (there is one stupid file for which there are no phrases within this distance threshold)
                    per_rec_words = ' '.join([com_id_2_uuid[p[0]] for p in per_rec_data[mbid_list[train_ind]]])
                else:
                    per_rec_words = ''
                docs_train.append(per_rec_words)
            
            #Computing term frequencies (training set) Our vocab is only learned from training examples
            X_train_counts = count_vect.fit_transform(docs_train)
            
            if feature_type == 'tf':
                features_train = X_train_counts.toarray()
            elif feature_type == 'tp':
                features_train = X_train_counts.toarray()
                features_train = np.where(features_train >= 1,1,features_train)
            elif feature_type == 'tf-idf':
                #computing features from term frequencies (training set)
                features_train = tfidf_transformer.fit_transform(X_train_counts)            
                features_train = features_train.toarray()
            else:
                print "Please specify a valid feature type"
                return False

            #checking the input classifier params
            if not isinstance(classifier[1], dict):
                classifier_params = {}
            else:
                classifier_params = classifier[1]
            #training the model with the obtained tf-idf features
            if not mlObj_var1.skl_classifiers.has_key(classifier[0]):
                print "Please provide a valid clsasifier name"
                return False
            
            clf = mlObj_var1.skl_classifiers[classifier[0]]['handle'](**classifier_params)

            if mlObj_var1.skl_classifiers[classifier[0]]['norm_feat_req']:
                features_train = sca.fit_transform(features_train.astype(np.float))
                
            clf.fit(features_train, label_list[train_inds])
            
            docs_test = []
            #preparing the tf-idf matrix for the testing data.
            for test_ind in test_inds:
                if per_rec_data.has_key(mbid_list[test_ind]):
                    per_rec_words = ' '.join([com_id_2_uuid[p[0]] for p in per_rec_data[mbid_list[test_ind]]])
                else:
                    per_rec_words = ''
                docs_test.append(per_rec_words)
            
            #Computing term frequencies (testing set)
            X_test_counts = count_vect.transform(docs_test)

            if feature_type == 'tf':
                features_test = X_test_counts.toarray()
            elif feature_type == 'tp':
                features_test = X_test_counts.toarray()
                features_test = np.where(features_test >= 1,1,features_test)
            elif feature_type == 'tf-idf':
                #computing features from term frequencies (training set)
                features_test = tfidf_transformer.transform(X_test_counts)            
                features_test = features_test.toarray()
            else:
                print "Please specify a valid feature type"
                return False        

            #performing prediction of labels using the trained model
            if mlObj_var1.skl_classifiers[classifier[0]]['norm_feat_req']:
                features_test = sca.transform(features_test.astype(np.float))            
            
            predicted = clf.predict(features_test)
            
            label_list_pred[test_inds] = predicted
            for ii, pred_val in enumerate(predicted):
                predicted_raga[test_inds[ii]] = map_raga[pred_val]
                
        cnt = 0
        for i in range(len(predicted_raga)):
            if raga_list[i] == predicted_raga[i]:
                cnt+=1
        print "You got %d number of ragas right for a total of %d number of recordings (Variant1)"%(cnt, len(predicted_raga))
        
        cMTC_var1 =  confusion_matrix(label_list, label_list_pred)

        accuracy_var1 = float(cnt)/float(len(predicted_raga))
    
        ########################## End of variant 1 of cross fold testing ##################################
    accuracy_var2 = -1
    if var2:
        ########################### Performing cross fold train testing Variant 2 ###################################
        #In this variant tf-idf vectors are computed for the entire dataset. The only affect this will have is in the computatino of idf term. 
        # This way is actually better because testing files will also contribute to the importance that is given to a work in the idf term.
        # Since computation of tf-idf is unsupervised (no raga label used) even this variant should be a valid experimental setup.
        count_vect = CountVectorizer(stop_words = stop_words)
        tfidf_transformer = TfidfTransformer(norm=norm_tfidf, smooth_idf=False)
        mlObj  = ml.experimenter()
        mlObj.setExperimentParams(nExp = n_expts, typeEval = type_eval, nInstPerClass = -1, classifier = classifier, balanceClasses=balance_classes)   #Note that only balanced classes option does stratified kfold exp
        #NOTE: balanceClasses will make sure each fold has equal number of samples from each class. and threre are equal number of feature instances per class
        
        docs_train = [] # this time we will use all the document in out dataset
        for ii in range(label_list.size):
                if per_rec_data.has_key(mbid_list[ii]):  #not every file has phrases found!! (there is one stupid file for which there are no phrases within this distance threshold)
                    per_rec_words = ' '.join([com_id_2_uuid[p[0]] for p in per_rec_data[mbid_list[ii]]])
                else:
                    per_rec_words = ''
                docs_train.append(per_rec_words)
                
        count_all = count_vect.fit_transform(docs_train)
        features = tfidf_transformer.fit_transform(count_all) 
        if False:
            dump = {'features': features.toarray(), 'labels': np.array(label_list)}
            pickle.dump(dump, open('features_dump_300_Hindustani.pkl','w'))
            return True
        mlObj.setFeaturesAndClassLabels(features.toarray(), np.array(raga_list))
        mlObj.runExperiment()

        accuracy_var2 = mlObj.overallAccuracy
        
        print "You got %d number of ragas right for a total of %d number of recordings (Variant2)"%(np.round(mlObj.overallAccuracy*len(label_list)), len(label_list))
        
        ########################## End of variant 2 of cross fold testing ##################################
    
    ##saving experimental results
    fid = open(os.path.join(out_dir,'experiment_results.pkl'),'w')
    results = {}
    if var1:
        results.update({'var1': {'cm': cMTC_var1, 'gt_label': label_list, 'pred_label': label_list_pred, 'mbids': mbid_list, 'mapping': map_raga, 'accuracy': float(cnt)/len(predicted_raga)}})
    if var2:
        results.update({'var2': {'cm': mlObj.cMTXExp, 'gt_label': mlObj.classLabelsInt, 'pred_label':mlObj.decArray, 'mbids': mbid_list, 'mapping':mlObj.cNames, 'accuracy': mlObj.overallAccuracy, 'pf_accuracy':mlObj.accuracy}})
    
    pickle.dump(results, fid)
    #also dumping the input params to this function
    params_input = {}
    for k in inspect.getargspec(raga_recognition_V2).args:
        params_input[k] = locals()[k]
    
    fid = open(os.path.join(out_dir,'experiment_params.json'),'w')
    json.dump(params_input, fid)

    return np.max([accuracy_var1, accuracy_var2])
    
        
def raga_recognition_V3():
    """
    This version combines V1 and V2. TF-IDF are computed using KNN words of the given words, since many documents have no words from the given vocabulary. 
    """
        
        
    
    
    