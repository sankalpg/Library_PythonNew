import sys
import os
import numpy as np
import snap
sys.path.append(os.path.join(os.path.dirname(__file__), '../batchProcessing'))
import networkx as nx
import matplotlib.pyplot as plt

import batchProcessing as BP
import psycopg2 as psy
import networkProcessing as netPro

FLT_MAX = np.finfo(np.float).max

myUser = 'sankalp'

ISMIR2015_10RAGASCarnatic = ['55', '10', '27', '8', '9', '137', '159', '20', '13', '210']
colors = {'55':'red', '10':'blue', '27':'black', '8':'green', '9':'yellow', '137':'pink', '159':'orange', '20':'magenta', '13':'cyan', '210':'purple'}

ThshldArray = np.array([0.0,
2082.47,
8329.86,
18742.19,
33319.45,
52061.64,
74968.76,
102040.82,
133277.80,
168679.72,
208246.56,
251978.34,
299875.05,
351936.69,
408163.27,
468554.77,
533111.20,
601832.57,
674718.87,
751770.10,
832986.26,
918367.35,
1007913.37,
1101624.32,
1199500.21,
1301541.02,
1407746.77,
1518117.45,
1632653.06,
1751353.60,
1874219.08,
2001249.48,
2132444.81,
2267805.08,
2407330.28,
2551020.41,
2698875.47,
2850895.46,
3007080.38,
3167430.24,
3331945.02,
3500624.74,
3673469.39,
3850478.97,
4031653.48,
4216992.92,
4406497.29,
4600166.60,
4798000.83,
5000000.00])

def constructNetworkSNAP_NodeLabels(fileListFile, outputNetworkFile, thresholdBin, pattDistExt, myDatabase=''):
    """
    This function uses Snap library to generate a labelled network. 
    :fileListFile: file which contains list of filenames which are to be used for the network generation
    :outputNetworkFile: output file in which we want to save our network
    :thresholdBin: bin index based on which a threshold will be applied on the distance between the patterns. This mapping between the bin and the distance threshold is same across every module, including the module for computing clustering coefficient.
    :pattDistExt: file extension of the file which has pattern distances 
    :myDatabase: the name of the PSQL database from which raga labels are fetched for each pattern
    """   
    
    #initialize an graph
    G = snap.PNEANet.New()
    
    #G.AddFltAttrE("dist", 0.0) #defining label on the edge
    
    # creating the network
    filenames = open(fileListFile).readlines()
    for filename in filenames:
        fname = filename.strip() + pattDistExt
        data = np.loadtxt(fname)
        #lines = open(fname).readlines()
        for ii in range(data.shape[0]):
            #id1 = int(line.split('\t')[0].strip())
            #id2 = int(line.split('\t')[1].strip())
            #dist = float(line.split('\t')[2].strip())
            id1 = int(data[ii][0])
            id2 = int(data[ii][1])
            dist = float(data[ii][2])
            #print id1, id2, dist
            
            if dist >= ThshldArray[thresholdBin-1]:
                continue
            if not G.IsNode(id1):
                G.AddNode(id1)
            if not G.IsNode(id2):
                G.AddNode(id2)
            if not G.IsEdge(id1,id2):
                G.AddEdge(id1,id2)
                eId = G.GetEId(id1, id2)
                #G.AddFltAttrDatE(eId,dist,'dist');
                
    
    #G.AddIntAttrN("RagaID", 0) #defining label on the node
    NIdColorH = snap.TIntStrH()
    NIdLabelH = snap.TIntStrH()
    
    #annotating netowrk nodes with raga labels and edges with distances
    cmd1 = "select raagaId from file where id = (select file_id from pattern where id =%d)"
    
    try:
        con = psy.connect(database=myDatabase, user=myUser) 
        cur = con.cursor()
        print "Successfully connected to the server"
        
        for NI in G.Nodes():
            nid = NI.GetId()
            cur.execute(cmd1%(nid))
            ragaId = cur.fetchone()[0]
            print ragaId
            NIdLabelH[nid]=str(ragaId)
            NIdColorH[nid]=colors[str(ragaId)]
            #G.AddIntAttrDatN(nid,int(ragaId),'RagaID');
            
    except psy.DatabaseError, e:
        print 'Error %s' % e
        if con:
            con.rollback()
            con.close()
        sys.exit(1)
    
    if con:
        con.close()
    
    snap.SavePajek_PNEANet(G, outputNetworkFile)
    #snap.SavePajek_PNEANet(G, outputNetworkFile, NIdColorH,NIdLabelH)
    
def constructNetworkSNAP(fileListFile, outputNetworkFile, thresholdBin, pattDistExt):
    """
    This function uses Snap library to generate a network. 
    :fileListFile: file which contains list of filenames which are to be used for the network generation
    :outputNetworkFile: output file in which we want to save our network
    :thresholdBin: bin index based on which a threshold will be applied on the distance between the patterns. This mapping between the bin and the distance threshold is same across every module, including the module for computing clustering coefficient.
    :pattDistExt: file extension of the file which has pattern distances 
    """    
    
    #initialize an graph
    G = snap.PNEANet.New()
    
    # creating the network
    filenames = open(fileListFile).readlines()
    for filename in filenames:
        fname = filename.strip() + pattDistExt
        data = np.loadtxt(fname)
        for ii in range(data.shape[0]):
            id1 = int(data[ii][0])
            id2 = int(data[ii][1])
            dist = float(data[ii][2])
            
            if dist >= ThshldArray[thresholdBin-1]:
                continue
            if not G.IsNode(id1):
                G.AddNode(id1)
            if not G.IsNode(id2):
                G.AddNode(id2)
            if not G.IsEdge(id1,id2):
                G.AddEdge(id1,id2)
    
    #saving network as a pajek file
    snap.SavePajek_PNEANet(G, outputNetworkFile)
    

    
def constructNetwork_Weighted_NetworkX(fileListFile, outputNetworkFile, thresholdBin, pattDistExt, wghtMethod, DISPARITY_FILTER_SIGNIF_LEVEL):
    """
    This function uses networkX library to generate a weighted network.
    
    :fileListFile:  file which contains list of filenames which are to be used for the network generation
    :outputNetworkFile: output file in which we want to save our network (as pajeck)
    :thresholdBin:  bin index based on which a distance threshold will be applied on the distance between 
                    the patterns. This mapping between the bin and the distance threshold is same across 
                    every module, including the module for computing clustering coefficient. 
                    The threshold values are found to be an optimal range for this problem.
    :pattDistExt:   file extension for the file where we store pattern distances 
    :wghtMethod:    0 -> wght = 1/distance, 1 -> wght = exp(-distance/mean(distances)), -1->unweighted
    :DISPARITY_FILTER_SIGNIF_LEVEL: significance level of the disparity filtering. If negative no filtering is performed
    """
    
    if thresholdBin > 0:
        distanceThsld = ThshldArray[thresholdBin-1]
    else:
        distanceThsld = FLT_MAX
        
    #initialize an graph
    G=nx.Graph()
    meanDist=1
    #estimating mean distance value in case wghtMethod=1
    if wghtMethod == 1:
        distances = 0
        n_distances = 0
        filenames = open(fileListFile).readlines()
        for ii, filename in enumerate(filenames):
            print "Processing file for average distance %d\n"%(ii+1)
            fname = filename.strip() + pattDistExt
            data = np.loadtxt(fname)
            ind_valid = np.where(data[:,2] < distanceThsld)[0]
            distances += np.sum(data[ind_valid,2])
            n_distances += len(ind_valid)
        meanDist = distances/float(n_distances)

    # creating the network
    filenames = open(fileListFile).readlines()
    for ii, filename in enumerate(filenames):
        print ii, filename
        print "processing file %d\n"%(ii+1)
        fname = filename.strip() + pattDistExt
        data = np.loadtxt(fname)
        ind_valid = np.where(data[:,2] < distanceThsld)[0]
        
        if wghtMethod == 0:
            edge_list = [(int(data[ii][0]), int(data[ii][1]), float(1.0/data[ii][2])) for ii in ind_valid]
        elif wghtMethod == 1:
            edge_list = [(int(data[ii][0]), int(data[ii][1]), float(np.exp(-data[ii][2]/meanDist))) for ii in ind_valid]
        elif wghtMethod == -1:
            edge_list = [(int(data[ii][0]), int(data[ii][1]), 1.0) for ii in ind_valid]
        
        G.add_weighted_edges_from(edge_list)
        
    if DISPARITY_FILTER_SIGNIF_LEVEL > 0:
        G = netPro.filter_graph_edges(G,DISPARITY_FILTER_SIGNIF_LEVEL)
        
    #saving weighted network/graph as a pajek file
    if not outputNetworkFile == None:
        nx.write_pajek(G, outputNetworkFile)
    else:
        return G
    

        
    
