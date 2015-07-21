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

def filter_graph_edges(G, DISPARITY_FILTER_SIGNIF_LEVEL, verbose=True, print_prefix='', field='weight'):
    '''
    A large number of complex systems find a natural abstraction in the form of weighted networks whose nodes represent
    the elements of the system and the weighted edges identify the presence of an interaction and its relative strength.
    In recent years, the study of an increasing number of large-scale networks has highlighted the statistical
    heterogeneity of their interaction pattern, with degree and weight distributions that vary over many orders of
    magnitude. These features, along with the large number of elements and links, make the extraction of the truly
    relevant connections forming the network's backbone a very challenging problem. More specifically, coarse-graining
    approaches and filtering techniques come into conflict with the multiscale nature of large-scale systems. Here, we
    define a filtering method that offers a practical procedure to extract the relevant connection backbone in complex
    multiscale networks, preserving the edges that represent statistically significant deviations with respect to a
    null model for the local assignment of weights to edges. An important aspect of the method is that it does not
    belittle small-scale interactions and operates at all scales defined by the weight distribution. We apply our
    method to real-world network instances and compare the obtained results with alternative backbone
    extraction techniques. (http://www.pnas.org/content/106/16/6483.abstract)
    '''

    if verbose:
        print '%sFiltering with ' % print_prefix + str(100*(1-DISPARITY_FILTER_SIGNIF_LEVEL))+'% confidence ...',

    # FOR DIRECTED
    #indegree = G.in_degree(weight=None)
    #outdegree = G.out_degree(weight=None)
    #instrength = G.in_degree(weight='weight')
    #outstrength = G.out_degree(weight='weight')

    # FOR UNDIRECTED
    degree = G.degree(weight=None)
    strength = G.degree(weight=field)

    edges = G.edges()
    for i, j in edges:
            # FOR DIRECTED
            #pij = float(G[i][j]['weight'])/float(outstrength[i])
            #pji = float(G[i][j]['weight'])/float(instrength[j])
            #aij = (1-pij)**(outdegree[i]-1)
            #aji = (1-pji)**(indegree[j]-1)
            #if aij < DISPARITY_FILTER_SIGNIF_LEVEL or aji < DISPARITY_FILTER_SIGNIF_LEVEL:
            #    continue

            # FOR UNDIRECTED
            pij = float(G[i][j][field])/float(strength[i])
            aij = (1-pij)**(degree[i]-1)
            if aij < DISPARITY_FILTER_SIGNIF_LEVEL:
                continue

            G.remove_edge(i, j)
    nodes = G.nodes()
    for n in nodes:
        if G.degree(n) < 1:
            G.remove_node(n)
    if verbose:
        print G.number_of_nodes(), 'nodes,', G.number_of_edges(), 'edges'

    return G

def filterAndSaveNetwork(networkFile, outputNetworkFile, DISPARITY_FILTER_SIGNIF_LEVEL):
    """
    This function filters the network and saves it in another file for further analysis
    
    """
    
    G = nx.read_pajek(networkFile)
    
    #converting to undirected network
    G = nx.Graph(G)
    
    if DISPARITY_FILTER_SIGNIF_LEVEL>0:
        G = filter_graph_edges(G,DISPARITY_FILTER_SIGNIF_LEVEL)
    
    #saving weighted network/graph as a pajek file
    nx.write_pajek(G, outputNetworkFile)
    

def detectCommunitiesInNewtworkNX(networkFile, outputFile):
    """
    This function detects communities in the network using "community" module for networkX.
    
    :networkFile: network file which is saved as a pajek file
    :outputFile: file in whicn community information is written
    """    
    G = nx.read_pajek(networkFile)
    
    #converting to undirected network
    G = nx.Graph(G)
    
    partition = community.best_partition(G)
    comm_data = {}
    for com in set(partition.values()):
        if not comm_data.has_key(int(com)):
            comm_data[int(com)] = []
        list_nodes = [nodes for nodes in partition.keys()
                                if partition[nodes] == com]
        for n in list_nodes:
            comm_data[int(com)].append({'nId': int(n), 'ragaId':'', 'mbid':'', 'compId':''})
    
    json.dump(comm_data, open(outputFile, 'w'))
    
    
def attachRagaLabelToNodes(communityFile, outputFile, myDatabase = ''):
    """
    Temporary function
    """
    #annotating netowrk nodes with raga labels and edges with distances
    cmd1 = "select raagaId from file where id = (select file_id from pattern where id =%d)"
    data = np.loadtxt(communityFile)
    fid = open(outputFile,'w')
    
    try:
        con = psy.connect(database=myDatabase, user=myUser) 
        cur = con.cursor()
        print "Successfully connected to the server"
        output = []
        for ii in range(data.shape[0]):
            nodeId = data[ii,0]
            cur.execute(cmd1%(int(nodeId)))
            ragaId = cur.fetchone()[0]
            fid.write("%d\t%d\t%s\n"%(nodeId, data[ii,1], ragaId))
            
    except psy.DatabaseError, e:
        print 'Error %s' % e
        if con:
            con.rollback()
            con.close()
        sys.exit(1)
    
    if con:
        con.close()
    
    fid.close()
        


def computeBetweenesCentrality(netFile, outFile):
    
    G = nx.read_pajek(netFile)
    G = nx.Graph(G)
    bc = nx.betweenness_centrality(G)
    json.dump(bc, open(outFile,'w'))
    
    
def computeDegreeDistribution(netFile, outFile):
    
    G = nx.read_pajek(netFile)
    G = nx.Graph(G)
    nodes = G.nodes()
    d = G.degree(nodes)
    json.dump(d, open(outFile,'w'))    


if __name__=="__main__":    
#unit test stuff
# community detection
    detectCommunitiesInNewtworkNX('unitTests/Weight_-1___DTshld_10___PostFilt___C.edges', 'unitTests/Weight_-1___DTshld_10___PostFilt___C.community')
    