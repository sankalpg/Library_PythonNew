import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import networkx as net

def plotClusteringCoff(root_dir, nFiles, plotName=-1):
    """
    This function plots clustering cofficient as a function of threshold using which the network was build. It also plots the CC corresponding to the randomized network.
    """
    
    baseClusterCoffFileName = '_ClusteringCoff'
    randomizationSuffix = '_RANDOM'
    baseNetworkPropFileName = '_NetworkInfo'
    CC = []
    CC_Rand = []    
    for ii in range(1, nFiles+1):

        try:
            cc = np.loadtxt(os.path.join(root_dir, str(ii)+baseClusterCoffFileName+'.txt'))
            cc_rand = np.loadtxt(os.path.join(root_dir,str(ii)+baseClusterCoffFileName+randomizationSuffix+'.txt'))
        except:
            cc = 0
            cc_rand = 0
         
        CC.append(cc)
        CC_Rand.append(cc_rand)
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.hold(True)
    fsize = 22
    fsize2 = 14
    font="Times New Roman"
    
    plt.xlabel("threshold bin", fontsize = fsize, fontname=font)
    plt.ylabel("Clustering Coff", fontsize = fsize, fontname=font, labelpad=fsize2)
    
    
    pLeg = []
    p, = plt.plot(CC, 'r', linewidth=2)
    pLeg.append(p)
    p, = plt.plot(CC_Rand, 'b--', linewidth=2)
    pLeg.append(p)
    
    xLim = ax.get_xlim()
    yLim = ax.get_ylim()
    
    ax.set_aspect((xLim[1]-xLim[0])/(2*float(yLim[1]-yLim[0])))
    plt.legend(pLeg, ['Original Network', 'Randomized Network'], loc ='lower right', ncol = 1, fontsize = fsize2, scatterpoints=1, frameon=True, borderaxespad=0.1)
    #plt.tick_params(axis='both', which='major', labelsize=fsize2)
    
    
    if isinstance(plotName, int):
        plt.show()
    elif isinstance(plotName, str):
        fig.savefig(plotName)
    
    
    


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
            pij = float(G[i][j][0][field])/float(strength[i])
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


def computeRaagaClusters(networkFile):
    
    #reading network
    G = networkx.read_pajek(networkFile)

    



    
    