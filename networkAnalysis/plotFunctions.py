import numpy as np
import os,sys
import matplotlib.pyplot as plt


baseClusterCoffFileName = '_ClusteringCoff'
randomizationSuffix = '_RANDOM'
baseNetworkPropFileName = '_NetworkInfo'
    
def readClusterinCoffCurve(root_dir, nFiles):
    """
    This function just reads multiple files and accumulate CC values and returns it
    """
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
    
    return CC, CC_Rand
    

def plotClusteringCoff(root_dir, nFiles, plotName=-1, legData = []):
    """
    This function plots clustering cofficient as a function of threshold using which the network was build. It also plots the CC corresponding to the randomized network.
    
    Example 
    plt.plotClusteringCoff(['networkData/pasaNorm/weight0/conf200.00/', 'networkData/pasaNorm/weight0/conf99.90/', 'networkData/pasaNorm/weight0/conf99.00/', 'networkData/pasaNorm/weight0/conf90.00/', 'networkData/pasaNorm/weight0/conf80.00/', 'networkData/pasaNorm/weight0/conf50.00/'], 25, legData=['NoFilt', 'Conf99.9', 'Conf99', 'Conf90', 'Conf80','Conf50'], plotName='networkData/pasaNorm/weight0/CC_Curves.pdf')
    """
    
    CC = []
    CC_Rand = []
    legendArr= []
    if isinstance(root_dir,str):
        cc, cc_rand = readClusterinCoffCurve(root_dir, nFiles)
        CC.append(cc)
        CC_Rand.append(cc_rand)
        legendArr.append('M%d_orig'%1)
        legendArr.append('M%d_rand'%1)
    elif isinstance(root_dir,list):
        for ii, r in enumerate(root_dir):
            cc, cc_rand = readClusterinCoffCurve(r, nFiles)
            CC.append(cc)
            CC_Rand.append(cc_rand)
            if len(legData)-1>=ii:
                legendArr.append('%s($G$)'%legData[ii])
                legendArr.append('%s($G_r$)'%legData[ii])
            else:
                legendArr.append('M%d'%ii)
                legendArr.append('M%d'%ii)
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.hold(True)
    fsize = 22
    fsize2 = 16
    font="Times New Roman"
    
    plt.xlabel("Threshold", fontsize = fsize, fontname=font)
    plt.ylabel("Clustering Coefficient", fontsize = fsize, fontname=font, labelpad=fsize2)
    
    
    pLeg = []
    
    markers = ['.', 'o', 's', '^', '<', '>', 'p']    
    colors = ['r', 'b', 'm', 'c', 'g', 'k']
    colors_dotted = ['r--', 'b--', 'm--', 'c--', 'g--', 'k--']
    for ii, cc in enumerate(CC):
        p, = plt.plot(cc, colors[ii], linewidth=2, marker = markers[ii])
        pLeg.append(p)
        p, = plt.plot(CC_Rand[ii], colors_dotted[ii], linewidth=2, marker = markers[ii])
        pLeg.append(p)
    
    ax.set_ylim([0,0.4])
    ax.set_xlim([0,25])
    
    xLim = ax.get_xlim()
    yLim = ax.get_ylim()
    
    ax.set_aspect((xLim[1]-xLim[0])/(2*float(yLim[1]-yLim[0])))
    plt.legend(pLeg, legendArr, ncol = 6, fontsize = 12, scatterpoints=1, frameon=True, borderaxespad=0.6, bbox_to_anchor=(-0.08, 1.02, 1., .102), loc=3, columnspacing=0.5, handletextpad=0.1)
    plt.tick_params(axis='both', which='major', labelsize=fsize2)
    
    
    if isinstance(plotName, int):
        plt.show()
    elif isinstance(plotName, str):
        fig.savefig(plotName, bbox_inches='tight')
