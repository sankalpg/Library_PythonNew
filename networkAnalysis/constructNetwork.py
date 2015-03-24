import sys
import os
import numpy as np
import snap
sys.path.append(os.path.join(os.path.dirname(__file__), '../../batchProcessing'))

import batchProcessing as BP
import psycopg2 as psy

myUser = 'sankalp'

ISMIR2015_10RAGASCarnatic = ['55', '10', '27', '8', '9', '137', '159', '20', '13', '210']
colors = {'55':'red', '10':'blue', '27':'black', '8':'green', '9':'yellow', '137':'pink', '159':'orange', '20':'magenta', '13':'cyan', '210':'purple'}

ThshldArray = np.array([0.0,
 208.247,
 832.986,
1874.219,
3331.945,
5206.164,
7496.876,
10204.082,
13327.780,
16867.972,
20824.656,
25197.834,
29987.505,
35193.669,
40816.327,
46855.477,
53311.120,
60183.257,
67471.887,
75177.010,
83298.626,
91836.735,
100791.337,
110162.432,
119950.021,
130154.102,
140774.677,
151811.745,
163265.306,
175135.360,
187421.908,
200124.948,
213244.481,
226780.508,
240733.028,
255102.041,
269887.547,
285089.546,
300708.038,
316743.024,
333194.502,
350062.474,
367346.939,
385047.897,
403165.348,
421699.292,
440649.729,
460016.660,
479800.083,
500001.000])

def constructNetworkCollection(fileListFile, outputNetworkFile, thresholdBin, pattDistExt = '.pattDistances', myDatabase=''):
    
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
            
            if dist >= ThshldArray[thresholdBin-1]*10:  # this *10 is temp TODO: remove it asap. This is to make the rsults bit exact with the C code
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
    
    
        
        
    
