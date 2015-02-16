import json, os, sys
import numpy as np
import compmusic
from compmusic import dunya as dn
from compmusic.dunya import hindustani as hn
from compmusic.dunya import docserver as ds

sys.path.append(os.path.join(os.path.dirname(__file__), '../batchProcessing'))

import batchProcessing as BP



def fetchHindustaniCollection(out_dir):
    
    dn.set_token("60312f59428916bb854adaa208f55eb35c3f2f07")
    
    #fetch all the recordings in the collection
    recordings = hn.get_recordings()
    
    
    for ii, rec in enumerate(recordings):
        print "processing %d out of %d\n"%(ii+1, len(recordings))
        
        recording = rec['mbid']
        recInfo = hn.get_recording(recording)
        
        recName = recInfo['title']
        releaseName = recInfo['release']['title']
        artist = recInfo['artists'][0]['name']
        
        dirName = os.path.join(out_dir, artist, releaseName)
        fileName = os.path.join(out_dir, artist, releaseName, recName)
        fileName_mp3 = fileName + '.mp3'
        fileName_tonic = fileName + '.tonic'
        fileName_pitch = fileName + '.pitch'
        
        
        if not os.path.exists(dirName):
            os.makedirs(dirName)
        
        #download mp3 file
        if not os.path.isfile(fileName_mp3):
            try:
                mp3File = ds.get_mp3(recording)
                fid = open(fileName_mp3, "w")
                fid.write(mp3File)
                fid.close()
            except:
                print "couldn't fetch mp3 file: %s\n"%fileName
        
        #download tonic information
        if not os.path.isfile(fileName_tonic):
            try:
                tonic = json.loads(ds.file_for_document(recording, 'ctonic', 'tonic'))
                np.savetxt(fileName_tonic, np.array([tonic]), fmt='%.3e', delimiter='\t')
            except:
                print "couldn't fetch tonic data: %s\n"%fileName

        #download pitch of the recording
        if not os.path.isfile(fileName_pitch):
            try:
                pitch = json.loads(ds.file_for_document(recording, 'pitch', 'pitch'))
                np.savetxt(fileName_pitch, pitch, fmt='%.3e', delimiter='\t')
            except:
                print "couldn't fetch pitch data: %s\n"%fileName       
                
                
def generateDummyTaniSegments(root_dir):
    
    files = BP.GetFileNamesInDir(root_dir, '.mp3')
    
    for f in files:
        fname, ext = os.path.splitext(f)
        fid = open(fname + '.tablaSec', "w")
        fid.close()
        
        
        
        
        
        
        
