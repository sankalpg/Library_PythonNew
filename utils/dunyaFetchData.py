from carnatic import models
from docserver import util
import os, sys, shutil
import numpy as np

def obtainRagaRecordingMapping():
    
    ragas = models.Raaga.objects.all()
    
    ragaRecordings={}
    
    for raga in ragas:
        ragaRecordings[raga.name]=[]
        
        recordings = raga.recordings()
        for recording in recordings:
            ragaRecordings[raga.name].append(recording.mbid)
    
    return ragaRecordings



def obtainPitchTonicForRaga(raganame, outputfolder):


        if not os.path.exists(outputfolder):
            os.makedirs(outputfolder)
        
        #get raga recording mapping
        ragaRecordings = obtainRagaRecordingMapping()
        print len(ragaRecordings.keys())
        
        for mbid in ragaRecordings[raganame]:
            pitchData = util.docserver_get_json(mbid,'pitch','pitch')
            tonic = float(util.docserver_get_contents(mbid,'ctonic','tonic'))
            HopSize = 196
            
            filename =util.docserver_get_filename(mbid,'mp3')
            
            shutil.copyfile(filename, outputfolder+ '/'+filename.split('/')[-1])
            
            filename = filename.split('/')[-1].split('.')[0]
            
            #TStamps = np.array(range(0,len(pitchData)))*np.float(HopSize)/44100.0
            
            #dump = np.array([TStamps, pitchData]).transpose()
            
            np.savetxt(outputfolder+ '/'+filename+'.pitch', pitchData, delimiter = "\t")
            
            np.savetxt(outputfolder+ '/'+filename+'.tonic', np.array([tonic]), delimiter = "\t")
            
            
            
            
