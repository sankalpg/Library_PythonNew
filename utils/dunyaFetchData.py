import os, sys, shutil
import numpy as np
from mutagen import easyid3
import requests



sys.path.append(os.path.join(os.path.dirname(__file__), '../batchProcessing'))

import batchProcessing as BP


def obtainPitchTonicInDB(root_dir):

        audiofiles = BP.GetFileNamesInDir(root_dir,'mp3')
        cnt=0
        for audiofile in audiofiles:
            fname, ext = os.path.splitext(audiofile)
            
            pithcFile = fname + '.pitch'
            tonicFile = fname + '.tonic'

            if not os.path.exists(pithcFile):
                
                try:
                    mbid = easyid3.ID3(audiofile)['UFID:http://musicbrainz.org'].data
                    pitchUrl = 'http://dunya.compmusic.upf.edu/document/by-id/' + mbid + '/pitch?v=0.4&subtype=pitch'
                    tonicUrl = 'http://dunya.compmusic.upf.edu/document/by-id/' + mbid + '/ctonic?v=0.2&subtype=tonic'
                    pitchData = np.array(eval(requests.get(pitchUrl).content))
                    tonicData = float(requests.get(tonicUrl).content)
                    np.savetxt(pithcFile, pitchData, delimiter = "\t", fmt='%0.2f')
                    np.savetxt(tonicFile, np.array([tonicData]), delimiter = "\t" , fmt='%0.2f')
                except:
                    cnt=cnt+1
                    print audiofile

        print cnt




