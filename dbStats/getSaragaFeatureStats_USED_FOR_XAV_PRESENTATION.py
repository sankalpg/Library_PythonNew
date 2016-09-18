import os, sys
from compmusic import dunya as dn
from compmusic.dunya import hindustani as hn
from compmusic.dunya import carnatic as ca
from compmusic.dunya import docserver as ds
from compmusic import musicbrainz


def GetFileNamesInDir(dir_name, filter=".wav"):
    names = []
    for (path, dirs, files) in os.walk(dir_name):
        for f in files:
            #if filter in f.lower():
            if filter.split('.')[-1].lower() == f.split('.')[-1].lower():
                #print(path+"/"+f)
                #print(path)
                #ftxt.write(path + "/" + f + "\n")
                names.append(path + "/" + f)
                
    return names


def getFileCounts(root_dir, exts = []):

	for ext in exts:
		filenames = GetFileNamesInDir(root_dir, ext)
		print "Files with ext: %s, are %d in numbers"%(ext, len(filenames))

def getFeatureCount(root_dir, exts = []):

	for ext in exts:
		filenames = GetFileNamesInDir(root_dir, ext)
		cnt = 0
		for filename in filenames:
			lines = open(filename, 'r').readlines()
			cnt+= len(lines)
		print "Files with ext: %s, have %d numbers of annotations"%(ext, cnt)		

	

# Read file was used to get stats for Xavier's presentation in Bombay, we just coundted stuff, and no seeding from musicbrainz was done!