import numpy as np
from scipy.io import loadmat
import string

NULL_TOKEN = None

# LEVELS = range(-500, 2701, 100)
# CODES = string.ascii_letters[:len(LEVELS)]
LEVELS = range(-800, 2501, 100)
NOTES = ['G','m','M','P','d','D','n','N','S','r','R','g']*3
CODES = string.ascii_letters[:len(NOTES)]


class Data:
    def __init__(self, pdata, thres=8, width=35, verbose=False):
        #self.songname = songname
        vals = []
        #with open('metadata.txt') as f:
            #for line in f:
                #vals = line.strip().split(',')
                #if vals[1] == songname:
                    #break
        #if len(vals) == 0:
            #raise Exception('Song not found in metadata.txt')
        #self.index = int(vals[0])
        #self.tonic = float(vals[2])
        #self.raga_index = int(vals[3])
        self.start = int(0)
        self.end = int(300)
        vals = 'RgmPdn'
        print vals
        self.ignore = ''.join([c for n,c in zip(NOTES,CODES) if n in vals])
        self.pdata = pdata
        t, s = self.get_timeSeries()
        self.times = t
        self.ts = s
        sym, st, en = self.get_symbols(thres=thres, width=width, verbose=verbose)
        self.string = sym
        self.st = st
        self.en = en
        note = self.get_notes(self.string)
        self.transcribed = note
        qts = self.get_quantized_ts()
        self.levels = qts
        #label, gst, gen = self.get_gt()
        #self.label = label
        #self.gst = gst
        #self.gen = gen

    def get_raw_ts(self):
        with open('data_pitch/'+self.songname+'.tpe') as f:
            for line in f:
                yield map(float, line.strip().split())[:2]
            
    def get_ts(self):
        data = loadmat('data_pitch/'+self.songname+'.mat')['output']
        st, en = self.start, self.end
        t, s = data[st*100:en*100], data[st*100:en*100]
        return t, s
    
    def get_timeSeries(self):
        #print data, data.shape
        st, en = self.start, self.end
        print st, en
        nSamp = np.ceil(en/float(self.pdata[2]))
        print nSamp
        #t,s = 0,0
        t, s = self.pdata[0][:nSamp], self.pdata[1][:nSamp]
        return t, s
        
    @staticmethod
    def _quantize(val, ignore='', width=35, verbose=False):
        if np.isnan(val):
            return None
        if verbose and (val < LEVELS[0] or val > LEVELS[-1]) :
            print 'Warning: Time series value %f out of quantization range' % (val,)
        ind = [l for l in LEVELS if abs(l-val) < width]
        if len(ind) == 0:
            return NULL_TOKEN
        val = CODES[LEVELS.index(ind[0])]
        if val in ignore:
            return None
        return val


    def get_symbols(self, thres=10, width=35, verbose=False):
        '''
        Returned values t_st and t_en contain relative indexes from the 
        time series corresponding to beginning of clipped song  
        instead of absolute time in seconds
        '''
        # print 'Notes %s being ignored for song %s' % (set(self.get_notes(self.ignore)), self.songname)
        symbs = map(lambda x:self._quantize(x,self.ignore, width, verbose), self.ts)
        t_st, t_en = [0], []
        prev, count = symbs[0], 0
        symbols = []
        for i, s in enumerate(symbs):
            if s == prev:
                count += 1
                continue
            if prev is None or count < thres:
                t_st[-1] = i
            else:
                symbols.append(prev)
                t_en.append(i-1)
                t_st.append(i)
            count = 1
            prev = s
        else:
            t_st.pop()
        # combine duplicate symbols in second pass
        prev = ''
        index = 0
        for i, s in enumerate(symbols):
            if s != prev:
                prev = s
                symbols[index] = s
                t_st[index] = t_st[i]
                t_en[index] = t_en[i]
                index += 1
            else:
                t_en[index-1] = t_en[i]
        symbols[index:] = []
        t_st[index:] = []
        t_en[index:] = []
        string = ''.join(symbols)
        return string, t_st, t_en
    
    
    
    def get_quantized_ts(self):
        qts = np.array([None]*(max(self.en)+1))
        for i in xrange(len(self.st)):
            qts[self.st[i]:self.en[i]+1] = LEVELS[CODES.index(self.string[i])]
        return qts
            
    def get_gt(self):
        '''
        Returns the start and end time of ground truth phrases in seconds
        along with the label representing the annotated name of the phrase
        '''
        with open('data_annotation/' + self.songname + '.txt') as f:
            label, count = f.readline().strip().split()
            count = int(count)
            st, en = [0.0]*count, [0.0]*count
            for i in xrange(count):
                s, e = map(float, f.readline().strip().split())
                st[i], en[i] = s, e
            return label, st, en
            
    def get_string_part(self, beg, fin):
        ''' beg and fin are the starting and ending time instants
        respectively. This function looks at the boundaries of all
        symbols and figures out which symbols lie within the range.
        This function is used to calculate the string representation
        of a section of time series.
        '''
        beg, fin = np.min(np.where(self.times>=beg)), np.min(np.where(self.times>=fin))
        try:
            s = self.st.index(max(x for x in self.st if x <= beg))
        except:
            s = 0
        try:
            e = self.en.index(min(x for x in self.en if x >= fin))
        except:
            e = len(self.en) - 1 
        return self.string[s:e+1], s, e 
    
    @staticmethod
    def get_notes(string):
        return ''.join(NOTES[CODES.index(c)] for c in string)
    
    def get_ts_index(self, t):
        '''Gives the index in the time series for a corresponding time instant'''
        for i, ti in enumerate(self.times):
            if ti > t: return i
        return len(self.times)-1
    


    def get_precall(self, predicted, overlap = 0.25):
        if len(predicted) == 0:
            return 0.0, 0.0, 0.0
        gst = map(self.get_ts_index, self.gst)
        gen = map(self.get_ts_index, self.gen)
        predicted = map(lambda (x,y):(self.st[x],self.en[y]), predicted)
        pst, pen = zip(*predicted)
    
        overlap_mat = np.zeros([len(pst),len(gst)])
        for i in xrange(len(pst)):
            for j in xrange(len(gst)):
                s1, e1, s2, e2 = pst[i], pen[i], gst[j], gen[j]
                overlap_mat[i,j] = 1.0*max(0, min(e1,e2) - max(s1,s2) + 1)/(e2 - s2 + 1)
        # print overlap_mat
        precision = np.mean(np.any(overlap_mat>=overlap, axis=1))
        recall = np.mean(np.any(overlap_mat>=overlap, axis=0))
        f_score = (2*precision*recall)/(precision+recall)
        return precision, recall, f_score
    
    def get_behavioral_seq(self, win=3):
        behav_dict = {}
        seq = [0]*(len(self.string) - win + 1)
        next = 0
        for i in xrange(len(self.string)-win+1):
            if self.string[i:i+win] not in behav_dict:
                behav_dict[self.string[i:i+win]] = next
                next += 1
            seq[i] = behav_dict[self.string[i:i+win]]
        return seq, behav_dict
                

if __name__=='__main__':
    pass
