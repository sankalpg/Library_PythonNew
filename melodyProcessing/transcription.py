import numpy as np
from scipy.io import loadmat
import string

NULL_TOKEN = None

# LEVELS = range(-500, 2701, 100)
# CODES = string.ascii_letters[:len(LEVELS)]
LEVELS = range(-800, 2501, 100)
NOTES = ['G','m','M','P','d','D','n','N','S','r','R','g']*3
CODES = string.ascii_letters[:len(NOTES)]


class SvarTranscription:
    
    def __init__(self, pdata):
        vals = []
        self.pdata = pdata
        self.start = int(0)
        self.end = int(len(self.pdata[1])*self.pdata[2])
        
        
    def perform_transcription(self, ignoreNotes, thres=8, width=35, verbose=False):
        self.ignore = ignoreNotes
        t, s = self.get_timeSeries()
        self.times = t
        self.ts = s
        sym, st, en = self.get_symbols(thres=thres, width=width, verbose=verbose)
        self.string = sym
        self.st = st
        self.en = en
        note = self.get_notes(self.string)
        self.transcribed = note
        level = self.convert_to_levels(self.string)
        self.svara = level
        #qts = self.get_quantized_ts()
        #self.levels = qts
        output = []
        for ii, start_time in enumerate(self.st):
            output.append([self.st[ii]*self.pdata[2], self.en[ii]*self.pdata[2], self.svara[ii]])
        return output
      
    
    def get_timeSeries(self):
        st, en = self.start, self.end
        print "Getting time series..."
        nSamp = np.ceil(en/float(self.pdata[2]))
        #print nSamp
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
        print "Transcribing melody..."
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
            #print s, prev, t_st[i], t_en[i], t_en[i] - t_st[i], t_st[i] - t_en[index-1], t_en[i] - t_en[index-1], t_en[index-1]
            if s != prev:
                prev = s
                symbols[index] = s
                t_st[index] = t_st[i]
                t_en[index] = t_en[i]
                index += 1                   
            elif t_st[i] - t_en[index-1] > 120:
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
    
    
    @staticmethod
    def convert_to_levels(string):
        level = []
        for i in xrange(len(string)):
            level.append(LEVELS[CODES.index(string[i])])
        return level    
    
    
    def get_ts_index(self, t):
        '''Gives the index in the time series for a corresponding time instant'''
        for i, ti in enumerate(self.times):
            if ti > t: return i
        return len(self.times)-1
    


if __name__=='__main__':
    pass
