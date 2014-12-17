import numpy as np
import operator

class piecewise_linear_segmenter:
	def __init__(self,x):
		self.reinit(x)
		return None
	# Init
	def reinit(self,x):
		# ---- Setting ----
		#self.fitting='regression'
		self.fitting='extremes'
		#self.errormeasure='max'
		self.errormeasure='mean'
		# -----------------
		# Copy data
		self.x=x.copy()
		self.xi=np.array(range(len(x)))
		self.merge_historic=[]
		# Init segments and costs
		self.segment=[]
		for i in range(0,len(self.x)-1,2):
			self.segment.append(self._create_segment(i,i+1))
		self.merge_cost=[]
		for i in range(len(self.segment)-1):
			self.merge_cost.append(self._error(self._merge(i,i+1)))
		# Estimate noise in x using the derivative
		# [[J. Rice (1984), Bandwidth choice for nonparametric regression, Ann. Stat 12, 1215-1230.]]
		#self.sigma_noise=np.sqrt(np.sum(np.diff(self.x)**2)/float(2*(len(x)-1)))
		#self.sigma_noise=np.std(self.x)
		#foo,self.sigma_noise=medmad(self.x)
		self.sigma_noise=np.mean(np.abs(np.diff(self.x)))
		# Done
		return None
	# Segmentation
	def do_segmentation(self,tol=3,absolute_max_error=None,force_min_len=None):
		if self.sigma_noise==0: return False
		# Init max_error and min_len (operation mode selection)
		max_error=float(tol)*self.sigma_noise
		if absolute_max_error!=None: max_error=absolute_max_error
		min_len=0
		if force_min_len!=None: min_len=force_min_len
		# Iterate
		p,min_cost=min(enumerate(self.merge_cost),key=operator.itemgetter(1))
		while min_cost<max_error or self._segment_min_len()<min_len:
			self.merge_historic.append([self.segment[p]['i'],self.segment[p]['j'],self.segment[p+1]['i'],self.segment[p+1]['j']])
			self.segment[p]=self._merge(p,p+1)
			del self.segment[p+1]
			if len(self.segment)==1: break
			if p+1<len(self.merge_cost):
				self.merge_cost[p]=self._error(self._merge(p,p+1))
				del self.merge_cost[p+1]
			else:
				del self.merge_cost[p]
			if p>0:
				self.merge_cost[p-1]=self._error(self._merge(p-1,p))
			p,min_cost=min(enumerate(self.merge_cost),key=operator.itemgetter(1))
		# Done
		return True
	# Return reconstruction
	def get_reconstruction(self,mark_transitions=None):
		xhat=self.x.copy()
		for segment in self.segment:
			xi=self.xi[segment['i']:segment['j']+1]
			xhat[xi]=segment['xhat']
			if mark_transitions!=None:
				xhat[xi[0]]=mark_transitions
				xhat[xi[-1]]=mark_transitions
		return xhat
	# Return segments
	def get_segments(self):
		return self.segment[:]
	# Return merge costs
	def get_merge_costs(self):
		return self.merge_cost[:]
	def get_merge_history(self):
		return self.merge_historic[:]
	# ------------------------------------
	# Current minimum segment length
	def _segment_min_len(self):
		mn=np.infty
		for segment in self.segment:
			if len(segment['xhat'])<mn: mn=len(segment['xhat'])
		return mn
	# Create segment
	def _create_segment(self,i,j,initing=False):
		segment={}
		segment['i']=i
		segment['j']=j
		xi=self.xi[segment['i']:segment['j']+1]
		x=self.x[segment['i']:segment['j']+1]
		if self.fitting=='regression':
			segment['coefs']=np.polyfit(xi,x,1)
		elif self.fitting=='extremes':
			slope=(x[-1]-x[0])/(xi[-1]-xi[0])
			segment['coefs']=np.array([slope,x[0]-slope*xi[0]])
		else: self._die()
		segment['xhat']=np.polyval(segment['coefs'],xi)
		return segment
	# Merge segment
	def _merge(self,i,j):
		return self._create_segment(self.segment[i]['i'],self.segment[j]['j'])
	# Error calculation function
	def _error(self,segment):
		xi=self.xi[segment['i']:segment['j']+1]
		if self.errormeasure=='max': return np.max(np.abs(self.x[xi]-segment['xhat']))
		elif self.errormeasure=='mean': return np.mean(np.abs(self.x[xi]-segment['xhat']))
		else: self._die()
		return None
	# Die
	def _die(self):
		print 'Error in piecewise linear segmentation.'
		sys.exit()
		return None
	# ------------------------------------

