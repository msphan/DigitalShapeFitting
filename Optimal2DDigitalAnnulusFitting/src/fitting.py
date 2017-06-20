import timeit
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

class DigitalAnnulus:
	""" 
	Digital annulus class
	"""
	def __init__(self,c,r,w):
		""" 
		initialize a digital annulus with center c, radius r and width w 
		"""
		self.c = c
		self.r = r
		self.w = w

	def build_data(self, I, J):
		""" 
		build fitting data function with
			- I : number of inliers
			- J : number of outliers
		"""
		cx, cy, r, w = self.c[0], self.c[1], self.r, self.w
		k, inliers = I, []
		while(k > 0):
			x = np.random.uniform(cx-r,cx+r,k)
			y = np.random.uniform(cy-r,cy+r,k)
			inliers = list(set(inliers+\
						[tuple(t for t in (i,j)) for (i,j) in zip(x,y) 
							if (i - cx)**2+(j-cy)**2 >= (r - w)**2 and 
						  		    (i - cx)**2+(j-cy)**2 <= r**2]))
			k = I - len(inliers)

		self.inliers = np.array([list(i) for i in inliers])

		k, outliers = J, []
		while(k > 0):
			x = np.random.uniform(0,1,k)
			y = np.random.uniform(0,1,k)	
			outliers = list(set(outliers+\
						[tuple(t for t in (i,j)) for (i,j) in zip(x,y) 
							if (i - cx)**2+(j-cy)**2 < (r - w)**2 or 
						  		    (i - cx)**2+(j-cy)**2 > r**2]))
			k = J - len(outliers)

		self.outliers = np.array([list(i) for i in outliers])
		
	def plot_data(self):
		"""
		plot inliers and outliers
		"""
		ax = plt.gca()
		ax.cla()
		ax.plot(self.outliers[:,0], self.outliers[:,1],'o',self.inliers[:,0], self.inliers[:,1],'o')
		exCircle = plt.Circle((self.c[0],self.c[1]),self.r,fill=False)
		inCircle = plt.Circle((self.c[0],self.c[1]),self.r - self.w,fill=False)
		fig = plt.gcf()
		
		fig.gca().add_artist(exCircle)
		fig.gca().add_artist(inCircle)

		plt.grid()
		plt.show()
		
	def plot_fit_data(self):
		"""
		plot inliers and outliers and the fitted digital annulus
		"""
		ax = plt.gca()
		ax.cla()
		ax.plot(self.outliers[:,0], self.outliers[:,1],'o',self.inliers[:,0], self.inliers[:,1],'o')
		exCircle = plt.Circle((self.c[0],self.c[1]),self.r,fill=False)
		inCircle = plt.Circle((self.c[0],self.c[1]),self.r - self.w,fill=False)
		fig = plt.gcf()
		
		fig.gca().add_artist(exCircle)
		fig.gca().add_artist(inCircle)

		fitdata = self.get_fit_para()
		fc = (fitdata[0][0], fitdata[1][0])
		fr = fitdata[2]

		fexCircle = plt.Circle((fc[0],fc[1]),fr,color='r',fill=False)
		finCircle = plt.Circle((fc[0],fc[1]),fr - self.w,color='r',fill=False)

		fig.gca().add_artist(fexCircle)
		fig.gca().add_artist(finCircle)

		plt.grid()
		plt.show()

	def fusion(self,seed):
		"""
		randomly mix outliers and inliers array with a seed number
		return a matrix N x 2, where N = inliers + outliers
		"""
		res = np.concatenate((self.inliers, self.outliers),axis=0)
		np.random.seed(seed)
		idx = np.random.permutation(res.shape[0])
		return res[idx,:]
		
	def get_fit_para(self):
		"""
		retrieve fitted parameters from dual space
		return a tuple with fitted center coordinates and radius
		"""
		itv = np.array([self.bestfit[0][0],self.bestfit[0][1]])
		p1, p2 = self.bestfit[0][2], self.bestfit[0][3]
		para = self.get_line(p1,p2)
		xo, yo, xv, yv = para[0,0], para[0,1], para[1,0], para[1,1]
		xc = xo + itv * xv
		yc = yo + itv * yv
		dc = self.__get_dis(para,itv[0],p1)

		return (xc,yc,dc)

	def fit(self):
		"""
		optimal consensus set fitting algorithm
		from paper of Rita Zrour et al. IWCIA 2015
		N.B : the paper guarantees the optimal solution, not unique solution !!!
		"""
		seed = 20 # initialize random number
		data = self.fusion(seed) # generate data from inliers and outliers
		N = data.shape[0] # get number of points
		self.maxval = 0 # maximal number of inliers
		self.bestfit = [] # best fitted intervals
		
		start_time = timeit.default_timer()

		for idx1 in np.arange(N-1):
			p1 = data[idx1,:] # first point			
			for idx2 in np.arange(idx1+1,N):
				p2 = data[idx2,:] # second point			
				para = self.get_line(p1,p2) # get line parameters
				self.interlist = [] # list contain all interval
				for idx3 in np.arange(N):
					if idx3 not in (idx1, idx2): # not check idx1, idx2
						p3 = data[idx3,:] # check point
						res = self.intersect(para,p1,p3) # get res in dual space
						# analyse interval of p3 such that p3 is inlier
						self.analyse_interval(para,res,p1,p3)
				
				dtype = [('inter',float),('f',int)]
				srtlst = np.array(self.interlist,dtype=dtype)
				srtlst.sort(order='inter')

				fval = 2
				for ival in np.arange(len(srtlst)-1):
					fval = fval + srtlst[ival][1]
					if fval > self.maxval:
						self.maxval = fval
						self.bestfit = [(srtlst[ival][0],srtlst[ival+1][0],p1,p2)]
					elif fval == self.maxval:
						self.bestfit.append((srtlst[ival][0],srtlst[ival+1][0],p1,p2))
		
		elapsed = timeit.default_timer() - start_time	
		self.runtime = elapsed

	def analyse_interval(self,para,res,p,q):
		"""
		search the interval in dual space for which q is inliers
		"""
		w = self.w
		if res[0]!=[]: # once
			if len(res[1])==1: # once
				dq1 = self.__get_dis(para,res[1][0]-1,q)
				dpw1 = self.__get_dis(para,res[1][0]-1,p-w)
				dq2 = self.__get_dis(para,res[1][0]+1,q)
				dpw2 = self.__get_dis(para,res[1][0]+1,p-w)
				if res[0][0] < res[1][0]:
					self.interlist.append((res[0][0],1))
					if dq2 < dpw2:
						self.interlist.append((res[1][0],-1))
					else:
						self.interlist.append((float('inf'),-1))		
				else:
					if dq1 < dpw1:
						self.interlist.append((res[1][0],1))
					else:
						self.interlist.append((-float('inf'),1))		
					self.interlist.append((res[0][0],-1))
			elif len(res[1])==2: # twice
				if res[0][0] < res[1][0]:
					if res[1][0]<res[1][1]:
						self.interlist.append((res[0][0],1))
						self.interlist.append((res[1][0],-1))
						self.interlist.append((res[1][1],1))
						self.interlist.append((float('inf'),-1))
					else:
						self.interlist.append((res[0][0],1))
						self.interlist.append((res[1][1],-1))
						self.interlist.append((res[1][0],1))
						self.interlist.append((float('inf'),-1))	
				else:
					if res[1][0]<res[1][1]:
						self.interlist.append((-float('inf'),1))	
						self.interlist.append((res[1][0],-1))
						self.interlist.append((res[1][1],1))
						self.interlist.append((res[0][0],-1))

					else:
						self.interlist.append((-float('inf'),1))	
						self.interlist.append((res[1][1],-1))
						self.interlist.append((res[1][0],1))
						self.interlist.append((res[0][0],-1))
			else: # zero #verify
				dp = self.__get_dis(para,res[0][0]+1,p)
				dq = self.__get_dis(para,res[0][0]+1,q)
				if dq < dp and dp - w >= 0:
					self.interlist.append((res[0][0],1))
					self.interlist.append((float('inf'),-1))	
				else:
					self.interlist.append((-float('inf'),1))	
					self.interlist.append((res[0][0],-1))
		else: # zero
			if len(res[1])==1: # once
				dp1 = self.__get_dis(para,res[1][0]-1,p-w)
				dq1 = self.__get_dis(para,res[1][0]-1,q)
				dp2 = self.__get_dis(para,res[1][0]+1,p-w)
				dq2 = self.__get_dis(para,res[1][0]+1,q)	
				if dq1 > dp1 and dq2 < dp2:
					self.interlist.append((-float('inf'),1))
					self.interlist.append((res[1][0],-1))
				elif dq1 < dp1 and dq2 > dp2:
					self.interlist.append((res[1][0],1))
					self.interlist.append((float('inf'),-1))
				elif dq1 > dp1 and dq2 > dp2:
					self.interlist.append((-float('inf'),1))	
					self.interlist.append((float('inf'),-1))
			elif len(res[1]==2): # twice
				if res[1][0] < res[1][1]:
					self.interlist.append((-float('inf'),1))
					self.interlist.append((res[1][0],-1))
					self.interlist.append((res[1][1],1))
					self.interlist.append((float('inf'),-1))
				else:
					self.interlist.append((-float('inf'),1))
					self.interlist.append((res[1][1],-1))
					self.interlist.append((res[1][0],1))
					self.interlist.append((float('inf'),-1))
			else: # zero
				dp = self.__get_dis(para,0,p)
				dpw = self.__get_dis(para,0,p-w)
				dq = self.__get_dis(para,0,q)
				if dq < dp and dq > dpw:
					self.interlist.append((-float('inf'),1))	
					self.interlist.append((float('inf'),-1))

	def __refine(self,para,p,q,t):
		"""
		refine result of intersection point of q and p when the the digital annulus is invalid
		precond : internal border intersects x axis
		return the refined point
		"""			
		xp, yp, xq, yq = p[0], p[1], q[0], q[1]
		xo, yo, xv, yv = para[0,0], para[0,1], para[1,0], para[1,1]
		w = self.w
		res = sp.solve((((xp-xo-t*xv)**2+(yp-yo-t*yv)**2)**0.5)-w,t) # review for speed ?
		dtp = self.__get_dis(para,float(res[0]),p)
		dtq = self.__get_dis(para,float(res[0]),q)
		if dtq < dtp:
			return float(res[0])
		else:
			return float(res[1])
	
	def intersect(self,para,p,q):
		"""
		get intersected result of points p,q and the digital annulus parameterized by para in dual space
		return a list of valided intersection points
		"""
		
		# start_time = timeit.default_timer()

		xp, yp, xq, yq = p[0], p[1], q[0], q[1]
		xo, yo, xv, yv = para[0,0], para[0,1], para[1,0], para[1,1]
		w = self.w
		t = sp.Symbol('t') #symbolic variable

		# intersect bt q and the external circle
		a = (xp**2-xq**2)+(yp**2-yq**2)-2*xo*(xp-xq)-2*yo*(yp-yq)
		b = -2*xv*(xp-xq)-2*yv*(yp-yq)

		# start_time = timeit.default_timer()
		t1 = sp.solve(a+b*t,t) # solve eq
		# elapsed = timeit.default_timer() - start_time
		# print 't1 = {}'.format(t1)
		# print 't1 time = {}'.format(elapsed)
		#... how to manipulate t1 ?
		if t1 != []:	
			if t1[0] is sp.EmptySet:
				t1 = []
			else:
				t1[0] = float(t1[0])
				if (self.__get_dis(para,t1[0],p) - w) < 0:
					t1 = [self.__refine(para,p,q,t)]

		# intersect bt q and the internal circle
		t2 = []
		# check condition 1
		a2 = xv**2+yv**2
		b2 = -2*xv*(xp-xo)-2*yv*(yp-yo)
		c2 = (xp-xo)**2+(yp-yo)**2-w**2
		# start_time = timeit.default_timer()
		cond1 = sp.solve_poly_inequality(sp.Poly(a2*t**2+b2*t+c2,t),">=")
		# elapsed = timeit.default_timer() - start_time
		# print 'cond1 = {}'.format(cond1)
		# print 'cond1 time = {}'.format(elapsed)
		# N.B. cond1 can return a list with 2 intervals as maximum
		# must verify cond1 == []

		# check condition 2
		a3 = (a + w**2)/(2*w)
		b3 = b/(2*w)
		# start_time = timeit.default_timer()
		# cond2 = sp.solve_poly_inequality(sp.Poly(a3+b3*t,t),">=")
		cond2 = sp.solve_univariate_inequality(a3+b3*t >= 0,t,relational=False)
		# elapsed = timeit.default_timer() - start_time
		# print 'cond2 = {}'.format(cond2)
		# print 'cond2 time = {}'.format(elapsed)

		# N.B. cond2 return a list with only 1 interval
		# should verify cond2 == []

		# start_time = timeit.default_timer()
		cond = []
		if cond1 != [] and cond2 != []:
			cond = [cd.intersect(cond2) for cd in cond1]
		# elapsed = timeit.default_timer() - start_time
		# print 'cond = {}'.format(cond)
		# print 'cond time = {}'.format(elapsed)	

		a4 = b3**2-xv**2-yv**2
		b4 = 2*a3*b3+2*xv*(xp-xo)+2*yv*(yp-yo)
		c4 = a3**2-(xp-xo)**2-(yp-yo)**2
		# start_time = timeit.default_timer()
		temp = sp.solve(a4*t**2+b4*t+c4,t) # solve eq
		# elapsed = timeit.default_timer() - start_time
		# print 'temp time = {}'.format(elapsed)	

		# N.B : imagery number, there is alway solution
		# start_time = timeit.default_timer()
		if cond != []:
			# print cond
			for itemp in temp:
				# print itemp
				if True in [cd.contains(itemp) for cd in cond]:
					t2 = t2 + [float(itemp)]
		# elapsed = timeit.default_timer() - start_time
		# print 'check temp time = {}'.format(elapsed)	

		# elapsed = timeit.default_timer() - start_time
		# print 'time = {}'.format(elapsed)	
		# print (t1,t2)

		return (t1,t2)
	
	def get_line(self, p1, p2):
		"""
		return line parameters in which the line is perpendicular to the line bt p1 and p2
		"""
		o = (p1+p2)/2
		v = p2 - p1
		return np.reshape(np.concatenate((o,np.array([-v[1],v[0]]))),(2,2))
					
	def __get_dis(self, para, t, p):
		"""
		return distance of point p to the center t in dual space
		"""
		# para = self.__get_line(p1,p2)
		xo, yo, xv, yv = para[0,0], para[0,1], para[1,0], para[1,1]
		xt, yt = xo + t * xv, yo + t * yv
		xp, yp = p[0], p[1]
		dt = np.sqrt(((xp - xt)**2 + (yp - yt)**2))
		return dt

	def plot_dis(self,p1,p2,p3,it):
		"""
		plot distance in dual space
		"""
		t = np.sort(np.random.uniform(it[0],it[1],5000))
		para = self.get_line(p1,p2)

		itct = self.intersect(para,p1,p3)
		ritct = []
		if itct[0] != []:
			if itct[0][0] is not sp.EmptySet:
				ritct = ritct + [float(itct[0][0])]

		if itct[1] != []:
				ritct = ritct + [float(val) for val in itct[1]]

		ritct = np.array(ritct)

		dtitct = self.__get_dis(para,ritct,p3)

		dt1 = self.__get_dis(para,t,p1)
		dt2 = dt1 - self.w
		dt3 = self.__get_dis(para,t,p3)

		ax = plt.gca()  # gca stands for 'get current axis'
		ax.cla()
		ax.spines['left'].set_position('zero')
		ax.spines['right'].set_color('none')
		ax.spines['bottom'].set_position('zero')
		ax.spines['top'].set_color('none')
		ax.spines['left'].set_smart_bounds(True)
		ax.spines['bottom'].set_smart_bounds(True)
		ax.xaxis.set_ticks_position('bottom')
		ax.yaxis.set_ticks_position('left')

		ax.plot(t,dt1,t,dt2,t,dt3,ritct,dtitct,'o')

		st = 'intersect : '
		for i in np.arange(ritct.shape[0]):
			st = st+'('+str(round(ritct[i],2))+','+str(round(dtitct[i],2))+')'
		plt.title(st)
		plt.grid()
		plt.show()

		
