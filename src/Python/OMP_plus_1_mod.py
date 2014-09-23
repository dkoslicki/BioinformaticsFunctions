#!/local/cluster/bin/python2.7
import os
import sys
from StringIO import StringIO
import scipy.optimize.nnls
import scipy.sparse
import numpy as np
from subprocess import *
import gzip
import itertools



 
def OMP_plus_1_mod(X,mu,nu,I):
	#X needs to be a numpy array
	#mu needs to be a numpy vector
	#nu needs to be a real
	#I needs to be an int
	

	#Finding the l_2 norms of columns of X
	M, N = np.shape(X)
	ColNorms = np.zeros((N,));
	for i in range(N):
		ColNorms[i]=np.linalg.norm(X[:,i])


	gamma = np.zeros((N,))
	r = mu
	T = list()
	it = 0
	#Iterations
	Xtrans=np.transpose(X) #Pre-compute the transpose
	
	while  np.absolute(np.sum(gamma)-1) > nu and  it < I:
		#Step 3 of the 'Iterations of Algorithm 1'
		e = Xtrans.dot(r) #e = Xtrans*r;
		for t in T:	#e(T)=0;
			e[t] = 0
			
		for ind in range(len(e)): #e(find(e <= 0))=0; e = e ./ ColNorms;
			if e[ind]<=0:
				e[ind] = 0
			e[ind]/ColNorms[ind]

		j0 = np.argmax(e) #[~,j0]=max(e);

		#Step 4 of the 'Iterations of Algorithm 1'
		T.append(j0) #T = [T j0];
    
		#Step 5 of the 'Iterations of Algorithm 1'
		Xp=np.zeros((M,len(T)))
		for ind in range(len(T)): #Xp = X(:,T);
			Xp[:,ind]=X[:,T[ind]]
		
		gamma_temp, rnorm = scipy.optimize.nnls(Xp, mu) #gamma(T)=lsqnonneg(Xp,mu);
		for ind in range(len(T)):
			gamma[T[ind]] = gamma_temp[ind]
		
		r = mu - Xp.dot(gamma_temp) #r = mu - Xp*gamma(T);
		it = it + 1
		
		#Increment of iterations count
		it=it+1

	gamma = np.divide(gamma,float(np.sum(gamma))) #gamma / sum(gamma)
	return gamma, it

#Test data
##X = np.array([[1,2,3],[4,5,6]])
#X=np.random.rand(4096,100000)
#mu = np.random.rand(4096)
##mu = np.array([1,2])
#nu = .01
#I = 1000

#gamma, it = OMP_plus_1_mod(X,mu,nu,I)
#print(gamma)
#print(it)

#if __name__ == "__main__":
#	main(sys.argv[1:])