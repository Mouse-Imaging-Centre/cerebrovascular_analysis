#!/usr/bin/env python
# -*- coding: utf-8 -*-

#### the logarithmic Naive Bayes classifier Training module

#### Gets all the training graphs as input, creates a feature Vector, use feature Vector to create the training result
# reads in the graph and put the features in a feature-vector of type numpy matrix
# reads in the graph labels and put them in target-vector of type numpy array
# reads in the graphs and create adjecancy matrix

# NOTE: in the training we can select best features (most discriminative)

#features (from graph edge_properties calculated by feature extraction ):
#label = int(0)
#0.curvature 	[>0]
#1.diameter  	[>0]
#2.length    	[>0]
#3.tortuosity
#4.midpointX
#5.midpointY
#6.midpointZ
#7.directionX  		[in +Z-plane]
#8.directionY		[in +Z-plane]
#9.directionZ		[in +Z-plane]
#10.directionX_neg  	[in -Z-plane]
#11.directionY_neg	[in -Z-plane]
#12.directionZ_neg	[in -Z-plane]
#13.angle with +X	[0< <pi]		=> can be converted to 0<..<pi/2; if pi/2<x<pi => pi-x
#14.angle with +Y	[0< <pi]		=> can be converted to 0<..<pi/2; if pi/2<x<pi => pi-x
#15.angle with +Z	[0< <pi]		=> can be converted to 0<..<pi/2; if pi/2<x<pi => pi-x
#16.proximity		[with 113 center_points of anatomical regions in MR atlas] 2D array [[MRregion,proximity],[MRregion,proximity],..]
#17.angleswref		[with 32 reference directions of 32 vascular labels] 2D array [[vascular_label,angleswref],[vascular_label,angleswref],..]
#18.rel_dir		[relative angle with 4 neighbors]
#19.rel_diameter 	[relative diameter with 4 neighbors]
# NOTE: input graphs should have been reduced to intermediaries and have diameter, length, intermediaries and label as edge_properties
#


#### option to perform PCA on all features_vectors and make them orthogonal
##1. suppose all Gaussian idd => learning =  a vector of mean value for each feature, a vector of std value for each feature 	& class priors
##2. suppose some Gaussian idd, some Gamma{, some multivariate Gaussain} => learning = output vector of parameters of each feature distribution 	& class priors
##3. suppose a multivariate gaussian on all features (GDA) => learning = mean vector and covariance matrix of the distribution 	& class priors


#
#
#  Created May 27, 2012
#  modified Oct 30, 2012  change in the kent distribution calculation to get rid of outliers
#  Last modified July 3, 2013  add features rel_diameter and rel_dir
#  Last modified April 8, 2015  line 235: mri_label_dist was modified
#  Sahar Ghanavati
 

from sys import argv
import math, os, shelve, string ,time
import commands
import copy
#import matplotlib.pyplot as plt
#from pylab import *		#needed for command find

from scipy import special
from scipy.optimize import fmin,fmin_ncg
import scipy.linalg

from vessel_tracking import graph_analysis
from numpy import *
import numpy
import numpy.matlib as Matlib
from optparse import OptionParser, Option, OptionValueError
from minc_util.progress import progress_report
from sys import argv
import math, os, shelve, string ,time, sys
import commands
import copy
#import matplotlib.pyplot as plt
#from pylab import *		#needed for command find
import pylab		#needed for command find
from scipy import special
import scipy
import random

#---------------------------------------------------------------------------------
#
def print_adj_mat (adjacency_matrix,labelNumerics):
	print  "\nlabelNumerics: ",labelNumerics 
	print "\n\nadjacency_matrix:"
	for i in range(adjacency_matrix.shape[0]):
		print "adj_mat row ", i, " : min = ",(adjacency_matrix[i,:].min())," max = ",(adjacency_matrix[i,:].max()), " sum = ",sum(adjacency_matrix[i,:]), "\n"
	
	#print  labelNumerics[0:15]
	if (len(labelNumerics)> 14):
		print "\t",
		for j in range(15):
			print labelNumerics[j],"\t",
		print " "	
		for i in range(adjacency_matrix.shape[0]):
			print "\n", labelNumerics[i],"\t",
			for j in range(15):
				print "%.4f"%adjacency_matrix[i,j],"\t",
		print "\n\n"
	else:
		print "\t",
		for j in range(len(labelNumerics)):
			print labelNumerics[j],"\t",
		print " "	
		for i in range(adjacency_matrix.shape[0]):
			print "\n", labelNumerics[i],"\t",
			for j in range(len(labelNumerics)):
				print "%.4f"%adjacency_matrix[i,j],"\t",
		print "\n\n"
		return 0
		
		
	if (len(labelNumerics)> 29):
		print "\t",
		for j in range(15,30):
			print labelNumerics[j],"\t",
		print " "	
		for i in range(adjacency_matrix.shape[0]):
			print "\n", labelNumerics[i],"\t",
			for j in range(15,30):
				print "%.4f"%adjacency_matrix[i,j],"\t",
		print "\n\n"
	else:
		print "\t",
		for j in range(15,len(labelNumerics)):
			print labelNumerics[j],"\t",
		print " "	
		for i in range(adjacency_matrix.shape[0]):
			print "\n", labelNumerics[i],"\t",
			for j in range(15,len(labelNumerics)):
				print "%.4f"%adjacency_matrix[i,j],"\t",
		print "\n\n"
		return 0
		
		
	if (len(labelNumerics)> 44):
		print "\t",
		for j in range(30,45):
			print labelNumerics[j],"\t",
		print " "	
		for i in range(adjacency_matrix.shape[0]):
			print "\n", labelNumerics[i],"\t",
			for j in range(30,45):
				print "%.4f"%adjacency_matrix[i,j],"\t",
		print "\n\n"
	else:
		print "\t",
		for j in range(30,len(labelNumerics)):
			print labelNumerics[j],"\t",
		print " "	
		for i in range(adjacency_matrix.shape[0]):
			print "\n", labelNumerics[i],"\t",
			for j in range(30,len(labelNumerics)):
				print "%.4f"%adjacency_matrix[i,j],"\t",
		print "\n\n"
		return 0
		
	
	#print  labelNumerics[15:]
	print "\t",
	for j in range(45,len(labelNumerics)):
		print labelNumerics[j],"\t",
	print "End"	
	for i in range(adjacency_matrix.shape[0]):
		print "\n", labelNumerics[i],"\t",
		for j in range(45,len(labelNumerics)+1):
			print "%.4f" %adjacency_matrix[i,j],"\t",
	print "\n\n"

	return 0
	

	
	

	
#
def find_key(dic, val, current_cnt):
#"""return the key of dictionary dic given the value"""
#if more than one key, then return the one that is closest to current_cnt
	found_keys = [k for k, v in dic.iteritems() if v == val]		#returns a list
	if len(found_keys)==1:
		#return [k for k, v in dic.iteritems() if v == val][0]
		return found_keys[0]
	elif len(found_keys)>1:
		closest= [abs(i-current_cnt) for i in found_keys]
		return found_keys[closest.index(min(closest))]
		
		
def graph2featureVect(training_files,featureNames):
	graphs=[]
	for training_file in training_files:
		try:
			graphs.append( graph_analysis.input_graph(training_file))  	# open graph data and copy contents to memory
			print ("Succefully read in %s\n" %training_file)	
		except:
			print("Error reading in %s\n" %training_file)
			
	#### initialize the feature_vector to be zeros(#datapoints=sum(len(g.edge_list()) w label!=0 , len(featureNames))
	##	initialize the target_vector to be zeros(#datapoints=sum(len(g.edge_list()) w label!=0 , 1)
	MRI_labels=[]
	if 'mri_label_dist' in featureNames:
		MRI_labels=[l[0] for l in graphs[0].edge_property(graphs[0].edge_list()[0],'mri_label_dist')]
		feature_vector = Matlib.empty((0,len(featureNames)+len(MRI_labels)-1),float)
	else:
		feature_vector = Matlib.empty((0,len(featureNames)),float)
	target_vector = []	#Matlib.empty((0,1),float)
	labelNumerics = []
	edge_w_indx = {}
	
	g_cnt=-1
	e_cnt = 0
	for g in graphs:
		g_cnt =g_cnt+1
		for e in g.edge_list():
			if 'label' not in g.edge_properties(e).keys():
				print ("ERROR! The edge (%d,%d) in graph %s doesn't have label property.\nAborted!\n" %(e[0],e[1],training_files[g_cnt]))
				exit (0)
			#else:
			if g.edge_property(e,'label')>-1:	#including the edges with label 0
				#target_vector= numpy.vstack([target_vector, array(g.edge_property(e,'label'))])
				target_vector.append(int(g.edge_property(e,'label')))
				if g.edge_property(e,'label') not in labelNumerics and g.edge_property(e,'label')>0:
					labelNumerics.append(int(g.edge_property(e,'label')))
				features=[]
				for f in featureNames:
					if f not in g.edge_properties(e).keys():
						print ("ERROR! The edge (%d,%d) in graph %s doesn't have feature %s property.\nAborted!\n" %(e[0],e[1],training_files[g_cnt],f))
						exit (0)
					if not f== 'mri_label_dist' and not f=='rel_dir' and not f=='rel_diameter':	
						features.append(g.edge_property(e,f))
					elif f== 'mri_label_dist':
						for mr_l in MRI_labels:
							#if type(g.edge_property(e,'mri_label_dist'))==list:
								#print mr_l,":",e, " list ", training_files[g_cnt]
								#ind=[l[0] for l in g.edge_property(e,'mri_label_dist')].index(mr_l)
							#else:
								#print mr_l,":",e, ", ", training_files[g_cnt]
								#ind=[l[0] for l in g.edge_property(e,'mri_label_dist').tolist()].index(mr_l)
							mri_label_dist_list = [int(l[0]) for l in g.edge_property(e,'mri_label_dist')]
							if mr_l in mri_label_dist_list:
								ind = mri_label_dist_list.index(mr_l)
								features.append(g.edge_property(e,'mri_label_dist')[ind][1])
					elif f=='rel_dir' or f=='rel_diameter':
						rel_f = 1.0
						for rel_i in range(min(len(g.edge_property(e,f)),4)):	#we only consider up to 4 adjacent edges, if there are less 1.0 will be used
							rel_f = rel_f * g.edge_property(e,f)[rel_i]
						features.append(rel_f)	
							
				feature_vector	= numpy.vstack([feature_vector, array(features)])		#hstack #a = matrix([[10,20,30]]); a=append(a,[[1,2,3]],axis=0); a=append(a,[[15],[15]],axis=1)
			edge_w_indx[e_cnt]= e 
			e_cnt = e_cnt+1
			
	
	#### !!find neighbouring edges and in iterations find their
	
	adjacency_matrix = (1.0/len(labelNumerics))*Matlib.ones((len(labelNumerics),len(labelNumerics)+1),float)		#### smoothing (for 0s in the adj_matrix)! to be minimum of 1 occurance for the nieghbourhood!
	edge_neighboring_indx = {}
	
	g_cnt=-1
	e_cnt = 0
	for g in graphs:
		g_cnt =g_cnt+1
		for e in g.edge_list():
			#if 'label' not in g.edge_properties(e).keys():
				#print ("ERROR! The edge (%d,%d) in graph %s doesn't have label property.\nAborted!\n" %(e[0],e[1],training_files[g_cnt]))
                                #exit (0)
                        #else:
			if g.edge_property(e,'label')>0:
				indx0= labelNumerics.index(g.edge_property(e,'label'))
					
				e1_neighbours= g.vertices[e[0]].edges
				#e1_neighbours.remove(e2)
				for v in e1_neighbours:
					if v!=e[1]:
						neighbor_edge=tuple((min(e[0],v),max(e[0],v)))
						if g.edge_property(neighbor_edge,'label')>0:
							indx1= labelNumerics.index(g.edge_property(neighbor_edge,'label'))
							adjacency_matrix[indx0,indx1]= adjacency_matrix[indx0,indx1]+1
				e2_neighbours= g.vertices[e[1]].edges
				#e2_neighbours.remove(e1)
				for v in e2_neighbours:
					if v!=e[0]:
						neighbor_edge=tuple((min(e[1],v),max(e[1],v)))
						if g.edge_property(neighbor_edge,'label')>0:
							indx1= labelNumerics.index(g.edge_property(neighbor_edge,'label'))
							adjacency_matrix[indx0,indx1]= adjacency_matrix[indx0,indx1]+1
					
				if len(g.vertices[e[0]].edges)==1:		#end point
					adjacency_matrix[indx0,len(labelNumerics)] +=1
				if len(g.vertices[e[1]].edges)==1:		#end point
					adjacency_matrix[indx0,len(labelNumerics)] +=1					
			#find edge neighbouring indeces
			neighbor_indx=[]
			e1_neighbours= g.vertices[e[0]].edges
			for v in e1_neighbours:
				if v!=e[1]:
					neighbor_edge=tuple((min(e[0],v),max(e[0],v)))
					neighbor_indx.append(find_key(edge_w_indx, neighbor_edge, e_cnt))
			e2_neighbours= g.vertices[e[1]].edges
			for v in e2_neighbours:
				if v!=e[0]:
					neighbor_edge=tuple((min(e[1],v),max(e[1],v)))
					neighbor_indx.append(find_key(edge_w_indx, neighbor_edge, e_cnt))
					
			#if (len(neighbor_indx)!=2 and len(neighbor_indx)!=4):
				#print ("\n\nERROR! edge (%d,%d) has %d neighbouring vertices!" %(e[0],e[1],len(neighbor_indx)))
				#if (len(neighbor_indx)!=0):
					#print "Aborted!\n"
					#exit(0)
			edge_neighboring_indx[e_cnt]=neighbor_indx
			e_cnt = e_cnt+1
			
	#### adjacency matrix normalization 
	for i in range(adjacency_matrix.shape[0]):
		adjacency_matrix[i,:]= adjacency_matrix[i,:]/sum(adjacency_matrix[i,:])
		
	
	print_adj_mat (adjacency_matrix,labelNumerics)
	
	return labelNumerics, feature_vector, target_vector, adjacency_matrix, edge_w_indx, edge_neighboring_indx 	#featureVect


def machineEpsilon(func=float):
	machine_epsilon = func(1)
	while func(1)+func(machine_epsilon) != func(1):
		machine_epsilon_last = machine_epsilon
		machine_epsilon = func(machine_epsilon) / func(2)
	return machine_epsilon_last

def PCA (feature_vector,CDF_thresh,Y_bar=-inf, normalizing_stds = -inf, eigenvects=-inf,eigenvals=-inf):
	#### 1. centering the data to mean=0
	if type(Y_bar)==float:
		Y_bar = numpy.mean(feature_vector,0)	#mean of each column of feature_vector = mean of each feature data points
	feature_vector_cntr = feature_vector - tile(Y_bar,(feature_vector.shape[0],1))
	#### 2. normalizing the data to std=1
	if type(normalizing_stds)==float:
		normalizing_stds = numpy.std(feature_vector_cntr,0)
	feature_vector_cntrN = feature_vector_cntr / tile(normalizing_stds,(feature_vector_cntr.shape[0],1))

	if type(eigenvects)==float or type(eigenvals)==float:
		#### test: cov_matrix symmetric positive semidefinite (PSD) matrix => has a nonnegative eigenvalue
		cov_matrix = numpy.cov(feature_vector_cntrN,None,0)	#cov(m, y=None, rowvar=1, bias=0, ddof=None) rowvar=0:each column represents a variable, while the rows contain observations
		#print type(cov_matrix), "\ncov_matrix:", cov_matrix
		D,V = numpy.linalg.eig(cov_matrix)		# a: array_like shape (M, M), D:eigenvalues ndarray shape (M,) NOT ordered, V: eignevectors ndarray shape (M, M) The normalized (unit “length”) eigenvectors, column V[:,i] is the eigenvector corresponding to the eigenvalue D[i]			 
		#i = numpy.argsort(D)	#Returns the indices that would sort an array, ascending
		i = numpy.argsort(numpy.array(D))	#D[::-1] #Returns the indices that would sort an array, descending
		inv_i = i[::-1]
		eigenvals = D[inv_i]
		eigenvects = V[:,inv_i]
		#print type(eigenvals)
		#print type(eigenvects)
		print "sorted eigen values are (should be >0) : " , eigenvals		#### test: eigvals should be >0 => orthogonal_feature_vector: real
		#print "sorted eigen vects are (should be >0) : " , eigenvects		#### test: eigvals should be >0 => orthogonal_feature_vector: real
	
	cumulative_eigval = numpy.cumsum(eigenvals)
	cumulative_eigval = cumulative_eigval/cumulative_eigval[-1]*100.0	# % CDF
	M_ind = min(pylab.find(cumulative_eigval >= CDF_thresh))
	W = eigenvects[:,0:M_ind]
	
	print "\n\nM_ind:",M_ind, " out of ", numpy.matrix(V).shape[0], " features"
	print (numpy.matrix(W).H).shape
	print numpy.matrix(feature_vector_cntrN).shape
	print  numpy.matrix(numpy.diag(1.0/numpy.sqrt(eigenvals))).shape
	print ((numpy.matrix(feature_vector_cntrN) * numpy.matrix(numpy.diag(1.0/numpy.sqrt(eigenvals)))).transpose()).shape
	
	#### orthogonal_feature_vector = feature_vector_cntr * V * D^(-.5)??   
	#based on karhunen-Leove theorem and data projection from http://en.wikipedia.org/wiki/Principal_component_analysis:
	# projected data = conjugate_transpose(eigenvects)* (feature_vects/Diag(sqrt(eigenvals)))
	# a.H conjugate transpose
	orthogonal_feature_vector = numpy.matrix(W).H * (numpy.matrix(feature_vector_cntrN) * numpy.matrix(numpy.diag(1.0/numpy.sqrt(eigenvals)))).transpose()
	orthogonal_feature_vector = orthogonal_feature_vector.transpose()
	#### test: correctness of PCA: cov(orthogonal_feature_vector) almost = I
	return [orthogonal_feature_vector, numpy.array(Y_bar), numpy.array(normalizing_stds), numpy.matrix(eigenvects), numpy.array(eigenvals), M_ind]


def plot_univariate_gaussian(f,color_f,l,feature_vect,mean_val, std_val,c,saveoption,directory):
	#pylab.figure(c)
	pylab.hist(feature_vect,bins=100, normed=True, color=color_f)#plotting things as a probability distribution,normed=1, cumulative=True)
	#pylab.show()
	pylab.hold(True)
	y=[]
	#x=arange(min(feature_vect),max(feature_vect),0.001)
	x=arange(mean_val-3.0*std_val,mean_val+3.0*std_val,0.001)
	for i in x:
		y.append(10*univaraite_gaussian ( i, mean_val , std_val))
	pylab.plot(x,y, color_f, label=labelNumeric2Name[l])
	#pylab.show()
	pylab.hold(True)



def plot_gamma(f,l,feature_vect,k,theta,c,saveoption,directory):
	t= ('Gamma distribution fit with k %f and theta %f to feature %s for label %d with total length %d' %(k,theta,f,l, len(feature_vect)))
	print t
	pylab.figure(c)
	pylab.hist(feature_vect,bins=100, normed=True)#plotting things as a probability distribution,normed=1, cumulative=True)
	#pylab.show()
	hold(True)
	y=[]
	x=arange(min(feature_vect),max(feature_vect),0.001)
	for i in x:
		y.append(gamma_distribution ( i, k , theta))
	pylab.plot(x,y,'r')
	#pylab.show()
	title(t)
	if saveoption:
		figname= ("%s/gamma_fit_feature_%s_label%s.png" %(directory,f,l))
		pylab.savefig(figname, format="png")
	
	
def gamma_distribution(x,k,theta):
	p = pow(x,k-1.0) * numpy.exp(-x/theta)/ (special.gamma(k) * pow(theta,k))
	return p
	
def univaraite_gaussian(x,mean_val,std_val):
	p = (1.0/numpy.sqrt(2.0*numpy.pi*pow(std_val,2))) * numpy.exp( pow(x-mean_val,2) / (-2.0*pow(std_val,2)) )
	return p	
	
def multivariate_gaussian(X,mean_vect, cov_mat):
	m = cov_mat.shape[0]
	n = cov_mat.shape[1]
	if (mean_vect.shape[1]!=n):
		print "ERROR! The length of mean vector is not equal to the size of covariance matrix of multivariate gaussian!\nAborted!"
		exit (0)
	if (n!=m) :
		print "ERROR! The covariance matrix of multivariate gaussian is not squared!\nAborted!"
		exit (0)

	exp_arg = numpy.matrix((X-mean_vect))* numpy.matrix(numpy.linalg.inv(cov_mat)) * numpy.matrix((X-mean_vect).transpose())	
	p = pow(1.0/(2*numpy.pi) , n/2.0) * pow(numpy.linalg.det(cov_mat) , -0.5) * numpy.exp(-0.5 * exp_arg)
	return p


#def Kent_calculation (L):	#receiving L= feature_vector[indx,axialfeatureidx], Li=[lx,ly,lz])
##f(x) = 1/c(k,beta) exp(k*gamma1*x + beta*((gamma2*x)^2 - (gamma3*x)^2))
##g(x) = exp(k*gamma1*x + beta*((gamma2*x)^2 - (gamma3*x)^2))
##LL = log likelihood = k*gamma1*x + beta*((gamma2*x)^2 - (gamma3*x)^2)
##dLL/dk = gamma1*x
##dLL/dbeta = (gamma2*x)^2 - (gamma3*x)^2

	##calculate gamma1,2,3
	#T= Matlib.zeros((3,3),float)
	#for i in range(L.shape[0]):
		#T += L[i,:]*L[i,:].transpose()  		#sigma Li(3,1)*Li'(1,3) => T(3,3)
	
	#D,V = numpy.linalg.eig(T)		# a: array_like shape (M, M), D:eigenvalues ndarray shape (M,) NOT ordered, V: eignevectors ndarray shape (M, M) The normalized (unit “length”) eigenvectors, column V[:,i] is the eigenvector corresponding to the eigenvalue D[i]			 
	##i = numpy.argsort(D)	#Returns the indices that would sort an array, ascending
	#i = numpy.argsort(D[::-1])	#Returns the indices that would sort an array, descending
	#Tau = D[i]	#sorted eigenvals
	#t = V[:,i]	#sorted a=eigenvects
	#Tau_bar = Tau/float(L.shape[0])
	#gamma1 = t[:,0]
	#gamma2 = t[:,1]
	#gamma3 = t[:,2]
	
	##### run EM to find k and beta
	###1.initialize
	#k0 = 10
	#beta0 = 0 
	#v0 = [k0,beta0]
	###2.Parametric function: 'v=[k,beta]' is the parameter vector, 'x=[L]' the independent variable		x is N-by-3 and gamma is 3-by-1 => fp is N-by-1 => sum =1-by-1 
	###negLL = neg log likelihood = -k*gamma1*x - beta*((gamma2*x)^2 - (gamma3*x)^2)
	#fp = lambda v, x: (-v[0]*x*gamma1 - v[1]*(pow(x*gamma2,2) - pow(x*gamma3,2))).sum()
	#dfp = lambda x: [x*gamma1, pow(x*gamma2,2) - pow(x*gamma3,2)]
	###3.maximize fp w.r.t v
	###simplex: fmin(func, x0, args=(), xtol=0.0001, ftol=0.0001, maxiter=None, maxfun=None, full_output=0, disp=1, retall=0, callback=None)
	## v = fmin(fp,v0, args=(L),maxiter=10000,maxfun=10000)
	###conjugate gradient: fmin_ncg(f, x0, fprime, fhess_p = None, fhess = None, args = (), avextol = 1.0000000000000001e-05, epsilon = _epsilon, maxiter = None, full_output = 0, disp = 1, retall = 0, callback = None)
	#v = fmin_ncg(fp,v0, fprime = dfp, args=(L), avextol=1e-8,maxiter=10000)

	
	#k=v[0]
	#beta=v[1]
	
	#ck= numpy.sqrt((k-2*beta)*(k+2*beta))/(2*numpy.pi*numpy.exp(k))
	#return [k,beta,ck,gamma1,gamma2,gamma3]

def fp (v,x,gamma1, gamma2,gamma3): 
#### negative log likelihood
	N = x.shape[0]  #x =[N-by-3]
	k = v[0]
	b = v[1]
	
	if (k==0):
		k=1
	elif (k<0):
		k=-k
	if (b<0):
		b=-b

	#we should have 2*abs(b) < abs(k)
	if (2*b >= k):
		b = (k-1)/2
	if (b <0):
		b = 0  

	
	el1 =  abs(array(x*gamma1))
	el2 =  array(x*gamma2)
	el3 =  array(x*gamma3)
	
	
	neg_LL = -N * (0.5*log((k-2*b)*(k+2*b))-k) - (k*el1 + b*(numpy.power(el2,2) - numpy.power(el3,2))).sum()
	return neg_LL
	
	
	
def dfp (v,x,gamma1, gamma2,gamma3):
#### derivative negative log likelihood w.r.t k and beta
	N = x.shape[0]  #x =[N-by-3]
	k = v[0]
	b = v[1]
	
	if (k==0):
		k=1
	elif (k<0):
		k=-k
	if (b<0):
		b=-b

	#we should have 2*abs(b) < abs(k)
	if (2*b >= k):
		b = (k-1)/2
	if (b <0):
		b = 0  


	el1 =  abs(array(x*gamma1))
	el2 =  array(x*gamma2)
	el3 =  array(x*gamma3)
	
	neg_dLL = [ -N*(k/((k-2*b)*(k+2*b))-1) - el1.sum() , 4*N*b/((k-2*b)*(k+2*b)) - (numpy.power(el2,2) - numpy.power(el3,2)).sum() ]

	return neg_dLL

def Kent_calculation (Lpos):	#receiving Lpos= feature_vector[indx,axialfeatureidx], Li=[lx,ly,lz]) of type numpy.matrix
	print "Kent_calculation:"
	for i in range(Lpos.shape[0]):
		Lpos[i] = Lpos[i]/numpy.linalg.norm(Lpos[i])

	

	Lneg=-Lpos
	L=bmat('Lpos; Lneg')	#?matrix(numpy.negative(array(Lpos)))	#[a b; c d] {MATLAB}= vstack([hstack([a,b]),hstack([c,d])]) {numpy.array}= bmat('a b; c d') {numpy.matrix}		
#f(x) = 1/c(k,beta) exp(k*gamma1*x + beta*((gamma2*x)^2 - (gamma3*x)^2))
#g(x) = exp(k*gamma1*x + beta*((gamma2*x)^2 - (gamma3*x)^2))
#LL = log likelihood = k*gamma1*x + beta*((gamma2*x)^2 - (gamma3*x)^2)
#dLL/dk = gamma1*x
#dLL/dbeta = (gamma2*x)^2 - (gamma3*x)^2

	T= (L.transpose())*L
	D,V = numpy.linalg.eig(T)		# a: array_like shape (M, M), D:eigenvalues ndarray shape (M,) NOT ordered, V: eignevectors ndarray shape (M, M) The normalized (unit “length”) eigenvectors, column V[:,i] is the eigenvector corresponding to the eigenvalue D[i]			 
	#i = numpy.argsort(D)	#Returns the indices that would sort an array, ascending
	i = numpy.argsort(D[::-1])	#Returns the indices that would sort an array, descending
	Tau = D[i]	#sorted eigenvals
	t = V[:,i]	#sorted a=eigenvects
	Tau_bar = Tau/float(L.shape[0])
	gamma1 = t[:,0]
	gamma2 = t[:,1]
	gamma3 = t[:,2]
	
	noOutlierL = Matlib.empty((0,L.shape[1]),float)
	for i in range(L.shape[0]):
		if ( arccos((dot(L[i,:],gamma1))/(linalg.norm(L[i,:])*linalg.norm(gamma1))) < (pi/4.0)):
			noOutlierL = numpy.vstack([noOutlierL, array(L[i,:])])
	
	
	#### run EM to find k and beta
	#### 1.initialize
	k0 = 10
	beta0 = 0 
	v0 = [k0,beta0]
	v = scipy.optimize.fmin(fp,v0, args=(noOutlierL,gamma1,gamma2,gamma3),maxiter=10000,maxfun=10000)
	
	k=v[0]
	beta=v[1]
	
	if (k==0):
		k=1
	elif (k<0):
		k=-k
	if (beta<0):
		beta=-beta
	#we should have 2*abs(b) < abs(k)
	if (2*beta >= k):
		beta = (k-1)/2
	if (beta <0):
		beta = 0  

	
	ck = sqrt((k-2*beta)*(k+2*beta))/(2*pi*exp(k))
	
	if math.isnan(ck) or math.isnan(k) or math.isnan(beta) or math.isinf(ck) or math.isinf(k) or math.isinf(beta) or (k <0) or (beta< 0):
		print "Error!\n"
		print "ck :", ck , "\nbeta : ", beta, "\nk :", k, "\ngamma :", gamma1, gamma2, gamma3
		#print "\nAborted!"
		#exit (0)
		#put unit kent distribution there!
		k = 0.1
		beta = 0
		
		

	#print type(gamma1), " ", gamma1
	return [k,beta,ck,gamma1,gamma2,gamma3]
	
############################################################################################################################################0
program_name = 'bayesTrainer.py'


if __name__ == '__main__':
        
	usage = "Usage: "+program_name+" [options] training_graph1.db training_graph2.db ... training_graphN.db training_output(distribution_parameters).db\n"+\
            "   or  "+program_name+" --help";
            
	parser = OptionParser(usage)
	
	parser.add_option("--clobber", action="store_true", dest="clobber",
                       default=0, help="overwrite output file")
	
	parser.add_option("--labelLU", type="string", dest="labelLU",
                      help="give the name of labelLU.config file")

	#parser.add_option("--diameter", action="store_true", dest="midpoint",
                       #default=0, help="reads features midpointX, midpointY and midpointZ")
	
	parser.add_option("--length", action="store_true", dest="length",
                       default=0, help="reads feature length")
	
	parser.add_option("--curvature", action="store_true", dest="curvature",
                       default=0, help="reads feature curvature")
	
	parser.add_option("--tortuosity", action="store_true", dest="tortuosity",
                       default=0, help="reads feature tortuosity")
	
	parser.add_option("--midpoint", action="store_true", dest="midpoint",
                       default=0, help="reads features midpointX, midpointY and midpointZ")
	
	parser.add_option("--direction", action="store_true", dest="direction",
                       default=0, help="reads features directionX, directionY and directionZ in +Z hemisphere of the Euclidean coordinate space.")

	parser.add_option("--direction_neg", action="store_true", dest="direction_neg",
                       default=0, help="reads features directionX_neg, directionY_neg and directionZ_neg in -Z hemisphere of the Euclidean coordinate space.")
	
	parser.add_option("--angle", action="store_true", dest="angle",
                       default=0, help="reads features angleX, angleY and angleZ with +X, +Y, +Z axis of the Euclidean coordinate space.")
	
	parser.add_option("--axial_direction", action="store_true", dest="axial_direction",
                       default=0, help="reads features direction cosine in +Z hemisphere and calculate Bingham (Kent) distribution on the unit sphere.")
	
	parser.add_option("--mri_labels", action="store_true", dest="mri_labels",
                       default=0, help="reads features vessel distance from each 40 mri labels.")
	
	parser.add_option("--rel_diameter", action="store_true", dest="rel_diameter",
                       default=0, help="reads features vessel diameter relative to its adjacent vessels.")
	
	parser.add_option("--rel_dir", action="store_true", dest="rel_dir",
                       default=0, help="reads features vessel direction relative to its adjacent vessels.")

	parser.add_option("--classifier_type", action="store", type="string", dest="classifier_type",
                       default="NB", help="specify the classifier type:			1.NB=iid Gaussian NB (default)			2.NB2=univariate Gaussian and Gamma NB			3.GDA=multivariate Gaussian")

	parser.add_option("--PCA", action="store_true", dest="pca",
                       default=0, help="run PCA on feature vectors to make it orthogonal")
	
	parser.add_option("--plot", action="store_true", dest="plot",
                       default=0, help="plot feature distributions")
	
	parser.add_option("--save_plot", action="store_true", dest="save_plot",
                       default=0, help="save plot of feature distributions")
                       
	parser.add_option("--laplace_smooth", action="store_true", dest="laplace_smooth",
                       default=0, help="add noisy data to feature vector for better feature distributions inference")
                       
	parser.add_option("--PCA_CDF_thresh", type="float", dest="PCA_CDF_thresh",
                       default=100.0, help="set threshold % for CDF of PCA eigenvalues, in order to perform feature selection.")
	
	options, args = parser.parse_args()
	

	if len(args) < 2:
		parser.error("incorrect number of arguments")
		
	training_files = args[0:len(args)-1]
	
	output_file = args[len(args)-1]

	#if len(args) != 2:
		#parser.error("incorrect number of arguments")
		
	#featureVectFile, output_file = args

	print '>>> %s [%s]: %s\n' % (time.ctime(time.time()), os.getcwd(), string.join(argv))
	sys.stdout.flush()
	
	if not options.clobber and os.path.exists(output_file):
		raise SystemExit, \
        	"The --clobber option is needed to overwrite an existing file."

	#if options.axial_direction and options.pca:
		#print "ERROR! can't have both axial_direction and PCA options at the same time!\nAborted!\n"
		#exit(0)

	if options.pca and options.axial_direction:
		raise SystemExit, \
			"axial_direction and PCA options cannot be used together!"

	featureNames= ['diameter']

	if options.curvature:
		featureNames.append('curvature')
	
	if options.tortuosity:
		featureNames.append('tortuosity')
	
	if options.length:
		featureNames.append('length')
	
	if options.rel_dir:
		featureNames.append('rel_dir')
	
	if options.rel_diameter:
		featureNames.append('rel_diameter')

	if options.midpoint:
		featureNames.append('midpointX')
		featureNames.append('midpointY')
		featureNames.append('midpointZ')
	
	if options.direction or options.axial_direction:
		featureNames.append('directionX')
		featureNames.append('directionY')
		featureNames.append('directionZ')
	
	if options.direction_neg:
		featureNames.append('directionX_neg')
		featureNames.append('directionY_neg')
		featureNames.append('directionZ_neg')
	
	if options.angle:
		featureNames.append('angleX')
		featureNames.append('angleY')
		featureNames.append('angleZ')
	
	#if options.axial_direction:
		#featureNames.append('directionX')
		#featureNames.append('directionY')
		#featureNames.append('directionZ')

	if options.mri_labels:
		featureNames.append('mri_label_dist')
	
	#if options.proximity:
		#featureNames.append('proximity')

	#if options.anglewref:
		#featureNames.append('anglewref')
		
	
	#### for all possible labels in cerebral vasculature
	if options.labelLU and not os.path.exists(options.labelLU):
		raise SystemExit, ("--labelLU specifed yet file %s does not exist." % options.labelLU)
	if not options.labelLU:
		labelNumeric2Name ={0:"No label",35:"Anterior Cerebral Artery", 191:"R. Middle Cerebral Artry", 190:"L. Middle Cerebral Artry", 2:"R. Intern Carotid Artery", 43:"L. Intern Carotid Artery", 200:"R. Posterior Comm. Artry", 9:"L. Posterior Comm. Artry", 8:"R. Posterior Cereb Artry", 5:"L. Posterior Cereb Artry", 68:"R. Superior Cereb Artery", 227:"L. Superior Cereb Artery", 46:"R. Ant. Inf. Cereb Artry", 12:"L. Ant. Inf. Cereb Artry", 196:"Basilar Artery", 7:"Vertebral Artery", 49:"R. Internal Audit Artery", 45:"L. Internal Audit Artery", 7: "Vertebral Artery" , 3: "R. Paraolivary Artry", 4: "L. Paraolivary Artry" , 11:"Superior Saggital Sinus", 6:"Great Cerbral Vein Galen", 30:"R. Transverse Sinus", 246:"L. Transverse Sinus", 192:"R. Caudal Rhinal Vein", 34:"L. Caudal Rhinal Vein", 20:"R. Rostral Rhinal Vein", 21:"L. Rostral Rhinal Vein", 101:"R. Sigmoid Sinus", 24:"L. Sigmoid Sinus", 58:"R. Longitud. Hippo. Vein", 57:"L. Longitud. Hippo. Vein", 56:"R. Thalamostriate Vein", 54:"L. Thalamostriate Vein", 1:"R. Medial Colicular Vein", 16:"L. Medial Colicular Vein", 170:"Unknown Sinus/Vein #01", 171:"R. Lateral collicular V.", 172:"L. Lateral collicular V.", 250:"L. Unknown Sinus/Vein #2", 251:"R. Unknown Sinus/Vein #2"}											
	else:
		f = open(options.labelLU, 'r')
		lines=[]
		for line in f:
			if ((not line[0]=='#') and (not line=='')):
				lines.append(line)
		labelNumeric2Name ={}
		for l in lines:
			i0=l.index(';')
			i1=l[i0+1:].index(';')+i0+1
			labelNumeric2Name[int(l[0:i0])]=l[i0+1:i1]
		f.close()	
	
	

	##get info from training set	
	#f_handler = shelve.open(featureVectFile,"r")
	#labelNumerics = f_handler['labelNumerics'] 
	#labelNumeric2Name = f_handler['labelNumeric2Name'] 
	#featureNames = f_handler['featureNames'] 
	#feature_vector = f_handler['feature_vector'] 
	#target_vector = f_handler['target_vector'] 
	#adjacency_matrix = f_handler['adjacency_matrix']
	#f_handler.close()
	
	labelNumerics, feature_vector, target_vector, adjacency_matrix, edge_w_indx, edge_neighboring_indx = graph2featureVect(training_files,featureNames)		#featureVect
	
	if options.laplace_smooth:
		for i in range(feature_vector.shape[0]):
			noisy_d = [random.gauss(0,10e-10) for j in range(feature_vector.shape[1])]
			feature_vector	= numpy.vstack([feature_vector, numpy.multiply(array(noisy_d),array(feature_vector[i]))+array(feature_vector[i])])
			target_vector.append(int(target_vector[i]))
			feature_vector	= numpy.vstack([feature_vector, array(feature_vector[i])-numpy.multiply(array(noisy_d),array(feature_vector[i]))])
			target_vector.append(int(target_vector[i]))
		
	axialfeatureidx=[]
	if options.axial_direction:
		axialfeatureNames=['directionX','directionY','directionZ']	
		axialfeatureidx=[featureNames.index('directionX'),featureNames.index('directionY'),featureNames.index('directionZ')]

	if options.pca and len(axialfeatureidx)==0:
		[feature_vector, Y_bar, normalizing_stds, eigenvects, eigenvals, feat_select_indx] = PCA(feature_vector, options.PCA_CDF_thresh )
	elif options.pca and len(axialfeatureidx)>0:	#leave axial direction out of PCA and add at the end
		tmp_indx = range(feature_vector.shape[1])
		for i in axialfeatureidx:
			tmp_indx.remove(i)
		tmp_feature_vector = feature_vector.take(tmp_indx,axis=1)
		[tmp_feature_vector, Y_bar, normalizing_stds, eigenvects, eigenvals, feat_select_indx] = PCA(tmp_feature_vector, options.PCA_CDF_thresh )
		for i in axialfeatureidx:
			tmp_feature_vector = numpy.hstack([tmp_feature_vector, array(feature_vector[:,i])])	
		feature_vector = tmp_feature_vector
		axialfeatureidx = [feature_vector.shape[1]-1,feature_vector.shape[1]-2 ,feature_vector.shape[1]-3 ]

	if options.classifier_type == "NB":
		#### 1. NaiveBayes iid Gaussian	
		print "\n\nNB:\n\n"
		### calculate class priors
		class_prior = Matlib.empty((len(labelNumerics),1),float)
		mean_val = Matlib.empty((len(labelNumerics),feature_vector.shape[1]),float)
		std_val = Matlib.empty((len(labelNumerics),feature_vector.shape[1]),float)
		k_axial = -9999*Matlib.ones((len(labelNumerics),1),float)
		beta_axial = -9999*Matlib.ones((len(labelNumerics),1),float)
		ck_axial = -9999*Matlib.ones((len(labelNumerics),1),float)
		gamma_axial = -9999*Matlib.ones((3*len(labelNumerics),3),float)
		for i in range(len(labelNumerics)):
			l = labelNumerics[i]
			indx = list(pylab.find (numpy.array(target_vector)==l))
			class_prior[i,0]= float(len(indx))/float(len(target_vector))
			if class_prior[i,0] <=0:
				print "ERROR: label ", l, " # ", float(len(indx)), " / ", float(len(target_vector)), " prior " , class_prior[i,0]
				exit(0)
			kent_flag = 0	#has it been calculated once for this label?
			for j in range(feature_vector.shape[1]):
				if j not in axialfeatureidx:
					mean_val[i,j] = numpy.mean (feature_vector[indx,j])
					std_val[i,j] = numpy.std (feature_vector[indx,j], None,None, None, 1)		#numpy.std(a, axis=None, dtype=None, out=None, ddof=0, skipna=False, keepdims=False). The divisor used in calculations is N - ddof, where N represents the number of elements.
				elif j in axialfeatureidx and kent_flag==0:
					kent_flag = 1	#it has been calculated for this label
					feature_vector_arg = Matlib.empty((len(indx),len(axialfeatureidx)),float)		###feature_vector[indx,axialfeatureidx]
					for tmp_j in range(len(axialfeatureidx)):
						feature_vector_arg[:,tmp_j]=feature_vector[indx,axialfeatureidx[tmp_j]]
					
					print "for label: ", l , " has ",len(indx), " observations!" 
					[kent_k,kent_beta,kent_ck,kent_gamma1,kent_gamma2,kent_gamma3]= Kent_calculation (feature_vector_arg)
					k_axial[i,0] = kent_k
					beta_axial[i,0] = kent_beta
					ck_axial[i,0] = kent_ck
					gamma_axial[3*i:3*i+3,0] = kent_gamma1
					gamma_axial[3*i:3*i+3,1] = kent_gamma2
					gamma_axial[3*i:3*i+3,2] = kent_gamma3
					
		f_handler = shelve.open(output_file)
		f_handler['classifier_type'] = "NB"
		if options.pca:
			f_handler['PCA'] = "Y"
			f_handler['feat_select_indx'] = feat_select_indx			
			f_handler['Y_bar'] = Y_bar
			f_handler['normalizing_stds'] = normalizing_stds
			f_handler['eigenvects'] = eigenvects
			f_handler['eigenvals'] = eigenvals
		else:
			f_handler['PCA'] = "N"
		f_handler['class_prior'] = class_prior
		f_handler['labelNumerics'] = labelNumerics
		f_handler['labelNumeric2Name'] = labelNumeric2Name
		f_handler['mean_val'] = mean_val
		f_handler['std_val'] = std_val
		f_handler['axialfeatureidx'] = axialfeatureidx
		if options.axial_direction:
			f_handler['k_axial'] = k_axial
			f_handler['beta_axial'] = beta_axial
			f_handler['ck_axial'] = ck_axial
			f_handler['gamma_axial'] = gamma_axial
	
		f_handler['featureNames'] = featureNames
		f_handler['feature_vector'] = feature_vector 
		f_handler['target_vector'] = target_vector 
		f_handler['adjacency_matrix'] = adjacency_matrix 
		f_handler['history'] = '>>> %s: %s' % (time.ctime(time.time()), string.join(argv))			
		f_handler.close()
		print ("Succesfully wrote %s" %output_file)
		print "class_prior" , class_prior

	elif options.classifier_type == "NB2":
		#### 2. NaiveBayes with univaraite Gaussian and Gamma distribution for non-zero features	
		gamma_feature_indx=[] 
		for f_indx in range(len(featureNames)):
			if min(feature_vector[:,f_indx])>0:
				gamma_feature_indx.append(f_indx)	#if k is relatively large, gamma distribution can be approximated by Gaussian distribution with mean k*theta and variance k*theta^2 (Ross 1998)
				for i in range(len(labelNumerics)):
					if f_indx in gamma_feature_indx:	#if f_index is not already removed
						l = labelNumerics[i]
						indx = pylab.find (numpy.array(target_vector)==l)
						N = float(len(indx))
						s = log( (1.0/N)*sum(feature_vector[indx,f_indx]) ) - ( (1.0/N)*sum(log(feature_vector[indx,f_indx])) )
						k = (3.0 - s +sqrt((s-3.0)**2 + 24.0*s))/(12.0*s)
						if (k> 20):
							gamma_feature_indx.remove(f_indx)		#will be gaussian
							
			
		### calculate class priors
		class_prior = Matlib.empty((len(labelNumerics),1),float)
		mean_val = -9999*Matlib.ones((len(labelNumerics),feature_vector.shape[1]),float)
		std_val = -9999*Matlib.ones((len(labelNumerics),feature_vector.shape[1]),float)
		k_val = -9999*Matlib.ones((len(labelNumerics),feature_vector.shape[1]),float)
		theta_val = -9999*Matlib.ones((len(labelNumerics),feature_vector.shape[1]),float)
		k_axial = -9999*Matlib.ones((len(labelNumerics),2),float)
		dk_axial = -9999*Matlib.ones((len(labelNumerics),1),float)
		mu_axial = -9999*Matlib.ones((3*len(labelNumerics),2),float)
		for i in range(len(labelNumerics)):
			l = labelNumerics[i]
			indx = pylab.find (numpy.array(target_vector)==l)
			class_prior[i,0]= float(len(indx))/float(len(target_vector))
			if class_prior[i,0] <=0:
				print "ERROR: label ", l, " # ", float(len(indx)), " / ", float(len(target_vector)), " prior " , class_prior[i,0]
				exit(0)
			kent_flag = 0	#has it been calculated once for this label?
			for j in range(feature_vector.shape[1]):
				if j not in gamma_feature_indx and j not in axialfeatureidx:
					mean_val[i,j] = numpy.mean (feature_vector[indx,j])
					std_val[i,j] = numpy.std (feature_vector[indx,j], None,None, None, 1)		#numpy.std(a, axis=None, dtype=None, out=None, ddof=0, skipna=False, keepdims=False). The divisor used in calculations is N - ddof, where N represents the number of elements.
				elif j not in axialfeatureidx:
					N = float(len(indx))
					if (N <= 0):
						print ("ERROR! In calculation of Gamma distribution for feature %s: The number of data points in feature vector with label %d is %d.\nAborted!\n" %(featureNames[j],l,len(indx)))
						exit (0)
					if numpy.isinf(log( (1.0/N)*sum(feature_vector[indx,j]))) or  numpy.isnan(log( (1.0/N)*sum(feature_vector[indx,j]))):
						print ("ERROR! In calculation of Gamma distribution for feature %s: The log of mean of data points in feature vector (%f) with label %d is nan/inf.\nAborted!\n" %(featureNames[j],(1.0/N)*sum(feature_vector[indx,j]),l))
						exit (0)
					if numpy.isinf(sum(log(feature_vector[indx,j]))) or  numpy.isnan(sum(log(feature_vector[indx,j]))):
						print ("ERROR! In calculation of Gamma distribution for feature %s: sum of log of data points in feature vector (%f) with label %d is nan/inf.\nAborted!\n" %(featureNames[j],sum(log(feature_vector[indx,j])),l))
						exit (0)
					s = log( (1.0/N)*sum(feature_vector[indx,j]) ) - ( (1.0/N)*sum(log(feature_vector[indx,j])) )
					k_val[i,j] = (3.0 - s +sqrt((s-3.0)**2 + 24.0*s))/(12.0*s)			#k approximation within 1.5% of the correct value.
					theta_val[i,j] = (1.0/ (k_val[i,j]*N)) * sum(feature_vector[indx,j])
				elif j in axialfeatureidx and kent_flag==0:
					kent_flag = 1	#it has been calculated for this label
					feature_vector_arg = Matlib.empty((len(indx),len(axialfeatureidx)),float)		###feature_vector[indx,axialfeatureidx]
					for tmp_j in range(len(axialfeatureidx)):
						feature_vector_arg[:,tmp_j]=feature_vector[indx,axialfeatureidx[tmp_j]]

					[kent_k,kent_beta,kent_ck,kent_gamma1,kent_gamma2,kent_gamma3]= Kent_calculation (feature_vector_arg)
					k_axial[i,0] = kent_k
					beta_axial[i,0] = kent_beta
					ck_axial[i,0] = kent_ck
					gamma_axial[3*i:3*i+3,0] = kent_gamma1
					gamma_axial[3*i:3*i+3,1] = kent_gamma2
					gamma_axial[3*i:3*i+3,2] = kent_gamma3
					
					
		f_handler = shelve.open(output_file)
		f_handler['classifier_type'] = "NB2"
		if options.pca:
			f_handler['PCA'] = "Y"
			f_handler['feat_select_indx'] = feat_select_indx			
			f_handler['Y_bar'] = Y_bar
			f_handler['normalizing_stds'] = normalizing_stds
			f_handler['eigenvects'] = eigenvects
			f_handler['eigenvals'] = eigenvals
		else:
			f_handler['PCA'] = "N"
		f_handler['class_prior'] = class_prior
		f_handler['labelNumerics'] = labelNumerics
		f_handler['labelNumeric2Name'] = labelNumeric2Name
		f_handler['mean_val'] = mean_val
		f_handler['std_val'] = std_val
		f_handler['k_val'] = k_val
		f_handler['theta_val'] = theta_val
		f_handler['axialfeatureidx'] = axialfeatureidx
		if options.axial_direction:
			f_handler['k_axial'] = k_axial
			f_handler['beta_axial'] = beta_axial
			f_handler['ck_axial'] = ck_axial
			f_handler['gamma_axial'] = gamma_axial
		
		f_handler['gamma_feature_indx'] = gamma_feature_indx
		f_handler['featureNames'] = featureNames
		f_handler['feature_vector'] = feature_vector 
		f_handler['target_vector'] = target_vector 
		f_handler['adjacency_matrix'] = adjacency_matrix 
		f_handler['history'] = '>>> %s: %s' % (time.ctime(time.time()), string.join(argv))			
		f_handler.close()
		print ("Succesfully wrote %s" %output_file)
		print "class_prior" , class_prior


	
	elif options.classifier_type == "GDA":	
		#### 3. GDA with Gaussian distribution	(generalized discriminate analysis (GDA))
		### calculate class priors
		class_prior = Matlib.empty((len(labelNumerics),1),float)
		mean_vector = Matlib.empty((len(labelNumerics),feature_vector.shape[1]),float)
		covariance_mat = Matlib.empty((len(labelNumerics)*feature_vector.shape[1],feature_vector.shape[1] ),float)
		###covariance_mat = covariance_mat.reshape(feature_vector.shape[1],feature_vector.shape[1],len(labelNumerics) )
		for i in range(len(labelNumerics)):
			l = labelNumerics[i]
			indx = pylab.find (numpy.array(target_vector)==l)
			class_prior[i,0]= float(len(indx))/float(len(target_vector))
			if class_prior[i,0] <=0:
				print "ERROR: label ", l, " # ", float(len(indx)), " / ", float(len(target_vector)), " prior " , class_prior[i,0]
				exit(0)
			for j in range(feature_vector.shape[1]):
				mean_vector[i,j] = numpy.mean (feature_vector[indx,j])
			covariance_mat[i*feature_vector.shape[1]:(i+1)*feature_vector.shape[1],:] = numpy.cov (feature_vector[indx,:],None,0)	#numpy.cov(m, y=None, rowvar=1, bias=0, ddof=None). rowvar=0 ach column represents a variable, while the rows contain observations. bias: Default normalization is by (N - 1), where N is the number of observations.   

		f_handler = shelve.open(output_file)
		f_handler['classifier_type'] = "GDA"
		if options.pca:
			f_handler['PCA'] = "Y"
			f_handler['feat_select_indx'] = feat_select_indx			
			f_handler['Y_bar'] = Y_bar
			f_handler['normalizing_stds'] = normalizing_stds
			f_handler['eigenvects'] = eigenvects
			f_handler['eigenvals'] = eigenvals
		else:
			f_handler['PCA'] = "N"
		f_handler['class_prior'] = class_prior
		f_handler['labelNumerics'] = labelNumerics
		f_handler['labelNumeric2Name'] = labelNumeric2Name
		f_handler['mean_vector'] = mean_vector
		f_handler['covariance_mat'] = covariance_mat
		f_handler['featureNames'] = featureNames
		f_handler['feature_vector'] = feature_vector 
		f_handler['target_vector'] = target_vector 
		f_handler['adjacency_matrix'] = adjacency_matrix 
		f_handler['history'] = '>>> %s: %s' % (time.ctime(time.time()), string.join(argv))			
		f_handler.close()
		print ("Succesfully wrote %s" %output_file)


		print "class_prior" , class_prior

#### in order to plot should enable plotting first by showing nonsense plot and close it, also give a feature_vect of a selected feature and label with its mean,std or k,theta to the correct distribution plot function
	if options.plot:
		directory = (os.path.abspath(os.path.dirname(output_file))+"/plots_"+os.path.basename(output_file)[:-3])
		if options.save_plot:
			if not os.path.exists (directory):
				os.makedirs (directory)

		#### plotting the calculated pdf for each of class for each feature in the training set along with individual points!
		#x=arange(0,2,0.01)
		#y=2*sin(2*pi*(x-1/4))
		#plot(x,y)
		#show()
		label_list = [35, 190, 191, 11, 246, 30]
		color_list = ['r','g', 'c', 'm', 'k', 'b','y']
		c=1
		#x=[float(i)/400 for i in range(1,100)]
		#x=arange(0,0.25,0.001)

		c= 1
		if options.classifier_type == "NB":
			for j in range(len(featureNames)):
				pylab.figure(c)
				pylab.hold(True)
				f = featureNames[j]
				for l_i in range(len(label_list)):
					l = label_list[l_i]
					if l in labelNumerics:
						indx = pylab.find (numpy.array(target_vector)==l)
						i = labelNumerics.index(l)
						plot_mean_value = mean_val[i,j]
						plot_std_value = std_val[i,j]		#numpy.std(a, axis=None, dtype=None, out=None, ddof=0, skipna=False, keepdims=False). The divisor used in calculations is N - ddof, where N represents the number of elements.
						plot_feat_vect = feature_vector[indx,j]
						plot_univariate_gaussian(f,color_list[l_i],l,plot_feat_vect,plot_mean_value, plot_std_value,c,options.save_plot,directory)
				c =c+1
				t= ('Univariate Gaussian fit for feature %s ' %(f))
				print t
				pylab.title(t)
				pylab.legend(loc=1)
				if options.save_plot:
					figname= ("%s/gaussian_feature_%s.png" %(directory,f))
					pylab.savefig(figname, format="png")

					#else:
						#print "target_vector=", target_vector,"\nlabelNumerics:", labelNumerics,"\nl:",l," ",indx

		c= 1
		if options.classifier_type == "NB2":
			for j in range(len(featureNames)):
				if j not in gamma_feature_indx:
					for l_i in range(len(label_list)):
						l = label_list[l_i]
						if l in labelNumerics:
							indx = pylab.find (numpy.array(target_vector)==l)
							i = labelNumerics.index(l)
							f = featureNames[j]
							plot_mean_value = mean_val[i,j]
							plot_std_value = std_val[i,j]		#numpy.std(a, axis=None, dtype=None, out=None, ddof=0, skipna=False, keepdims=False). The divisor used in calculations is N - ddof, where N represents the number of elements.
							plot_feat_vect = feature_vector[indx,j]
							pylab.figure(c)
							pylab.hold(True)
							plot_univariate_gaussian(f,color_list[l_i], l, plot_feat_vect, plot_mean_value, plot_std_value,c,options.save_plot,directory)
							c = c+1

				else:
					f = featureNames[j]
					plot_k_value = k_val[i,j]
					plot_theta_vect = theta_val[i,j]		#numpy.std(a, axis=None, dtype=None, out=None, ddof=0, skipna=False, keepdims=False). The divisor used in calculations is N - ddof, where N represents the number of elements.
					plot_feat_vect = feature_vector[indx,j]
					pylab.figure(c)
					pylab.hold(True)
					plot_gamma(f, l, plot_feat_vect, plot_k_value, plot_theta_vect,c,options.save_plot,directory)
					c = c+1
					
		#show()
		
###convert this code to be a function
###need to perform smoothing on the label and features ??? if some observations are 0?

	###### plotting the calculated pdf for each of class for each feature in the training set along with individual points!
	##x=arange(0,2,0.01)
	##y=2*sin(2*pi*(x-1/4))
	##plot(x,y)
	##show()
	##label_list = [35, 190, 191, 11, 246, 30]
	##c=1
	###x=[float(i)/400 for i in range(1,100)]
	##x=arange(0,0.25,0.001)
##for l in label_list:
	##for f in featureNameList:
			##if f not in feature_listtype and f not in feature_dicttype.keys() and f in feature_negatives:
					##t= ('feature %s for label %d' %(f,l))
					##figure(c)
					##c+=1
					##hist(featureVector[l][f][:-2],bins=100,normed=1)
					##hold(True)
					##y=[]
					##x=arange(min(featureVector[l][f][:-2]),max(featureVector[l][f][:-2]),0.001)
					##for i in x:
						##y.append(normal_pdf ( i, featureVector[l][f][-2]['mean'] , featureVector[l][f][-1]['std']))
					##plt.plot(x,y,'r')
					##title(t)
			##elif f not in feature_listtype and f not in feature_dicttype.keys() and f not in feature_negatives:
					##t= ('feature %s for label %d' %(f,l)) 
					##figure(c)
					##c+=1
					##hist(featureVector[l][f][:-2],bins=100,normed=1)
					##hold(True)
					##y=[]
					##x=arange(min(featureVector[l][f][:-2]),max(featureVector[l][f][:-2]),0.001)
					##for i in x:
						##y.append(gamma_distribution ( i, featureVector[l][f][-2]['gamma_k'] , featureVector[l][f][-1]['gamma_theta']))
					##plt.plot(x,y,'r')
					##title(t)	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
