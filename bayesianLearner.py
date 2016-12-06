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
#  Created May 18, 2015 
#  Modified June 12, 2015 added the STUN optimization
#  Modified June 30, 2015 update total system energy after each label swap in GD, SA and STUN. This should both improve the accuracy and speed of the calculations. In GD, we also added a condition to swap labels only if deltaE < 0
#  Modified Aug 26, 2015 added estimated_label_potentials property to write_posterior_prob_map, make the estimate_label log scale in write_posterior_prob_map, add write_prior_prob_map
#  Modified Oct 21, 2015 added the estimated_target_vector_prob to the estimate_label edge_property of GD, SA and STUN outputs, corrected wrong indent in calc_likelihood and calc_log_likelihood. PCA feature reduction use np.dot for matrix product  
#  Sahar Ghanavati


from sys import argv
from scipy import special
from scipy.optimize import fmin,fmin_ncg
import scipy.linalg
from vessel_tracking import graph_analysis
from numpy import *
import numpy as np
import numpy.matlib as Matlib
from optparse import OptionParser, Option, OptionValueError
from minc_util.progress import progress_report
import math, os, shelve, string ,time, sys
import commands
import copy
from operator import itemgetter
#import matplotlib.pyplot as plt
#from pylab import *		#needed for command find
import pylab		#needed for command find
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
                #target_vector= np.vstack([target_vector, array(g.edge_property(e,'label'))])
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

                feature_vector	= np.vstack([feature_vector, array(features)])		#hstack #a = matrix([[10,20,30]]); a=append(a,[[1,2,3]],axis=0); a=append(a,[[15],[15]],axis=1)
            edge_w_indx[e_cnt]= e 
            e_cnt = e_cnt+1
    return labelNumerics, feature_vector, target_vector, edge_w_indx

def make_adjmatrix(training_files,labelNumerics,edge_w_indx):
    #### !!find neighbouring edges and in iterations find their
    graphs=[]
    for training_file in training_files:
        try:
            graphs.append( graph_analysis.input_graph(training_file))  	# open graph data and copy contents to memory
            print ("Succefully read in %s\n" %training_file)	
        except:
            print("Error reading in %s\n" %training_file)

    adjacency_matrix = (1.0/len(labelNumerics))*Matlib.ones((len(labelNumerics),len(labelNumerics)+1),float)		#### smoothing (for 0s in the adj_matrix)! to be minimum of 1 occurance for the nieghbourhood!	
    for g in graphs:
        for e in g.edge_list():
            if g.edge_property(e,'label')>0 and g.edge_property(e,'label') in labelNumerics:
                indx0= labelNumerics.index(g.edge_property(e,'label'))					
                e1_neighbours= g.vertices[e[0]].edges
                for v in e1_neighbours:
                    if v!=e[1]:
                        neighbor_edge=tuple((min(e[0],v),max(e[0],v)))
                        if g.edge_property(neighbor_edge,'label')>0 and g.edge_property(neighbor_edge,'label') in labelNumerics:
                            indx1= labelNumerics.index(g.edge_property(neighbor_edge,'label'))
                            adjacency_matrix[indx0,indx1]= adjacency_matrix[indx0,indx1]+1
                e2_neighbours= g.vertices[e[1]].edges
                for v in e2_neighbours:
                    if v!=e[0]:
                        neighbor_edge=tuple((min(e[1],v),max(e[1],v)))
                        if g.edge_property(neighbor_edge,'label')>0 and g.edge_property(neighbor_edge,'label') in labelNumerics:
                            indx1= labelNumerics.index(g.edge_property(neighbor_edge,'label'))
                            adjacency_matrix[indx0,indx1]= adjacency_matrix[indx0,indx1]+1					
                if len(g.vertices[e[0]].edges)==1:		#end point
                    adjacency_matrix[indx0,len(labelNumerics)] +=1
                if len(g.vertices[e[1]].edges)==1:		#end point
                    adjacency_matrix[indx0,len(labelNumerics)] +=1					

    #### adjacency matrix normalization 
    for i in range(adjacency_matrix.shape[0]):
        adjacency_matrix[i,:]= adjacency_matrix[i,:]/sum(adjacency_matrix[i,:])			
    #print_adj_mat (adjacency_matrix,labelNumerics)	
    return adjacency_matrix	

def find_neighbors (g, edge_w_indx):
    edge_neighboring_indx = {}
    e_cnt = 0
    for e in g.edge_list():
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
        edge_neighboring_indx[e_cnt]=neighbor_indx
        e_cnt = e_cnt+1
    return edge_neighboring_indx

def machineEpsilon(func=float):
    machine_epsilon = func(1)
    while func(1)+func(machine_epsilon) != func(1):
        machine_epsilon_last = machine_epsilon
        machine_epsilon = func(machine_epsilon) / func(2)
    return machine_epsilon_last

def PCA (feature_vector,CDF_thresh,Y_bar=-inf, normalizing_stds = -inf, eigenvects=-inf,eigenvals=-inf,M_ind=-1):
    #### 1. centering the data to mean=0
    if type(Y_bar)==float:
        Y_bar = np.mean(feature_vector,0)	#mean of each column of feature_vector = mean of each feature data points
    feature_vector_cntr = feature_vector - tile(Y_bar,(feature_vector.shape[0],1))
    #### 2. normalizing the data to std=1
    if type(normalizing_stds)==float:
        normalizing_stds = np.std(feature_vector_cntr,0)
    feature_vector_cntrN = feature_vector_cntr / tile(normalizing_stds,(feature_vector_cntr.shape[0],1))

    #""" SVD:http://docs.scipy.org/doc/numpy/reference/generated/numpy.linalg.svd.html
        # X.H is conjugate transpose matrix of X
        # if matrix X is orthogonal: 1.invertible X^T (transpose) = X^-1 (inverse) 2.unitary X^-1 (inv) =X* (X.H conjugate transpose) 3.normal X.H X = X X.H
        # feature_vector_cntrN = u * np.diag(s) * v.H
        # U, s, V = np.linalg.svd(feature_vector_cntrN, full_matrices=False) : U=u, V=v.H
        # The rows of V are the eigenvectors of A.H A which is of size of features => it is what we are looking for
        #SO EIGENVECTORS ARE COLUMNS OF V.transpose
        # The columns of U are the eigenvectors of A A.H which is of size of datapoints
        # For row i in V and column i in U, the corresponding eigenvalue is s[i]**2.
        # Decomposition(data mapping) is T = A eigenvects = A V.transpose  http://en.wikipedia.org/wiki/Principal_component_analysis#Further_components and http://en.wikipedia.org/wiki/Principal_component_analysis#Dimensionality_reduction
    #"""
    if type(eigenvects)==float or type(eigenvals)==float:
        #### test: cov_matrix symmetric positive semidefinite (PSD) matrix => has a nonnegative eigenvalue
        cov_matrix = np.cov(feature_vector_cntrN,None,0)    #cov(m, y=None, rowvar=1, bias=0, ddof=None) rowvar=0:each column represents a variable, while the rows contain observations
        #print type(cov_matrix), "\ncov_matrix:", cov_matrix
       #'''
        #### ERROR: This is wrong because we shouldn't calculate the SVD of feature_vector_cntrN but the cov_matrix!!!
        #U, s, V = np.linalg.svd(feature_vector_cntrN, full_matrices=False)   #feature_vector_cntrN= U * np.diag(s) * V
        ##The rows of v are the eigenvectors of a.H a. The columns of u are the eigenvectors of a a.H.
        ##if feature_vector_cntrN is Nxm # U is descending ordered eigenvectors of Nxm
        #sqrt_eigenvals = s    # s is descending ordered list of eigenvals of len m
        #eigenvals = [e*e for e in s]      
        #eigenvects = V.H #which is equal to V.transpose = inv(V)  
        ##eigenvects = U #U contains vectors of singular components  
        #eigenvalsmat = np.diag(s)*np.diag(s)    #is diagonal matrix of mxm where eigenvals are on the main diagonal and the rest is 0
        #'''
        
        D,V = np.linalg.eig(cov_matrix)		# a: array_like shape (M, M), D:eigenvalues ndarray shape (M,) NOT ordered, V: eignevectors ndarray shape (M, M) The normalized (unit “length”) eigenvectors, column V[:,i] is the eigenvector corresponding to the eigenvalue D[i]			 
        #i = np.argsort(D)	#Returns the indices that would sort an array, ascending
        i = np.argsort(np.array(D))	#D[::-1] #Returns the indices that would sort an array, descending
        inv_i = i[::-1]
        eigenvals = D[inv_i]
        eigenvects = V[:,inv_i]
        #print type(eigenvals)
        #print type(eigenvects)
        #print "V:", V 
        #print "sorted eigen values are (should be >0) : " , eigenvals		#### test: eigvals should be >0 => orthogonal_feature_vector: real
        #print "sorted eigen vects are (should be >0) : " , eigenvects		#### test: eigvals should be >0 => orthogonal_feature_vector: real

    if M_ind==-1:
        cumulative_eigval = np.cumsum(eigenvals)
        cumulative_eigval = cumulative_eigval/cumulative_eigval[-1]*100.0	# % CDF
        M_ind = min(pylab.find(cumulative_eigval >= CDF_thresh))
    W = eigenvects[:,0:M_ind]
    orthogonal_feature_vector = np.dot(np.matrix(feature_vector_cntrN) , W)   #if W is eigenvectors T = X W maps a data vector x(i) from an original space of p variables to a new space of p variables which are uncorrelated over the dataset.

    print "\n\nM_ind:",M_ind, " out of ", np.matrix(eigenvects).shape[0], " features and feature vect shape of ", np.matrix(feature_vector_cntrN).shape
    print "significant eigenvals:", eigenvals[0:M_ind]
    #print (np.matrix(W).H).shape
    #print np.matrix(feature_vector_cntrN).shape
    #print np.matrix(np.diag(1.0/np.sqrt(eigenvals))).shape
    #print ((np.matrix(feature_vector_cntrN) * np.matrix(np.diag(1.0/np.sqrt(eigenvals)))).transpose()).shape

    #### orthogonal_feature_vector = feature_vector_cntr * V * D^(-.5)??   
    ##based on karhunen-Leove theorem and data projection from http://en.wikipedia.org/wiki/Principal_component_analysis:
    ## projected data = conjugate_transpose(eigenvects)* (feature_vects/Diag(sqrt(eigenvals)))
    ## a.H conjugate transpose
    #orthogonal_feature_vector = np.matrix(W).H * (np.matrix(feature_vector_cntrN) * np.matrix(np.diag(1.0/np.sqrt(eigenvals)))).transpose()
    #orthogonal_feature_vector = orthogonal_feature_vector.transpose()
    #### test: correctness of PCA: cov(orthogonal_feature_vector) almost = I
    return [orthogonal_feature_vector, np.array(Y_bar), np.array(normalizing_stds), np.matrix(eigenvects), np.array(eigenvals), M_ind]


def univaraite_gaussian(x,mean_val,std_val):
    p = (1.0/np.sqrt(2.0*np.pi*pow(std_val,2))) * np.exp( pow(x-mean_val,2) / (-2.0*pow(std_val,2)) )
    return p	


def error_calculation(target_vector, estimated_target_vector, iteration, g,edge_w_indx, E_current,T=0):
    ####only calculate the error if the ground truth label is available!
    availabel_labels = list(np.unique(np.array(target_vector)))
    labeled_num = 0
    error_num = 0
    second_corr=0
    third_corr=0
    fourth_corr=0
    total_volume=0
    error_volume=0
    if len(availabel_labels)>1:	#otherwise test set was not labeled and all target_vector is 0
        for i in range(len(target_vector)):
            if (target_vector[i] > 0) :		#not label= 0 or label= -1 (for label propagation!!!)
                labeled_num += 1
                if target_vector[i] != estimated_target_vector[i] : #[i] for ith edge
                    error_num += 1
                    if (target_vector[i] ==  local_potentials[i][1][0]) :		#not label= 0 or label= -1 (for label propagation!!!)
                        second_corr=second_corr+1
                    elif (target_vector[i] ==  local_potentials[i][2][0]) :		#not label= 0 or label= -1 (for label propagation!!!)
                        third_corr=third_corr+1
                    elif (target_vector[i] ==  local_potentials[i][2][0]) :		#not label= 0 or label= -1 (for label propagation!!!)
                        fourth_corr=fourth_corr+1
        #if ((labeled_num > 0) and (iteration%1==0)):	#otherwise test set was not labeled
            #print (" \n\nIteration %d: Out of %d edges that were labeled, %d were labelled correctly with Naive-Bayes. The error is %f%%. Recognition Rate is %f%%." %(iteration,labeled_num,labeled_num-error_num, 100*float(error_num)/float(labeled_num), 100-100*float(error_num)/float(labeled_num) ) )
            #print ("Out of %d labelled edges error_num %d=%f%%, second_corr %d=%f%%, third_corr %d=%f%%, fourth_corr %d=%f%%.\n\n" %(labeled_num, error_num, 100.0*float(error_num)/float(labeled_num), second_corr, 100.0*float(second_corr)/float(labeled_num), third_corr, 100.0*float(third_corr)/float(labeled_num), fourth_corr, 100.0*float(fourth_corr)/float(labeled_num)))  
    
    
        #### error by volume
        for i in range(len(g.edge_list())):
            e=edge_w_indx[i]
            if not e in g.edge_list() or 'cyl_radius' not in g.edge_properties(e).keys() or 'cyl_height' not in g.edge_properties(e).keys():
                print e,' ', g.edge_properties(e).keys()
            else:
            #if e in g.edge_list() and 'cyl_radius' in g.edge_properties(e).keys() and 'cyl_height' in g.edge_properties(e).keys():
                e_radius = g.edge_property(e,'cyl_radius')
                e_height = g.edge_property(e,'cyl_height')
                e_volume = 0
                for j in range(len(e_radius)):
                    e_volume += cyl_volume_calculation(e_radius[j],e_height[j])
                total_volume += e_volume
                if target_vector[i] != estimated_target_vector[i]:
                    error_volume += e_volume
        
        if ((labeled_num > 0) and (total_volume > 0) and (iteration%1==0)):	#otherwise test set was not labeled
            print ((" \n\nIteration %d" %iteration)+"(T=" +str(T)+(",E=%f): Out of %d edges that were labeled, %d were labelled correctly with Naive-Bayes. The error is %f%% (volumetric error %f%%). Recognition Rate is %f%% (volumetric RR %f%%)." %(E_current,labeled_num,labeled_num-error_num, 100*float(error_num)/float(labeled_num), 100*float(error_volume)/float(total_volume), 100-100*float(error_num)/float(labeled_num), 100-100*float(error_volume)/float(total_volume) )) )
            print ("Out of %d labelled edges error_num %d=%f%%, second_corr %d=%f%%, third_corr %d=%f%%, fourth_corr %d=%f%%.\n\n" %(labeled_num, error_num, 100.0*float(error_num)/float(labeled_num), second_corr, 100.0*float(second_corr)/float(labeled_num), third_corr, 100.0*float(third_corr)/float(labeled_num), fourth_corr, 100.0*float(fourth_corr)/float(labeled_num)))  
                
    return [error_num,labeled_num,total_volume,error_volume]

def error_calculation_by_label(RR_by_label,target_vector, estimated_target_vector, local_potentials,feature_vector, iteration,edge_neighboring_indx, g,edge_w_indx):
    ####only calculate the error if the ground truth label is available!
    availabel_labels = list(np.unique(np.array(target_vector)))
    if 0 in availabel_labels:
        availabel_labels.remove(0)
    each_label_err={}
    each_label_num={}
    each_label_err90={}		#label error of major edges with diameter in [90-100] percentile range
    each_label_num90={}
    each_label_mindiameter={}
    each_label_maxdiameter={}
    for l in availabel_labels:
        each_label_err[l] = 0
        each_label_num[l] = 0
        each_label_err90[l] = 0
        each_label_num90[l] = 0
        indx = pylab.find (np.array(target_vector==l))
        each_label_mindiameter[l]= np.min (feature_vector[indx,0])
        each_label_maxdiameter[l]= np.max (feature_vector[indx,0])
    for i in range(len(target_vector)):
        if target_vector[i] in availabel_labels:		#not label= 0 or label= -1 (for label propagation!!!)
            l=target_vector[i]
            each_label_num[l] += 1
            if target_vector[i] != estimated_target_vector[i] :
                each_label_err[l] += 1
            if (feature_vector[i,0]> (0.75*(each_label_maxdiameter[l]-each_label_mindiameter[l]) + each_label_mindiameter[l]) ):
                each_label_num90[l] += 1
                if target_vector[i] != estimated_target_vector[i] :
                    each_label_err90[l] += 1	
    
    #### error by volume
    each_label_total_volume={}
    each_label_error_volume={}
    for l in availabel_labels:
        each_label_total_volume[l]=0
        each_label_error_volume[l]=0
    for i in range(len(g.edge_list())):
        e=edge_w_indx[i]
        if e in g.edge_list() and 'cyl_radius' in g.edge_properties(e).keys() and 'cyl_height' in g.edge_properties(e).keys():
            e_radius = g.edge_property(e,'cyl_radius')
            e_height = g.edge_property(e,'cyl_height')
            e_volume = 0
            for j in range(len(e_radius)):
                e_volume += cyl_volume_calculation(e_radius[j],e_height[j])
            l=target_vector[i]		#g.edge_property(g.edge_list()[e],'label')
            if l>0:
                each_label_total_volume[l] += e_volume
                if target_vector[i] != estimated_target_vector[i]:
                    each_label_error_volume[l] += e_volume

    for l in availabel_labels:
        if l not in RR_by_label.keys():
            RR_by_label[l] = [each_label_num[l],each_label_total_volume[l]]
        print  ("\nOut of %d edges (volume %f) with label %d, %d were labelled correctly (correct volume %f). Error rate is %f%% (volumetric error %f%%). Recognition Rate is %f%% (volumetric RR %f%%)." %(each_label_num[l],each_label_total_volume[l],l,each_label_num[l]-each_label_err[l],each_label_total_volume[l]-each_label_error_volume[l],100*float(each_label_err[l])/float(each_label_num[l]),100*float(each_label_error_volume[l])/float(each_label_total_volume[l]),100-100*float(each_label_err[l])/float(each_label_num[l]),100-100*float(each_label_error_volume[l])/float(each_label_total_volume[l]) ))                                                 
        RR_by_label[l].append(100-100*float(each_label_err[l])/float(each_label_num[l]))
        RR_by_label[l].append(100-100*float(each_label_error_volume[l])/float(each_label_total_volume[l]))
        RR_by_label[l].append(int(each_label_num[l]-each_label_err[l]))		
        RR_by_label[l].append(each_label_total_volume[l]-each_label_error_volume[l])
    return RR_by_label

def cyl_volume_calculation(r,h):
    vol = h*(pi*r*r)
    return vol
    
def confusion_matrix_calculation(target_vector, estimated_target_vector,labelNumerics):
    print "\n\n*****************************************************\nconfusion_matrix:\n"
    availabel_labels = list(np.unique(np.array(target_vector)))
    if 0 in availabel_labels:
        availabel_labels.remove(0)
    
    confusion_matrix = Matlib.zeros((len(availabel_labels),len(labelNumerics)),float)		
    for i in range(len(target_vector)):
        l=target_vector[i]
        if l>0:
            l2=estimated_target_vector[i]
            #print availabel_labels
            #print l
            #print l2
            i= list(availabel_labels).index(l)
            j= list(labelNumerics).index(l2)
            confusion_matrix[i,j]+=1
        
    #### adjacency matrix normalization 
    for i in range(confusion_matrix.shape[0]):
        confusion_matrix[i,:]= confusion_matrix[i,:]/sum(confusion_matrix[i,:])
    
    if (len(labelNumerics)> 14):
        print "\t",
        for j in range(15):
            print int(labelNumerics[j]),"\t",
        print " "	
        for i in range(confusion_matrix.shape[0]):
            print "\n", int(availabel_labels[i]),"\t",
            for j in range(15):
                print "%.3f"%confusion_matrix[i,j],"\t",
        print "\n\n"
    else:
        print "\t",
        for j in range(len(labelNumerics)):
            print int(labelNumerics[j]),"\t",
        print " "	
        for i in range(confusion_matrix.shape[0]):
            print "\n", int(availabel_labels[i]),"\t",
            for j in range(len(labelNumerics)):
                print "%.3f" %confusion_matrix[i,j],"\t",
        print "\n\n"
        return 0
    
    if (len(labelNumerics)> 29):
        print "\t",
        for j in range(15,30):
            print int(labelNumerics[j]),"\t",
        print " "	
        for i in range(confusion_matrix.shape[0]):
            print "\n", int(availabel_labels[i]),"\t",
            for j in range(15,30):
                print "%.3f"%confusion_matrix[i,j],"\t",
        print "\n\n"
    else:
        print "\t",
        for j in range(15,len(labelNumerics)):
            print int(labelNumerics[j]),"\t",
        print " "	
        for i in range(confusion_matrix.shape[0]):
            print "\n", int(availabel_labels[i]),"\t",
            for j in range(15,len(labelNumerics)):
                print "%.3f" %confusion_matrix[i,j],"\t",
        print "\n\n"
        return 0
        
    if (len(labelNumerics)> 44):
        print "\t",
        for j in range(30,45):
            print int(labelNumerics[j]),"\t",
        print " "	
        for i in range(confusion_matrix.shape[0]):
            print "\n", int(availabel_labels[i]),"\t",
            for j in range(30,45):
                print "%.3f"%confusion_matrix[i,j],"\t",
        print "\n\n"
    else:	
        print "\t",
        for j in range(30,len(labelNumerics)):
            print int(labelNumerics[j]),"\t",
        print " "	
        for i in range(confusion_matrix.shape[0]):
            print "\n", int(availabel_labels[i]),"\t",
            for j in range(30,len(labelNumerics)):
                print "%.3f" %confusion_matrix[i,j],"\t",
        print "\n\n"
        return 0
        

    print "\t",
    for j in range(45,len(labelNumerics)):
        print int(labelNumerics[j]),"\t",
    print " "	
    for i in range(confusion_matrix.shape[0]):
        print "\n", int(availabel_labels[i]),"\t",
        for j in range(45,len(labelNumerics)):
            print "%.3f" %confusion_matrix[i,j],"\t",
    print "\n\n"

    #print  "\nlabelNumerics: ",labelNumerics 
    return 0
    
    
def calc_globalE(estimated_target_vector,local_potentials,adjacency_matrix,edge_neighboring_indx,labelNumerics):
    ##local_potentials contains potential of all labels for each edge based on local features! => P(l/features_of_edge)
    ##adjacency_matrix rows = P(l_neighbor/l_this_edge) = P(l_this_edge/l_neighbor)*P(l_neighbor/neighbor's_features)/P(l_this_edge)
    ##for each edge E = P(l_this_edge/features_of_edge)*P(l_this_edge/l_neighbors)*P(l_neighbor/neighbor's_features) = P(l_this_edge/features_of_edge)*P(l_neighbors/l_this_edge)	### here I am neglecting the term  P(l_this_edge) in the bayes rule of adj_matrix term!!!					
    neg_E = 0
    for i in range(len(estimated_target_vector)):
        l_e = estimated_target_vector[i]
        if l_e>0 and l_e in labelNumerics:
            local_labels = [local_potentials[i][j][0] for j in range(len(local_potentials[i]))]
            #if (-1/(machineEpsilon(float)))!=local_potentials[i][local_labels.index(l_e)][1]:
            neg_E = neg_E + local_potentials[i][local_labels.index(l_e)][1]
            neighb_idx = edge_neighboring_indx[i]
            mat_idx0 = labelNumerics.index(l_e)
            for j in neighb_idx:
                l_neighb = estimated_target_vector[j]
                if l_neighb>0 and l_neighb in labelNumerics:
                    mat_idx1 = labelNumerics.index(l_neighb)
                    neg_E = neg_E + np.log10(adjacency_matrix[mat_idx0,mat_idx1])
            if (len(neighb_idx)==2):		# Y edge that has an end-point 
                neg_E = neg_E + np.log10(adjacency_matrix[mat_idx0,len(labelNumerics)])		#neighbor=End_point
                
    E = -neg_E
    if np.isnan(E) or np.isinf(E):
        E = 1.0/machineEpsilon()
    return E


def calc_deltaE(current_e,current_l,new_l,estimated_target_vector,local_potentials,adjacency_matrix,edge_neighboring_indx,labelNumerics):
    if current_l<=0 or new_l<=0:
        print "\nERROR! label <=0, current_l:",current_l, " new_l:",new_l,"\nAborted!"
        exit(0)
    
    U_current = 0
    U_new = 0

    local_labels = [local_potentials[current_e][j][0] for j in range(len(local_potentials[current_e]))]
    neighb_idx = edge_neighboring_indx[current_e]
    
    U_current += local_potentials[current_e][local_labels.index(current_l)][1]
    mat_idx0 = labelNumerics.index(current_l)
    for j in neighb_idx:
        l_neighb = estimated_target_vector[j]
        if l_neighb>0:
            mat_idx1 = labelNumerics.index(l_neighb)
            U_current += np.log10(adjacency_matrix[mat_idx0,mat_idx1])   #once for the current edge and once for its neighbor
            U_current += np.log10(adjacency_matrix[mat_idx1,mat_idx0])
    if (len(neighb_idx)==2):
        U_current += np.log10(adjacency_matrix[mat_idx0,len(labelNumerics)])		#neighbor=End_point
    if np.isnan(U_current) or np.isinf(U_current):
        U_current = -1.0/machineEpsilon()
           
    U_new += local_potentials[current_e][local_labels.index(new_l)][1]
    mat_idx0 = labelNumerics.index(new_l)
    for j in neighb_idx:
        l_neighb = estimated_target_vector[j]
        if l_neighb>0:
            mat_idx1 = labelNumerics.index(l_neighb)
            U_new += np.log10(adjacency_matrix[mat_idx0,mat_idx1])
            U_new += np.log10(adjacency_matrix[mat_idx1,mat_idx0])
    if (len(neighb_idx)==2):
        U_new += np.log10(adjacency_matrix[mat_idx0,len(labelNumerics)])		#neighbor=End_point
    if np.isnan(U_new) or np.isinf(U_new):
        U_new = -1.0/machineEpsilon()
                
    delta_E = U_current-U_new		#these U_new and U_current are actually negative of energies so the negative of differences will be delta_U	
    return delta_E

def calc_STUN_E(estimated_target_vector,local_potentials,adjacency_matrix,edge_neighboring_indx,labelNumerics,Umin,tunnelparam):
    U = calc_globalE(estimated_target_vector,local_potentials,adjacency_matrix,edge_neighboring_indx,labelNumerics)
    STUN = 1 - np.exp(-tunnelparam*(U-Umin))
    return STUN,U 

def calc_STUN_E_quick(estimated_target_vector,local_potentials,adjacency_matrix,edge_neighboring_indx,labelNumerics,Umin,tunnelparam,U_current,current_e,current_l,new_l):
    if current_l==new_l:
        STUN = 1-np.exp(-tunnelparam*(U_current-Umin))
        return STUN, U_current
    # U = Ucurrent - Ucurrent_l + Unew_l = Ucurrent + deltaU
    U = U_current + calc_deltaE(current_e,current_l,new_l,estimated_target_vector,local_potentials,adjacency_matrix,edge_neighboring_indx,labelNumerics)
    STUN = 1 - np.exp(-tunnelparam*(U-Umin))
    return STUN,U 

#def calc_deltaE_STUN():    

def calc_E_1edge (current_e,current_l,estimated_target_vector,local_features_potentials,adjacency_matrix,edge_neighboring_indx,labelNumerics):
    if current_l<=0 :
        print "\nERROR! label <=0, current_l:",current_l,"\nAborted!"
        exit(0)
            
    neg_U_current=0
    local_labels = [local_features_potentials[current_e][j][0] for j in range(len(local_features_potentials[current_e]))]
    neighb_idx = edge_neighboring_indx[current_e]
    neg_U_current += local_features_potentials[current_e][local_labels.index(current_l)][1]
    mat_idx0 = labelNumerics.index(current_l)
    for j in neighb_idx:
        l_neighb = estimated_target_vector[j]
        if l_neighb>0:
            mat_idx1 = labelNumerics.index(l_neighb)
            tmp = adjacency_matrix[mat_idx0,mat_idx1]
            if tmp<=0 or numpy.isinf(tmp) or  numpy.isnan(tmp):
                neg_U_current += numpy.log(machineEpsilon(float)/10e10)          #log = -61    #10e306
            else:   
                neg_U_current += numpy.log(tmp)
            tmp = adjacency_matrix[mat_idx1,mat_idx0]
            if tmp<=0 or numpy.isinf(tmp) or  numpy.isnan(tmp):
                neg_U_current += numpy.log(machineEpsilon(float)/10e10)          #log = -61    #10e306
            else:   
                neg_U_current += numpy.log(tmp)

    if (len(neighb_idx)==2):
    #for k in range(4-len(neighb_idx)):
        tmp = adjacency_matrix[mat_idx0,len(labelNumerics)]		#neighbor=End_point
        if tmp<=0 or numpy.isinf(tmp) or  numpy.isnan(tmp):
            neg_U_current += numpy.log(machineEpsilon(float)/10e10)          #log = -61    #10e306
        else:   
            neg_U_current += numpy.log(tmp)
        
    U_current = -neg_U_current
    return U_current
       
def calc_deltaE_STUN (current_e,current_l,new_l,E_min, E_current,tunparam, estimated_target_vector,local_features_potentials,adjacency_matrix,edge_neighboring_indx,labelNumerics, normalizer = 1.0):
    if current_l<=0 or new_l<=0:
        print "\nERROR! label <=0, current_l:",current_l, " new_l:",new_l,"\nAborted!"
        exit(0)
            
    This_e_U_current = calc_E_1edge (current_e,current_l,estimated_target_vector,local_features_potentials,adjacency_matrix,edge_neighboring_indx,labelNumerics)
    This_e_U_new = calc_E_1edge (current_e,new_l,estimated_target_vector,local_features_potentials,adjacency_matrix,edge_neighboring_indx,labelNumerics)
    
    neighb_idx = edge_neighboring_indx[current_e]
    tmp_estimated_target_vector = estimated_target_vector
    tmp_estimated_target_vector[current_e] = new_l

    for j in neighb_idx:
            This_e_U_current += calc_E_1edge (j,estimated_target_vector[j],estimated_target_vector,local_features_potentials,adjacency_matrix,edge_neighboring_indx,labelNumerics)
            This_e_U_new += calc_E_1edge (j,estimated_target_vector[j],tmp_estimated_target_vector,local_features_potentials,adjacency_matrix,edge_neighboring_indx,labelNumerics)

    E_exclude_this_e = E_current - This_e_U_current
    delta_E = calc_E_stun (E_exclude_this_e+This_e_U_new, E_min, tunparam,normalizer) - calc_E_stun (E_current, E_min, tunparam,normalizer)
    return delta_E
    
    
def draw_from_dist(pdf,num):
    ##draw num samples from a distribution with pdf
    samples = []
    N = len(pdf)
    cnt=0
    while (cnt<num):
        x = int(round(random.uniform(0, N-1)))		#random.uniform(a, b) : random floating point number N such that a <= N <= b
        f = pdf[x]
        if f> random.uniform(0, 1):
            samples.append(x)
            cnt+=1		
    return samples
    
    
def write_output_result(g,edge_w_indx,local_potentials,estimated_target_vector,target_vector, output_name, estimated_target_vector_prob = 0):
    history = '>>> %s: auto labelling' % (time.ctime(time.time()))	
    if attributes.has_key("history"):
        history = attributes["history"] + "\n" + history
        del attributes['history']
        
    h=copy.deepcopy(g)
    for i in range(len(g.edge_list())):
        e=edge_w_indx[i]
        h.set_edge_property(e,'estimated_label_potentials',[list(j) for j in local_potentials[i]])	#local_potentials is a list for all edges, for each edge = local_potentials[i] is list of tuples=[(label,highest_prob),(l,prob),...,(l,lowest_prob)]
        if estimated_target_vector_prob == 0:
            h.set_edge_property(e,'estimated_label',[[estimated_target_vector[i],1.0]])	#local_potentials is a list for all edges, for each edge = local_potentials[i] is list of tuples=[(label,highest_prob),(l,prob),...,(l,lowest_prob)]
        else:
            h.set_edge_property(e,'estimated_label',[[estimated_target_vector[i],estimated_target_vector_prob[i]]])
        if (target_vector[i]>0 and (target_vector[i] == estimated_target_vector[i])):
            h.set_edge_property(e,'error_label',1)	#correct labelling	
        elif (target_vector[i]>0 and (target_vector[i] != estimated_target_vector[i])):
            h.set_edge_property(e,'error_label',2)	#error labelling	
        else:
            h.set_edge_property(e,'error_label',0)
    graph_analysis.output_graph(output_name, h, history, attributes)   
    
    print ("Successfully wrote the %s\n" %output_name)
    
    cmd=("graph2cylinder.py %s %s --use_estimated_label --clobber " %(output_name,output_name[:-3]+"_cyl.db"))	
    os.system(cmd)

    cmd=("graph2cylinder.py %s %s --use_error_label --clobber " %(output_name,output_name[:-3]+"_error_cyl.db"))	
    os.system(cmd)
    

def write_posterior_prob_map(g,edge_w_indx, posterior_prob,target_vector, posterior_probas, output_name):
    history = '\n>>> %s: posteriror probability map' % (time.ctime(time.time()))	
    if attributes.has_key("history"):
        history = attributes["history"] + "\n" + history
        del attributes['history']

    h=copy.deepcopy(g)
    print len(g.edge_list())
    for i in range(len(g.edge_list())):
        e=edge_w_indx[i]
        print i, posterior_prob[i]
        if math.isnan(posterior_prob[i][1]) or math.isinf( posterior_prob[i][1]) :
            h.set_edge_property(e,'estimated_label',[[-100,0]])
        elif math.isnan(np.log10(posterior_prob[i][1])) or math.isinf( np.log10(posterior_prob[i][1])):
            h.set_edge_property(e,'estimated_label',[[-100,0]])
        else:    
            h.set_edge_property(e,'estimated_label',[[posterior_prob[i][0],posterior_prob[i][1]]])	#255 is for mapping the color to [0..255] for viewing in brain-view2
        posterior_probas_sorted = sorted(posterior_probas[i],key=lambda l:l[1], reverse=True)
        h.set_edge_property(e,'estimated_label_potentials',posterior_probas_sorted)	#local_potentials is a list for all edges, for each edge = local_potentials[i] is list of tuples=[(label,highest_prob),(l,prob),...,(l,lowest_prob)]
    graph_analysis.output_graph(output_name, h, history, attributes)   
    print ("Succefully wrote the %s\n" %output_name)

    cmd=("graph2cylinder.py %s %s --use_estimated_label --clobber " %(output_name,output_name[:-3]+"_cyl.db"))	
    #print(cmd)
    os.system(cmd)

def write_prior_prob_map(g,edge_w_indx, prior_prob,target_vector,labelNumerics, output_name):
    history = '\n>>> %s: prior probability map' % (time.ctime(time.time()))	
    if attributes.has_key("history"):
        history = attributes["history"] + "\n" + history
        del attributes['history']

    h=copy.deepcopy(g)
    for i in range(len(g.edge_list())):
        e=edge_w_indx[i]
        l = target_vector[i]
        if l in labelNumerics:
            j = labelNumerics.index(l)
            if math.isnan(prior_prob[j,0]) or math.isinf( prior_prob[j,0]):
                h.set_edge_property(e,'estimated_label',[[0,0]])
            else:    
                h.set_edge_property(e,'estimated_label',[[np.log10(prior_prob[j,0]),prior_prob[j,0]]])	#255 is for mapping the color to [0..255] for viewing in brain-view2
        else:
            h.set_edge_property(e,'estimated_label',[[0,0]])
    graph_analysis.output_graph(output_name, h, history, attributes)   
    print ("Successfully wrote the %s\n" %output_name)

    cmd=("graph2cylinder.py %s %s --use_estimated_label --clobber " %(output_name,output_name[:-3]+"_cyl.db"))	
    os.system(cmd)
   
    
def calc_local_likelihood(feature_vector,labelNumerics,class_prior, mean_val, std_val):
    posteriorprobs = []		#each element = local_potentials[i] = potential of all possible labels for edge[i] sorted from high probability to low = {label:probability}
    for i in range(feature_vector.shape[0]):
        label_prob = {} #p(feature1/l)*p(feature2/l)*..*p(l) save as edge_property  estimated_label in a list for all labels
        for k in range(len(labelNumerics)):
            l = labelNumerics[k]
            p = class_prior[k,0]	# p(l)
            for j in range(feature_vector.shape[1]):
                tmp = univaraite_gaussian(feature_vector[i,j],mean_val[k,j],std_val[k,j]) #this should be >=0
                #if np.isinf(tmp) or  np.isnan(tmp):
                    #print ("ERROR!log is nan/inf!\n In calculation of log univaraite_gaussian for feature %s with true label %d for label %d: feature_vector[i,j] = %f,mean_val[k,j] = %f,std_val[k,j] = %f, gaussian= %f .\n" %(featureNames[j],target_vector[i,0],l,feature_vector[i,j],mean_val[k,j],std_val[k,j],univaraite_gaussian(feature_vector[i,j],mean_val[k,j],std_val[k,j])))
                    #p+= np.log(machineEpsilon(float)/10e10)	#10e306
                if tmp < 0:
                    tmp = 0
                    print ("ERROR! univaraite_gaussian for feature %s with true label %d for label %d: feature_vector[i,j] = %f,mean_val[k,j] = %f,std_val[k,j] = %f, gaussian= %f .\n" %(featureNames[j],target_vector[i],l,feature_vector[i,j],mean_val[k,j],std_val[k,j],univaraite_gaussian(feature_vector[i,j],mean_val[k,j],std_val[k,j])))
                #elif tmp >= 1:
                    #tmp = 1
                    #print ("ERROR! univaraite_gaussian for feature %s with true label %d for label %d: feature_vector[i,j] = %f,mean_val[k,j] = %f,std_val[k,j] = %f, gaussian= %f .\n" %(featureNames[j],target_vector[i],l,feature_vector[i,j],mean_val[k,j],std_val[k,j],univaraite_gaussian(feature_vector[i,j],mean_val[k,j],std_val[k,j])))
                p = p*tmp
            label_prob[l] = p
            if math.isnan(p) or math.isinf(p):
                label_prob[l] = 0 #machineEpsilon(float)  #-1/(machineEpsilon(float))#np.log(machineEpsilon(float)/10e306)#
        label_prob_sorted = sorted(label_prob.iteritems(), key=lambda (k,v): (v,k), reverse=True)		#sort from high to low probability, put in tuple [(key,value),(key,value),...]		
        #label_prob_sorted = sorted(label_prob.items(), key=itemgetter(1), reverse=True)		#sort from high to low probability and return a list instead of original dict => sorted_list=[(label1,prob1),(label2,prob2),...]		
        label_prob_sorted_list = [list(l) for l in label_prob_sorted] 
        posteriorprobs.append(label_prob_sorted_list)
    return posteriorprobs
    
def calc_log_local_likelihood(feature_vector,labelNumerics,class_prior, mean_val, std_val):
    posteriorprobs = []		#each element = local_potentials[i] = potential of all possible labels for edge[i] sorted from high probability to low = {label:probability}
    for i in range(feature_vector.shape[0]):
        label_prob = {} #p(feature1/l)*p(feature2/l)*..*p(l) save as edge_property  estimated_label in a list for all labels
        for k in range(len(labelNumerics)):
            l = labelNumerics[k]
            p = np.log10(class_prior[k,0])	# p(l)
            for j in range(feature_vector.shape[1]):
                tmp = univaraite_gaussian(feature_vector[i,j],mean_val[k,j],std_val[k,j]) #this should be >=0
                #if np.isinf(tmp) or  np.isnan(tmp):
                    #print ("ERROR!log is nan/inf!\n In calculation of log univaraite_gaussian for feature %s with true label %d for label %d: feature_vector[i,j] = %f,mean_val[k,j] = %f,std_val[k,j] = %f, gaussian= %f .\n" %(featureNames[j],target_vector[i,0],l,feature_vector[i,j],mean_val[k,j],std_val[k,j],univaraite_gaussian(feature_vector[i,j],mean_val[k,j],std_val[k,j])))
                    #p+= np.log(machineEpsilon(float)/10e10)	#10e306
                if tmp < 0:
                    tmp = 0 
                    print ("ERROR! univaraite_gaussian for feature %s with true label %d for label %d: feature_vector[i,j] = %f,mean_val[k,j] = %f,std_val[k,j] = %f, gaussian= %f .\n" %(featureNames[j],target_vector[i],l,feature_vector[i,j],mean_val[k,j],std_val[k,j],univaraite_gaussian(feature_vector[i,j],mean_val[k,j],std_val[k,j])))
                #elif tmp >= 1:
                    #tmp = 1
                    #print ("ERROR! univaraite_gaussian for feature %s with true label %d for label %d: feature_vector[i,j] = %f,mean_val[k,j] = %f,std_val[k,j] = %f, gaussian= %f .\n" %(featureNames[j],target_vector[i],l,feature_vector[i,j],mean_val[k,j],std_val[k,j],univaraite_gaussian(feature_vector[i,j],mean_val[k,j],std_val[k,j])))
                #if not np.isinf(np.log10(tmp)): # so log10 is not -inf
                p += np.log10(tmp)
                #elif tmp>1:
                    #p += 1000
                #elif tmp < 1:
                    #p += -1000
                #else:
                    #print ("ERROR! univaraite_gaussian for feature %s with true label %d for label %d: feature_vector[i,j] = %f,mean_val[k,j] = %f,std_val[k,j] = %f, gaussian= %f .\n" %(featureNames[j],target_vector[i],l,feature_vector[i,j],mean_val[k,j],std_val[k,j],univaraite_gaussian(feature_vector[i,j],mean_val[k,j],std_val[k,j])))
            label_prob[l] = p
            if math.isnan(p) or math.isinf(p):
                label_prob[l] = np.log10(machineEpsilon(float)) #machineEpsilon(float)  #-1/(machineEpsilon(float))#np.log(machineEpsilon(float)/10e306)#
        label_prob_sorted = sorted(label_prob.iteritems(), key=lambda (k,v): (v,k), reverse=True)		#sort from high to low probability, put in tuple [(key,value),(key,value),...]		
        #label_prob_sorted = sorted(label_prob.items(), key=itemgetter(1), reverse=True)		#sort from high to low probability and return a list instead of original dict => sorted_list=[(label1,prob1),(label2,prob2),...]		
        label_prob_sorted_list = [list(l) for l in label_prob_sorted] 
        posteriorprobs.append(label_prob_sorted_list)
    return posteriorprobs


def calc_likelihood(feature_vector,target_vector,labelNumerics,class_prior, mean_val, std_val, true_label=0):
    posteriorprobs = calc_local_likelihood(feature_vector,labelNumerics,class_prior, mean_val, std_val)
    local_potentials = copy.copy(posteriorprobs)
    for i in range(feature_vector.shape[0]):
        for j in range(len(posteriorprobs[i])):
            l_e =  posteriorprobs[i][j][0] # the label for this edge    (we calculate the posterior probability first for every possible label, given that neighbours have the true labels)       
            neighb_idx = edge_neighboring_indx[i]
            if l_e in labelNumerics:
                mat_idx0 = labelNumerics.index(l_e)
                for idx in neighb_idx:
                    if not true_label:
                        l_neighb = local_potentials[idx][0][0]   #current best label for neighbour edge
                    else:   
                        l_neighb = target_vector[idx]   #true label for neighbour edge (assume true label for the nieghbours)
                    if l_neighb>0 and l_neighb in labelNumerics:
                        mat_idx1 = labelNumerics.index(l_neighb)
                        posteriorprobs[i][j][1] *= adjacency_matrix[mat_idx0,mat_idx1]
                if (len(neighb_idx)==2):		# Y edge that has an end-point 
                    posteriorprobs[i][j][1] *= adjacency_matrix[mat_idx0,len(labelNumerics)]		#neighbor=End_point
    return posteriorprobs

def calc_log_likelihood(feature_vector,target_vector,labelNumerics,class_prior, mean_val, std_val):
    posteriorprobs = calc_log_local_likelihood(feature_vector,labelNumerics,class_prior, mean_val, std_val)
    local_potentials = copy.copy(posteriorprobs)
    for i in range(feature_vector.shape[0]):
        for j in range(len(posteriorprobs[i])):
            l_e =  posteriorprobs[i][j][0] # the label for this edge           
            neighb_idx = edge_neighboring_indx[i]
            if l_e in labelNumerics:
                mat_idx0 = labelNumerics.index(l_e)
                for idx in neighb_idx:
                    l_neighb = local_potentials[idx][0][0]   #current best label for neighbour edge
                    if l_neighb>0 and l_neighb in labelNumerics:
                        mat_idx1 = labelNumerics.index(l_neighb)
                        posteriorprobs[i][j][1] += np.log10(adjacency_matrix[mat_idx0,mat_idx1])
                if (len(neighb_idx)==2):		# Y edge that has an end-point 
                    posteriorprobs[i][j][1] += np.log10(adjacency_matrix[mat_idx0,len(labelNumerics)])		#neighbor=End_point
    return posteriorprobs

def calc_posteriorprob(feature_vector,target_vector,labelNumerics,class_prior, mean_val, std_val, true_label=0):
    posteriorprobs = calc_likelihood(feature_vector,target_vector,labelNumerics,class_prior, mean_val, std_val,true_label)
    posteriorprobs_truelabel = [[target_vector[i],0] for i in range(feature_vector.shape[0])]
    for i in range(feature_vector.shape[0]):
        denominator = 0
        for j in range(len(posteriorprobs[i])):
            denominator += posteriorprobs[i][j][1]
        for j in range(len(posteriorprobs[i])):
            posteriorprobs[i][j][1] = posteriorprobs[i][j][1]/denominator  #(denominator+machineEpsilon(float))           
    for i in range(feature_vector.shape[0]):
        for j in range(len(posteriorprobs[i])):
            l_e =  posteriorprobs[i][j][0] # the label for this edge
            if l_e == target_vector[i]:    #if the label is the true label for this edge
                posteriorprobs_truelabel[i]=[posteriorprobs[i][j][0], posteriorprobs[i][j][1]]        
    return posteriorprobs, posteriorprobs_truelabel

  
def current_state (temperatures, global_U,global_err, global_vol_err,minErr, minErrvol, minErr_it, minErr_l, minErr_energy, minErr_T, minErr_optimizer,minvol_Err, minvolErr_vol, minvolErr_it, minvolErr_l, minvolErr_energy, minvolErr_T, minvolErr_optimizer,minE_optimizer,minE , minE_it ,minE_l ,minE_err ,minE_T, target_vector, estimated_target_vector,labeling_potentials,global_cnt,g,edge_w_indx,T, tunparam,optimizer="GD",normalizer = 1.0):
    #### optimization step 1. calculate current state
    temperatures.append(T)
    U_current = calc_globalE_weightlabel(estimated_target_vector,local_features_potentials,adjacency_matrix,edge_neighboring_indx,labelNumerics)	#with current labelling of all edges that come from labeling_potentials[:][0][0]
    if optimizer == "STUN":
        U_current_stun = calc_E_stun (U_current, minE, tunparam,normalizer)
        print ("F_STUN = %e U_current = %e, minE = %e" %(U_current_stun,U_current,minE))  
    if numpy.isnan(U_current) or numpy.isinf(U_current):
        print "Error! U_current undefined!\nAborted"
        exit(0)
    global_U.append(U_current)
    [error_num,labeled_num,total_volume,error_volume]= error_calculation(target_vector, estimated_target_vector,labeling_potentials,global_cnt,g,edge_w_indx,U_current,T)
    if (labeled_num > 0):	##otherwise test set was not labeled and we can't calculate error of labeling
        global_err.append(100*float(error_num)/float(labeled_num+machineEpsilon(float)))
        global_vol_err.append(100*float(error_volume)/float(total_volume+machineEpsilon(float)))			
        if (100*float(error_num)/float(labeled_num+machineEpsilon(float))) < minErr:		## saving best error configuration
            minErr_it = global_cnt
            minErr = 100*float(error_num)/float(labeled_num+machineEpsilon(float))
            minErrvol = 100*float(error_volume)/float(total_volume+machineEpsilon(float))
            minErr_l = estimated_target_vector
            minErr_energy = U_current
            minErr_T = T	
            minErr_optimizer = optimizer
        if (100*float(error_volume)/float(total_volume+machineEpsilon(float))) < minvolErr_vol:		## saving best error configuration
            minvolErr_it = global_cnt
            minvol_Err = 100*float(error_num)/float(labeled_num+machineEpsilon(float))
            minvolErr_vol = 100*float(error_volume)/float(total_volume+machineEpsilon(float))
            minvolErr_l = estimated_target_vector
            minvolErr_energy = U_current
            minvolErr_T = T	
            minvolErr_optimizer = optimizer
    if U_current<minE:		## saving best energy configuration	
        minE = U_current
        minE_it = global_cnt
        minE_l = estimated_target_vector
        minE_err = 100*float(error_num)/float(labeled_num+machineEpsilon(float))
        minE_err = 100*float(error_volume)/float(total_volume+machineEpsilon(float))
        minE_T = T
        minE_optimizer = optimizer
    return temperatures, global_U,global_err, global_vol_err,minErr, minErrvol, minErr_it, minErr_l, minErr_energy, minErr_T,minErr_optimizer,minvol_Err, minvolErr_vol, minvolErr_it, minvolErr_l, minvolErr_energy, minvolErr_T, minvolErr_optimizer,minE_optimizer, minE , minE_it ,minE_l ,minE_err ,minE_T, U_current

def calc_transition_pdf (target_vector, estimated_target_vector, local_features_potentials, adjacency_matrix, edge_neighboring_indx, labelNumerics, minE, tunparam, optimizer = "GD",normalizer = 1.0):
    #### optimization step 2. calculate deltaE and pdf of transition to new labels for each node	
    labeling_potentials = [0 for i in range(len(target_vector))]		#each element = labeling_potentials[i] = potential of all possible labels for edge[i] sorted from high probability to low = {label:probability}
    pos_delE_cnt=0 #for SA

    #print "optimizer is ", optimizer,
    
    for i in edgelist_visit:					#for each node: X = feature_vector[i,:], current_l = labeling_potentials[i][0][0]
        l_current = estimated_target_vector[i]					
        delta_U = []							#energy variation of transition to other labels for current edge[i]
        total_ltransition = 0
        this_edge_label_prob = {} 				#p(l/feature1)*p(l/feature2)*..*p(l) save as edge_property  estimated_label in a list for all labels
        for k in range(len(labelNumerics)):		#calculate energy variation for transition to l_new
            l_new = labelNumerics[k]
            if optimizer == "GD":
                current_delta_U = calc_deltaE(i, l_current, l_new,estimated_target_vector, local_features_potentials, adjacency_matrix,edge_neighboring_indx,labelNumerics)
                this_edge_label_prob[l_new] = current_delta_U
            elif optimizer == "SA":
                current_delta_U = calc_deltaE(i, l_current, l_new,estimated_target_vector, local_features_potentials, adjacency_matrix,edge_neighboring_indx,labelNumerics)
            elif optimizer == "STUN":
                current_delta_U = calc_deltaE_STUN (i,l_current,l_new,minE, U_current,tunparam, estimated_target_vector,local_features_potentials,adjacency_matrix,edge_neighboring_indx,labelNumerics,normalizer)
            else:
                current_delta_U = calc_deltaE_generic (i, l_current, l_new, minE,estimated_target_vector, local_features_potentials, adjacency_matrix,edge_neighboring_indx,labelNumerics, tunparam, optimizer, normalizer)
            
            delta_U.append(current_delta_U)
        if optimizer == "GD":
            this_edge_label_prob_sorted = sorted(this_edge_label_prob.iteritems(), key=lambda (k,v): (v,k), reverse=False)		#sort from low to high probability and return a list instead of original dict => sorted_list=[(label1,prob1),(label2,prob2),...]							 
            labeling_potentials[i] = this_edge_label_prob_sorted  #### we can't just append because we do not iterate edges in order, the order is shuffled
            if not estimated_target_vector[i] == labelNumerics[delta_U.index(min(delta_U))]:
                estimated_target_vector[i] = labelNumerics[delta_U.index(min(delta_U))]
        elif not optimizer == "GD":
            pdf_ltransition=[]	#probability density function of transition to other labels for current edge[i]
            normal_delta_U = [(x /(max(abs(max(delta_U)),abs(min(delta_U)))+machineEpsilon(float))) for x in delta_U]	## normalize delta_U to be in range [-1...1]: (so map delta_U =[0...max] -> [0 ..1] and negative delta_U will be interpolated linearly	
            for k in range(len(labelNumerics)):	#calculate energy variation for transition to l_new
                if exp(-normal_delta_U[k]/(T+machineEpsilon(float)))==inf:
                    exp_val= 1/machineEpsilon(float)
                else:
                    exp_val=exp(-normal_delta_U[k]/(T+machineEpsilon(float)))
                total_ltransition = total_ltransition + exp_val
                pdf_ltransition.append(exp_val)
                this_edge_label_prob[labelNumerics[k]] = exp_val
            pdf_ltransition = [(x/(total_ltransition+machineEpsilon(float))) for x in pdf_ltransition]		##labeling_potentials contains potential of all labels for each edge based on local features!
            ltransition_indx = draw_from_dist(pdf_ltransition,1)[0]		## actual transition by drawing from pdf_ltransition
            estimated_target_vector[i] = labelNumerics[ltransition_indx]
            this_edge_label_prob_sorted = sorted(this_edge_label_prob.iteritems(), key=lambda (k,v): (v,k), reverse=True)		#sort from high to low probability and return a list instead of original dict => sorted_list=[(label1,prob1),(label2,prob2),...]							 
            labeling_potentials[i] = [(sr[0],sr[1]/(total_ltransition+machineEpsilon(float))) for sr in this_edge_label_prob_sorted]
            if normal_delta_U[ltransition_indx]>0:
                pos_delE_cnt+=1
            percent_pos_delE.append(float(pos_delE_cnt)/float(len(labelNumerics)))
            
    return estimated_target_vector, labeling_potentials, percent_pos_delE

def optimizer_T_update (global_cnt,global_err, T, eta,estimated_target_vector, edge_neighboring_indx, global_U , global_U_threshold, max_global_iteration,feature_vector, anneal_change, nonlin_annealing, minE, tunparam, optimizer = "GD", normalizer = 1.0, fthresh = 0.03):
    #### 3. update T and check number of iterations
    global_cnt = global_cnt + 1
    while_param = True
    if optimizer == "GD":
        if (global_cnt > max_global_iteration):
            while_param = False
        elif (global_cnt > 5):	#run at least 5 iterations	
            tmp_cnt = 0
            for converge_i in range(1,6):#if for 5 consecutive iterations the global energy has not changed => energy converged
                tmp_converge = abs(global_U[global_cnt-converge_i-1]-global_U[global_cnt-converge_i])
                if (tmp_converge < global_U_threshold):
                    tmp_cnt = tmp_cnt + 1
            if (tmp_cnt > 4):
                while_param = False
                
    elif optimizer == "STUN":
        if (global_cnt > max_global_iteration):
            while_param = False
        elif (global_cnt > 250):	#run at least 100 iterations	
            tmp_cnt = 0
            for converge_i in range(1,19):#if for 5 consecutive iterations the global energy has not changed => energy converged
                tmp_converge = abs(global_U[global_cnt-converge_i-1]-global_U[global_cnt-converge_i])
                if (tmp_converge < global_U_threshold):
                    tmp_cnt = tmp_cnt + 1
            if (tmp_cnt > 10):
                while_param = False
            tmp_cnt = 0
            for converge_i in range(1,19):#if for 19 last iterations on average f_stun > fthresh => T increase
                if abs(calc_E_stun (global_U[global_cnt-converge_i-1], minE, tunparam,normalizer)-calc_E_stun (global_U[global_cnt-converge_i], minE, tunparam,normalizer))< fthresh:
                    tmp_cnt = tmp_cnt + 1
            if (tmp_cnt > 9 and T < 1e-10):
                T = T/eta		#on average f_stun > fthresh => T increase (so it gets out of the barrier by tunneling)
        else:
            T = eta*T		#on average f_stun < fthresh => T decrease (so it searches)
    else:
        if (global_cnt > 250):	#run at least 250 iterations	
            tmp_cnt = 0
            for converge_i in range(1,19):#if for 5 consecutive iterations the global energy has not changed => energy converged
                tmp_converge = abs(global_U[global_cnt-converge_i-1]-global_U[global_cnt-converge_i])
                if (tmp_converge < global_U_threshold):
                    tmp_cnt = tmp_cnt + 1
            if (tmp_cnt > 10):
                while_param = False
        T = eta*T
        if nonlin_annealing:
            eta = 1.0- exp(-(global_cnt+2500.0)/800.0)		# eta = 1-exp(-(it+2500)/800)
        if (global_cnt > max_global_iteration):
            while_param = False
        elif (global_cnt > 750) and std(global_err[global_cnt-100:global_cnt-1])<2.5:		## if for the last 100 iterations globalerr has not changed much
            while_param = False	
        elif (global_cnt > 20):	#run at least 20 iterations	
            tmp_cnt = 0
            for converge_i in range(1,11):#if for 10 consecutive iterations the global energy has not changed => energy converged
                tmp_converge = abs(global_U[global_cnt-converge_i-1]-global_U[global_cnt-converge_i])
                if (tmp_converge < global_U_threshold):
                    tmp_cnt = tmp_cnt + 1
            if (tmp_cnt > 9):
                while_param = False
        if not nonlin_annealing:
            if global_cnt>10 and global_err<50 and not anneal_change:
                anneal_change=1
                eta = 0.9+eta/10.0		#eta+(1-eta)*0.9		if eta<1 => eta/10<0.1 => 0.9+eta/10.0<1
        if global_cnt%10==0:
            #### special iteration: switch each edge's label with the current label of a random adjacent edge 
            special_it_vector = []
            for i in range(feature_vector.shape[0]):
                n_e_idx_list = edge_neighboring_indx[i]
                if len(n_e_idx_list)>0 :
                    n_e_idx = random.sample(n_e_idx_list,1)[0]
                else:
                    n_e_idx = i
                special_it_vector.append(estimated_target_vector[n_e_idx])
                
            for i in edgelist_visit:
                estimated_target_vector[i] = special_it_vector [i]
    
    return	while_param, global_cnt, T, eta, anneal_change, estimated_target_vector		

############################################################################################################################################0
program_name = 'bayesLearner.py'


if __name__ == '__main__':

    usage = "Usage: "+program_name+" [options] test_graph.db [output_filename(test_labeled).db]\n"+\
        "   or  "+program_name+" --help";

    parser = OptionParser(usage)

    parser.add_option("--trainer", type="string", dest="trainer",
            help="The training input with prior, adjacency matrix and feature distributions calculated.")
    
    parser.add_option("--training_inputs", type="string", dest="training_inputs",
            default=0, help="The comma-seperated list of files to do the training and to calculate prior, adjacency matrix and feature distributions.\n[training_graph1.db,training_graph2.db,...,training_graphN.db]")
    
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

    #parser.add_option("--angle", action="store_true", dest="angle",
    #        default=0, help="reads features angleX, angleY and angleZ with +X, +Y, +Z axis of the Euclidean coordinate space.")

    #parser.add_option("--axial_direction", action="store_true", dest="axial_direction",
    #                   default=0, help="reads features direction cosine in +Z hemisphere and calculate Bingham (Kent) distribution on the unit sphere.")

    parser.add_option("--mri_labels", action="store_true", dest="mri_labels",
            default=0, help="reads features vessel distance from each 40 mri labels.")

    parser.add_option("--rel_diameter", action="store_true", dest="rel_diameter",
            default=0, help="reads features vessel diameter relative to its adjacent vessels.")

    parser.add_option("--rel_dir", action="store_true", dest="rel_dir",
            default=0, help="reads features vessel direction relative to its adjacent vessels.")

    #parser.add_option("--classifier_type", action="store", type="string", dest="classifier_type",
    #                   default="NB", help="specify the classifier type:			1.NB=iid Gaussian NB (default)			2.NB2=univariate Gaussian and Gamma NB			3.GDA=multivariate Gaussian")

    parser.add_option("--PCA", action="store_true", dest="pca",
            default=0, help="run PCA on feature vectors to make it orthogonal")

    #parser.add_option("--plot", action="store_true", dest="plot",
    #                   default=0, help="plot feature distributions")

    #parser.add_option("--save_plot", action="store_true", dest="save_plot",
    #                   default=0, help="save plot of feature distributions")

    #parser.add_option("--laplace_smooth", action="store_true", dest="laplace_smooth",
    #        default=0, help="add noisy data to feature vector for better feature distributions inference")

    parser.add_option("--PCA_CDF_thresh", type="float", dest="PCA_CDF_thresh", default=100.0, help="set threshold % for CDF of PCA eigenvalues, in order to perform feature selection(Default 100).")

    parser.add_option("--training_output", type="string", dest="training_output", help="training_output(distribution_parameters).db to save the output of training stage")

    parser.add_option("--gradient_iteration_num",type="int", dest="gradient_iteration_num", default=0, help="the number of gradient iterations before performing Gibbs sampling to refine initial labels [ Use >5]")                  

    parser.add_option("--iteration_num",type="int", dest="iteration_num", default=0, help="the number of global iterations of Gibbs sampling to refine initial labels [Use >1000]")                  

    parser.add_option("--temprature",type="float", dest="temprature", default=1.8e-15, help="the initial temperature for the relaxation. T=1.8e-15=1e-12*(.9**60) is in transition (around iteration 60 when T0=1e-12 and annealing .9)")   

    parser.add_option("--annealing", type="float", dest="anneal_param", default= 0.99, help = "temperature is multiplied by annealing parameter at the end of each global iteration (Default is 0.99)")

    parser.add_option("--nonlin_annealing", action="store_true", dest="nonlin_annealing", default=0, help="nonlinear annealing of 1-exp(-x) format")

    parser.add_option("--STUN_iteration_num",type="int", dest="STUN_iteration_num",
                    default=2000, help="the number of global iterations of stochastic tunneling to refine initial labels")  
    
    parser.add_option("--tunneling", type="float", dest="tunneling_param",
                    default= 0.05, help = "tunneling parameter which is the cut off for steepness of barriers: try 0.05 or 0.1")
    '''
    parser.add_option("--tunneling_threshold", type="float", dest="tunneling_threshold",
                    default= 0.03, help = "fthreshold on the E-Emin in stochastic tunneling")
    '''
    parser.add_option("--plot", action="store_true", dest="plot",default=0, help="plot the output labelled graph and error graph")

    parser.add_option("--verbose", action="store_true", dest="verbose", default=0, help="spit out detailed execution report to shell")	

    parser.add_option("--animation", action="store_true", dest="animation", default=0, help="Option to save intermediate h5s of labelling results for the purpose of labelling progression movie.")  

    parser.add_option("--posteriorprob", action="store_true", dest="posteriorprob", default=0, help="Option to save the posterior probability map for true labels.")

    parser.add_option("--decimate", type="int", dest="decimate",
            metavar="decimation_factor", default = 20, 
            help="decimate the iteration steps for animation making (default 20)")


    options, args = parser.parse_args()


    if len(args)==2:
        test_file, output_file = args
    elif len(args)==1:
        test_file = args[0]
        output_file = (test_file[:-3]+ "_autolabel.h5")
    else:    
        parser.error("incorrect number of arguments")


    print '>>> %s [%s]: %s\n' % (time.ctime(time.time()), os.getcwd(), string.join(argv))
    sys.stdout.flush()

    if not options.trainer and not options.training_inputs:
        raise SystemExit, \
            "Either --trainer or --training_inputs options is needed to get the training information."
            
    if options.trainer and options.training_inputs:
        raise SystemExit, \
            "Only one of --trainer or --training_inputs options should be used to get the training information."

    if not options.clobber and os.path.exists(output_file):
        raise SystemExit, \
            "The --clobber option is needed to overwrite an existing file."
            
        
    if options.trainer:
        #get the training data set
        f_handler = shelve.open(options.trainer,"r")
        classifier_type = f_handler['classifier_type']
        pca = f_handler['PCA']
        if pca=="Y":
            feat_select_indx = f_handler['feat_select_indx'] 		
            Y_bar = f_handler['Y_bar']
            normalizing_stds = f_handler['normalizing_stds']
            eigenvects = f_handler['eigenvects']
            eigenvals = f_handler['eigenvals'] 
        class_prior = f_handler['class_prior'] 
        labelNumerics = f_handler['labelNumerics']
        labelNumeric2Name = f_handler['labelNumeric2Name'] 
        featureNames = f_handler['featureNames']  
        adjacency_matrix = f_handler['adjacency_matrix']
        mean_val = f_handler['mean_val'] 
        std_val = f_handler['std_val'] 
        
        f_handler.close()
        print  ("Succesfully read the training file %s\n" %options.trainer)
        print "labelNumerics: ", labelNumerics
     
                
    if options.training_inputs:
        training_files = options.training_inputs.split(',')
        print "The list of training files are:"
        print training_files
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
        """    
        if options.angle:
            featureNames.append('angleX')
            featureNames.append('angleY')
            featureNames.append('angleZ')
        """
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
            labelNumeric2Name ={0:"No label",35:"Anterior Cerebral Artery", 191:"R. Middle Cerebral Artry", 190:"L. Middle Cerebral Artry", 2:"R. Intern Carotid Artery", 43:"L. Intern Carotid Artery", 200:"R. Posterior Comm. Artry", 9:"L. Posterior Comm. Artry", 8:"R. Posterior Cereb Artry", 5:"L. Posterior Cereb Artry", 68:"R. Superior Cereb Artery", 227:"L. Superior Cereb Artery", 46:"R. Ant. Inf. Cereb Artry", 12:"L. Ant. Inf. Cereb Artry", 196:"Basilar Artery", 7:"Vertebral Artery", 49:"R. Internal Audit Artery", 45:"L. Internal Audit Artery", 7: "Vertebral Artery" , 3: "R. Paraolivary Artry", 4: "L. Paraolivary Artry" , 11:"Superior Saggital Sinus", 6:"Great Cerbral Vein Galen", 30:"R. Transverse Sinus", 246:"L. Transverse Sinus", 192:"R. Caudal Rhinal Vein", 34:"L. Caudal Rhinal Vein", 20:"R. Rostral Rhinal Vein", 21:"L. Rostral Rhinal Vein", 101:"R. Sigmoid Sinus", 24:"L. Sigmoid Sinus", 58:"R. Longitud. Hippo. Vein", 57:"L. Longitud. Hippo. Vein", 56:"R. Thalamostriate Vein", \
            54:"L. Thalamostriate Vein", 1:"R. Medial Colicular Vein", 16:"L. Medial Colicular Vein", 170:"Unknown Sinus/Vein #01", 171:"R. Lateral collicular V.", 172:"L. Lateral collicular V.", 250:"L. Unknown Sinus/Vein #2", 251:"R. Unknown Sinus/Vein #2"}											
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
        labelNumerics, feature_vector, target_vector, edge_w_indx = graph2featureVect(training_files,featureNames)		#featureVect
        adjacency_matrix = make_adjmatrix(training_files,labelNumerics,edge_w_indx)
        """
        if options.laplace_smooth:
            for i in range(feature_vector.shape[0]):
                noisy_d = [random.gauss(0,10e-10) for j in range(feature_vector.shape[1])]
                feature_vector	= np.vstack([feature_vector, np.multiply(array(noisy_d),array(feature_vector[i]))+array(feature_vector[i])])
                target_vector.append(int(target_vector[i]))
                feature_vector = np.vstack([feature_vector, array(feature_vector[i])-np.multiply(array(noisy_d),array(feature_vector[i]))])
                target_vector.append(int(target_vector[i]))
        """        
    
        if options.pca:
            [feature_vector, Y_bar, normalizing_stds, eigenvects, eigenvals, feat_select_indx] = PCA(feature_vector, options.PCA_CDF_thresh )
    
        #### calculate class priors 
        class_prior = Matlib.empty((len(labelNumerics),1),float)
        mean_val = Matlib.empty((len(labelNumerics),feature_vector.shape[1]),float)
        std_val = Matlib.empty((len(labelNumerics),feature_vector.shape[1]),float)
        for i in range(len(labelNumerics)):
            l = labelNumerics[i]
            indx = list(pylab.find (np.array(target_vector)==l))
            class_prior[i,0]= float(len(indx))/float(len(target_vector))
            if class_prior[i,0] <=0:
                print "ERROR: label ", l, " # ", float(len(indx)), " / ", float(len(target_vector)), " prior " , class_prior[i,0]
                exit(0)
            for j in range(feature_vector.shape[1]):
                mean_val[i,j] = np.mean (feature_vector[indx,j])
                std_val[i,j] = np.std (feature_vector[indx,j], None,None, None, 1)		#np.std(a, axis=None, dtype=None, out=None, ddof=0, skipna=False, keepdims=False). The divisor used in calculations is N - ddof, where N represents the number of elements.
    
        if options.training_output:		
            f_handler = shelve.open(options.training_output)
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
            f_handler['featureNames'] = featureNames
            f_handler['feature_vector'] = feature_vector 
            f_handler['target_vector'] = target_vector 
            f_handler['adjacency_matrix'] = adjacency_matrix 
            f_handler['history'] = '>>> %s: %s' % (time.ctime(time.time()), string.join(argv))			
            f_handler.close()
            print ("Succesfully wrote %s" %options.training_output)
        #print "class_prior" , class_prior
    
        """
        #  We should check mean_val and std_val not to be Nan
        # We also need the class_prior, labelNumerics, featureNames and adjacency_matrix
        """
 
    ##get info from test set	
    try:
        g, attributes= graph_analysis.input_graph(test_file, ["history", "vertex_offsets"])  	# open graph data and copy contents to memory
    #print ("Succefully read in %s\n" %test_file)	
    except:
        print("Error reading in %s\n" %test_file)

    inp_vect=[]
    inp_vect.append(test_file)
    unimportantlabellist,feature_vector, target_vector, edge_w_indx = graph2featureVect(inp_vect,featureNames)		#featureVect
    edge_neighboring_indx = find_neighbors (g, edge_w_indx)
    

    if options.pca:
        [feature_vector, Y_bar, normalizing_stds, eigenvects, eigenvals,feat_select_indx] = PCA (feature_vector,options.PCA_CDF_thresh,Y_bar, normalizing_stds, eigenvects,eigenvals,feat_select_indx)


    print  ("Succesfully read %s\n" %test_file)

    """
    # We should check    feature_vector,  mean_val and std_val not to be Nan
    """
    ##################################################################
    #### 0. calculate the posterior probability for the true labels
    ##################################################################    
    write_prior_prob_map(g,edge_w_indx, class_prior,target_vector,labelNumerics, output_file[:-3]+"_priorprob.db")
    #local_potentials = calc_likelihood(feature_vector,target_vector,labelNumerics,class_prior, mean_val, std_val)
    #best_local_potentials = [local_potentials[i][0][0] for i in range(len(local_potentials))]		#for each edge: the label with highest probability  #[i] for ith edge,[0] the highest label probability, [0] label name
    #write_posterior_prob_map(g,edge_w_indx, best_local_potentials,output_file[:-3]+"_localposteriorprob.db")
    posterior_probas, posteriorprobs_truelabel = calc_posteriorprob(feature_vector,target_vector,labelNumerics,class_prior, mean_val, std_val, 1)
    write_posterior_prob_map(g,edge_w_indx, posteriorprobs_truelabel,target_vector, posterior_probas, output_file[:-3]+"_posteriorprob.db")
   
    ##################################################################
    #### 1. initialize the labels
    ##################################################################    
    RR_by_label = {}
    local_potentials = calc_log_local_likelihood(feature_vector,labelNumerics,class_prior, mean_val, std_val)
    """
    local_potentials = []		#each element = local_potentials[i] = potential of all possible labels for edge[i] sorted from high probability to low = {label:probability}
    for i in range(feature_vector.shape[0]):
        label_prob = {} #p(feature1/l)*p(feature2/l)*..*p(l) save as edge_property  estimated_label in a list for all labels
        for k in range(len(labelNumerics)):
            l = labelNumerics[k]
            p = np.log(class_prior[k,0])	# p(l)
            for j in range(feature_vector.shape[1]):
                tmp=univaraite_gaussian(feature_vector[i,j],mean_val[k,j],std_val[k,j])
                #if np.isinf(tmp) or  np.isnan(tmp):
                    #print ("ERROR!log is nan/inf!\n In calculation of log univaraite_gaussian for feature %s with true label %d for label %d: feature_vector[i,j] = %f,mean_val[k,j] = %f,std_val[k,j] = %f, gaussian= %f .\n" %(featureNames[j],target_vector[i,0],l,feature_vector[i,j],mean_val[k,j],std_val[k,j],univaraite_gaussian(feature_vector[i,j],mean_val[k,j],std_val[k,j])))
                    #p+= np.log(machineEpsilon(float)/10e10)	#10e306
                if tmp <=0:
                    tmp = 0
                    print ("ERROR! univaraite_gaussian for feature %s with true label %d for label %d: feature_vector[i,j] = %f,mean_val[k,j] = %f,std_val[k,j] = %f, gaussian= %f .\n" %(featureNames[j],target_vector[i,0],l,feature_vector[i,j],mean_val[k,j],std_val[k,j],univaraite_gaussian(feature_vector[i,j])))
                elif tmp > = 1:
                    tmp = 1
                    print ("ERROR! univaraite_gaussian for feature %s with true label %d for label %d: feature_vector[i,j] = %f,mean_val[k,j] = %f,std_val[k,j] = %f, gaussian= %f .\n" %(featureNames[j],target_vector[i,0],l,feature_vector[i,j],mean_val[k,j],std_val[k,j],univaraite_gaussian(feature_vector[i,j])))
                p+= np.log(tmp)
            label_prob[l] = p
            if math.isnan(p):
                label_prob[l] = -1/(machineEpsilon(float))#np.log(machineEpsilon(float)/10e306)#
        label_prob_sorted = sorted(label_prob.iteritems(), key=lambda (k,v): (v,k), reverse=True)		#sort from high to low probability, put in tuple [(key,value),(key,value),...]		
        #label_prob_sorted = sorted(label_prob.items(), key=itemgetter(1), reverse=True)		#sort from high to low probability and return a list instead of original dict => sorted_list=[(label1,prob1),(label2,prob2),...]		
        local_potentials.append(label_prob_sorted)
    """

    ##################################################################
    #### 1.1. initialization error calculation
    ##################################################################
    #print "Initial ERRORS:" 
    estimated_target_vector = [local_potentials[i][0][0] for i in range(len(local_potentials))]		#for each edge: the label with highest probability  #[i] for ith edge,[0] the highest label probability, [0] label name
    estimated_target_vector_prob = [local_potentials[i][0][1] for i in range(len(local_potentials))]     #for each edge: the label with highest probability  #[i] for ith edge,[0] the highest label probability, [0] label name
    [error_num,labeled_num,total_volume,error_volume] = error_calculation (target_vector, estimated_target_vector,0,g,edge_w_indx,0)
    RR_by_label = error_calculation_by_label(RR_by_label,target_vector, estimated_target_vector, local_potentials,feature_vector, 0,edge_neighboring_indx,g,edge_w_indx)
    if options.verbose:
        confusion_matrix_calculation(target_vector, estimated_target_vector,labelNumerics)

    ###################################################################
    ##### 1.2. write labeled output
    ###################################################################
    if options.plot:	
        write_output_result(g,edge_w_indx,local_potentials,estimated_target_vector,target_vector, output_file[:-3]+"_initial.db", estimated_target_vector_prob)
    #initial status 
    U_current = calc_globalE(estimated_target_vector,local_potentials,adjacency_matrix,edge_neighboring_indx,labelNumerics)	#with current labelling of all edges that come from local_potentials[:][0][0]
    if np.isnan(U_current) or np.isinf(U_current):
        print "Error! U_current undefined!\nAborted"
        exit(0)

    print "\nthe initial global energy = ", U_current
    init_E = U_current	
    RR_initial = 100 - 100*float(error_num)/float(labeled_num+machineEpsilon(float))
    RRvol_initial = 100 -100*float(error_volume)/float(total_volume+machineEpsilon(float))
    initial_num = int(labeled_num - error_num)
    initial_vol = total_volume - error_volume

    U_optimum = calc_globalE(target_vector,local_potentials,adjacency_matrix,edge_neighboring_indx,labelNumerics)	#with current labelling of all edges that come from local_potentials[:][0][0]
    if np.isnan(U_optimum) or np.isinf(U_optimum):
        print "Error! U_optimum undefined!\nAborted"
        exit(0)

    """
    ####################################################################
    ###### 1.3. calculated posterior probability of true label and write the output with the posterior probabilities as the labels
    ####################################################################
   if options.posteriorprob:
   #1. posterior probability without adj_matrix (local_potentials)
       posterior_prob = [0.0 for i in range(len(estimated_target_vector))]
       for i in range(len(estimated_target_vector)):
           l_e = target_vector[i]
           if l_e>0 and l_e in labelNumerics:
               local_labels = [local_potentials[i][j][0] for j in range(len(local_potentials[i]))]
           posterior_prob[i] = local_potentials[i][local_labels.index(l_e)][1]
       write_posterior_prob_map(g,edge_w_indx, posterior_prob, output_file[:-3]+"_localpotentials.db")

           #2. whole posterior probability with adj_matrix
       posterior_prob = [0.0 for i in range(len(estimated_target_vector))]
       for i in range(len(estimated_target_vector)):
           l_e = target_vector[i]
           if l_e>0 and l_e in labelNumerics:
               local_labels = [local_potentials[i][j][0] for j in range(len(local_potentials[i]))]
           posterior_prob[i] = local_potentials[i][local_labels.index(l_e)][1]
           neighb_idx = edge_neighboring_indx[i]
           mat_idx0 = labelNumerics.index(l_e)
           for j in neighb_idx:
               l_neighb = target_vector[j]
               if l_neighb>0 and l_neighb in labelNumerics:
                   mat_idx1 = labelNumerics.index(l_neighb)
               posterior_prob[i] += log(adjacency_matrix[mat_idx0,mat_idx1])
               if (len(neighb_idx)==2):		# Y edge that has an end-point 
                   posterior_prob[i] += log(adjacency_matrix[mat_idx0,len(labelNumerics)])		#neighbor=End_point
       write_posterior_prob_map(g,edge_w_indx, posterior_prob,output_file[:-3]+"_posteriorprob.db")
    """
    ####################################################################
    ###### 2. iterate to refine the labels based on their adjacency with gradient_descent
    ####################################################################
    if options.gradient_iteration_num==0:
        sys.exit("No gradient iteration is set")
    max_gradient_iteration= options.gradient_iteration_num
    #if (max_gradient_iteration < 3):
    #max_gradient_iteration = 3

    print "\n*************************** Start Gradient Descent autolabelling ********************************\n"
    print "\nthe current global energy = ", U_current,", the optimum global energy with true labels = " ,U_optimum
    edgelist_visit = range(feature_vector.shape[0])
    ##### initialize T, eta and max_global_iteration => they come from options
    global_U_threshold = 1e-200	#threshold for global energy variations in each global iteration
    global_U = []			#global energy
    global_err = []		#global labelling error
    global_vol_err = []		#global labelling volumetric error
    percent_pos_delE = []
    temperatures=[]

    minE=U_current
    minE_it=0
    minE_l=estimated_target_vector
    minE_err=100*float(error_num)/float(labeled_num+machineEpsilon(float))
    minE_errvol=100*float(error_volume)/float(total_volume+machineEpsilon(float))
    minErr=100*float(error_num)/float(labeled_num+machineEpsilon(float))
    minErrvol=100*float(error_volume)/float(total_volume+machineEpsilon(float))
    minErr_it=0
    minErr_l=estimated_target_vector
    minErr_energy=U_current


    #### run gradient descent for 20 iterations then do stochastic while loop until global iterations reach maximum or global energy converges:
    U_current = calc_globalE(estimated_target_vector,local_potentials,adjacency_matrix,edge_neighboring_indx,labelNumerics)	#with current labelling of all edges that come from local_potentials[:][0][0]

    global_cnt = 0	
    while_param = True
    while (while_param):
        random.shuffle(edgelist_visit)	#### browse nodes in random order
        #U_current = calc_globalE(estimated_target_vector,local_potentials,adjacency_matrix,edge_neighboring_indx,labelNumerics)	#with current labelling of all edges that come from local_potentials[:][0][0]
        if np.isnan(U_current) or np.isinf(U_current):
            print "Error! U_current undefined!\nAborted"
            exit(0)
        global_U.append(U_current)
        [error_num,labeled_num,total_volume,error_volume]= error_calculation(target_vector, estimated_target_vector,global_cnt,g,edge_w_indx,U_current)
        if (labeled_num > 0):	#otherwise test set was not labeled and we can't calculate error of labeling
            global_err.append(100*float(error_num)/float(labeled_num+machineEpsilon(float)))
            global_vol_err.append(100*float(error_volume)/float(total_volume+machineEpsilon(float)))
            #saving best error configuration
            if global_err[global_cnt]<minErr:
                minErr = global_err[global_cnt]
                minErrvol=global_vol_err[global_cnt]
                minErr_it = global_cnt
                minErr_l=estimated_target_vector
                minErr_energy=U_current
                minErr_T = options.temprature
    
        #saving best energy configuration		
        if U_current<minE:
            minE = U_current
            minE_it = global_cnt
            minE_l=estimated_target_vector
            minE_err=100*float(error_num)/float(labeled_num+machineEpsilon(float))
            minE_errvol=100*float(error_volume)/float(total_volume+machineEpsilon(float))
            minE_T = options.temprature
    
        for i in edgelist_visit:			#for each node: X = feature_vector[i,:], current_l = local_potentials[i][0][0]
            l_current = estimated_target_vector[i]  
            #### calculate pdf of transition to new labels for current node			
            delta_U = [] #energy variation of transition to other labels for current edge[i]
            for k in range(len(labelNumerics)):	#calculate energy variation for transition to l_new
                l_new = labelNumerics[k]
                current_delta_U = calc_deltaE(i,l_current,l_new,estimated_target_vector,local_potentials,adjacency_matrix,edge_neighboring_indx,labelNumerics)
                delta_U.append(current_delta_U)    
            #### options.optimization_nature=0	#gradient_descent optimization 
            ### for now let's only choose the label with minimum delta_U
            if min(delta_U) < 0:    #if it's minimizing the energy function    
                estimated_target_vector[i] = labelNumerics[delta_U.index(min(delta_U))]
                U_current = U_current + min(delta_U)    #### update energy function with this label swap
                #print "relabeled to ",labelNumerics[delta_U.index(min(delta_U))]    
    
        #### at the end of each global iteration
        global_cnt = global_cnt + 1
        if (global_cnt > max_gradient_iteration):
            while_param = False
        elif (global_cnt > 10):	#run at least 10 iterations	
            tmp_cnt = 0
            for converge_i in range(1,11):#if for 10 consecutive iterations the global energy has not changed => energy converged
                tmp_converge = global_U[global_cnt-converge_i-1]-global_U[global_cnt-converge_i]
            if (tmp_converge < global_U_threshold):
                tmp_cnt = tmp_cnt + 1
                if (tmp_cnt > 9):
                    while_param = False

    ###################################################################
    ##### 2.1. write labeled output
    ###################################################################
    [error_num,labeled_num,total_volume,error_volume]=error_calculation (target_vector, estimated_target_vector,global_cnt,g,edge_w_indx,U_current)
    RR_by_label = error_calculation_by_label(RR_by_label,target_vector, estimated_target_vector, local_potentials,feature_vector, global_cnt,edge_neighboring_indx,g,edge_w_indx)
    if options.verbose:
        confusion_matrix_calculation(target_vector, estimated_target_vector,labelNumerics)
    if options.plot and options.gradient_iteration_num!=0:
        posterior_probas, posteriorprobs_currentlabel = calc_posteriorprob(feature_vector,estimated_target_vector,labelNumerics,class_prior, mean_val, std_val)
        estimated_target_vector_prob = [lable_prob[1] for lable_prob in posteriorprobs_currentlabel]
        write_output_result(g,edge_w_indx,local_potentials,estimated_target_vector,target_vector, output_file[:-3]+"_GD.db",estimated_target_vector_prob)
    #GD state
    E_GD = U_current
    RR_GD = 100 - 100*float(error_num)/float(labeled_num+machineEpsilon(float))
    RRvol_GD = 100 -100*float(error_volume)/float(total_volume+machineEpsilon(float))
    GD_num = int(labeled_num-error_num)
    GD_vol = total_volume-error_volume
    GD_it = global_cnt

    ####################################################################
    ###### 3. iterate to refine the labels based on their adjacency with simulated annealing
    ####################################################################
    if options.iteration_num==0:
        sys.exit("No stochastic simulated annealing iteration is set")
    print "\n*************************** Start Simulated Annealing autolabelling ********************************\n"
    print "\nthe current global energy = ", U_current,", the optimum global energy with true labels = " ,U_optimum,", the initial temperature = " ,options.temprature
    max_global_iteration= options.iteration_num
    T = options.temprature
    eta = options.anneal_param

    if options.nonlin_annealing:
        #T = 1e-10
        eta = 1.0- exp(-(0+2500.0)/800.0)		# eta = 1-exp(-(it+2500)/800)
    
    init_T = T
    global_cnt=0
    while_param = True
    anneal_change=0
    
    #randomize the labelling state
    for i in edgelist_visit:
        estimated_target_vector[i] = labelNumerics[int(round(random.uniform(0, len(labelNumerics)-1)))]

    U_current = calc_globalE(estimated_target_vector,local_potentials,adjacency_matrix,edge_neighboring_indx,labelNumerics) #with current labelling of all edges that come from local_potentials[:][0][0]
    tic = time.time()
    while (while_param):
        temperatures.append(T)
        ##print "\niteration:",global_cnt
        ##it_tic= time.time()
        pos_delE_cnt=0
        ##random.seed(422)
        random.shuffle(edgelist_visit)  #### browse nodes in random order
        #U_current = calc_globalE(estimated_target_vector,local_potentials,adjacency_matrix,edge_neighboring_indx,labelNumerics) #with current labelling of all edges that come from local_potentials[:][0][0]
        global_U.append(U_current)
        [error_num,labeled_num,total_volume,error_volume]= error_calculation(target_vector, estimated_target_vector,global_cnt,g,edge_w_indx,U_current,T)
        if (labeled_num > 0):   #otherwise test set was not labeled and we can't calculate error of labeling
            global_err.append(100*float(error_num)/float(labeled_num+machineEpsilon(float)))
            global_vol_err.append(100*float(error_volume)/float(total_volume+machineEpsilon(float)))
            #saving best error configuration
            if global_err[global_cnt]<minErr:
                minErr = global_err[global_cnt]
                minErrvol=global_vol_err[global_cnt]
                minErr_it = global_cnt
                minErr_l=estimated_target_vector
                minErr_energy=U_current
                minErr_T = T
    
        #saving best energy configuration       
        if U_current<minE:
            minE = U_current
            minE_it = global_cnt
            minE_l = estimated_target_vector
            minE_err = 100*float(error_num)/float(labeled_num+machineEpsilon(float))
            minE_err = 100*float(error_volume)/float(total_volume+machineEpsilon(float))
            minE_T = T
            
        for i in edgelist_visit:            #for each node: X = feature_vector[i,:], current_l = local_potentials[i][0][0]
            l_current = estimated_target_vector[i]
            #print "\ni:",i, " current l=", l_current, " true l=",int(target_vector[i])  #, "\nlabels visited:"
            ##U_current = calc_globalE(estimated_target_vector,local_potentials,adjacency_matrix,edge_neighboring_indx,labelNumerics)   #with current labelling of all edges that come from local_potentials[:][0][0]
            
            #### calculate pdf of transition to new labels for current node         
            #pdf_ltransition =[]    #probability density function of transition to other labels for current edge[i]
            delta_U = [] #energy variation of transition to other labels for current edge[i]
            total_ltransition = 0
            for k in range(len(labelNumerics)): #calculate energy variation for transition to l_new
                l_new = labelNumerics[k]
                #estimated_target_vector_new = copy.deepcopy(estimated_target_vector)
                #estimated_target_vector_new[i] = l_new
                #U_new = calc_globalE(estimated_target_vector_new,local_potentials,adjacency_matrix,edge_neighboring_indx,labelNumerics)   #with current labelling of all edges that come from local_potentials[:][0][0] and replacing local_potentials[i][0][0] to l_new
                #current_delta_U = U_new-U_current
                current_delta_U = calc_deltaE(i,l_current,l_new,estimated_target_vector,local_potentials,adjacency_matrix,edge_neighboring_indx,labelNumerics)
                delta_U.append(current_delta_U)
            #### stochastic optimization
            #### normalized pdf
            pdf_ltransition=[]
            #### normalize delta_U to be in range [-1...1]: (so map delta_U =[0...max] -> [0 ..1] and negative delta_U will be interpolated linearly
            ##normal_delta_U = [(2.0*(x - min(delta_U))/(max(delta_U)-min(delta_U))-1) for x in delta_U]
            normal_delta_U = [(x /(max(abs(max(delta_U)),abs(min(delta_U)))+machineEpsilon(float))) for x in delta_U]           #### ??? pick a constant value for max(delta_U) for all edges and all datasets?
            #print "normal_delta_U:",normal_delta_U 
            for k in range(len(labelNumerics)): #calculate energy variation for transition to l_new
                exp_val=exp(-normal_delta_U[k]/(T+machineEpsilon(float)))
                if np.isinf(exp_val) or np.isnan(exp_val):
                    exp_val= 1/machineEpsilon(float)
                total_ltransition = total_ltransition + exp_val
                pdf_ltransition.append(exp_val)
            #print "total_ltransition:",total_ltransition
            pdf_ltransition = [(x/(total_ltransition+machineEpsilon(float))) for x in pdf_ltransition]
            #print "pdf_ltransition:",pdf_ltransition
            #### actual transition by drawing from pdf_ltransition
            # #local_potentials contains potential of all labels for each edge based on local features!
            # we change the local_potentials[i][0][0], better to keep probablities of local features and have another vector for current estimated labelings!!!
            ltransition_indx = draw_from_dist(pdf_ltransition,1)[0]
            estimated_target_vector[i] = labelNumerics[ltransition_indx]
            #### Update energy with this label swap
            U_current = U_current + calc_deltaE(i,l_current,estimated_target_vector[i],estimated_target_vector,local_potentials,adjacency_matrix,edge_neighboring_indx,labelNumerics)
    #### update energy function with this label swap
            #print "relabeled to ",labelNumerics[ltransition_indx]
            if normal_delta_U[ltransition_indx]>0:
                pos_delE_cnt+=1
        
        percent_pos_delE.append(100.0*float(pos_delE_cnt)/float(len(edgelist_visit)))
    
        #### at the end of each global iteration
        #print ("Temprature is %f " %T)
        if options.animation and global_cnt%options.decimate==0:
            write_output_result(g,edge_w_indx,local_potentials,estimated_target_vector,target_vector, os.path.dirname(os.path.abspath(output_file))+"/it_"+str(global_cnt)+".db")

        global_cnt = global_cnt + 1
        T = eta*T
        if options.nonlin_annealing:
            eta = 1.0- exp(-(global_cnt+2500.0)/800.0)      # eta = 1-exp(-(it+2500)/800)
            
        if (global_cnt > max_global_iteration):
            while_param = False
        elif (global_cnt > 50): #run at least 50 iterations 
            tmp_cnt = 0
            for converge_i in range(1,11):#if for 10 consecutive iterations the global energy has not changed => energy converged
                tmp_converge = global_U[global_cnt-converge_i-1]-global_U[global_cnt-converge_i]
                if (tmp_converge < global_U_threshold):
                    tmp_cnt = tmp_cnt + 1
            if (tmp_cnt > 9):
                while_param = False
                
        if not options.nonlin_annealing:
            if global_cnt>10 and global_err<50 and not anneal_change:
                anneal_change=1
                eta = 0.9+eta/10.0      #eta+(1-eta)*0.9        if eta<1 => eta/10<0.1 => 0.9+eta/10.0<1
                            
    print "\n\nLabelling DONE!"             
    final_T = T
    ####################################################################
    ###### 3.1 write labeled output
    ####################################################################
    [error_num,labeled_num,total_volume,error_volume]=error_calculation (target_vector, estimated_target_vector,global_cnt,g,edge_w_indx,U_current,T)
    RR_by_label = error_calculation_by_label(RR_by_label,target_vector, estimated_target_vector, local_potentials,feature_vector, global_cnt,edge_neighboring_indx,g,edge_w_indx)
    if options.verbose:
        confusion_matrix_calculation(target_vector, estimated_target_vector,labelNumerics)
    if options.plot and options.iteration_num!=0:
        posterior_probas, posteriorprobs_currentlabel = calc_posteriorprob(feature_vector,estimated_target_vector,labelNumerics,class_prior, mean_val, std_val)
        estimated_target_vector_prob = [lable_prob[1] for lable_prob in posteriorprobs_currentlabel]
        write_output_result(g,edge_w_indx,local_potentials,estimated_target_vector,target_vector, output_file[:-3]+"_S.db",estimated_target_vector_prob)
    #SA state
    E_SA = U_current
    RR_SA = 100 - 100*float(error_num)/float(labeled_num+machineEpsilon(float))
    RRvol_SA = 100 -100*float(error_volume)/float(total_volume+machineEpsilon(float))
    SA_num = int(labeled_num-error_num)
    SA_vol = total_volume-error_volume
    SA_it = global_cnt
    SA_T = T

    ####################################################################
    ###### 4. iterate to refine the labels based on their adjacency with stochastic tunneling
    ####################################################################
    if options.STUN_iteration_num==0:
        sys.exit("No stochastic tunneling iteration is set")
    print "\n*************************** Start Stochastic tunneling autolabelling ********************************\n"
    print "\nthe current global energy = ", U_current,", the optimum global energy with true labels = " ,U_optimum,", the initial temperature = " ,T
    max_global_iteration= options.STUN_iteration_num
    T = options.temprature
    eta = options.anneal_param
    tunparam = options.tunneling_param
    #fthresh = options.tunneling_threshold

    if options.nonlin_annealing:
        #T = 1e-10
        eta = 1.0- exp(-(0+2500.0)/800.0)		# eta = 1-exp(-(it+2500)/800)
    
    init_T = T
    global_cnt=0
    while_param = True
    anneal_change=0
    
    #randomize the labelling state
    for i in edgelist_visit:
        estimated_target_vector[i] = labelNumerics[int(round(random.uniform(0, len(labelNumerics)-1)))]
        
    E0 = calc_globalE(estimated_target_vector,local_potentials,adjacency_matrix,edge_neighboring_indx,labelNumerics) #with current labelling of all edges that come from local_potentials[:][0][0]
    U_current = E0
                
    tic = time.time()
    while (while_param):
        temperatures.append(T)
        pos_delE_cnt=0
        random.shuffle(edgelist_visit)  #### browse nodes in random order
        #U_current = calc_globalE(estimated_target_vector,local_potentials,adjacency_matrix,edge_neighboring_indx,labelNumerics) #with current labelling of all edges that come from local_potentials[:][0][0]
        #### U_current was initialized before the while loop, and then is updated at every label swap        
        global_U.append(U_current)
        [error_num,labeled_num,total_volume,error_volume]= error_calculation(target_vector, estimated_target_vector,global_cnt,g,edge_w_indx,U_current,T)
        if (labeled_num > 0):   #otherwise test set was not labeled and we can't calculate error of labeling
            global_err.append(100*float(error_num)/float(labeled_num+machineEpsilon(float)))
            global_vol_err.append(100*float(error_volume)/float(total_volume+machineEpsilon(float)))
            #saving best error configuration
            if global_err[global_cnt]<minErr:
                minErr = global_err[global_cnt]
                minErrvol=global_vol_err[global_cnt]
                minErr_it = global_cnt
                minErr_l=estimated_target_vector
                minErr_energy=U_current
                minErr_T = T
    
        #saving best energy configuration       
        if U_current<minE:
            minE = U_current
            minE_it = global_cnt
            minE_l = estimated_target_vector
            minE_err = 100*float(error_num)/float(labeled_num+machineEpsilon(float))
            minE_err = 100*float(error_volume)/float(total_volume+machineEpsilon(float))
            minE_T = T
            
        for i in edgelist_visit:            #for each node: X = feature_vector[i,:], current_l = local_potentials[i][0][0]
            l_current = estimated_target_vector[i]     
            #USTUNcurr, U_current = calc_STUN_E(estimated_target_vector,local_potentials,adjacency_matrix,edge_neighboring_indx,labelNumerics,E0,tunparam)
            #USTUNcurr, Ucurr = calc_STUN_E_quick(estimated_target_vector,local_potentials,adjacency_matrix,edge_neighboring_indx,labelNumerics,E0,tunparam,U_current,edge_indx,l_current,l_current):
            USTUNcurr = 1-np.exp(-tunparam*(U_current-E0))
            
            delta_U = [] #energy variation of transition to other labels for current edge[i]
            total_ltransition = 0
            for k in range(len(labelNumerics)): #calculate energy variation for transition to l_new
                l_new = labelNumerics[k]
                estimated_target_vector_new = copy.deepcopy(estimated_target_vector)
                estimated_target_vector_new[i] = l_new
                #USTUNnew, Unew = calc_STUN_E(estimated_target_vector_new,local_potentials,adjacency_matrix,edge_neighboring_indx,labelNumerics,E0,tunparam)
                USTUNnew, Unew = calc_STUN_E_quick(estimated_target_vector_new,local_potentials,adjacency_matrix,edge_neighboring_indx,labelNumerics,E0,tunparam,U_current,i,l_current,l_new)
                current_delta_U = USTUNnew - USTUNcurr
                delta_U.append(current_delta_U)
            #### stochastic optimization
            #### normalized pdf
            pdf_ltransition=[]
            #### normalize delta_U to be in range [-1...1]: (so map delta_U =[0...max] -> [0 ..1] and negative delta_U will be interpolated linearly
            ##normal_delta_U = [(2.0*(x - min(delta_U))/(max(delta_U)-min(delta_U))-1) for x in delta_U]
            normal_delta_U = [(x /(max(abs(max(delta_U)),abs(min(delta_U)))+machineEpsilon(float))) for x in delta_U]           #### ??? pick a constant value for max(delta_U) for all edges and all datasets?
            #print "normal_delta_U:",normal_delta_U 
            for k in range(len(labelNumerics)): #calculate energy variation for transition to l_new
                exp_val=exp(-normal_delta_U[k]/(T+machineEpsilon(float)))
                if np.isinf(exp_val) or np.isnan(exp_val):
                    exp_val= 1/machineEpsilon(float)
                total_ltransition = total_ltransition + exp_val
                pdf_ltransition.append(exp_val)
            #print "total_ltransition:",total_ltransition
            pdf_ltransition = [(x/(total_ltransition+machineEpsilon(float))) for x in pdf_ltransition]
            #print "pdf_ltransition:",pdf_ltransition
            #### actual transition by drawing from pdf_ltransition
            # #local_potentials contains potential of all labels for each edge based on local features!
            # we change the local_potentials[i][0][0], better to keep probablities of local features and have another vector for current estimated labelings!!!
            ltransition_indx = draw_from_dist(pdf_ltransition,1)[0]
            estimated_target_vector[i] = labelNumerics[ltransition_indx]
            #### calculate U_current with this label swap
            USTUNnew, U_current = calc_STUN_E_quick(estimated_target_vector_new,local_potentials,adjacency_matrix,edge_neighboring_indx,labelNumerics,E0,tunparam,U_current,i,l_current,estimated_target_vector[i])
            if U_current < E0:
                E0 = U_current

            #print "relabeled to ",labelNumerics[ltransition_indx]
            if normal_delta_U[ltransition_indx]>0:
                pos_delE_cnt+=1
        
        percent_pos_delE.append(100.0*float(pos_delE_cnt)/float(len(edgelist_visit)))
    
        #### at the end of each global iteration
        #print ("Temprature is %f " %T)
        if options.animation and global_cnt%options.decimate==0:
            write_output_result(g,edge_w_indx,local_potentials,estimated_target_vector,target_vector, os.path.dirname(os.path.abspath(output_file))+"/it_"+str(global_cnt)+".db")

        global_cnt = global_cnt + 1
        T = eta*T
        if options.nonlin_annealing:
            eta = 1.0- exp(-(global_cnt+2500.0)/800.0)      # eta = 1-exp(-(it+2500)/800)
            
        if (global_cnt > max_global_iteration):
            while_param = False
        elif (global_cnt > 50): #run at least 50 iterations 
            tmp_cnt = 0
            for converge_i in range(1,11):#if for 10 consecutive iterations the global energy has not changed => energy converged
                tmp_converge = global_U[global_cnt-converge_i-1]-global_U[global_cnt-converge_i]
                if (tmp_converge < global_U_threshold):
                    tmp_cnt = tmp_cnt + 1
            if (tmp_cnt > 9):
                while_param = False
                
        if not options.nonlin_annealing:
            if global_cnt>10 and global_err<50 and not anneal_change:
                anneal_change=1
                eta = 0.9+eta/10.0      #eta+(1-eta)*0.9        if eta<1 => eta/10<0.1 => 0.9+eta/10.0<1
                            
    print "\n\nLabelling DONE!"             
    final_T = T
    
    #### 1. calculate current state
    #[temperatures, global_U,global_err, global_vol_err,minErr, minErrvol, minErr_it, minErr_l, minErr_energy, minErr_T,minErr_optimizer,minvol_Err, minvolErr_vol, minvolErr_it, minvolErr_l, minvolErr_energy, minvolErr_T, minvolErr_optimizer,minE_optimizer,minE , minE_it ,minE_l ,minE_err ,minE_T,U_current] = current_state (temperatures, global_U,global_err, global_vol_err,minErr, minErrvol, minErr_it, minErr_l, minErr_energy, minErr_T,minErr_optimizer,minvol_Err, minvolErr_vol, minvolErr_it, minvolErr_l, minvolErr_energy, minvolErr_T, minvolErr_optimizer,minE_optimizer, minE , minE_it ,minE_l ,minE_err ,minE_T, target_vector, estimated_target_vector,labeling_potentials,global_cnt,g,edge_w_indx,T, tunparam,optimizer,U_optimum)		#,U_optimum
    #### 2. calculate deltaE and pdf of transition to new labels for each node	
    [estimated_target_vector, labeling_potentials, percent_pos_delE] = calc_transition_pdf (target_vector, estimated_target_vector, local_features_potentials, adjacency_matrix, edge_neighboring_indx, labelNumerics, minE, tunparam, optimizer,U_optimum)		#,U_optimum
    #### 3. update T and check number of iterations
    [while_param, global_cnt, T, eta, anneal_change, estimated_target_vector] = optimizer_T_update (global_cnt,global_err, T, eta,estimated_target_vector, edge_neighboring_indx, global_U , global_U_threshold, max_stun_iteration,feature_vector, anneal_change, options.nonlin_annealing, minE, tunparam, optimizer,U_optimum,fthresh)		#,U_optimum
    
    U_current = calc_globalE(estimated_target_vector,local_potentials,adjacency_matrix,edge_neighboring_indx,labelNumerics) #with current labelling of all edges that come from local_potentials[:][0][0]
    global_U.append(U_current)
    [error_num,labeled_num,total_volume,error_volume]= error_calculation(target_vector, estimated_target_vector,global_cnt,g,edge_w_indx,U_current,T)
    if (labeled_num > 0):   #otherwise test set was not labeled and we can't calculate error of labeling
        global_err.append(100*float(error_num)/float(labeled_num+machineEpsilon(float)))
        global_vol_err.append(100*float(error_volume)/float(total_volume+machineEpsilon(float)))
        #saving best error configuration
        if global_err[global_cnt]<minErr:
            minErr = global_err[global_cnt]
            minErrvol=global_vol_err[global_cnt]
            minErr_it = global_cnt
            minErr_l=estimated_target_vector
            minErr_energy=U_current
            minErr_T = T

    #saving best energy configuration       
    if U_current<minE:
        minE = U_current
        minE_it = global_cnt
        minE_l = estimated_target_vector
        minE_err = 100*float(error_num)/float(labeled_num+machineEpsilon(float))
        minE_err = 100*float(error_volume)/float(total_volume+machineEpsilon(float))
        minE_T = T
        
    #### at the end of each global iteration
    #print ("Temprature is %f " %T)
    if options.animation and global_cnt%options.decimate==0:
        write_output_result(g,edge_w_indx,local_potentials,estimated_target_vector,target_vector, os.path.dirname(os.path.abspath(output_file))+"/it_STUN_"+str(global_cnt)+".db")
    
    global_cnt = global_cnt + 1
    T = eta*T
    if options.nonlin_annealing:
        eta = 1.0- exp(-(global_cnt+2500.0)/800.0)      # eta = 1-exp(-(it+2500)/800)
        
    if (global_cnt > max_global_iteration):
        while_param = False
    elif (global_cnt > 50): #run at least 50 iterations 
        tmp_cnt = 0
        for converge_i in range(1,11):#if for 10 consecutive iterations the global energy has not changed => energy converged
            tmp_converge = global_U[global_cnt-converge_i-1]-global_U[global_cnt-converge_i]
            if (tmp_converge < global_U_threshold):
                tmp_cnt = tmp_cnt + 1
        if (tmp_cnt > 9):
            while_param = False
            
    if not options.nonlin_annealing:
        if global_cnt>10 and global_err<50 and not anneal_change:
            anneal_change=1
            eta = 0.9+eta/10.0      #eta+(1-eta)*0.9        if eta<1 => eta/10<0.1 => 0.9+eta/10.0<1
    print "\n\nLabelling DONE!"             
    toc = time.time()-tic
    hr = int(toc)/3600
    minut = int(toc - (hr*3600))/60
    sec = toc - (hr*3600) - (minut*60)
    print ("\n\n\nThe STUN optimization took %d hours and %d minutes and %d seconds\nEND!" %(hr,minut,sec)) 	
    final_T_stun = T
    ####################################################################
    ###### 4.1 write labeled output
    ####################################################################
    [error_num,labeled_num,total_volume,error_volume]=error_calculation (target_vector, estimated_target_vector,global_cnt,g,edge_w_indx,U_current,T)
    RR_by_label = error_calculation_by_label(RR_by_label,target_vector, estimated_target_vector, local_potentials,feature_vector, global_cnt,edge_neighboring_indx,g,edge_w_indx)
    if options.verbose:
        confusion_matrix_calculation(target_vector, estimated_target_vector,labelNumerics)
    if options.plot and options.STUN_iteration_num!=0:
        posterior_probas, posteriorprobs_currentlabel = calc_posteriorprob(feature_vector,estimated_target_vector,labelNumerics,class_prior, mean_val, std_val)
        estimated_target_vector_prob = [lable_prob[1] for lable_prob in posteriorprobs_currentlabel]
        write_output_result(g,edge_w_indx,local_potentials,estimated_target_vector,target_vector, output_file[:-3]+"_STUN.db",estimated_target_vector_prob)

    #STUN state
    E_STUN = U_current
    RR_STUN = 100 - 100*float(error_num)/float(labeled_num+machineEpsilon(float))
    RRvol_STUN = 100 -100*float(error_volume)/float(total_volume+machineEpsilon(float))
    STUN_num = int(labeled_num-error_num)
    STUN_vol = total_volume-error_volume
    STUN_it = global_cnt
    STUN_T = T

    ####################################################################
    ###### 5.1 write best case scenario labeled output
    ####################################################################
    if options.verbose:
        write_output_result(g,edge_w_indx,labeling_potentials,minE_l,target_vector, output_file[:-3]+"_minE.db")
        write_output_result(g,edge_w_indx,labeling_potentials,minErr_l,target_vector, output_file[:-3]+"_minErr.db")
        write_output_result(g,edge_w_indx,labeling_potentials,minvolErr_l,target_vector, output_file[:-3]+"_minvolErr.db")

    ####################################################################
    ###### 5. Final report!
    ####################################################################
    print( "global_Energywiterations = "),global_U
    print( "\nglobal_Errorwiterations = "),global_err
    print( "\nglobal_volErrorwiterations = "),global_vol_err
    #print( "\npercent_pos_delEwiterations = "),percent_pos_delE
    #print( "\ntempreturewiterations = "),temperatures


    print "iteration ",minErr_it, " had the minimum error ", minErr, " and its energy is ", minErr_energy	
    print "iteration ",minE_it, " had the minimum Energy ", minE, " and its error is", minE_err
    print "\nthe current global energy = ", U_current,", the optimum global energy with true labels = " ,U_optimum
    print (" \n\nFinal Iteration(T="+str(T)+("E=%f): Out of %d edges that were labeled, %d were labelled correctly with Naive-Bayes. The error is %f%% (volumetric error %f%%). Recognition Rate is %f%% (volumetric RR %f%%)." %(U_current,labeled_num,labeled_num-error_num, 100*float(error_num)/float(labeled_num), 100*float(error_volume)/float(total_volume), 100-100*float(error_num)/float(labeled_num), 100-100*float(error_volume)/float(total_volume) )) )

    toc = time.time()-tic

    hr = int(toc)/3600
    minut = int(toc - (hr*3600))/60
    sec = toc - (hr*3600) - (minut*60)
    print ("\n\n\nThe labeling took %d hours and %d minutes and %d seconds\nEND!" %(hr,minut,sec)) 	

    txtfile = output_file[:-3]+"_summary.txt"
    f = open(txtfile, 'w')
    txt = ("time %d hrs %d min %d sec\n"%(hr,minut,sec)) 
    f.write(txt)
    txt = ("initial temperature "+str(init_T)+"\n") 
    f.write(txt)
    txt = ("ground truth labeling energy %e\n"%(U_optimum))
    f.write(txt)
    txt = ("final iteration %d \n"%global_cnt) 
    f.write(txt)
    txt = ("final temperature "+str(final_T)+"\n") 
    f.write(txt)
    txt = ("final energy "+str(U_current)+"\n") 
    f.write(txt)
    txt = ("initial energy "+str(init_E)+"\n") 
    f.write(txt)
    txt = ("GD energy "+str(E_GD)+"\n") 
    f.write(txt)
    txt = ("true_label energy "+str(U_optimum)+"\n") 
    f.write(txt)

    txt = ("local MAP state: Iteration 0, Temperature %e, Energy %e, RR %f, RRvol %f\n"%(init_T,init_E,RR_initial, RRvol_initial))
    f.write(txt)
    txt = ("GD state: Iteration %d, Temperature %e, Energy %e, RR %f, RRvol %f\n"%(GD_it,init_T,E_GD,RR_GD,RRvol_GD))
    f.write(txt)
    txt = ("SA state: Iteration %d, Temperature %e, Energy %e, RR %f, RRvol %f\n"%(SA_it, SA_T, E_SA,RR_SA,RRvol_SA))
    f.write(txt)
    txt = ("STUN state: Iteration %d, Temperature %e, Energy %e, RR %f, RRvol %f\n"%(STUN_it, STUN_T,E_STUN,RR_STUN,RRvol_STUN))
    f.write(txt)
    #
    #txt = ("highest RR state: Energy %e, RR %f, RRvol %f \n"%(minErr_energy,(100-minErr),(100-minErrvol)))
    #f.write(txt)
    #txt = ("highest RRvol state: Energy %e, RR %f, RRvol %f\n"%(minvolErr_energy,(100-minvol_Err),(100-minvolErr_vol)))  
    #f.write(txt)
    #txt = ("lowest energy state: Energy %e, RR %f, RRvol %f\n"%(minE,(100-minE_err),(100-minE_errvol)))
    #f.write(txt)
    #txt = ("final Iteration %d, Temperature %e, Energy %e \n"%(global_cnt,final_T,U_current)) 
    #f.write(txt)

    #txt = ("highest RR state: Iteration %d, Temperature %e, Energy %e, RR %f, RRvol %f, Optimizer %s\n"%(minErr_it,minErr_T,minErr_energy,(100-minErr),(100-minErrvol),minErr_optimizer))
    #f.write(txt)
    #txt = ("highest RRvol state: Iteration %d, Temperature %e, Energy %e, RR %f, RRvol %f, Optimizer %s\n"%(minvolErr_it,minvolErr_T,minvolErr_energy,(100-minvol_Err),(100-minvolErr_vol),minvolErr_optimizer))  
    #f.write(txt)
    #txt = ("lowest energy state: Iteration %d, Temperature %e, Energy %e, RR %f, RRvol %f, Optimizer %s\n"%(minE_it,minE_T,minE,(100-minE_err),(100-minE_errvol),minE_optimizer))
    #f.write(txt)
    #txt = ("final Iteration %d, Temperature %e, Energy %e \n"%(global_cnt,final_T,U_current)) 
    #f.write(txt)

    #txt = ("initial RR "+str(RR_initial)+"%\tvolumetricRR " + str(RRvol_initial)+"\n") 
    #f.write(txt)
    #txt = ("GD RR "+str(RR_GD)+"%\tvolumetricRR " + str(RRvol_GD)+"\n") 
    #f.write(txt)
    #txt = ("SR RR "+str(100-100*float(error_num)/float(labeled_num))+"%\tvolumetricRR " + str(100-100*float(error_volume)/float(total_volume))+"\n\n") 
    #f.write(txt)
    #txt = ("highest RR "+str(100-minErr)+"% (volRR "+str(100-minErrvol)+"%)\titeration " + str(minErr_it)+"\ttemperature "+str(minErr_T)+"\tenergy "+str(minErr_energy)+"\n") 
    #f.write(txt)
    #txt = ("lowest energy RR "+str(100-minE_err)+"% (volRR "+str(100-minE_errvol)+"%)\titeration " + str(minE_it)+"\ttemperature "+str(minE_T)+"\tenergy "+str(minE)+"\n\n") 
    #f.write(txt)
    #
    txt = ("label,\t#vessels,\tvolume,\t\tRR%_initial,\t\t\tRRvol_initial,\t\t#vessels_initial,\tvolume_initial,\t\tRR%_GD,\t\tRRvol_GD,\t\t#vessels_GD,\t\tvolume_GD,\t\t\tRR%_SR,\t\t\tRRvol_SR,\t\t#vessels_SR,\t\tvolume_SR \n") 
    f.write(txt)
    txt = ("total,\t"+str(int(labeled_num))+",\t\t\t"+"%.4f,\t\t"%(total_volume)+"%.2f,\t\t\t\t"%(RR_initial)+"%.2f,\t\t\t\t"%(RRvol_initial)+"%d,\t\t\t\t\t"%(int(initial_num))+"%.4f,\t\t\t\t"%(initial_vol)+"%.2f,\t\t\t\t"%(RR_GD)+"%.2f,\t\t\t\t"%(RRvol_GD)+"%d,\t\t\t\t\t"%(int(GD_num))+"%.4f,\t\t\t\t"%(GD_vol)+"%.2f,\t\t\t\t"%(100-100*float(error_num)/float(labeled_num))+"%.2f,\t\t\t\t"%(100-100*float(error_volume)/float(total_volume))+"%d,\t\t\t\t\t"%(int(labeled_num-error_num))+"%.4f\n"%(total_volume-error_volume))
    f.write(txt)
    for l in RR_by_label.keys():
        f.write(str(int(l))+",\t\t")
        f.write(str(int(RR_by_label[l][0]))+",\t\t\t")
        f.write("%.4f,\t\t"%(RR_by_label[l][1]))
        f.write("%.2f,\t\t\t\t"%(RR_by_label[l][2]))
        f.write("%.2f,\t\t\t\t"%(RR_by_label[l][3]))
        f.write("%d,\t\t\t\t\t"%(int(RR_by_label[l][4])))		
        f.write("%.4f,\t\t\t\t"%(RR_by_label[l][5]))
        f.write("%.2f,\t\t"%(RR_by_label[l][6]))	
        f.write("%.2f,\t\t\t"%(RR_by_label[l][7]))
        f.write("%d"%(int(RR_by_label[l][8]))+",\t\t\t\t")
        f.write("%.4f,\t\t\t\t"%(RR_by_label[l][9]))	
        f.write("%.2f,\t\t\t"%(RR_by_label[l][10]))
        f.write("%.2f,\t\t\t"%(RR_by_label[l][11]))
        f.write("%d"%(int(RR_by_label[l][12]))+",\t\t\t\t")	
        f.write("%.4f"%(RR_by_label[l][13])+"\n")
    f.close()












