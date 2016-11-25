#!/usr/bin/env python
# -*- coding: utf-8 -*-

# get the feature_vector of test dataset and the result of training stage (distribution parameters and the adjacency_matrix and classifier type) and perform classification
# speed up delta_U calc compare to v0_1.py
# direction kent distribution out of PCA!
#
# Added option --plot and --verbose to bayesClassifier_v4
#
#
#
#  Created May 27, 2012
#  modified July 25, 2013  add features rel_diameter and rel_dir & change output to summary.txt
#  modified April 14, 2015  in error calculation check if edge has cyl_height and cyl_radius
#  Sahar Ghanavati


from sys import argv
import math, os, shelve, string ,time, sys
import commands
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
#import matplotlib.pyplot as plt
from pylab import *		#needed for command find
from operator import itemgetter
import random
import copy
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
            print labelNumerics[j],"\t\t",
        print " "	
        for i in range(adjacency_matrix.shape[0]):
            print "\n", labelNumerics[i],"\t",
            for j in range(15):
                print "%.4f"%adjacency_matrix[i,j],"\t",
        print "\n\n"
    else:
        print "\t",
        for j in range(len(labelNumerics)):
            print labelNumerics[j],"\t\t",
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
            print labelNumerics[j],"\t\t",
        print " "	
        for i in range(adjacency_matrix.shape[0]):
            print "\n", labelNumerics[i],"\t",
            for j in range(15,30):
                print "%.4f"%adjacency_matrix[i,j],"\t",
        print "\n\n"
    else:
        print "\t",
        for j in range(15,len(labelNumerics)):
            print labelNumerics[j],"\t\t",
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
            print labelNumerics[j],"\t\t",
        print " "	
        for i in range(adjacency_matrix.shape[0]):
            print "\n", labelNumerics[i],"\t",
            for j in range(30,45):
                print "%.4f"%adjacency_matrix[i,j],"\t",
        print "\n\n"
    else:
        print "\t",
        for j in range(30,len(labelNumerics)):
            print labelNumerics[j],"\t\t",
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
        print labelNumerics[j],"\t\t",
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
        

def graph2featureVect(g,featureNames):		#featureVect
    #### initialize the feature_vector to be zeros(#datapoints=sum(len(g.edge_list()) w label!=0 , len(featureNames))
    ##	initialize the target_vector to be zeros(#datapoints=sum(len(g.edge_list()) w label!=0 , 1)
    MRI_labels=[]
    if 'mri_label_dist' in featureNames:
        MRI_labels = [l[0] for l in g.edge_property(g.edge_list()[0],'mri_label_dist')]
        feature_vector = Matlib.empty((len(g.edge_list()),len(featureNames)+len(MRI_labels)-1),float)
    else:
        feature_vector = Matlib.empty((len(g.edge_list()),len(featureNames)),float)
    #feature_vector = Matlib.empty((len(g.edge_list()),len(featureNames)),float)
    target_vector = []	#Matlib.empty((len(g.edge_list()),1),float)
    labelNumerics = []
    edge_w_indx = {}
    
    
    
    for i in range(len(g.edge_list())):
        e=g.edge_list()[i]
        if 'label' not in g.edge_properties(e).keys():
            print ("ERROR! The edge (%d,%d) in test graph doesn't have label property.\nAborted!\n" %(e[0],e[1]))
            exit (0)
        if g.edge_property(e,'label')>-1:	#including the edges with label 0
            #target_vector= numpy.vstack([target_vector, array(g.edge_property(e,'label'))])
            #target_vector[i] =  array(g.edge_property(e,'label'))
            target_vector.append(int(g.edge_property(e,'label')))
            if g.edge_property(e,'label') not in labelNumerics and g.edge_property(e,'label')>0:
                labelNumerics.append(g.edge_property(e,'label'))
            features=[]
            for f in featureNames:
                #print f, features
                if f not in g.edge_properties(e).keys():
                    print ("ERROR! The edge (%d,%d) in test graph doesn't have feature %s property.\nAborted!\n" %(e[0],e[1],f))
                    exit (0)
                if not f== 'mri_label_dist' and not f=='rel_dir' and not f=='rel_diameter':	
                    features.append(g.edge_property(e,f))
                elif f== 'mri_label_dist':
                    for mr_l in MRI_labels:
                        ind=[l[0] for l in g.edge_property(e,'mri_label_dist')].index(mr_l)
                        features.append(g.edge_property(e,'mri_label_dist')[ind][1])
                elif f=='rel_dir' or f=='rel_diameter':
                    rel_f = 1.0
                    for rel_i in range(min(len(g.edge_property(e,f)),4)):	#we only consider up to 4 adjacent edges, if there are less 1.0 will be used
                        rel_f = rel_f * g.edge_property(e,f)[rel_i]
                    features.append(rel_f)	
                #features.append(g.edge_property(e,f))
            #feature_vector	= numpy.vstack([feature_vector, array(features)])		#hstack #a = matrix([[10,20,30]]); a=append(a,[[1,2,3]],axis=0); a=append(a,[[15],[15]],axis=1)
            #print len(features),features
            feature_vector[i] = array(features)
        edge_w_indx[i]= e 
            
    
    #### !!find neighbouring edges and in iterations find their
    #adjacency_matrix = (1.0/len(labelNumerics))*Matlib.ones((len(labelNumerics),len(labelNumerics)+1),float)		#### smoothing (for 0s in the adj_matrix)! to be minimum of 1 occurance for the nieghbourhood!
    edge_neighboring_indx = {}
    
    for i in range(len(g.edge_list())):
        e=edge_w_indx[i]		#g.edge_list()[i]

        #find edge neighbouring indeces
        neighbor_indx=[]
        e1_neighbours= g.vertices[e[0]].edges
        for v in e1_neighbours:
            if v!=e[1]:
                neighbor_edge=tuple((min(e[0],v),max(e[0],v)))
                neighbor_indx.append(find_key(edge_w_indx, neighbor_edge, i))
        e2_neighbours= g.vertices[e[1]].edges
        for v in e2_neighbours:
            if v!=e[0]:
                neighbor_edge=tuple((min(e[1],v),max(e[1],v)))
                neighbor_indx.append(find_key(edge_w_indx, neighbor_edge, i))
                
        edge_neighboring_indx[i]=neighbor_indx
    
    return feature_vector, target_vector, edge_w_indx, edge_neighboring_indx
    
    

def machineEpsilon(func=float):
    machine_epsilon = func(1)
    while func(1)+func(machine_epsilon) != func(1):
        machine_epsilon_last = machine_epsilon
        machine_epsilon = func(machine_epsilon) / func(2)
    return machine_epsilon_last

def PCA (feature_vector,M_indx,Y_bar=-inf, normalizing_stds = -inf, eigenvects=-inf,eigenvals=-inf):
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
        D,V = numpy.linalg.eig(cov_matrix)		# a: array_like shape (M, M), D:eigenvalues ndarray shape (M,) NOT ordered, V: eignevectors ndarray shape (M, M) The normalized (unit “length”) eigenvectors, column V[:,i] is the eigenvector corresponding to the eigenvalue D[i]			 
        #i = numpy.argsort(D)	#Returns the indices that would sort an array, ascending
        i = numpy.argsort(numpy.array(D))	#D[::-1] #Returns the indices that would sort an array, descending
        inv_i = i[::-1]
        eigenvals = D[inv_i]
        eigenvects = V[:,inv_i]
        
    W = eigenvects[:,0:M_indx]
    
    #### orthogonal_feature_vector = feature_vector_cntr * V * D^(-.5)??   
    #based on karhunen-Leove theorem and data projection from http://en.wikipedia.org/wiki/Principal_component_analysis:
    # projected data = conjugate_transpose(eigenvects)* (feature_vects/Diag(sqrt(eigenvals)))
    # a.H conjugate transpose
    orthogonal_feature_vector = numpy.matrix(W).H * (numpy.matrix(feature_vector_cntrN) * numpy.matrix(numpy.diag(1.0/numpy.sqrt(eigenvals)))).transpose()
    orthogonal_feature_vector = orthogonal_feature_vector.transpose()
        
    
    #### test: correctness of PCA: cov(orthogonal_feature_vector) almost = I
    return [orthogonal_feature_vector, numpy.array(Y_bar), numpy.array(normalizing_stds), numpy.matrix(eigenvects), numpy.array(eigenvals)]
    
    

    
def kent_distribution(L,k,beta,ck,gamma):		#kent_distribution(feature_vector[i,axialfeatureidx],k_axial[k,0],beta_axial[k,0],ck_axial[k,0],gamma_axial[3*k:3*k+3,:])
#f(x) = 1/c(k,beta) exp(k*gamma1*x + beta*((gamma2*x)^2 - (gamma3*x)^2))
    #L is a row
    gamma1 = gamma[:,0]	#a column
    gamma2 = gamma[:,1]	#a column
    gamma3 = gamma[:,2]	#a column

    el1 =  abs(array(L*gamma1))
    el2 =  array(L*gamma2)
    el3 =  array(L*gamma3)

    #bingham: p = (1.0/(4.0*numpy.pi*dk)) * numpy.exp( k1*pow(L*mu1,2) + k2*pow(L*mu2,2) )
    #p = (1.0/(ck)) * numpy.exp( k*(L*gamma1) + beta*(pow(L*gamma2,2) - pow(L*gamma3,2)) )
    p = ck * numpy.exp( k*el1 + beta*(numpy.power(el2,2) - numpy.power(el3,2)) ) 
    return p
    
    
    
def gamma_distribution(x,k,theta):
    p = pow(x,k-1.0) * numpy.exp(-x/theta)/ (special.gamma(k) * pow(theta,k))		# + machineEpsilon()
    return p
    
def univaraite_gaussian(x,mean_val,std_val):
    p = (1.0/(numpy.sqrt(2.0*numpy.pi*pow(std_val,2)))) * numpy.exp( pow(x-mean_val,2) / (-2.0*pow(std_val,2)) )		# + machineEpsilon() ### this is not a good idea, because if the std_val=0 means there was only 1 sample in training set => probability should be 0 (small)for this label and not very large when divided by machineEpsilon
    return p	
    
def multivariate_gaussian(X,mean_vect, cov_mat):
    m = cov_mat.shape[0]
    n = cov_mat.shape[1]
    if (mean_vect.shape[1]!=n):
        print( "ERROR! The length of mean vector is not equal to the size of covariance matrix of multivariate gaussian!\nAborted!\n")
        #log_f.flush()
        exit (0)
    if (n!=m) :
        print( "ERROR! The covariance matrix of multivariate gaussian is not squared!\nAborted!\n")
        #log_f.flush()
        exit (0)

    exp_arg = numpy.matrix((X-mean_vect))* numpy.matrix(numpy.linalg.inv(cov_mat)) * numpy.matrix((X-mean_vect).transpose())	
    p = pow(1.0/(2*numpy.pi) , n/2.0) * pow(numpy.linalg.det(cov_mat) , -0.5) * numpy.exp(-0.5 * exp_arg)
    return p

def error_calculation(target_vector, estimated_target_vector, iteration, g,edge_w_indx, E_current,T=0):
    ####only calculate the error if the ground truth label is available!
    availabel_labels = list(numpy.unique(numpy.array(target_vector)))
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
    availabel_labels = list(numpy.unique(numpy.array(target_vector)))
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
        indx = find (numpy.array(target_vector==l))
        each_label_mindiameter[l]= numpy.min (feature_vector[indx,0])
        each_label_maxdiameter[l]= numpy.max (feature_vector[indx,0])
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
        #if not e in g.edge_list() or 'cyl_radius' not in g.edge_properties(e).keys() or 'cyl_height' not in g.edge_properties(e).keys():
            #print e,' ', g.edge_properties(e).keys()
        #else:
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
    availabel_labels = list(numpy.unique(numpy.array(target_vector)))
    if 0 in availabel_labels:
        availabel_labels.remove(0)

    #availabel_labels2 = numpy.unique(numpy.array(local_potentials[:][0][0]))
    #if 0 in availabel_labels2:
        #availabel_labels2.remove(0)
    
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
    neg_E=0
    for i in range(len(estimated_target_vector)):
        l_e = estimated_target_vector[i]
        if l_e>0 and l_e in labelNumerics:
            local_labels = [local_potentials[i][j][0] for j in range(len(local_potentials[i]))]
            if (-1/(machineEpsilon(float)))!=local_potentials[i][local_labels.index(l_e)][1]:
                neg_E = neg_E + local_potentials[i][local_labels.index(l_e)][1]
                neighb_idx = edge_neighboring_indx[i]
                mat_idx0 = labelNumerics.index(l_e)
                for j in neighb_idx:
                    l_neighb = estimated_target_vector[j]
                    if l_neighb>0 and l_neighb in labelNumerics:
                        mat_idx1 = labelNumerics.index(l_neighb)
                        neg_E = neg_E + log(adjacency_matrix[mat_idx0,mat_idx1])
                #if (len(neighb_idx)!=4):
                if (len(neighb_idx)==2):		# Y edge that has an end-point 
                #for k in range(4-len(neighb_idx)):
                    neg_E = neg_E + log(adjacency_matrix[mat_idx0,len(labelNumerics)])		#neighbor=End_point
                    #else:
                        #print("\n\nERROR! edge %d has %d neighbouring vertices!\nAborted!\n" %(i,len(neighb_idx)))
                        ##log_f.write(write_str)
                        ##log_f.flush()
                        #exit(0)
    E= -neg_E 			
    return E

def calc_deltaE(current_e,current_l,new_l,estimated_target_vector,local_potentials,adjacency_matrix,edge_neighboring_indx,labelNumerics):
    ##delta_U = U_newl - U_currentl = - local_potentials[current_e][newl] + local_potentials[current_e][currentl] \
    ##						    - sum (log(adjacency_matrix[newl,neighbors_l])) + sum (log(adjacency_matrix[currentl,neighbors_l])) \
    ##						    -log(adjacency_matrix[neighbors_l,newl]) + log(adjacency_matrix[neighbors_l,currentl]) :for all neighbors
    
    if current_l<=0 or new_l<=0:
        print "\nERROR! label <=0, current_l:",current_l, " new_l:",new_l,"\nAborted!"
        exit(0)
    
    ##print "\nestimated_target_vector:",estimated_target_vector,
    #delta_E=0
    #if (current_l != new_l):
        #U_current= calc_globalE(estimated_target_vector,local_potentials,adjacency_matrix,edge_neighboring_indx,labelNumerics)
        
        #tmp_target_vect2=[x for x in estimated_target_vector]
        #tmp_target_vect2[current_e]=new_l
        ##print "\ttmp_target_vect2:",tmp_target_vect2,
        #U_new=calc_globalE(tmp_target_vect2,local_potentials,adjacency_matrix,edge_neighboring_indx,labelNumerics)
        
        #delta_E = U_new-U_current
        
        
        
    U_current=0
    U_new=0

    local_labels = [local_potentials[current_e][j][0] for j in range(len(local_potentials[current_e]))]
    neighb_idx = edge_neighboring_indx[current_e]
    
    U_current += local_potentials[current_e][local_labels.index(current_l)][1]
    mat_idx0 = labelNumerics.index(current_l)
    for j in neighb_idx:
        l_neighb = estimated_target_vector[j]
        if l_neighb>0:
            mat_idx1 = labelNumerics.index(l_neighb)
            U_current += log(adjacency_matrix[mat_idx0,mat_idx1])
            U_current += log(adjacency_matrix[mat_idx1,mat_idx0])
    #if (len(neighb_idx)!=4):
    if (len(neighb_idx)==2):
    #for k in range(4-len(neighb_idx)):
        U_current += log(adjacency_matrix[mat_idx0,len(labelNumerics)])		#neighbor=End_point
        #else:
            #print  ("\n\nERROR! edge %d has %d neighbouring vertices!\nAborted!\n" %(current_e,len(neighb_idx)))
            ##log_f.write(write_str)
            ##log_f.flush()
            #exit(0)
            
    U_new += local_potentials[current_e][local_labels.index(new_l)][1]
    mat_idx0 = labelNumerics.index(new_l)
    for j in neighb_idx:
        l_neighb = estimated_target_vector[j]
        if l_neighb>0:
            mat_idx1 = labelNumerics.index(l_neighb)
            U_new += log(adjacency_matrix[mat_idx0,mat_idx1])
            U_new += log(adjacency_matrix[mat_idx1,mat_idx0])
    #if (len(neighb_idx)!=4):
    if (len(neighb_idx)==2):
    #for k in range(4-len(neighb_idx)):
        U_new += log(adjacency_matrix[mat_idx0,len(labelNumerics)])		#neighbor=End_point
        #else:
            #print  ("\n\nERROR! edge %d has %d neighbouring vertices!\nAborted!\n" %(current_e,len(neighb_idx)))
            ##log_f.write(write_str)
            ##log_f.flush()
            #exit(0)
                
    delta_E = U_current-U_new		#these U_new and U_current are actually negative of energies so the negative of differences will be delta_U	
    #print " U_current:",U_current, "U_new:", U_new 	
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
    
    
def write_output_result(g,edge_w_indx,local_potentials,estimated_target_vector,target_vector, output_name):
    history = '>>> %s: auto labelling' % (time.ctime(time.time()))	
    if attributes.has_key("history"):
        history = attributes["history"] + "\n" + history
        del attributes['history']
    
    h=copy.deepcopy(g)
    for i in range(len(g.edge_list())):
        e=edge_w_indx[i]
        h.set_edge_property(e,'estimated_label_potentials',[list(j) for j in local_potentials[i]])	#local_potentials is a list for all edges, for each edge = local_potentials[i] is list of tuples=[(label,highest_prob),(l,prob),...,(l,lowest_prob)]
        h.set_edge_property(e,'estimated_label',[[estimated_target_vector[i],1.0]])	#local_potentials is a list for all edges, for each edge = local_potentials[i] is list of tuples=[(label,highest_prob),(l,prob),...,(l,lowest_prob)]
        if (target_vector[i]>0 and (target_vector[i] == estimated_target_vector[i])):
            h.set_edge_property(e,'error_label',1)	#correct labelling	
        elif (target_vector[i]>0 and (target_vector[i] != estimated_target_vector[i])):
            h.set_edge_property(e,'error_label',2)	#error labelling	
        else:
            h.set_edge_property(e,'error_label',0)
    graph_analysis.output_graph(output_name, h, history, attributes)   
    
    print ("Succefully wrote the %s\n" %output_name)
    
    cmd=("\npython /projects/souris/sghanavati/src/scripts/graph2cylinder.py %s %s --use_estimated_label --clobber " %(output_name,output_name[:-3]+"_cyl.db"))	#python /micehome/jgsled/bin/
    #print(cmd)
    os.system(cmd)

    cmd=("\npython /projects/souris/sghanavati/src/scripts/graph2cylinder.py %s %s --use_error_label --clobber " %(output_name,output_name[:-3]+"_error_cyl.db"))	#python /micehome/jgsled/bin/
    #print(cmd)
    os.system(cmd)
    
    #cmnd=("rm %s" %(output_name[:-3]+"*_cyl.db"))
    #print(cmd)
    #os.system(cmd)
############################################################################################################################################0
program_name = 'bayesClassifier_v5.py'


if __name__ == '__main__':
        
    usage = "Usage: "+program_name+" [options] training_output.db test_graph.db output_filename(test_labeled).db\n"+\
            "   or  "+program_name+" --help";
            
    parser = OptionParser(usage)
    
    parser.add_option("--gradient_iteration_num",type="int", dest="gradient_iteration_num",
                    default=10, help="the number of gradient iterations before performing Gibbs sampling to refine initial labels")                  
    
    parser.add_option("--iteration_num",type="int", dest="iteration_num",
                    default=1000, help="the number of global iterations of Gibbs sampling to refine initial labels")                  
    
    parser.add_option("--temprature",type="float", dest="temprature",
                    default=1.8e-15, help="the initial temperature for the relaxation. T=1.8e-15=1e-12*(.9**60) is in transition (around iteration 60 when T0=1e-12 and annealing .9)")   
                    
    parser.add_option("--annealing", type="float", dest="anneal_param", 
                    default= 0.99, help = "temperature is multiplied by annealing parameter at the end of each global iteration")

    parser.add_option("--nonlin_annealing", action="store_true", dest="nonlin_annealing",
                    default=0, help="nonlinear annealing of 1-exp(-x) format")
    
    parser.add_option("--plot", action="store_true", dest="plot",default=0, help="plot the output labelled graph and error graph")
    
    parser.add_option("--verbose", action="store_true", dest="verbose", default=0, help="spit out detailed execution report to shell")	

    parser.add_option("--animation", action="store_true", dest="animation", default=0, help="Option to save intermediate h5s of labelling results for the purpose of labelling progression movie.")  
    
    parser.add_option("--decimate", type="int", dest="decimate",
                      metavar="decimation_factor", default = 20, 
                      help="decimate the iteration steps for animation making (default 20)")

    options, args = parser.parse_args()
    

    if len(args) != 3:
        parser.error("incorrect number of arguments")
        
    print '\n\n>>> %s [%s]: %s\n' % (time.ctime(time.time()), os.getcwd(), string.join(argv))
    print "\n\n",options, "\n\n", args , "\n\n"
    sys.stdout.flush()
    
    tic= time.time()
    
    training_file, test_file, output_file = args
    


    
    
    #get the training data set
    f_handler = shelve.open(training_file,"r")
    classifier_type = f_handler['classifier_type']
    pca = f_handler['PCA']
    if pca=="Y":
        feat_select_indx = f_handler['feat_select_indx'] 		
        Y_bar = f_handler['Y_bar']
        normalizing_stds = f_handler['normalizing_stds']
        eigenvects = f_handler['eigenvects']
        eigenvals = f_handler['eigenvals']

    class_prior = f_handler['class_prior'] 
    #print "class_prior" , class_prior
    labelNumerics = f_handler['labelNumerics']
    #for i in range(len(labelNumerics)): 
        #print labelNumerics[i] ," : ", class_prior[i,0]
    labelNumeric2Name = f_handler['labelNumeric2Name'] 
    featureNames = f_handler['featureNames']  
    adjacency_matrix = f_handler['adjacency_matrix']
    
    if classifier_type=="NB":
        mean_val = f_handler['mean_val'] 
        std_val = f_handler['std_val'] 
        axialfeatureidx = f_handler['axialfeatureidx']
        if len(axialfeatureidx)>0:
            k_axial = f_handler['k_axial']
            beta_axial = f_handler['beta_axial']
            ck_axial = f_handler['ck_axial']
            gamma_axial = f_handler['gamma_axial']
            
    elif classifier_type == "NB2":
        mean_val = f_handler['mean_val'] 
        std_val = f_handler['std_val']
        k_val = f_handler['k_val'] 
        theta_val = f_handler['theta_val'] 
        axialfeatureidx = f_handler['axialfeatureidx']
        if len(axialfeatureidx)>0:
            k_axial = f_handler['k_axial']
            beta_axial = f_handler['beta_axial']
            ck_axial = f_handler['ck_axial']
            gamma_axial = f_handler['gamma_axial']
            
        gamma_feature_indx = f_handler['gamma_feature_indx'] 
    
    elif classifier_type == "GDA":
        mean_vector = f_handler['mean_vector'] 
        covariance_mat = f_handler['covariance_mat']    

    else:
        f_handler.close()
        print ( "ERROR! the classifier_type %s is not supported!\n" %classifier_type)
        #log_f.write(write_str)
        exit(0)
    
    #log_f.flush()
    f_handler.close()
    print  ("Succesfully read %s\n" %training_file)
    print "labelNumerics: ", labelNumerics
    #log_f.write(write_str)
    #log_f.write( "\n********************************\nlabelNumerics : "), log_f.write(str(labelNumerics)), log_f.write("\n\n")
    #print "\n********************************\nlabelNumerics : ", labelNumerics
    
    
    ##get info from test set	
    try:
        g, attributes= graph_analysis.input_graph(test_file, ["history", "vertex_offsets"])  	# open graph data and copy contents to memory
        #print ("Succefully read in %s\n" %test_file)	
    except:
        print("Error reading in %s\n" %test_file)
        
    feature_vector, target_vector, edge_w_indx, edge_neighboring_indx = graph2featureVect(g,featureNames)		#featureVect
    if pca=="Y" and len(axialfeatureidx)==0:
        [feature_vector, Y_bar, normalizing_stds, eigenvects, eigenvals] = PCA (feature_vector,feat_select_indx,Y_bar, normalizing_stds, eigenvects,eigenvals)
    elif pca=="Y" and len(axialfeatureidx)>0:	#leave axial direction out of PCA and add at the end
        tmp_indx = range(feature_vector.shape[1])
        for i in axialfeatureidx:
            tmp_indx.remove(i)
        tmp_feature_vector = feature_vector.take(tmp_indx,axis=1)
        [tmp_feature_vector, Y_bar, normalizing_stds, eigenvects, eigenvals] = PCA(tmp_feature_vector,feat_select_indx,Y_bar, normalizing_stds, eigenvects,eigenvals)
        for i in axialfeatureidx:
            tmp_feature_vector = numpy.hstack([tmp_feature_vector, array(feature_vector[:,i])])	
        feature_vector = tmp_feature_vector
        axialfeatureidx = [feature_vector.shape[1]-1,feature_vector.shape[1]-2 ,feature_vector.shape[1]-3 ]

    
    print  ("Succesfully read %s\n" %test_file)
    
    #f_handler = shelve.open(test_file,"r")
    ###labelNumerics = f_handler['labelNumerics'] 
    ###labelNumeric2Name = f_handler['labelNumeric2Name'] 
    ##featureNames = f_handler['featureNames'] 
    #feature_vector = f_handler['feature_vector'] 
    #target_vector = f_handler['target_vector'] 
    ###adjacency_matrix = f_handler['adjacency_matrix']
    #edge_w_indx= f_handler ['edge_w_indx'] 				##dictionary of {edge_indx : (e1,e2)}  where edge_indx is the i in feature_vector[i][:] and (e1,e2) is the edge tuple in the graph
    #edge_neighboring_indx = f_handler['edge_neighboring_indx']  		##dictionary of {edge_indx:[neighbour_edge1_index,neighbour_edge2_index,...]}
    #f_handler.close()
        

    
    
##################################################################
#### 1. initialize the labels
##################################################################
    RR_by_label = {}
    local_potentials = []		#each element = local_potentials[i] = potential of all possible labels for edge[i] sorted from high probability to low = {label:probability}
    print "classifier_type " , classifier_type
    
    if classifier_type=="NB":	
        for i in range(feature_vector.shape[0]):
            label_prob = {} #p(feature1/l)*p(feature2/l)*..*p(l) save as edge_property  estimated_label in a list for all labels
            for k in range(len(labelNumerics)):
                l = labelNumerics[k]
                p = numpy.log(class_prior[k,0])	# p(l)
                kent_flag = 0	#has it been calculated once for this label?
                for j in range(feature_vector.shape[1]):
                    if j not in axialfeatureidx:
                        tmp=univaraite_gaussian(feature_vector[i,j],mean_val[k,j],std_val[k,j])
                        #if numpy.isinf(tmp) or  numpy.isnan(tmp):
                        if tmp <=0:
                            #print ("ERROR!log is nan/inf!\n In calculation of log univaraite_gaussian for feature %s with true label %d for label %d: feature_vector[i,j] = %f,mean_val[k,j] = %f,std_val[k,j] = %f, gaussian= %f .\n" %(featureNames[j],target_vector[i,0],l,feature_vector[i,j],mean_val[k,j],std_val[k,j],univaraite_gaussian(feature_vector[i,j],mean_val[k,j],std_val[k,j])))
                            p+= numpy.log(machineEpsilon(float)/10e10)	#10e306
                        else:
                            p+= numpy.log(tmp)
                    elif j in axialfeatureidx and kent_flag==0:
                        kent_flag = 1	#it has been calculated for this label
                        p+= numpy.log(kent_distribution(feature_vector[i,axialfeatureidx],k_axial[k,0],beta_axial[k,0],ck_axial[k,0],gamma_axial[3*k:3*k+3,:]))
                #if not math.isnan(p):		
                label_prob[l] = p
                if math.isnan(p):
                    label_prob[l] = -1/(machineEpsilon(float))#numpy.log(machineEpsilon(float)/10e306)#
            label_prob_sorted = sorted(label_prob.iteritems(), key=lambda (k,v): (v,k), reverse=True)		#sort from high to low probability, put in tuple [(key,value),(key,value),...]		
            #label_prob_sorted = sorted(label_prob.items(), key=itemgetter(1), reverse=True)		#sort from high to low probability and return a list instead of original dict => sorted_list=[(label1,prob1),(label2,prob2),...]		
            local_potentials.append(label_prob_sorted)


    elif classifier_type == "NB2":
        for i in range(feature_vector.shape[0]):
            label_prob = {} #p(l/feature1)*p(l/feature2)*..*p(l) save as edge_property  estimated_label in a list for all labels
            for k in range(len(labelNumerics)):
                l = labelNumerics[k]
                p = numpy.log(class_prior[k,0])	# p(l)
                kent_flag = 0	#has it been calculated once for this label?
                for j in range(feature_vector.shape[1]):
                    if j not in gamma_feature_indx and j not in axialfeatureidx:
                        tmp=univaraite_gaussian(feature_vector[i,j],mean_val[k,j],std_val[k,j])
                        #if numpy.isinf(tmp) or  numpy.isnan(tmp):
                        if tmp<=0:
                            #print ("ERROR!log is nan/inf!\n In calculation of log univaraite_gaussian for feature %s with true label %d for label %d: feature_vector[i,j] = %f,mean_val[k,j] = %f,std_val[k,j] = %f, gaussian= %f .\n" %(featureNames[j],target_vector[i,0],l,feature_vector[i,j],mean_val[k,j],std_val[k,j],univaraite_gaussian(feature_vector[i,j],mean_val[k,j],std_val[k,j])))
                            p+= numpy.log(machineEpsilon(float)/10e10)	#10e306
                        else:
                            p+= numpy.log(tmp)
                    elif j not in axialfeatureidx:
                        p+= numpy.log(gamma_distribution(feature_vector[i,j],k_val[k,j],theta_val[k,j]))
                    elif j in axialfeatureidx and kent_flag==0:
                        kent_flag = 1	#it has been calculated for this label
                        p+= numpy.log(kent_distribution(feature_vector[i,axialfeatureidx],k_axial[k,0],beta_axial[k,0],ck_axial[k,0],gamma_axial[3*k:3*k+3,:]))
                                                
                #if not math.isnan(p):		
                label_prob[l] = p
                if math.isnan(p):
                    label_prob[l] = -1/(machineEpsilon(float))#numpy.log(machineEpsilon(float)/10e306)#
            label_prob_sorted = sorted(label_prob.iteritems(), key=lambda (k,v): (v,k), reverse=True)		#sort from high to low probability and return a list instead of original dict => sorted_list=[(label1,prob1),(label2,prob2),...]							 
            local_potentials.append(label_prob_sorted)

            
    elif classifier_type == "GDA":
        for i in range(feature_vector.shape[0]):
            label_prob = {} #p(l/feature1)*p(l/feature2)*..*p(l) save as edge_property  estimated_label in a list for all labels
            for k in range(len(labelNumerics)):
                l = labelNumerics[k]
                p = numpy.log(class_prior[k,0])	# p(l)
                X = feature_vector[i,:]
                tmp=multivariate_gaussian(X,mean_vector[k,:], covariance_mat[k*feature_vector.shape[1]:(k+1)*feature_vector.shape[1],:])
                if tmp<=0:
                    #print ("ERROR!log is nan/inf!\n In calculation of log univaraite_gaussian for feature %s with true label %d for label %d: feature_vector[i,j] = %f,mean_val[k,j] = %f,std_val[k,j] = %f, gaussian= %f .\n" %(featureNames[j],target_vector[i,0],l,feature_vector[i,j],mean_val[k,j],std_val[k,j],univaraite_gaussian(feature_vector[i,j],mean_val[k,j],std_val[k,j])))
                    p+= numpy.log(machineEpsilon(float)/10e10)		#10e306
                else:
                    p+= numpy.log(tmp)
                #p+= numpy.log(multivariate_gaussian(X,mean_vector[k,:], covariance_mat[k*feature_vector.shape[1]:(k+1)*feature_vector.shape[1],:]))
                #if not math.isnan(p):		
                label_prob[l] = p
                if math.isnan(p):
                    label_prob[l] = -1/(machineEpsilon(float))#numpy.log(machineEpsilon(float)/10e306)#
            label_prob_sorted = sorted(label_prob.iteritems(), key=lambda (k,v): (v,k), reverse=True)		#sort from high to low probability and return a list instead of original dict => sorted_list=[(label1,prob1),(label2,prob2),...]							 
            local_potentials.append(label_prob_sorted)
        
        
        

##################################################################
#### 1.1. initialization error calculation
##################################################################
    #print "Initial ERRORS:" 
    estimated_target_vector = [local_potentials[i][0][0] for i in range(len(local_potentials))]		#for each edge: the label with highest probability  #[i] for ith edge,[0] the highest label probability, [0] label name
    
    [error_num,labeled_num,total_volume,error_volume]=error_calculation (target_vector, estimated_target_vector,0,g,edge_w_indx,0)
    RR_by_label = error_calculation_by_label(RR_by_label,target_vector, estimated_target_vector, local_potentials,feature_vector, 0,edge_neighboring_indx,g,edge_w_indx)
    if options.verbose:
        confusion_matrix_calculation(target_vector, estimated_target_vector,labelNumerics)

###################################################################
##### 1.2. write labeled output
###################################################################
    if options.plot:	
        write_output_result(g,edge_w_indx,local_potentials,estimated_target_vector,target_vector, output_file[:-3]+"_initial.db")
    
####################################################################
###### 2. iterate to refine the labels based on their adjacency with gradient_descent
####################################################################
    max_global_iteration= options.iteration_num
    max_gradient_iteration= options.gradient_iteration_num
    if (max_gradient_iteration < 3):
        max_gradient_iteration = 3
    T = options.temprature
    eta = options.anneal_param
    
    if options.nonlin_annealing:
        T = 1e-10
        eta = 1.0- exp(-(0+2500.0)/800.0)		# eta = 1-exp(-(it+2500)/800)
        
    
    U_current = calc_globalE(estimated_target_vector,local_potentials,adjacency_matrix,edge_neighboring_indx,labelNumerics)	#with current labelling of all edges that come from local_potentials[:][0][0]
    if numpy.isnan(U_current) or numpy.isinf(U_current):
        print "Error! U_current undefined!\nAborted"
        exit(0)
    
    print "\nthe initial global energy = ", U_current
    
    
    init_T = T
    init_E = U_current	
    RR_initial = 100 - 100*float(error_num)/float(labeled_num+machineEpsilon(float))
    RRvol_initial = 100 -100*float(error_volume)/float(total_volume+machineEpsilon(float))
    initial_num = int(labeled_num - error_num)
    initial_vol = total_volume - error_volume

    U_optimum = calc_globalE(target_vector,local_potentials,adjacency_matrix,edge_neighboring_indx,labelNumerics)	#with current labelling of all edges that come from local_potentials[:][0][0]
    if numpy.isnan(U_optimum) or numpy.isinf(U_optimum):
        print "Error! U_optimum undefined!\nAborted"
        exit(0)
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
    
    global_cnt = 0	
    while_param = True
    while (while_param):
        random.shuffle(edgelist_visit)	#### browse nodes in random order
        U_current = calc_globalE(estimated_target_vector,local_potentials,adjacency_matrix,edge_neighboring_indx,labelNumerics)	#with current labelling of all edges that come from local_potentials[:][0][0]
        if numpy.isnan(U_current) or numpy.isinf(U_current):
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
                minErr_T = T
    
        #saving best energy configuration		
        if U_current<minE:
            minE = U_current
            minE_it = global_cnt
            minE_l=estimated_target_vector
            minE_err=100*float(error_num)/float(labeled_num+machineEpsilon(float))
            minE_errvol=100*float(error_volume)/float(total_volume+machineEpsilon(float))
            minE_T = T
            
        for i in edgelist_visit:			#for each node: X = feature_vector[i,:], current_l = local_potentials[i][0][0]
            l_current = estimated_target_vector[i]
           
            #### calculate pdf of transition to new labels for current node			
            delta_U = [] #energy variation of transition to other labels for current edge[i]
            for k in range(len(labelNumerics)):	#calculate energy variation for transition to l_new
                l_new = labelNumerics[k]
                current_delta_U = calc_deltaE(i,l_current,l_new,estimated_target_vector,local_potentials,adjacency_matrix,edge_neighboring_indx,labelNumerics)
                delta_U.append(current_delta_U)
            #print "delta_U:",delta_U

            #### options.optimization_nature=0	#gradient_descent optimization 
            ### for now let's only choose the label with minimum delta_U
            estimated_target_vector[i] = labelNumerics[delta_U.index(min(delta_U))]
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
    if options.plot:
        write_output_result(g,edge_w_indx,local_potentials,estimated_target_vector,target_vector, output_file[:-3]+"_GD.db")
    
    E_GD = U_current
    RR_GD = 100 - 100*float(error_num)/float(labeled_num+machineEpsilon(float))
    RRvol_GD = 100 -100*float(error_volume)/float(total_volume+machineEpsilon(float))
    GD_num = int(labeled_num-error_num)
    GD_vol = total_volume-error_volume
####################################################################
###### 3. iterate to refine the labels based on their adjacency with simulated annealing
####################################################################
    print "\n*************************** Start Simulated Annealing autolabelling ********************************\n"
    print "\nthe current global energy = ", U_current,", the optimum global energy with true labels = " ,U_optimum,", the initial temperature = " ,T
    global_cnt=0
    while_param = True
    anneal_change=0
    while (while_param):
        temperatures.append(T)
        ##print "\niteration:",global_cnt
        ##it_tic= time.time()
        pos_delE_cnt=0
        ##random.seed(422)
        random.shuffle(edgelist_visit)	#### browse nodes in random order
        U_current = calc_globalE(estimated_target_vector,local_potentials,adjacency_matrix,edge_neighboring_indx,labelNumerics)	#with current labelling of all edges that come from local_potentials[:][0][0]
        global_U.append(U_current)
        [error_num,labeled_num,total_volume,error_volume]= error_calculation(target_vector, estimated_target_vector,global_cnt,g,edge_w_indx,U_current,T)
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
                minErr_T = T
    
        #saving best energy configuration		
        if U_current<minE:
            minE = U_current
            minE_it = global_cnt
            minE_l=estimated_target_vector
            minE_err=100*float(error_num)/float(labeled_num+machineEpsilon(float))
            minE_err=100*float(error_volume)/float(total_volume+machineEpsilon(float))
            minE_T = T
            
        for i in edgelist_visit:			#for each node: X = feature_vector[i,:], current_l = local_potentials[i][0][0]
            l_current = estimated_target_vector[i]
            #print "\ni:",i, " current l=", l_current, " true l=",int(target_vector[i])  #, "\nlabels visited:"
            ##U_current = calc_globalE(estimated_target_vector,local_potentials,adjacency_matrix,edge_neighboring_indx,labelNumerics)	#with current labelling of all edges that come from local_potentials[:][0][0]
            
            #### calculate pdf of transition to new labels for current node			
            #pdf_ltransition =[]	#probability density function of transition to other labels for current edge[i]
            delta_U = [] #energy variation of transition to other labels for current edge[i]
            total_ltransition = 0
            for k in range(len(labelNumerics)):	#calculate energy variation for transition to l_new
                l_new = labelNumerics[k]
                ##print l_new,",",
                ##tmp_estimated_target_vector = copy.deepcopy(estimated_target_vector)
                ##tmp_estimated_target_vector[i] = l_new
                ##U_new = calc_globalE(tmp_estimated_target_vector,local_potentials,adjacency_matrix,edge_neighboring_indx,labelNumerics)	#with current labelling of all edges that come from local_potentials[:][0][0] and replacing local_potentials[i][0][0] to l_new
                ##delta_U.append(U_new-U_current)
                current_delta_U = calc_deltaE(i,l_current,l_new,estimated_target_vector,local_potentials,adjacency_matrix,edge_neighboring_indx,labelNumerics)
                delta_U.append(current_delta_U)
            #print "delta_U:",delta_U
            
            
            
            
            #### stochastic optimization
            #### normalized pdf
            pdf_ltransition=[]
            #### normalize delta_U to be in range [-1...1]: (so map delta_U =[0...max] -> [0 ..1] and negative delta_U will be interpolated linearly
            ##normal_delta_U = [(2.0*(x - min(delta_U))/(max(delta_U)-min(delta_U))-1) for x in delta_U]
            normal_delta_U = [(x /(max(abs(max(delta_U)),abs(min(delta_U)))+machineEpsilon(float))) for x in delta_U]			#### ??? pick a constant value for max(delta_U) for all edges and all datasets?
            #print "normal_delta_U:",normal_delta_U 
            for k in range(len(labelNumerics)):	#calculate energy variation for transition to l_new
                if exp(-normal_delta_U[k]/(T+machineEpsilon(float)))==inf:
                    exp_val= 1/machineEpsilon(float)
                else:
                    exp_val=exp(-normal_delta_U[k]/(T+machineEpsilon(float)))
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
            #print "relabeled to ",labelNumerics[ltransition_indx]
            if normal_delta_U[ltransition_indx]>0:
                pos_delE_cnt+=1
        
        percent_pos_delE.append(float(pos_delE_cnt)/float(len(labelNumerics)))
    
        #### at the end of each global iteration
        #print ("Temprature is %f " %T)
        if options.animation and global_cnt%options.decimate==0:
            write_output_result(g,edge_w_indx,local_potentials,estimated_target_vector,target_vector, os.path.dirname(os.path.abspath(output_file))+"/it_"+str(global_cnt)+".db")

        global_cnt = global_cnt + 1
        T = eta*T
        if options.nonlin_annealing:
            eta = 1.0- exp(-(global_cnt+2500.0)/800.0)		# eta = 1-exp(-(it+2500)/800)
            
        if (global_cnt > max_global_iteration):
            while_param = False
        elif (global_cnt > 50):	#run at least 50 iterations	
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
                eta = 0.9+eta/10.0		#eta+(1-eta)*0.9		if eta<1 => eta/10<0.1 => 0.9+eta/10.0<1
            

                
    print "\n\nLabelling DONE!"			
    
    final_T = T

####################################################################
###### 3.1 write labeled output
####################################################################
    if options.plot:
        write_output_result(g,edge_w_indx,local_potentials,estimated_target_vector,target_vector, output_file[:-3]+"_S.db")

    
    #if options.optimization_nature == 1:
        ##saving best error configuration
        #h=copy.deepcopy(g)
        #for i in range(len(g.edge_list())):
            #e=edge_w_indx[i]
            #h.set_edge_property(e,'estimated_label',[list(j) for j in local_potentials[i]])	#local_potentials is a list for all edges, for each edge = local_potentials[i] is list of tuples=[(label,highest_prob),(l,prob),...,(l,lowest_prob)]
            #if (target_vector[i]>0 and (target_vector[i] == minErr_l[i])):
                #h.set_edge_property(e,'error_label',1)	#correct labelling	
            #elif (target_vector[i]>0 and (target_vector[i] != minErr_l[i])):
                #h.set_edge_property(e,'error_label',2)	#error labelling	
            #else:
                #h.set_edge_property(e,'error_label',0)
        #graph_analysis.output_graph(output_file[:-3]+"_minErr.db", h, history, attributes)   
        
        #print ("Succefully wrote the %s\n" %output_file[:-3]+"_minErr.db")
        
        #cmd=("\npython /projects/souris/sghanavati/src/scripts/graph2cylinder.py %s %s --use_estimated_label --clobber " %(output_file[:-3]+"_minErr.db",output_file[:-3]+"_minErr_cyl.db"))	#python /micehome/jgsled/bin/
        #print(cmd)
        #os.system(cmd)

        #cmd=("\npython /projects/souris/sghanavati/src/scripts/graph2cylinder.py %s %s --use_error_label --clobber " %(output_file[:-3]+"_minErr.db",output_file[:-3]+"_minErr_error_cyl.db"))	#python /micehome/jgsled/bin/
        #print(cmd)
        #os.system(cmd)
        
        #cmnd=("rm %s" %(output_file[:-3]+"_minErr*_cyl.db"))
        #print(cmd)
        #os.system(cmd)	
        
        
        ##saving best Energy configuration
        #h=copy.deepcopy(g)
        #for i in range(len(g.edge_list())):
            #e=edge_w_indx[i]
            #h.set_edge_property(e,'estimated_label',[list(j) for j in local_potentials[i]])	#local_potentials is a list for all edges, for each edge = local_potentials[i] is list of tuples=[(label,highest_prob),(l,prob),...,(l,lowest_prob)]
            #if (target_vector[i]>0 and (target_vector[i] == minE_l[i])):
                #h.set_edge_property(e,'error_label',1)	#correct labelling	
            #elif (target_vector[i]>0 and (target_vector[i] != minE_l[i])):
                #h.set_edge_property(e,'error_label',2)	#error labelling	
            #else:
                #h.set_edge_property(e,'error_label',0)
        #graph_analysis.output_graph(output_file[:-3]+"_minE.db", h, history, attributes)   
        
        #print ("Succefully wrote the %s\n" %output_file[:-3]+"_minE.db")
        
        #cmd=("\npython /projects/souris/sghanavati/src/scripts/graph2cylinder.py %s %s --use_estimated_label --clobber " %(output_file[:-3]+"_minE.db",output_file[:-3]+"_minE_cyl.db"))	#python /micehome/jgsled/bin/
        #print(cmd)
        #os.system(cmd)

        #cmd=("\npython /projects/souris/sghanavati/src/scripts/graph2cylinder.py %s %s --use_error_label --clobber " %(output_file[:-3]+"_minE.db",output_file[:-3]+"_minE_error_cyl.db"))	#python /micehome/jgsled/bin/
        #print(cmd)
        #os.system(cmd)
        
        #cmnd=("rm %s" %(output_file[:-3]+"_minE*_cyl.db"))
        #print(cmd)
        #os.system(cmd)	
        

####################################################################
###### 4. Error Calculations!
####################################################################
    [error_num,labeled_num,total_volume,error_volume]=error_calculation (target_vector, estimated_target_vector,global_cnt,g,edge_w_indx,U_current,T)
    RR_by_label = error_calculation_by_label(RR_by_label,target_vector, estimated_target_vector, local_potentials,feature_vector, global_cnt,edge_neighboring_indx,g,edge_w_indx)
    confusion_matrix_calculation(target_vector, estimated_target_vector,labelNumerics)

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
    
    txt = ("initial RR "+str(RR_initial)+"%\tvolumetricRR " + str(RRvol_initial)+"\n") 
    f.write(txt)
    txt = ("GD RR "+str(RR_GD)+"%\tvolumetricRR " + str(RRvol_GD)+"\n") 
    f.write(txt)
    txt = ("SR RR "+str(100-100*float(error_num)/float(labeled_num))+"%\tvolumetricRR " + str(100-100*float(error_volume)/float(total_volume))+"\n\n") 
    f.write(txt)
    
    txt = ("highest RR "+str(100-minErr)+"% (volRR "+str(100-minErrvol)+"%)\titeration " + str(minErr_it)+"\ttemperature "+str(minErr_T)+"\tenergy "+str(minErr_energy)+"\n") 
    f.write(txt)
    txt = ("lowest energy RR "+str(100-minE_err)+"% (volRR "+str(100-minE_errvol)+"%)\titeration " + str(minE_it)+"\ttemperature "+str(minE_T)+"\tenergy "+str(minE)+"\n\n") 
    f.write(txt)
    
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
    
