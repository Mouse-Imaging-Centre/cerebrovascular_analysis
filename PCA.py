#!/usr/bin/env python
# -*- coding: utf-8 -*-


#### option to perform PCA on all features_vectors and make them orthogonal


#
#
#  Created May 18, 2015 
#  Modified Oct 21, 2015  PCA feature reduction use np.dot for matrix product  
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
#from pylab import *        #needed for command find
import pylab        #needed for command find
import scipy
import random

def PCA (input_matrix,CDF_thresh=1.0,Y_bar=-inf, normalizing_stds = -inf, eigenvects=-inf,eigenvals=-inf,M_ind=-1):
    #### 1. centering the data to mean=0
    if type(Y_bar)==float:
        Y_bar = np.mean(input_matrix,0)   #mean of each column of input_matrix = mean of each feature data points
    input_matrix_cntr = input_matrix - tile(Y_bar,(input_matrix.shape[0],1))
    #### 2. normalizing the data to std=1
    if type(normalizing_stds)==float:
        normalizing_stds = np.std(input_matrix_cntr,0)
    input_matrix_cntrN = input_matrix_cntr / tile(normalizing_stds,(input_matrix_cntr.shape[0],1))

    #""" SVD:http://docs.scipy.org/doc/numpy/reference/generated/numpy.linalg.svd.html
        # X.H is conjugate transpose matrix of X
        # if matrix X is orthogonal: 1.invertible X^T (transpose) = X^-1 (inverse) 2.unitary X^-1 (inv) =X* (X.H conjugate transpose) 3.normal X.H X = X X.H
        # input_matrix_cntrN = u * np.diag(s) * v.H
        # U, s, V = np.linalg.svd(input_matrix_cntrN, full_matrices=False) : U=u, V=v.H
        # The rows of V are the eigenvectors of A.H A which is of size of features => it is what we are looking for
        #SO EIGENVECTORS ARE COLUMNS OF V.transpose
        # The columns of U are the eigenvectors of A A.H which is of size of datapoints
        # For row i in V and column i in U, the corresponding eigenvalue is s[i]**2.
        # Decomposition(data mapping) is T = A eigenvects = A V.transpose  http://en.wikipedia.org/wiki/Principal_component_analysis#Further_components and http://en.wikipedia.org/wiki/Principal_component_analysis#Dimensionality_reduction
    #"""
    if type(eigenvects)==float or type(eigenvals)==float:
        #### test: cov_matrix symmetric positive semidefinite (PSD) matrix => has a nonnegative eigenvalue
        cov_matrix = np.cov(input_matrix_cntrN,None,0)    #cov(m, y=None, rowvar=1, bias=0, ddof=None) rowvar=0:each column represents a variable, while the rows contain observations
        #print type(cov_matrix), "\ncov_matrix:", cov_matrix
       #'''
        #### ERROR: This is wrong because we shouldn't calculate the SVD of input_matrix_cntrN but the cov_matrix!!!
        #U, s, V = np.linalg.svd(input_matrix_cntrN, full_matrices=False)   #input_matrix_cntrN= U * np.diag(s) * V
        ##The rows of v are the eigenvectors of a.H a. The columns of u are the eigenvectors of a a.H.
        ##if input_matrix_cntrN is Nxm # U is descending ordered eigenvectors of Nxm
        #sqrt_eigenvals = s    # s is descending ordered list of eigenvals of len m
        #eigenvals = [e*e for e in s]      
        #eigenvects = V.H #which is equal to V.transpose = inv(V)  
        ##eigenvects = U #U contains vectors of singular components  
        #eigenvalsmat = np.diag(s)*np.diag(s)    #is diagonal matrix of mxm where eigenvals are on the main diagonal and the rest is 0
        #'''
        
        D,V = np.linalg.eig(cov_matrix)     # a: array_like shape (M, M), D:eigenvalues ndarray shape (M,) NOT ordered, V: eignevectors ndarray shape (M, M) The normalized (unit “length”) eigenvectors, column V[:,i] is the eigenvector corresponding to the eigenvalue D[i]          
        #i = np.argsort(D)  #Returns the indices that would sort an array, ascending
        i = np.argsort(np.array(D)) #D[::-1] #Returns the indices that would sort an array, descending
        inv_i = i[::-1]
        eigenvals = D[inv_i]
        eigenvects = V[:,inv_i]
        #print type(eigenvals)
        #print type(eigenvects)
        #print "V:", V 
        #print "sorted eigen values are (should be >0) : " , eigenvals      #### test: eigvals should be >0 => orthogonal_input_matrix: real
        #print "sorted eigen vects are (should be >0) : " , eigenvects      #### test: eigvals should be >0 => orthogonal_input_matrix: real

    if M_ind==-1:
        cumulative_eigval = np.cumsum(eigenvals)
        cumulative_eigval = cumulative_eigval/cumulative_eigval[-1]*100.0   # % CDF
        M_ind = min(pylab.find(cumulative_eigval >= CDF_thresh))
    W = eigenvects[:,0:M_ind]
    orthogonal_input_matrix = np.dot(np.matrix(input_matrix_cntrN) , W)   #if W is eigenvectors T = X W maps a data vector x(i) from an original space of p variables to a new space of p variables which are uncorrelated over the dataset.

    print "\n\nM_ind:",M_ind, " out of ", np.matrix(eigenvects).shape[0], " features and feature vect shape of ", np.matrix(input_matrix_cntrN).shape
    print "significant eigenvals:", eigenvals[0:M_ind]
    #print (np.matrix(W).H).shape
    #print np.matrix(input_matrix_cntrN).shape
    #print np.matrix(np.diag(1.0/np.sqrt(eigenvals))).shape
    #print ((np.matrix(input_matrix_cntrN) * np.matrix(np.diag(1.0/np.sqrt(eigenvals)))).transpose()).shape

    #### orthogonal_input_matrix = input_matrix_cntr * V * D^(-.5)??   
    ##based on karhunen-Leove theorem and data projection from http://en.wikipedia.org/wiki/Principal_component_analysis:
    ## projected data = conjugate_transpose(eigenvects)* (feature_vects/Diag(sqrt(eigenvals)))
    ## a.H conjugate transpose
    #orthogonal_input_matrix = np.matrix(W).H * (np.matrix(input_matrix_cntrN) * np.matrix(np.diag(1.0/np.sqrt(eigenvals)))).transpose()
    #orthogonal_input_matrix = orthogonal_input_matrix.transpose()
    #### test: correctness of PCA: cov(orthogonal_input_matrix) almost = I
    return [orthogonal_input_matrix, np.array(Y_bar), np.array(normalizing_stds), np.matrix(eigenvects), np.array(eigenvals), M_ind]


