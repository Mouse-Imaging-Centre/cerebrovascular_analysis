# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 16:28:04 2015

@author: sghanavati
"""

#!/usr/bin/env python

# Get manuallabel_graph.db  and corresponding CT.mnc
# initialize the CT.mnc to 0 everywhere
# for each vertex in the graph give the label to the mnc voxels in a sphere around the graph vertex of vertex radius 
# it will save the labeled CT.mnc
#  




from vessel_tracking import graph_analysis, distribution_analysis
#from graph_labeling import *
from optparse import OptionParser, Option, OptionValueError
from sys import argv
import os, commands, time, string
#import py_minc 
import pyminc.volumes.factory as f
import pyminc.volumes.volumes as v
import numpy as np
import scipy
import copy
import operator

#---------------------------------------------------------------------------------
#
#def mean_list(numberList):
    ##floatNums = [float(x) for x in numberList]
    #return sum(numberList) / float(len(numberList))	
    
#def create_dist_voxels(N):
#    #if x,y,z goes from [-N..N], create a list of all possible coordinates and sort them by distance:
#    indxs =[]
#    for i in range(-N,N+1):
#        for j in range(-N,N+1):
#            for k in range(-N,N+1):
#                dist=calc_dist(i,j,k)
#                indxs.append([i,j,k,dist])
#    indx_l = sorted(indxs,key=lambda indxs:indxs[3])
#    return indx_l
    
def calc_dist(x,y,z):
    d = np.sqrt(x*x+y*y+z*z)
    return d

def mnc_labeling(g,CT_file,output_file, multiplier=1.0):
    #labeled_ct = f.volumeFromFile(CT_file)
    labeled_ct = f.volumeLikeFile(CT_file,output_file)
    #labeled_ct.openFile()
    size = labeled_ct.getSizes()
    print "CT size is ", size
    print "first and last voxel intensities are:" ,labeled_ct.data[0, 0, 0] 
    print labeled_ct.data[ size[0]-1, size[1]-1, size[2]-1]
    spacing = labeled_ct.getSeparations()
    start = labeled_ct.getStarts()
    #set the volume to 0
    for i in range(size[0]):
        for j in range(size[1]):
            for k in range(size[2]):
                labeled_ct.data[ i, j, k] = 0
    #set the labels            
    for e in g.edge_list():
        l = g.edge_property(e, 'label')
        vlist = [e[0]]+[e[1]]+  g.edge_property(e, 'intermediaries')  
        for v in vlist:
            w = g.vertices[v].centre
            r = multiplier*g.vertices[v].radius
            wtovfloat = labeled_ct.convertWorldToVoxel(w)
            wtov = np.array([int(wtovfloat[0]),int(wtovfloat[1]),int(wtovfloat[2])])
            if wtov[0]<size[0] and wtov[1]<size[1] and wtov[2]<size[2] and wtov[0]>-1 and wtov[1]>-1 and wtov[2]>-1: 
                labeled_ct.data[int(wtov[0]),int(wtov[1]),int(wtov[2])] = l
            for i in range(max(0,int(wtov[0]-(r/spacing[0]))), min(int(wtov[0]+(r/spacing[0])),size[0])):
                for j in range(max(0,int(wtov[1]-(r/spacing[1]))), min(int(wtov[1]+(r/spacing[1])),size[1])):
                    for k in range(max(0,int(wtov[2]-(r/spacing[2]))), min(int(wtov[2]+(r/spacing[2])),size[2])):
                        currw = labeled_ct.convertWorldToVoxel(np.array([i, j , k]))
                        if calc_dist(w[0]-currw[0],w[1]-currw[1],w[2]-currw[2]) < r:
                            labeled_ct.data[int(i), int(j), int(k)] = l
    
    labeled_ct.writeFile()
    labeled_ct.closeVolume()
    return 0
    
    




program_name = 'graph_label2mnc.py'


if __name__ == '__main__':
        
    usage = "Usage: "+program_name+" [options] CT.mnc labelled_graph.db [output_labeled.mnc] \n"+\
            "   or  "+program_name+" --help";
            
    parser = OptionParser(usage)
    
    parser.add_option("--clobber", action="store_true", dest="clobber",
                    default=0, help="overwrite output file")
    parser.add_option("--multiplier", type="float", dest="multiplier",
                    default=10.0, help="multiple the radius by this value")
    #parser.add_option("--remove_interms",action="store_true",default=0, dest="remove_interms",help="remove the intermediate mnc files that were created for calculating seeds in iterations")
        
        
    tic = time.time()
    #ct_labeling = ct_labeling()		#initialize the class ct_labeling, make an instance called ct_labeling
                
    options, args = parser.parse_args()
    print '\n\n>>> %s [%s]: %s\n' % (time.ctime(time.time()), os.getcwd(), string.join(argv))
    print "\n\n",options, "\n\n", args , "\n\n"
   

        
    if len(args) == 2:
        CT_file, input_file = args
        output_file = CT_file[:-4]+"_graphlabel.mnc"
    elif len(args) == 3:
        CT_file, input_file, output_file = args
    else:
        parser.error("incorrect number of arguments")
        
    
    if not options.clobber and os.path.exists(output_file):
        raise SystemExit, \
            "The --clobber option is needed to overwrite an existing file."
            
        
    # reads in the graph to be labeled
    try:
        g, attributes = graph_analysis.input_graph(input_file, ["history", "vertex_offsets"])
        print ("Succefully read in the %s\n" %input_file)	
    except:
        print("Error reading in the %s\n" %input_file)
        
    mnc_labeling(g,CT_file,output_file, options.multiplier)

#    # read in the labeled CT : used to make labling training dataset
#    #labeled_ct=py_minc.ArrayVolume(labeled_CT_file,dim_names=py_minc.FILE_ORDER_DIMENSION_NAMES) #which is in [x,y,z]
#    labeled_ct = f.volumeFromFile(CT_file)	# without this always [z,y,x], with this :dim_names=py_minc.FILE_ORDER_DIMENSION_NAMES order as the image	
#    
#    labeled_ct = mnc_labeling(g,labeled_ct)
#    
#    labeled_ct.writeFile()
#    labeled_ct.closeVolume()
    print ("Succefully wrote the %s\n" %output_file)
    

    toc = time.time()
    print ("total elapsed time: %f" %(toc - tic))
