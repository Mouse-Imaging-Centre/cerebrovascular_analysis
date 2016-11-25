#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Collection of libraries for analysis and comparison of segmented, labelled vasculature
#
#
#
#
#  Created Nov 11, 2013
#  modified Nov 11, 2013
#  Last modified April 17, 2015  added a masked version of the perfusion_region function 
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
import numpy as np
import numpy.matlib as Matlib
from optparse import OptionParser, Option, OptionValueError
from minc_util.progress import progress_report
#import matplotlib.pyplot as plt
from pylab import *		#needed for command find
import random
import copy
import py_minc 
import pyminc.volumes.factory as pymincf
import pyminc.volumes.volumes as pymincv
import operator  #for sorting dictionary by key or value
#from operator import itemgetter

#---------------------------------------------------------------------------------
def distance(x,y):
    d = sqrt((x[0]-y[0])*(x[0]-y[0]) + (x[1]-y[1])*(x[1]-y[1])+(x[2]-y[2])*(x[2]-y[2]))
    return d

def estimate_edge_distlengths(g, property='dist_length'):
    """Compute the length of each vessel segment using the specified length function.

    g                    : simplified skeletal graph with intermediaries
    length_fcn(g, edge)  : vessel length function
    property             : name of length property to set for each edge

    Each edge of h is given the specified property.
    """
    for e in g.edge_list():
        dist_length = distance(g.vertices[e[0]].centre,g.vertices[e[1]].centre)
        g.set_edge_property(e, property, float(dist_length))

def vessel_avg_diameter(g, edge):
    #"estimate vessel diameter by averaging the diameters at intermediate points"
    diameters = [g.vertices[i].radius*2 for i in \
                 [edge[0]] + g.edge_property(edge, 'intermediaries') \
                  + [edge[1]]]
    return mean(diameters)


def estimate_edge_nooutlier_diameters(g, property='nooutlier_diameter'):
    """Compute the average diameter of each vessel segment using the specified diameter function.

    g                      : simplified skeletal graph with intermediaries
    diameter_fcn(g, edge)  : vessel diameter function
    property               : name of diameter property to set for each edge

    Each edge of h is given the specified property.
    """
    
    for e in g.edge_list():
        edge_vs = [e[0]] + g.edge_property(e, 'intermediaries') + [e[1]]
        diameter = mean([g.vertices[v].radius*2 for v in edge_vs])             #vessel_avg_diameter(g, edge)
        s = std([g.vertices[v].radius*2 for v in edge_vs])
        nooutlier_vs = [v for v in edge_vs if (diameter - 2*s < g.vertices[v].radius*2 < diameter + 2*s)]
        nooutlier_diameter = mean([g.vertices[v].radius*2 for v in nooutlier_vs])
        g.set_edge_property(e, property, float(nooutlier_diameter))

    
def find_neighbor_edges_asymmetrical(g,e, c=0):
    #find edge neighbouring indeces
    #neighbor_indx=[]
    neighbor_edges=[]
    
    e1_neighbours= g.vertices[e[c]].edges
    for v in e1_neighbours:
        neighbor_edge=tuple((min(e[c],v),max(e[c],v)))
        if neighbor_edge!=e:
            neighbor_edges.append(neighbor_edge)
            
    return neighbor_edges

    
def find_neighbor_edges(g,e):
    #find edge neighbouring indeces
    #neighbor_indx=[]
    neighbor_edges=[]
    
    e1_neighbours= g.vertices[e[0]].edges
    for v in e1_neighbours:
        if v!=e[1]:
            neighbor_edge=tuple((min(e[0],v),max(e[0],v)))
            neighbor_edges.append(neighbor_edge)
    e2_neighbours= g.vertices[e[1]].edges
    for v in e2_neighbours:
        if v!=e[0]:
            neighbor_edge=tuple((min(e[1],v),max(e[1],v)))
            neighbor_edges.append(neighbor_edge)
            
    return neighbor_edges
    
def find_edge_indx (g,e):
    edge = tuple([min(e),max(e)])
    return g.edge_list().index(edge)
    #for i in range(len(g.edge_list())):
        #if edge==g.edge_list()[i]:
            #return i

def difference_update(original_list, list2remove):
    for i in list2remove:
        if i in original_list:
            original_list.remove(i)
    return original_list
    

def connected_components(g, edgelist_copy = []):
    if len(edgelist_copy) == 0:
        edgelist_copy = [e for e in g.edge_list()]
    edgelist = copy.deepcopy(edgelist_copy)
    connected_components = []
    
    while len(edgelist)>0:
        e = edgelist.pop()
        # This set will contain the next group of nodes connected to each other.
        group = [e]

        # Build a queue with this node in it.
        queue = [e]

        # Iterate the queue.
        # When it's empty, we finished visiting a group of connected nodes.
        while queue:

            # Consume the next item from the queue.
            e = queue.pop(0)

            # Fetch the neighbors.
            neighbors = find_neighbor_edges(g,e)

            # Remove the edges we already visited (edges in group) from neighbors.
            difference_update(neighbors,group)

            # Remove the neighbors from the edgelist.
            difference_update(edgelist,neighbors)

            # Add them to the group of connected nodes.
            group = group+neighbors

            # Add them to the queue, so we visit them in the next iterations.
            queue = queue+neighbors

        # Add the group to the list of groups.
        connected_components.append(group)

    return connected_components
    
def number_of_connected_components(g, edgelist_copy = []):
    return len(connected_components(g, edgelist_copy ))
    
def size_of_connected_components(g, edgelist_copy = []):
    cc= connected_components(g, edgelist_copy )
    size_list = [len(i) for i in cc]    
    return size_list
    
    
    
def retain_loops(g,output_file):
    #loop is defined: if an edge is broken and the connected graph breaks into 2 components, that edge is not part of a loop
    while graph_analysis.prune_leaves(g) > 0:
        pass

    N = number_of_connected_components(g)
    edge2remove = []	#edge not part of a loop
    ne=0
    h = copy.deepcopy(g)
    for e in g.edge_list():
        properties = g.edge_properties(e)
        h.remove_edge(e)
        if number_of_connected_components(h)>N:
            edge2remove.append(e)
        if number_of_connected_components(h)!=N:	
            ne+=1
        h.add_edge(e, properties)	
    print ne, " ", len(edge2remove)
    h = copy.deepcopy(g)
    for e in edge2remove:
        h.remove_edge(e)
    graph_analysis.output_graph(output_file, h)	
    return h
        
    

    
    
def number_of_loops(g):
    #loop is defined: if an edge is broken and the connected graph breaks into 2 components, that edge is not part of a loop
    while graph_analysis.prune_leaves(g) > 0:
        pass
    N = number_of_connected_components(g)
    edge2remove = []	#edge not part of a loop
    ne=0
    for e in g.edge_list():
        h = copy.deepcopy(g)
        h.remove_edge(e)
        if number_of_connected_components(h)>N:
            edge2remove.append(e)
        if number_of_connected_components(h)!=N:	
            ne+=1
    #print ne, " ", len(edge2remove)
    h = copy.deepcopy(g)
    for e in edge2remove:
        h.remove_edge(e)
    return ne



def perfusion_region_noremap(g,mr_atlas_file,labels,perfusion_file,Perfusionradius=1e10):
# 	mr_atlas_file = "/projects/souris/sghanavati/data/c57_male_female_mri_atlas/c57_brain_segmentation_mnc2_merged_reg2vascular.mnc"
#	mr_atlas_file = "/projects/souris/sghanavati/data/c57_male_female_mri_atlas/c57_brain_segmentation_mnc2_merged_reg2vascular_downsample_merge.mnc"
    endpoint = []
    vpoint=[]
    for e in g.edge_list():
        #if len(find_neighbor_edges(g,e))==2 and g.edge_property(e,'label') in labels:
        if g.edge_property(e,'label') in labels:
            if len(g.vertices[e[0]].edges)==1:
                vpoint.append(e[0])
                endpoint.append(e)
            elif len(g.vertices[e[1]].edges)==1:
                vpoint.append(e[1])
                endpoint.append(e)
    #### mri_distance feature:
    ## nc_data_type=py_minc.NC_BYTE is writing a byte (floating point) but we want the voxel intensities to be int labels
    ## nc_data_type=py_minc.NC_SHORT is writing a short integer data and we want the voxel intensities to be int labels
    #mri_atlas = py_minc.ArrayVolume(mr_atlas_file, nc_data_type=py_minc.NC_BYTE )# without this always [z,y,x], with ,dim_names=py_minc.FILE_ORDER_DIMENSION_NAMES order as the image	
    #perfusion_pattern = py_minc.ArrayVolume(mr_atlas_file, nc_data_type=py_minc.NC_BYTE )# without this always [z,y,x], with ,dim_names=py_minc.FILE_ORDER_DIMENSION_NAMES order as the image	
    #sizes = mri_atlas.get_sizes()
    mri_atlas = pymincf.volumeFromFile(mr_atlas_file, dtype='ubyte')   #ushort
    mri_atlas.openFile()
    sizes = mri_atlas.getSizes()
    spacing = mri_atlas.getSeparations()
    perfusion_pattern = pymincf.volumeFromInstance(mri_atlas, perfusion_file, dtype="ubyte") #ushort


    #make all voxels 0, give label for vpoints
    for z_vox in range(sizes[0]):
        for y_vox in range(sizes[1]):
            for x_vox in range(sizes[2]):
                #perfusion_pattern.array[z_vox,y_vox,x_vox]=0
                perfusion_pattern.data[z_vox,y_vox,x_vox]=0
    for i in range(len(vpoint)):
        v=vpoint[i]
        world = g.vertices[v].centre
        #voxel = floor(mri_atlas.convert_world_to_voxel(array([world[0],world[1],world[2]])))
        voxel = floor(mri_atlas.convertWorldToVoxel(np.array([world[0],world[1],world[2]])))
        if voxel[0] in range(sizes[0]) and voxel[1] in range(sizes[1]) and voxel[2] in range(sizes[2]) :
            perfusion_pattern.data[int(voxel[0]),int(voxel[1]),int(voxel[2])] = g.edge_property(endpoint[i],'label')
            #perfusion_pattern.array[int(voxel[0]),int(voxel[1]),int(voxel[2])] = g.edge_property(endpoint[i],'label')
        
                
    for z_vox in range(sizes[0]):
        for y_vox in range(sizes[1]):
            for x_vox in range(sizes[2]):
                #world=mri_atlas.convert_voxel_to_world(array([z_vox,y_vox,x_vox]))
                world = mri_atlas.convertVoxelToWorld(np.array([z_vox,y_vox,x_vox]))

                min_dist=1e10
                min_label = 0
                for i in range(len(vpoint)):
                    v=vpoint[i]
                    if distance(world,g.vertices[v].centre)<Perfusionradius and distance(world,g.vertices[v].centre)<min_dist:
                        min_dist = distance(world,g.vertices[v].centre)
                        min_label = g.edge_property(endpoint[i],'label')
                #perfusion_pattern.array[z_vox,y_vox,x_vox]=	int(min_label)
                perfusion_pattern.data[z_vox,y_vox,x_vox]= int(min_label)
                
    #perfusion_pattern.output(perfusion_file)			
    perfusion_pattern.writeFile()
    perfusion_pattern.closeVolume()
    mri_atlas.closeVolume()

   
#### masked version of perfusion_region   
def perfusion_region(g, mr_atlas_file, labels, mask_file, label_remap, perfusion_file, Perfusionradius=1e10):
#   mr_atlas_file = "/projects/souris/sghanavati/data/c57_male_female_mri_atlas/c57_brain_segmentation_mnc2_merged_reg2vascular.mnc"
#   mr_atlas_file = "/projects/souris/sghanavati/data/c57_male_female_mri_atlas/c57_brain_segmentation_mnc2_merged_reg2vascular_downsample_merge.mnc"

    if label_remap:
        lpairs = label_remap.split(',')
        labelpairs = []
        for pair in lpairs:
            plabels = pair.split(':')
            plabels = [int(l) for l in plabels]
            labelpairs.append(plabels)
        tomap_l = [int(l[0]) for l in labelpairs]
        print labelpairs

        for e in g.edge_list():
            if g.edge_property(e,'label') in tomap_l:
                g.set_edge_property(e,'label',int(labelpairs[tomap_l.index(g.edge_property(e,'label'))][1]))


    endpoint = []
    vpoint=[]
    for e in g.edge_list():
        #if len(find_neighbor_edges(g,e))==2 and g.edge_property(e,'label') in labels:
        if g.edge_property(e,'label') in labels:
            if len(g.vertices[e[0]].edges)==1:
                vpoint.append(e[0])
                endpoint.append(e)
            elif len(g.vertices[e[1]].edges)==1:
                vpoint.append(e[1])
                endpoint.append(e)
                
    if not len(vpoint)==len(endpoint):
        print "ERROR: The length of endpoints and leaf edges should match but they don't!"
        exit(0)
    
    #### mri_distance feature:
    ## nc_data_type=py_minc.NC_BYTE is writing a byte (floating point) but we want the voxel intensities to be int labels
    ## nc_data_type=py_minc.NC_SHORT is writing a short integer data and we want the voxel intensities to be int labels
    #mri_atlas = py_minc.ArrayVolume(mr_atlas_file, nc_data_type=py_minc.NC_BYTE )# without this always [z,y,x], with ,dim_names=py_minc.FILE_ORDER_DIMENSION_NAMES order as the image    
    #perfusion_pattern = py_minc.ArrayVolume(mr_atlas_file, nc_data_type=py_minc.NC_BYTE )# without this always [z,y,x], with ,dim_names=py_minc.FILE_ORDER_DIMENSION_NAMES order as the image    
    #mask_pattern = py_minc.ArrayVolume(mask_file, nc_data_type=py_minc.NC_BYTE )# without this always [z,y,x], with ,dim_names=py_minc.FILE_ORDER_DIMENSION_NAMES order as the image    
    ##make all voxels 0, give label for vpoints
    #sizes = mri_atlas.get_sizes()
    #mask_sizes = mask_pattern.get_sizes()

    mri_atlas = pymincf.volumeFromFile(mr_atlas_file, dtype='ubyte')   #ushort
    mri_atlas.openFile()
    sizes = mri_atlas.getSizes()
    spacing = mri_atlas.getSeparations()
    perfusion_pattern = pymincf.volumeFromInstance(mri_atlas, perfusion_file, dtype="ubyte")
    mask_pattern = pymincf.volumeFromFile(mask_file, dtype='ubyte')
    mask_pattern.openFile()
    mask_sizes = mask_pattern.getSizes()
    
    if not numpy.array_equal(sizes,mask_sizes):
        print "ERROR: The size of the mri_atlas and the mask should match but they don't!"
        exit(0)
        
    for z_vox in range(sizes[0]):
        for y_vox in range(sizes[1]):
            for x_vox in range(sizes[2]):
                #if mask_pattern.array[z_vox,y_vox,x_vox] < 0.5:
                    #perfusion_pattern.array[z_vox,y_vox,x_vox] = 0
                if mask_pattern[z_vox,y_vox,x_vox] < 0.5:
                    perfusion_pattern.data[z_vox,y_vox,x_vox] = 0
                else:
                    label_dist = {l:Perfusionradius for l in labels}
                    #world = mri_atlas.convert_voxel_to_world(array([z_vox,y_vox,x_vox]))
                    world = mri_atlas.convertVoxelToWorld(np.array([z_vox,y_vox,x_vox]))
                    for i in range(len(vpoint)):
                        v = vpoint[i]
                        l = g.edge_property(endpoint[i],'label')
                        d = distance(world,g.vertices[v].centre)
                        if d < label_dist[l]:
                            label_dist[l] = d
                    label_dist_sorted = sorted(label_dist.items(), key=operator.itemgetter(1))      #### sort dictionary by value: result list of tuples   #### sorted(label_dist.items(), key=operator.itemgetter(0)) ##sort dictionary by key
                    #perfusion_pattern.array[z_vox,y_vox,x_vox] =  int(round(label_dist_sorted[0][0]))
                    perfusion_pattern.data[z_vox,y_vox,x_vox] =  int(round(label_dist_sorted[0][0]))
    #for i in range(len(vpoint)):
        #v=vpoint[i]
        #world = g.vertices[v].centre
        #voxel = floor(mri_atlas.convert_world_to_voxel(array([world[0],world[1],world[2]])))
        #if voxel[0] in range(sizes[0]) and voxel[1] in range(sizes[1]) and voxel[2] in range(sizes[2]) :
            #perfusion_pattern.array[int(voxel[0]),int(voxel[1]),int(voxel[2])] = g.edge_property(endpoint[i],'label')
        
                
    #for z_vox in range(sizes[0]):
        #for y_vox in range(sizes[1]):
            #for x_vox in range(sizes[2]):
                #world=mri_atlas.convert_voxel_to_world(array([z_vox,y_vox,x_vox]))
                #min_dist=1e10
                #min_label = 0
                #for i in range(len(vpoint)):
                    #v=vpoint[i]
                    #if distance(world,g.vertices[v].centre)<Perfusionradius and distance(world,g.vertices[v].centre)<min_dist:
                        #min_dist = distance(world,g.vertices[v].centre)
                        #min_label = g.edge_property(endpoint[i],'label')
                #perfusion_pattern.array[z_vox,y_vox,x_vox]= int(min_label)
                
    #perfusion_pattern.output(perfusion_file)           
    perfusion_pattern.writeFile()
    perfusion_pattern.closeVolume()
    mri_atlas.closeVolume()
    mask_pattern.closeVolume()
    
#### masked version of perfusion_region   
def perfusion_region_nomask(g, mr_atlas_file, labels, label_remap, perfusion_file, Perfusionradius=1e10):
#   mr_atlas_file = "/projects/souris/sghanavati/data/c57_male_female_mri_atlas/c57_brain_segmentation_mnc2_merged_reg2vascular.mnc"
#   mr_atlas_file = "/projects/souris/sghanavati/data/c57_male_female_mri_atlas/c57_brain_segmentation_mnc2_merged_reg2vascular_downsample_merge.mnc"

    if label_remap:
        lpairs = label_remap.split(',')
        labelpairs = []
        for pair in lpairs:
            plabels = pair.split(':')
            plabels = [int(l) for l in plabels]
            labelpairs.append(plabels)
        tomap_l = [int(l[0]) for l in labelpairs]
        print labelpairs

        for e in g.edge_list():
            if g.edge_property(e,'label') in tomap_l:
                g.set_edge_property(e,'label',int(labelpairs[tomap_l.index(g.edge_property(e,'label'))][1]))


    endpoint = []
    vpoint=[]
    for e in g.edge_list():
        #if len(find_neighbor_edges(g,e))==2 and g.edge_property(e,'label') in labels:
        if g.edge_property(e,'label') in labels:
            if len(g.vertices[e[0]].edges)==1:
                vpoint.append(e[0])
                endpoint.append(e)
            elif len(g.vertices[e[1]].edges)==1:
                vpoint.append(e[1])
                endpoint.append(e)
                
    if not len(vpoint)==len(endpoint):
        print "ERROR: The length of endpoints and leaf edges should match but they don't!"
        exit(0)
    
    #### mri_distance feature:
    ## nc_data_type=py_minc.NC_BYTE is writing a byte (floating point) but we want the voxel intensities to be int labels
    ## nc_data_type=py_minc.NC_SHORT is writing a short integer data and we want the voxel intensities to be int labels
    #mri_atlas = py_minc.ArrayVolume(mr_atlas_file, nc_data_type=py_minc.NC_BYTE )# without this always [z,y,x], with ,dim_names=py_minc.FILE_ORDER_DIMENSION_NAMES order as the image    
    #perfusion_pattern = py_minc.ArrayVolume(mr_atlas_file, nc_data_type=py_minc.NC_BYTE )# without this always [z,y,x], with ,dim_names=py_minc.FILE_ORDER_DIMENSION_NAMES order as the image        
    #sizes = mri_atlas.get_sizes()
    mri_atlas = pymincf.volumeFromFile(mr_atlas_file, dtype='ubyte')   #ushort
    mri_atlas.openFile()
    sizes = mri_atlas.getSizes()
    spacing = mri_atlas.getSeparations()
    perfusion_pattern = pymincf.volumeFromInstance(mri_atlas, perfusion_file, dtype="ubyte")
    
    #make all voxels 0, give label for vpoints                
    for z_vox in range(sizes[0]):
        for y_vox in range(sizes[1]):
            for x_vox in range(sizes[2]):
                label_dist = {l:Perfusionradius for l in labels}
                #world = mri_atlas.convert_voxel_to_world(array([z_vox,y_vox,x_vox]))
                world = mri_atlas.convertVoxelToWorld(np.array([z_vox,y_vox,x_vox]))
                for i in range(len(vpoint)):
                    v = vpoint[i]
                    l = g.edge_property(endpoint[i],'label')
                    d = distance(world,g.vertices[v].centre)
                    if d < label_dist[l]:
                        label_dist[l] = d
                label_dist_sorted = sorted(label_dist.items(), key=operator.itemgetter(1))      #### sort dictionary by value: result list of tuples   #### sorted(label_dist.items(), key=operator.itemgetter(0)) ##sort dictionary by key
                #perfusion_pattern.array[z_vox,y_vox,x_vox] =  int(round(label_dist_sorted[0][0]))
                perfusion_pattern.data[z_vox,y_vox,x_vox] =  int(round(label_dist_sorted[0][0]))
    
    #perfusion_pattern.output(perfusion_file)            
    perfusion_pattern.writeFile()
    perfusion_pattern.closeVolume()
    mri_atlas.closeVolume()
    


def perfusion_volumes(perfusion_file,l):
    cnt =0
    ## nc_data_type=py_minc.NC_BYTE is writing a byte (floating point) but we want the voxel intensities to be int labels
    ## nc_data_type=py_minc.NC_SHORT is writing a short integer data and we want the voxel intensities to be int labels
    #perfusion_pattern = py_minc.ArrayVolume(perfusion_file, nc_data_type=py_minc.NC_BYTE )# without this always [z,y,x], with ,dim_names=py_minc.FILE_ORDER_DIMENSION_NAMES order as the image	
    #sizes = perfusion_pattern.get_sizes()
    #print sizes[0], sizes[1], sizes[2]
    #spacing = perfusion_pattern.get_separations()
    #print spacing[0], spacing[1], spacing[2]
    
    perfusion_pattern = pymincf.volumeFromFile(perfusion_file, dtype='ubyte')  #ushort
    perfusion_pattern.openFile()
    sizes = perfusion_pattern.getSizes()
    spacing = perfusion_pattern.getSeparations()

    for z_vox in range(sizes[0]):
        for y_vox in range(sizes[1]):
            for x_vox in range(sizes[2]):
                #if int(round(perfusion_pattern.array[z_vox,y_vox,x_vox]))==l:
                if perfusion_pattern[z_vox,y_vox,x_vox]==l:
                    cnt+=1
    #print cnt
    perfusion_pattern.closeVolume()
    return cnt*spacing[0]*spacing[1]*spacing[2]
                    
    
    
    
#####################################
#### by label graph calculations ####### 
#####################################
    
#fl can be 'label' (ground_truth) or 'cyl_label' (automatic labelling)
def num_segment_by_label(g,l, fl='label'):
    if not( fl=='label' or fl=='cyl_label'):
        print "ERROR! Feature can be either label for groundtruth labelling or cyl_label for automatic labelling"
        return -1
    n = 0
    for e in g.edge_list():
        if g.edge_property(e,fl)==l:
            n +=1
    return n	

def mean_mri_label_dist(g,l,mri_labels,f='mri_label_dist', fl='cyl_label'):
    if not( fl=='label' or fl=='cyl_label'):
        print "ERROR! Feature can be either label for groundtruth labelling or cyl_label for automatic labelling"
        return -1
    val = len(mri_labels)*[0]
    n = 0
    for e in g.edge_list():
        if g.edge_property(e,fl)==l:
            mrils = [i[0] for i in g.edge_property(e,f)]
            dists = [i[1] for i in g.edge_property(e,f)]
            n +=1
            for i in range(len(mri_labels)):
                ml = mri_labels[i]
                indx = mrils.index(ml)
                val[i] += dists[indx]
    if n>0: 
        val = [i/float(n) for i in val]
        return val   
    else:
        return -1

#fl can be 'label' (ground_truth) or 'cyl_label' (automatic labelling)
def mean_feature_by_label(g,l, f='diameter', fl='label'):
    if not( fl=='label' or fl=='cyl_label'):
        print "ERROR! Feature can be either label for groundtruth labelling or cyl_label for automatic labelling"
        return -1
    val = 0
    n = 0
    for e in g.edge_list():
        if f=='estimated_label':
            if g.edge_property(e,f)[0]==l:
                n +=1
                val += g.edge_property(e,f)[1]       
        else:    
            if g.edge_property(e,fl)==l:
                n +=1
                val += g.edge_property(e,f)
    if n>0:        
        return (val/float(n))	
    else:
        return -1

def sd_feature_by_label(g,l, f='diameter', fl='label'):
    if not( fl=='label' or fl=='cyl_label'):
        print "ERROR! Feature can be either label for groundtruth labelling or cyl_label for automatic labelling"
        return -1
    val = []
    for e in g.edge_list():
        if g.edge_property(e,fl)==l:
            val.append( g.edge_property(e,f))
    if len(val)>0:        
        return (std(val))	
    else:
        return -1

def min_feature_by_label(g,l, f='diameter', fl='label'):
    if not( fl=='label' or fl=='cyl_label'):
        print "ERROR! Feature can be either label for groundtruth labelling or cyl_label for automatic labelling"
        return -1
    val = []
    for e in g.edge_list():
        if g.edge_property(e,fl)==l:
            val.append( g.edge_property(e,f))
    if len(val)>0:        
        return (min(val))
    else:
        return -1

def max_feature_by_label(g,l, f='diameter', fl='label'):
    if not( fl=='label' or fl=='cyl_label'):
        print "ERROR! Feature can be either label for groundtruth labelling or cyl_label for automatic labelling"
        return -1
    val = []
    for e in g.edge_list():
        if g.edge_property(e,fl)==l:
            val.append( g.edge_property(e,f))
    if len(val)>0:        
        return (max(val))	
    else:
        return -1


#fl can be 'label' (ground_truth) or 'cyl_label' (automatic labelling)
def weighted_mean_feature_by_label(g,l, f='diameter',w='length', fl='label'):
    if not( fl=='label' or fl=='cyl_label'):
        print "ERROR! Feature can be either label for groundtruth labelling or cyl_label for automatic labelling"
        return -1
    val = 0
    n = 0
    for e in g.edge_list():
        if f=='estimated_label':
            if g.edge_property(e,f)[0]==l:
                if w=='volume':
                    n += g.edge_property(e,'diameter')*g.edge_property(e,'diameter')*g.edge_property(e,'length')
                    val += g.edge_property(e,'diameter')*g.edge_property(e,'diameter')*g.edge_property(e,'length')*g.edge_property(e,f)[1] 
                else:
                    n += g.edge_property(e,w)
                    val += g.edge_property(e,w)*g.edge_property(e,f)

        else:
            if g.edge_property(e,fl)==l:
                if w=='volume':
                    n += g.edge_property(e,'diameter')*g.edge_property(e,'diameter')*g.edge_property(e,'length')
                    val += g.edge_property(e,'diameter')*g.edge_property(e,'diameter')*g.edge_property(e,'length')*g.edge_property(e,f) 
                else:
                    n += g.edge_property(e,w)
                    val += g.edge_property(e,w)*g.edge_property(e,f)
    if n>0:        
        return (val/float(n))   
    else:
        return -1	


def mean_estimated_label_potentials_by_label(g,l, pl=[]):   #pl is paired labels
    val = 0
    n = 0
    for e in g.edge_list():
        if g.edge_property(e,'label')==l:
            n +=1
            estimated_list ={}
            for item in g.edge_property(e,'estimated_label_potentials'):
                estimated_list[item[0]] = item[1]
            val += estimated_list[l]  
            for otherlabel in pl:
                val += estimated_list[otherlabel] 
    if n>0:        
        return (val/float(n))   
    else:
        return -1

def weightedmean_estimated_label_potentials_by_label(g,l, w='volume', pl=[]):    #pl is paired labels
    val = 0
    n = 0
    for e in g.edge_list():
        if g.edge_property(e,'label')==l:
            if w=='volume':
                n += g.edge_property(e,'diameter')*g.edge_property(e,'diameter')*g.edge_property(e,'length')
                estimated_list ={}
                for item in g.edge_property(e,'estimated_label_potentials'):
                    estimated_list[item[0]] = item[1]
                val += g.edge_property(e,'diameter')*g.edge_property(e,'diameter')*g.edge_property(e,'length')*estimated_list[l]  
                for otherlabel in pl:
                    val += g.edge_property(e,'diameter')*g.edge_property(e,'diameter')*g.edge_property(e,'length')*estimated_list[otherlabel] 
            else:
                n += g.edge_property(e,w)
                estimated_list ={}
                for item in g.edge_property(e,'estimated_label_potentials'):
                    estimated_list[item[0]] = item[1]
                val += g.edge_property(e,w)*estimated_list[l]  
                for otherlabel in pl:
                    val += g.edge_property(e,w)*estimated_list[otherlabel] 
    if n>0:        
        return (val/float(n))   
    else:
        return -1


#fl can be 'label' (ground_truth) or 'cyl_label' (automatic labelling)
def total_feature_by_label(g,l, f='length', fl='label'):
    if not( fl=='label' or fl=='cyl_label'):
        print "ERROR! Feature can be either label for groundtruth labelling or cyl_label for automatic labelling"
        return -1
    val = 0
    for e in g.edge_list():
        if g.edge_property(e,fl)==l:
            val += g.edge_property(e,f)
    return val	
    
        


def volume_by_label(g,l, fl='label'):
    if not( fl=='label' or fl=='cyl_label'):
        print "ERROR! Feature can be either label for groundtruth labelling or cyl_label for automatic labelling"
        return -1
    val = 0
    for e in g.edge_list():
        if g.edge_property(e,fl)==l:
            val += (g.edge_property(e,'length')*numpy.pi*(g.edge_property(e,'diameter')*g.edge_property(e,'diameter')/4.0))
    
    return val

def corrected_volume_by_label(g,l, fl='label'):
    if not( fl=='label' or fl=='cyl_label'):
        print "ERROR! Feature can be either label for groundtruth labelling or cyl_label for automatic labelling"
        return -1
    val = 0
    for e in g.edge_list():
        if g.edge_property(e,fl)==l:
            val += (g.edge_property(e,'dist_length')*numpy.pi*(g.edge_property(e,'nooutlier_diameter')*g.edge_property(e,'nooutlier_diameter')/4.0))
    
    return val
#####################################
#### by label graph calculations ####### 
#####################################
    


#####################################
#### Total graph calculations ####### 
#####################################
def total_number_of_segments(g):
    return len(g.edge_list())
    
def total_number_of_segments(g,dthresh = 0.0):
    c = 0
    for e in g.edge_list():
        if g.edge_property(e,'diameter')>dthresh:
            c +=1
    return c

def total_volume_of_segments(g):
    vol = 0
    for e in g.edge_list():
        vol += (g.edge_property(e,'length')*numpy.pi*(g.edge_property(e,'diameter')*g.edge_property(e,'diameter')/4.0))   
    return vol

    
def total_volume_of_segments(g,dthresh = 0.0):
    vol = 0
    for e in g.edge_list():
        if g.edge_property(e,'diameter')>dthresh:
            vol += (g.edge_property(e,'length')*numpy.pi*(g.edge_property(e,'diameter')*g.edge_property(e,'diameter')/4.0))   
    return vol

def total_corrected_volume_of_segments(g,dthresh = 0.0):
    vol = 0
    for e in g.edge_list():
        if g.edge_property(e,'nooutlier_diameter')>dthresh:
            vol += (g.edge_property(e,'dist_length')*numpy.pi*(g.edge_property(e,'nooutlier_diameter')*g.edge_property(e,'nooutlier_diameter')/4.0))   
    return vol
    
def total_feature_of_segments(g, fl='length'):
    l = 0
    for e in g.edge_list():
        l += g.edge_property(e,fl)   
    return l
 
def total_feature_of_segments(g, fl='length',dthresh = 0.0):
    l = 0
    for e in g.edge_list():
        if g.edge_property(e,'diameter')>dthresh:
            l += g.edge_property(e,fl)   
    return l
 
def max_feature_of_segments(g, fl='length'):  
    l = []
    for e in g.edge_list():
        l.append(g.edge_property(e,fl))   
    return max(l)
    
def max_feature_of_segments(g, fl='length',dthresh = 0.0):  
    l = []
    for e in g.edge_list():
        if g.edge_property(e,'diameter')>dthresh:
            l.append(g.edge_property(e,fl))   
    return max(l)

def min_feature_of_segments(g, fl='length'):  
    l = []
    for e in g.edge_list():
        l.append(g.edge_property(e,fl))   
    return min(l)

def min_feature_of_segments(g, fl='length',dthresh = 0.0):  
    l = []
    for e in g.edge_list():
        if g.edge_property(e,'diameter')>dthresh:
            l.append(g.edge_property(e,fl))   
    return min(l)
#####################################
#### Total graph calculations ####### 
#####################################
    
    
def modified_simplify_graph(g,label = 'label'): #can also be cyl_label

    graph_analysis.remove_self_loops(g)

    for edge in g.edge_list():
        if 'intermediaries' not in g.edge_properties(edge).keys():
            g.set_edge_property(edge, 'intermediaries', [])
        
        
    vertex_list = numpy.unique([e[0] for e in g.edge_list()] + [e[1] for e in g.edge_list()]) 

    for i in vertex_list:
        if (len(g.vertices[i].edges) == 2):
            if not g.has_edge(tuple(g.vertices[i].edges)):
                new_e = g.order(tuple(g.vertices[i].edges))
                e1 = g.order(tuple((new_e[0], i)))
                e2 = g.order(tuple((i, new_e[1])))
                if len(g.edge_properties(e2).keys())> len(g.edge_properties(e1).keys()):
                    tmp=e2
                    e2=e1
                    e1=tmp
                #print i, " ", new_e, " ", e1, " ", e2
                
                properties={}
                if (label in g.edge_properties(e1).keys()) and (label in g.edge_properties(e2).keys()):
                    g.set_edge_property(e1,label,int(g.edge_property(e1,label)))    #make sure label is int
                    g.set_edge_property(e2,label,int(g.edge_property(e2,label)))
                    #print g.edge_property(e1,'cyl_label') ," ", g.edge_property(e2,'cyl_label')
                    if (g.edge_property(e1,label) == g.edge_property(e2,label)):
                        properties={label: g.edge_property(e1,label)}
                    else:   #if the 2 adjacent edges that should be merged into 1 edge have different cyl_labels then get the label of the longer edge_length
                        if (this_vessel_path_length(g, e1) > this_vessel_path_length(g, e2)):
                            properties={label: g.edge_property(e1,label)}
                        else:
                            properties={label: g.edge_property(e2,label)}
                            
                if ('intermediaries' in g.edge_properties(e1).keys()) and ('intermediaries' in g.edge_properties(e2).keys()):
                    first = list(copy.copy(g.edge_property(e1, 'intermediaries')))
                    second = list(copy.copy(g.edge_property(e2, 'intermediaries')))
                    if i < new_e[0]:
                        first.reverse()
                    if i > new_e[1]:
                        second.reverse()
                    intermediaries = first + [i] + second                
                    properties['intermediaries']=intermediaries
                    
                #print g.edge_properties(e1).keys()
                #print g.edge_properties(e2).keys()
                #for n in g.edge_properties(e2).keys():
                    #print n ,":", g.edge_property(e2,n)
                    #print n ,":", g.edge_property(e1,n)
                #print "\n\n"   
                for name in g.edge_properties(e2).keys():
                    if label == 'label' and name=='cyl_label':
                        properties['cyl_label'] = properties['label']
                    elif  label == 'cyl_label' and name=='label':
                        properties['label'] = properties['cyl_label']
                    elif not name=='intermediaries' and not name=='label' and not name=='cyl_label' and not name=='diameter' and not name=='length':
                        #if type(g.edge_property(e1,name))==list or type(g.edge_property(e1,name))==array:
                            #print name ,":", g.edge_property(e1,name), " , " , g.edge_property(e2, name)
                        if type(g.edge_property(e1,name))==list: 
                            properties[name] = []
                            properties[name].append(g.edge_property(e1, name))
                            properties[name].append(g.edge_property(e2, name))  
                        elif type(g.edge_property(e1,name))==array:
                            properties[name] = []
                            properties[name].append(list(g.edge_property(e1, name)))
                            properties[name].append(list(g.edge_property(e2, name)))  
                            properties[name]= array(properties[name])
                        else:
                            if (this_vessel_path_length(g, e1) > this_vessel_path_length(g, e2)):
                                properties[name] = g.edge_property(e1, name)
                            else:
                                properties[name] = g.edge_property(e2, name)
                g.add_edge(new_e, properties)
                g.disconnect_vertex(i)

    #print len(g.edge_list())       
    return g

def this_vessel_path_length(g, edge):
    #"estimate vessel lengths by summing the euclidian distance between intermediate points"
    edge = g.order(edge)
    edge_length = lambda g, edge: \
        linalg.norm(g.vertices[edge[0]].centre - g.vertices[edge[1]].centre)
    points = [edge[0]] + list(g.edge_property(edge, 'intermediaries')) + [edge[1]]
    segments = [(points[i], points[i+1]) for i in range(len(points)-1)]
    return sum([edge_length(g,segment) for segment in segments])


def this_vessel_avg_diameter(g, edge):
    #"estimate vessel diameter by averaging the diameters at intermediate points"
    diameters = [g.vertices[i].radius*2 for i in \
                [edge[0]] + list(g.edge_property(edge, 'intermediaries')) \
                + [edge[1]]]
    return mean(diameters)


def this_estimate_edge_lengths(g, length_fcn=this_vessel_path_length, property='length'):
    for edge in g.edge_list():
        #print edge, " ",
        length = length_fcn(g, edge)
        #print length
        g.set_edge_property(edge, property, numpy.float32(length))


def this_estimate_edge_diameters(g, diameter_fcn=this_vessel_avg_diameter, property='diameter'):
    for edge in g.edge_list():
        #print edge," ",max([edge[0]] + list(g.edge_property(edge, 'intermediaries')) + [edge[1]]), " / ",
        diameter = diameter_fcn(g, edge)
        #print diameter
        g.set_edge_property(edge, property, numpy.float32(diameter))




    