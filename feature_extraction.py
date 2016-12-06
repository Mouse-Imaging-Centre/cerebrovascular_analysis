#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Calculate all the following features for each edge
# NOTE: in the training we can select best features (most discriminative)

#features:
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


#  Created July 28, 2011
#  Last modified May 28, 2013
#  Sahar Ghanavati


from vessel_tracking import graph_analysis
from numpy import *
import math
from optparse import OptionParser, Option, OptionValueError
from minc_util.progress import progress_report
from sys import argv
import os, shelve, string,sys
import commands
import copy
import time
import py_minc 
import numpy 

#---------------------------------------------------------------------------------
#

def calculate_dist(u, v ):
    dist = math.sqrt( pow(u[0]-v[0],2) + pow(u[1]-v[1],2) + pow(u[2]-v[2],2))
    return dist


def angle_between(u,v):
    alpha = 0
    alpha = arccos((dot(u,v))/(linalg.norm(u)*linalg.norm(v)))
    #if alpha > (pi/2):
        ##print alpha
        #alpha = pi- alpha
    return alpha
    
def vector_dir(g,edge):
    e1=edge[0]
    e2=edge[1]
    direction = g.vertices[e2].centre - g.vertices[e1].centre
    direction = direction/linalg.norm(direction)
    return direction







program_name = 'feature_extraction.py'


if __name__ == '__main__':
        
    usage = "Usage: "+program_name+" [options] graph_simplified.db mr_atlas_file.mnc mr_atlas_centroids.db cba_direction_reference.db output_filename.db\n"+\
            "   or  "+program_name+" --help";
            
    parser = OptionParser(usage)
    
    parser.add_option("--new_edges", action="store_true", dest="new_edges",
                    default=0, help="Only recalculate feature properties for the new_edges added in brain-view")
    options, args = parser.parse_args()
    
    print '\n\n>>> %s [%s]: %s\n' % (time.ctime(time.time()), os.getcwd(), string.join(argv))
    print "\n\n",options, "\n\n", args , "\n\n"
    sys.stdout.flush()
    
    if len(args) != 5:
        parser.error("incorrect number of arguments")
        
    
    graph_file, mr_atlas_file, mr_centroids_file, dir_reference_file, output_file = args


    #### mri_distance feature:
    cmd=("feature_MRI_labels.py %s %s %s --clobber" %(graph_file,mr_atlas_file, output_file))	
    print(cmd)
    sys.stdout.flush()
    os.system(cmd)	
    

    try:
        g, attributes = graph_analysis.input_graph(output_file, ["history", "vertex_offsets"])  	# open graph data and copy contents to memory
        print ("Succefully read in the %s\n" %graph_file)	
        sys.stdout.flush()
    except:
        print("Error reading in the %s \n" %graph_file)
        sys.stdout.flush()



    
    
    mr = shelve.open(mr_centroids_file, "r")
    mr_centroids = mr['mr_centroids']
    mr.close()
    
    #### find curvature and midpoint and turtoisity of each edge
    print "total edges :", len(g.edge_list())
    sys.stdout.flush()
    i=0
    for e in g.edge_list():
        print i
        sys.stdout.flush()
        i+=1
        if "label" not in g.edge_properties(e).keys():
            g.set_edge_property(e,"label", int(0))
        if "estimated_label" not in g.edge_properties(e).keys():
            g.set_edge_property(e,"estimated_label", [[0,0.0]])
        
        interm=[]
        interm.append(e[0])
        if 'intermediaries' in g.edge_properties(e).keys():
            for intermed in g.edge_property(e,'intermediaries'):
                interm.append(intermed)
        interm.append(e[1])
        
        if len(interm)%2 == 1:
            indx=int(math.floor(len(interm)/2))
            midpoint = g.vertices[interm[indx]].centre
        else:
            indx1=int(len(interm)/2 -1)
            indx2=int(len(interm)/2) 
            midpoint = (g.vertices[interm[indx1]].centre + g.vertices[interm[indx2]].centre)/2
        g.set_edge_property(e,'midpointX', numpy.float32(midpoint[0]))	
        g.set_edge_property(e,'midpointY', numpy.float32(midpoint[1]))	
        g.set_edge_property(e,'midpointZ', numpy.float32(midpoint[2]))	
        median_point = (g.vertices[e[0]].centre + g.vertices[e[1]].centre)/2
        g.set_edge_property(e,'curvature',numpy.float32(calculate_dist(midpoint, median_point)/calculate_dist(g.vertices[e[0]].centre, g.vertices[e[1]].centre)))	#curvature calculate by distance of edge midpoint to median_point    
            
        #e_tortuosity = sqrt( (g.vertices[e[0]].centre[0] - g.vertices[e[1]].centre[0])*(g.vertices[e[0]].centre[0] - g.vertices[e[1]].centre[0]) + (g.vertices[e[0]].centre[1] - g.vertices[e[1]].centre[1])*(g.vertices[e[0]].centre[1] - g.vertices[e[1]].centre[1]) + (g.vertices[e[0]].centre[2] - g.vertices[e[1]].centre[2])*(g.vertices[e[0]].centre[2]) - g.vertices[e[1]].centre[2]) / g.edge_property(e,'length')
        g.set_edge_property(e,'tortuosity',numpy.float32(g.edge_property(e,'length')/calculate_dist(g.vertices[e[0]].centre, g.vertices[e[1]].centre)))	#tortuosity calculate by distance of edge first and last point divided by path length  =  ratio of the length of the curve (L) to the distance between the ends of it (C)
            
            
        #### calculate the proximity feature vector	
        proximity=[]
        for l in mr_centroids.keys():
            dist= calculate_dist(midpoint, mr_centroids[l] )
            #proximity[l]=dist
            proximity.append([l,dist])
        g.set_edge_property(e,'proximity', proximity)




    #### calculate the angles with references		
    try:
        ref, ref_attributes = graph_analysis.input_graph(dir_reference_file, ["history", "vertex_offsets"])  	# open graph data and copy contents to memory
        print ("Succefully read in the cba_vasculature_dir_reference_file.db\n")
        sys.stdout.flush()
    except:
        print("Error reading in the cba_vasculature_dir_reference_file.db\n")
        sys.stdout.flush()
    
    label_namelist=[]
    for refedge in ref.edge_list():
        if ref.edge_property(refedge,'label')!=0:
            if ref.edge_property(refedge,'label') not in label_namelist:
                label_namelist.append(ref.edge_property(refedge,'label'))
            else:
                print ("Error more than one reference direction for label %f\n" %(ref.edge_property(refedge,'label')))
                print("Aborted!")
                sys.stdout.flush()
                exit(0)
        
    
    for edge in g.edge_list():
        direction_feat=[]
        for refedge in ref.edge_list():
            if ref.edge_property(refedge,'label') in label_namelist:
                alpha = 0 #initialize it in case arccos gives error
                edge_dir = vector_dir(g,edge)
                ref_dir = vector_dir(ref,refedge)
                alpha = angle_between(edge_dir, ref_dir)
                #direction_feat[ref.edge_property(refedge,'label')] = alpha
                direction_feat.append([ref.edge_property(refedge,'label'),alpha])
        g.set_edge_property(edge,"anglewref", direction_feat)	
        
    
    
    #### directionX,Y,Z and negatives, angles with +X,+Y, +Z
    x_pos = (1.0, 0.0, 0.0)
    y_pos = (0.0, 1.0, 0.0)
    z_pos = (0.0 ,0.0 ,1.0)	
        
    for edge in g.edge_list():
        alpha = 0 #initialize it in case arccos gives error
        edge_dir = vector_dir(g,edge)
        alpha = angle_between(edge_dir, x_pos)
        g.set_edge_property(edge, 'angleX' , numpy.float32(alpha))
        alpha = angle_between(edge_dir, y_pos)
        g.set_edge_property(edge, 'angleY' , numpy.float32(alpha))
        alpha = angle_between(edge_dir, z_pos)
        g.set_edge_property(edge, 'angleZ' , numpy.float32(alpha))
        if alpha > (pi/2):				# the edge_dir is in z_neg
            #g.set_edge_property(edge, 'direction_neg' , [edge_dir[0],edge_dir[1],edge_dir[2]])
            #g.set_edge_property(edge, 'direction_pos' , [-edge_dir[0],-edge_dir[1],-edge_dir[2]])
            g.set_edge_property(edge,'directionX_neg', numpy.float32(edge_dir[0]))
            g.set_edge_property(edge,'directionY_neg', numpy.float32(edge_dir[1]))
            g.set_edge_property(edge,'directionZ_neg', numpy.float32(edge_dir[2]))
            g.set_edge_property(edge,'directionX', numpy.float32(-edge_dir[0]))
            g.set_edge_property(edge,'directionY', numpy.float32(-edge_dir[1]))
            g.set_edge_property(edge,'directionZ', numpy.float32(-edge_dir[2]))
        else:				# the edge_dir is in z_pos
            #g.set_edge_property(edge, 'direction_pos' , [edge_dir[0],edge_dir[1],edge_dir[2]])
            #g.set_edge_property(edge, 'direction_neg' , [-edge_dir[0],-edge_dir[1],-edge_dir[2]])
            g.set_edge_property(edge,'directionX_neg', numpy.float32(-edge_dir[0]))
            g.set_edge_property(edge,'directionY_neg', numpy.float32(-edge_dir[1]))
            g.set_edge_property(edge,'directionZ_neg', numpy.float32(-edge_dir[2]))
            g.set_edge_property(edge,'directionX', numpy.float32(edge_dir[0]))
            g.set_edge_property(edge,'directionY', numpy.float32(edge_dir[1]))
            g.set_edge_property(edge,'directionZ', numpy.float32(edge_dir[2]))
    
    #relative edge diameter and angle
    for e in g.edge_list():
        rel_diameter=[]
        rel_dir=[]
        e1_neighbours= g.vertices[e[0]].edges
        for v in e1_neighbours:                      
            if v!=e[1]:
                neighbor_edge=tuple((min(e[0],v),max(e[0],v)))
                rel_diameter.append(numpy.float32(g.edge_property(neighbor_edge,'diameter')/g.edge_property(e,'diameter')))
                rel_dir.append(numpy.float32(angle_between(vector_dir(g,neighbor_edge), vector_dir(g,e))))
        e1_neighbours= g.vertices[e[1]].edges
        for v in e1_neighbours:                      
            if v!=e[0]:
                neighbor_edge=tuple((min(e[1],v),max(e[1],v)))
                rel_diameter.append(g.edge_property(neighbor_edge,'diameter')/g.edge_property(e,'diameter'))
                rel_dir.append(numpy.float32(angle_between(vector_dir(g,neighbor_edge), vector_dir(g,e))))
        g.set_edge_property(e,'rel_diameter',rel_diameter)
        g.set_edge_property(e,'rel_dir',rel_dir)

    #print "\n\nmri_label_dist:\nnumber of g.vertices = ", len(g.vertices)
    #for edge in g.edge_list():
        ##vertex_list=[edge[0]]+g.edge_property(edge,'intermediaries')+[edge[1]]
        #interm=[]
        #interm.append(edge[0])
        #if 'intermediaries' in g.edge_properties(edge).keys():
            #for intermed in g.edge_property(edge,'intermediaries'):
                #interm.append(intermed)
        #interm.append(edge[1])

        #mr_label_prop = []
        #mr_label_prop_dict={}
        #for l in mr_labels[0][1].keys():
            #mr_label_prop_dict[l] = 0
        #for v in interm:
            #print v," ",
            #world = g.vertices[v].centre
            #voxel = floor(mri_atlas.convert_world_to_voxel(array([world[0],world[1],world[2]])))	#voxel[z,y,x], world[x,y,z]
            ##find voxel in mr_labels[:,0]
            #indx = [list(l[0]) for l in mr_labels].index([int(voxel[0]),int(voxel[1]),int(voxel[2])])		#mr_labels[indx][1] is {l:d, l:d,...}
            ##save in a dictionary and then average and save in edge properties
            #for l in mr_labels[0][1].keys():
                #mr_label_prop_dict[l]+=mr_labels[indx][1][l]
        #for l in mr_labels[0][1].keys():		
            #mr_label_prop.append([l, mr_label_prop_dict[l]/float(len(interm))])	
        #g.set_edge_property(edge,'mri_label_dist', mr_label_prop)	

    
    history = '>>> %s: %s' % (time.ctime(time.time()), string.join(argv))		##argv = ['code.py', 'input.db', 'output']
    if attributes.has_key("history"):
        history = attributes["history"] + "\n" + history
        del attributes['history']
        
    graph_analysis.output_graph(output_file, g, history, attributes)   
    
    print ("Succefully wrote the %s\n" %output_file)
    sys.stdout.flush()
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    
    
