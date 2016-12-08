#!/usr/bin/env python
#  Created July 16, 2012
#  Last modified Nov 16, 2012 
#  Sahar Ghanavati


from vessel_tracking import graph_analysis
from numpy import *
import math
from optparse import OptionParser, Option, OptionValueError
from minc_util.progress import progress_report
from sys import argv
import os, shelve, string
import commands
import copy
import time, sys
import py_minc 

#---------------------------------------------------------------------------------
#

def calculate_dist(u, v ):
    dist = math.sqrt( pow(u[0]-v[0],2) + pow(u[1]-v[1],2) + pow(u[2]-v[2],2))
    return dist

def calculate_mean(g, edge):
    #"estimate vessel diameter by averaging the diameters at intermediate points"
    diameters = [g.vertices[i].radius*2 for i in \
                [edge[0]] + list(g.edge_property(edge, 'intermediaries')) \
                + [edge[1]]]
    return mean(diameters)

program_name = 'feature_MRI_labels.py'
if __name__ == '__main__':
        
    usage = "Usage: "+program_name+" [options] graph_simplified.db mr_atlas_file.mnc output_filename.db\n"+\
            "   or  "+program_name+" --help";
            
    parser = OptionParser(usage)
    
    parser.add_option("--clobber", action="store_true", dest="clobber",
                    default=0, help="overwrite output file")
    options, args = parser.parse_args()
    print '\n\n>>> %s [%s]: %s\n' % (time.ctime(time.time()), os.getcwd(), string.join(argv))
    print "\n\n",options, "\n\n", args , "\n\n"
    

    if len(args) != 3:
        parser.error("incorrect number of arguments")
        
    
    graph_file, mr_atlas_file, output_file = args

    if not options.clobber and os.path.exists(output_file):
        raise SystemExit, \
            "The --clobber option is needed to overwrite an existing file."


    try:
        g, attributes = graph_analysis.input_graph(graph_file, ["history", "vertex_offsets"])  	# open graph data and copy contents to memory
        print ("Successfully read in the %s\n" %graph_file)	
        sys.stdout.flush()
    except:
        print("Error reading in the %s \n" %graph_file)
        sys.stdout.flush()



    #### mri_distance feature:
    mri_atlas = py_minc.ArrayVolume(mr_atlas_file, nc_data_type=py_minc.NC_BYTE)# without this always [z,y,x], with ,dim_names=py_minc.FILE_ORDER_DIMENSION_NAMES order as the image	
    sizes = mri_atlas.get_sizes()
    mri_labels =[]
    for z_vox in range(sizes[0]):
        for y_vox in range(sizes[1]):
            for x_vox in range(sizes[2]):
                if mri_atlas.array[z_vox,y_vox,x_vox]>0 and mri_atlas.array[z_vox,y_vox,x_vox] not in mri_labels:
                    mri_labels.append(mri_atlas.array[z_vox,y_vox,x_vox])
    mri_labels=map(lambda x:int(round(x)), mri_labels)
    
    
    for l in mri_labels:
        #print "mri_label",l
        sys.stdout.flush()
        mri_atlas = py_minc.ArrayVolume(mr_atlas_file[:-4]+"_%d_dt.mnc" %l, nc_data_type=py_minc.NC_BYTE)# without this always voxels[z,y,x], with ,dim_names=py_minc.FILE_ORDER_DIMENSION_NAMES order as the image	
        for e in g.edge_list():
            interm=[]
            interm.append(e[0])
            if 'intermediaries' in g.edge_properties(e).keys():
                for intermed in g.edge_property(e,'intermediaries'):
                    interm.append(intermed)
            interm.append(e[1])
            dt_val=0	
            for v in interm:
                world = g.vertices[v].centre
                voxel = floor(mri_atlas.convert_world_to_voxel(array([world[0],world[1],world[2]])))
                dt_val+=round(mri_atlas.array[min(int(voxel[0]),sizes[0]-1),min(int(voxel[1]),sizes[1]-1),min(int(voxel[2]),sizes[2]-1)])
            if 'mri_label_dist' not in g.edge_properties(e).keys():
                g.set_edge_property(e,'mri_label_dist',[[float(l),float(dt_val/float(len(interm)))]])
            else:
                mrlist = list(g.edge_property(e,'mri_label_dist'))
                mrlist.append([ float(l) , float(dt_val/float(len(interm))) ])
                g.set_edge_property(e,'mri_label_dist', mrlist)
                

    del(mri_atlas)

    
    history = '>>> %s: %s' % (time.ctime(time.time()), string.join(argv))		##argv = ['code.py', 'input.db', 'output']
    if attributes.has_key("history"):
        history = attributes["history"] + "\n" + history
        del attributes['history']
        
    graph_analysis.output_graph(output_file, g, history, attributes)   
    
    print ("Successfully wrote the %s\n" %output_file)
    
