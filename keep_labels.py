# -*- coding: utf-8 -*-
"""
Created on Wed Nov 17 21:29:09 2015

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
import os, commands, time, string, math
import numpy as np
import scipy
import copy
import operator



program_name = 'keep_labels.py'


if __name__ == '__main__':
        
    usage = "Usage: "+program_name+" [options] input_graph.db [output_graph.db] \n"+\
            "   or  "+program_name+" --help";
            
    parser = OptionParser(usage)
    
    parser.add_option("--clobber", action="store_true", dest="clobber",
                    default=0, help="overwrite output file")
    parser.add_option("--labels", dest="labels",
                    help="Label numbers to keep in the output_file, separated by commas, i.e. --labels 3,11,21.",
                    type="string")
        
        
                
    options, args = parser.parse_args()
    print '\n\n>>> %s [%s]: %s\n' % (time.ctime(time.time()), os.getcwd(), string.join(argv))
    print "\n\n",options, "\n\n", args , "\n\n"
   

        
    if len(args) == 1:
        input_file = args[0]
        output_file = input_file[:-3]+"_keeplabels.db"
    elif len(args) == 2:
        input_file, output_file = args
    else:
        parser.error("incorrect number of arguments")

    if not options.clobber and os.path.exists(output_file):
        raise SystemExit, \
              "The --clobber option is needed to overwrite an existing file."
        
    ##1. borrow from delete_internal_leaves to find segments that are inside
    
    ##2. borrow from smooth_graph to decide which one is parent and which is daughter branch that needs to be deleted     
    
    ##3. borrow from adjust_junction to find the nearest connection point with preserving the topology 
    
    
    
    
    history = '>>> %s: %s' % (time.ctime(time.time()), string.join(argv))

    ls = []
    if options.labels:
        ls = options.labels.split(',')
        ls = [int(l) for l in ls]
    
    
    g, attributes = graph_analysis.input_graph(input_file, ["history", "vertex_offsets"])
    h = copy.deepcopy(g)
    
    for e in h.edge_list():
        if not h.edge_property(e,'label') in ls:
            g.remove_edge(e) 
        
    
    
    #NOTE: what we should actually do is to check if all the vertices of an edge are inside another edge => remove that edge
    # if not, then only disconnect_vertex of those inside, then reconnect the rest of the vertices to the encompassing_e

    if attributes.has_key("history"):
        history = attributes["history"] + "\n" + history
        del attributes['history']

    graph_analysis.output_graph(output_file, g, history, attributes)
    
    cmd=("graph2cylinder.py --use_label %s" %(output_file)) 
    os.system(cmd)

    cmd=("graph2obj.py %s %s %s " %(output_file,output_file[:-3]+"obj", options.clobber)) 
    os.system(cmd)


    
    
    
    
    
    
    
    
    
    
    
