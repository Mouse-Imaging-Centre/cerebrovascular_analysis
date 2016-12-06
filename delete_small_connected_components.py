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
from cerebrovascular_analysis import vessel_analysis 


program_name = 'delete_small_connected_components.py'


if __name__ == '__main__':
        
    usage = "Usage: "+program_name+" [options] input_graph.db [output_graph.db] \n"+\
            "   or  "+program_name+" --help";
            
    parser = OptionParser(usage)
    
    parser.add_option("--clobber", action="store_true", dest="clobber",
                    default=0, help="overwrite output file")
        
    parser.add_option("--threshold", type="int", dest="threshold", default = 1.0,
                       help="threshold for number of vessels in a connected component to be deleted (default 1)")
        
                
    options, args = parser.parse_args()
    print '\n\n>>> %s [%s]: %s\n' % (time.ctime(time.time()), os.getcwd(), string.join(argv))
    print "\n\n",options, "\n\n", args , "\n\n"
   

        
    if len(args) == 1:
        input_file = args[0]
        output_file = input_file[:-3]+"_delcc.db"
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

    g, attributes = graph_analysis.input_graph(input_file, ["history", "vertex_offsets"])
    
    connected_components = vessel_analysis.connected_components(g)
    
    for cc in connected_components:
        if len(cc) < options.threshold:
            for e in cc:
                g.remove_edge(e)
    

    if attributes.has_key("history"):
        history = attributes["history"] + "\n" + history
        del attributes['history']

    graph_analysis.output_graph(output_file, g, history, attributes)
    
    cmd=("graph2obj.py %s %s %s " %(output_file,output_file[:-3]+"obj", options.clobber))
    os.system(cmd)

    cmd=("graph2cylinder.py --use_label %s %s %s " %(output_file,output_file[:-3]+".h5", options.clobber)) 
    os.system(cmd)

    
    
