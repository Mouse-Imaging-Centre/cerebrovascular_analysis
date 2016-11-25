#!/usr/bin/env python
# -*- coding: utf-8 -*-

#  Created July 28, 2011
#  Last modified Oct 6, 2011
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

program_name = 'modify_db_intermediaries.py'


if __name__ == '__main__':
        
    usage = "Usage: "+program_name+" [options] graph_from_cylinder.db ref_graph.db [output.db]\n"+\
            "   or  "+program_name+" --help";
            
    parser = OptionParser(usage)
    
    parser.add_option("--new_edges", action="store_true", dest="new_edges",
                    default=0, help="Only recalculate feature properties for the new_edges added in brain-view")
    options, args = parser.parse_args()
    
    print '\n\n>>> %s [%s]: %s\n' % (time.ctime(time.time()), os.getcwd(), string.join(argv))
    print "\n\n",options, "\n\n", args , "\n\n"
    sys.stdout.flush()
    

    if len(args) == 3:
        graph_file, ref_graph, output_file = args
    elif len(args)==2:
        graph_file, ref_graph = args
        output_file = graph_file[:-2]+"db"
    else:
        parser.error("incorrect number of arguments")


    try:
        g, attributes = graph_analysis.input_graph(graph_file, ["history", "vertex_offsets"])  	# open graph data and copy contents to memory
        print ("Succefully read in the %s\n" %graph_file)	
        sys.stdout.flush()
    except:
        print("Error reading in the %s \n" %graph_file)
        sys.stdout.flush()

    try:
        ref_g = graph_analysis.input_graph(ref_graph)
        print ("Succefully read in the %s\n" %ref_graph)	
        sys.stdout.flush()
    except:
        print("Error reading in the %s \n" %ref_graph)
        sys.stdout.flush()


    for e in g.edge_list():
        if e in ref_g.edge_list():
            g.set_edge_property(e,'intermediaries', ref_g.edge_property(e,'intermediaries'))
        else:
            g.set_edge_property(e,'intermediaries',[])
            
            
    history = '>>> %s: %s' % (time.ctime(time.time()), string.join(argv))           ##argv = ['code.py', 'input.db', 'output']
    if attributes.has_key("history"):
            history = attributes["history"] + "\n" + history
            del attributes['history']
    
    graph_analysis.output_graph(output_file, g, history, attributes)   
    
    print ("Succefully wrote the %s\n" %output_file)
    sys.stdout.flush()
            
        


