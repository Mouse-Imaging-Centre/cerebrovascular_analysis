#!/usr/bin/env python

#delete short leaves of *.db file and save it in *_delshort.db
#
#
#  Created Apr 13, 2011
#  Sahar Ghanavati

#connect to bianca (system with hardy os 
#export PYTHONPATH=/home/jgsled/lib64_ubuntu_lucid/python:$PYTHONPATH
#python delshortleaves.py --help
#python delshortleaves.py ../tree_mask_blur_delleaves.db ../outputname --clobber				

from vessel_tracking import graph_analysis, path_io, path_analysis, curvature, filters,\
    vessel_tracking, medial_atoms, path_visualize, tubes
from morphology import graph, cluster_skeleton, object_io
from numpy import *
from optparse import OptionParser, Option, OptionValueError
from minc_util.progress import progress_report
from sys import argv
import os, shelve, string ,time
import commands
import copy
#---------------------------------------------------------------------------------
#

program_name = 'delete_short_leaves.py'


if __name__ == '__main__':
        
    usage = "Usage: "+program_name+" [options] unsimplified_graph.db [output_filename.db]\n"+\
        "Usage: --radius_threshold should be set > 0 to remove any leaves, then either use --length_threshold and/or --intermed_threshold\n"+\
            "   or  "+program_name+" --help";
            
    parser = OptionParser(usage)
    
    parser.add_option("--clobber", action="store_true", dest="clobber",
                    default=0, help="overwrite output file")
                    
    parser.add_option("--radius_threshold",type="float",dest="radius_threshold",default = 0.06,help="The threshold radius in mm for masking vessels. Default value is 0.06mm")
    parser.add_option("--length_threshold",type="float",dest="length_threshold",default = 0.2,help="The threshold length in mm for masking vessels,Default is 0.2mm")
    parser.add_option("--intermed_threshold",type="int",dest="intermed_threshold",default = 0,help="The threshold on number of intermediaries for masking shorter vessels from an endpoint. Recommended value is 100.")
                
    options, args = parser.parse_args()
    print '\n\n>>> %s [%s]: %s\n' % (time.ctime(time.time()), os.getcwd(), string.join(argv))
    print "\n\n",options, "\n\n", args , "\n\n"
    
    if len(args) == 1:
        input_file = args[0]	 #eval: convert str to list	#input_file, output_file= args => input_file is list
    
        output_file= input_file[:-3]+"_delshort.db"

    elif len(args) == 2:	
        input_file, output_file = args
    else:
        parser.error("incorrect number of arguments")
    
    if not options.clobber and os.path.exists(output_file):
        raise SystemExit, \
            "The --clobber option is needed to overwrite an existing file."

        
    radius_threshold = options.radius_threshold 
    length_threshold = options.length_threshold 
    intermed_threshold = options.intermed_threshold	

    if not options.clobber and os.path.exists(output_file):
        raise SystemExit, \
            "The --clobber option is needed to overwrite an existing file."


    try:
        g, attributes = graph_analysis.input_graph(input_file, ["history", "vertex_offsets"])  	# open graph data and copy contents to memory
        #g = graph_analysis.input_graph("tree_mask_blur_delleaves.db")
        print ("Succefully read in the input.db\n")	
    except:
        print("Error reading in the input.db\n")
        
    h=copy.deepcopy(g)
    # simplify h by reducing to vessel segments
    graph_analysis.simplify_graph_retaining_intermediaries(h)
    # estimate lengths and average diameters of each segment
    graph_analysis.estimate_edge_diameters(h)
    graph_analysis.estimate_edge_lengths(h)
    
    remove_edge_num=0
    for e in h.edge_list():
        interm=[]
        ##radius thresholding:
        if len(g.vertices[e[0]].edges)==1 and g.vertices[e[0]].radius< radius_threshold :  	##an edge with an end point (leave) noraml range of 0.01<radius<0.25mm   
            #print("edge (%d,%d) is leaf!" %(e[0],e[1]) )
            interm=[v for v in h.edge_property(e, 'intermediaries')]
            interm.append(e[0])
            remove_edge_num +=1
        if len(g.vertices[e[1]].edges)==1 and g.vertices[e[0]].radius< radius_threshold:	###noraml range of 0.01<radius<0.25mm   
            #print("edge (%d,%d) is leaf!" %(e[0],e[1]) )
            interm=[v for v in h.edge_property(e, 'intermediaries')]
            interm.append(e[1])
            remove_edge_num +=1
        if len(g.vertices[e[0]].edges)==1 and len(g.vertices[e[1]].edges)==1:
            #print ("edge (%d,%d) is isolated!" %(e[0],e[1]) )
            interm=[v for v in h.edge_property(e, 'intermediaries')]
            interm.append(e[0])
            interm.append(e[1])
        interm=list(set(interm)) 	#get rid of repetition	
        ##length or number of intermediaries thresholding:
        if len(interm)<intermed_threshold or h.edge_property(e,'length')< length_threshold:		##find short leaves
            #print(len(interm))
            for v in interm:
                g.disconnect_vertex(v)
                    
    history = '>>> %s: %s' % (time.ctime(time.time()), string.join(argv))		##argv = ['code.py', 'input.db', 'output']
    if attributes.has_key("history"):
        history = attributes["history"] + "\n" + history
        del attributes['history']
        
    graph_analysis.output_graph(output_file, g, history, attributes)   
    
    print ("%d short leaves out of %s edges were removed!" %(remove_edge_num,len(h.edge_list())))	
    print ("Succefully wrote the %s\n" %output_file)				

