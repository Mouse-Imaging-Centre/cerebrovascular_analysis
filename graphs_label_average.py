# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 21:29:09 2015

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
import numpy as np
import scipy
import copy
import operator


program_name = 'graphs_label_average.py'


if __name__ == '__main__':
        
    usage = "Usage: "+program_name+" [options] labelled_graph1.db labelled_graph2.db ...  labelled_graphN.db output_labeled.db \n"+\
            "   or  "+program_name+" --help";
            
    parser = OptionParser(usage)
    
    parser.add_option("--clobber", action="store_true", dest="clobber",
                    default=0, help="overwrite output file")
    parser.add_option("--posterior_prob", action="store_true", dest="posterior_prob",
                    default=0, help="use the probability instead of the label")
    parser.add_option("--labelprop", type="string", dest="labelprop",
                    default="cyl_label", help="which label property to be used for averaging (Default: 'cyl_label'")
    parser.add_option("--reference_graph", type="string", dest="reference_graph",
                    default="", help="Write the output averaged labels on what graph")
    parser.add_option("--weighted", action="store_true", dest="weighted",
                    default=0, help="weighted average by vessel segment volume (diameter^2 * length)")
    #parser.add_option("--remove_interms",action="store_true",default=0, dest="remove_interms",help="remove the intermediate mnc files that were created for calculating seeds in iterations")
        
        
    tic = time.time()
    #ct_labeling = ct_labeling()		#initialize the class ct_labeling, make an instance called ct_labeling
                
    options, args = parser.parse_args()
    print '\n\n>>> %s [%s]: %s\n' % (time.ctime(time.time()), os.getcwd(), string.join(argv))
    print "\n\n",options, "\n\n", args , "\n\n"
   

        
    if len(args) < 2:
        parser.error("incorrect number of arguments")
    else:
        input_files = args[0:-1]
        output_file = args[-1]
        
    
    if not options.clobber and os.path.exists(output_file):
        raise SystemExit, \
            "The --clobber option is needed to overwrite an existing file."
            
    label_vals = {}
    weight_vals = {}
    weighted_label_vals = {}
    
    for  input_file in input_files:
        # reads in the graph to be labeled
        try:
            g, attributes = graph_analysis.input_graph(input_file, ["history", "vertex_offsets"])
            print ("Succefully read in the %s\n" %input_file)	
        except:
            print("Error reading in the %s\n" %input_file)
        
        for e in g.edge_list():
            l = g.edge_property(e,'label')
            if l not in weight_vals.keys():
                weight_vals[l] = 0
            if l not in weighted_label_vals.keys():
                weighted_label_vals[l] = 0
            if options.labelprop=='estimated_label':
                weight_vals[l] += g.edge_property(e,'diameter')*g.edge_property(e,'diameter')*g.edge_property(e,'length')       #cylinder volume = pi*r^2*length
                if options.posterior_prob:
                    #if not l==g.edge_property(e,'estimated_label')[0][0]:
                        #print "ERROR in file:", input_file, " edge ", e, " the true label", l, " is not equal the label in posterior_prob estimated_label", g.edge_property(e,'estimated_label')[0][0], " ! \n Aborted!"
                        #exit(0)
                    #else:
                    val = g.edge_property(e,'estimated_label')[0][1]
                    weighted_label_vals[l] += g.edge_property(e,'diameter')*g.edge_property(e,'diameter')*g.edge_property(e,'length')*g.edge_property(e,'estimated_label')[0][1]
                else:
                    val = g.edge_property(e,'estimated_label')[0][0]
                    weighted_label_vals[l] += g.edge_property(e,'diameter')*g.edge_property(e,'diameter')*g.edge_property(e,'length')*g.edge_property(e,'estimated_label')[0][0]
            else:
                val = g.edge_property(e,options.labelprop)
            if l not in label_vals.keys():
                label_vals[l] = [val]
            else:
                label_vals[l].append(val)
  
    if options.reference_graph!="":
        print "write the average onto ", options.reference_graph
        ref_graph, attributes = graph_analysis.input_graph(options.reference_graph, ["history", "vertex_offsets"])
    else:
        print "write the average onto ", input_files[0]
        ref_graph, attributes = graph_analysis.input_graph(input_files[0], ["history", "vertex_offsets"])
      
    for e in ref_graph.edge_list():
        l = ref_graph.edge_property(e,'label')
        if l in label_vals.keys():
            if options.weighted:
                ref_graph.set_edge_property(e,'cyl_label', int(np.mean(label_vals[l])))  #"No label                \0" #add label dataset to edge properties and initialize it to None with null-terminted character '\0'
                if weight_vals[l]==0:
                    ref_graph.set_edge_property(e,'estimated_label', [l,0])  #"No label                \0" #add label dataset to edge properties and initialize it to None with null-terminted character '\0'            
                else:    
                    ref_graph.set_edge_property(e,'estimated_label', [l,weighted_label_vals[l]/weight_vals[l]])  #"No label                \0" #add label dataset to edge properties and initialize it to None with null-terminted character '\0'            
            else:
                ref_graph.set_edge_property(e,'cyl_label', int(np.mean(label_vals[l])))  #"No label                \0" #add label dataset to edge properties and initialize it to None with null-terminted character '\0'
                ref_graph.set_edge_property(e,'estimated_label', [l,np.mean(label_vals[l])])  #"No label                \0" #add label dataset to edge properties and initialize it to None with null-terminted character '\0'
        else:
            ref_graph.set_edge_property(e,'cyl_label', int(0))
            ref_graph.set_edge_property(e,'estimated_label', [l,0])  #"No label                \0" #add label dataset to edge properties and initialize it to None with null-terminted character '\0'


    history = '>>> %s: %s' % (time.ctime(time.time()), string.join(argv))		##argv = ['code.py', 'input.db', 'output']
    if attributes.has_key("history"):
        history = attributes["history"] + "\n" + history
        del attributes['history']
        
    graph_analysis.output_graph(output_file[:-2]+"db", ref_graph, history, attributes)   
    #graph_analysis.output_graph("simplifiedgraph.db", h, history, attributes)

    print ("Succefully wrote the %s\n" %output_file[:-2]+"db")

    ## save the final cylinder file into h5 file for the Brainview
    #output_final=output_file.replace(".db", ".h5")
    cmd=("graph2graph.py %s %s --clobber " %(output_file[:-2]+"db",output_file[:-2]+"h5"))	#python /micehome/jgsled/bin/
    print(cmd)
    os.system(cmd)
