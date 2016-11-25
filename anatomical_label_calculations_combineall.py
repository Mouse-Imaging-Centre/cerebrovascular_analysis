from vessel_tracking import graph_analysis
from optparse import OptionParser, Option, OptionValueError
from sys import argv
import os, commands, time, string, math
import numpy as np
#from time import time, ctime
from sys import argv
import shelve, os, string
from copy import deepcopy
import vessel_analysis

        
        
program_name = 'anatomical_label_calculations_combineall.py'


if __name__ == '__main__':
        
    usage = "Usage: "+program_name+" [options]  input1.db input2.db ... inputN.db \n"+\
            "   or  "+program_name+" --help";
            
    parser = OptionParser(usage)
    
    parser.add_option("--clobber", action="store_true", dest="clobber",
                    default=0, help="overwrite output file")
    parser.add_option("--reference", type="string", dest="reference", default="",
            help="give the reference database name ")
    parser.add_option("--output_name", type="string", dest="output_name", default="",
            help="give the output name for the output db file")
        
    parser.add_option("--reconnect", action="store_true", dest="reconnect",
                    default=0, help="connect all components of each label")
    
                
    options, args = parser.parse_args()
    print '\n\n>>> %s [%s]: %s\n' % (time.ctime(time.time()), os.getcwd(), string.join(argv))
    print "\n\n",options, "\n\n", args , "\n\n"
   
    input_files = args[:]
    
    
    history = '\n>>> %s: %s' % (time.ctime(time.time()), string.join(argv))
    
    ref , ref_attributes = graph_analysis.input_graph(options.reference, ["history", "vertex_offsets"]) 
    #h = copy.deepcopy(ref)
    
    keep_edges = []
    
    for input_file in input_files:
        print ">>> ",input_file
        if os.path.exists(input_file):
            g, attributes = graph_analysis.input_graph(input_file, ["history", "vertex_offsets"])    
            for e in g.edge_list():
                keep_edges.append(e)
            
    
    for e in ref.edge_list():
        if e not in keep_edges:
            if ref.edge_property(e,'cyl_label') not in [9,200]: #don't remove pcommAs
                ref.remove_edge(e)
        
    if ref_attributes.has_key("history"):
        #history = attributes["history"] + "\n>>>Main labels retained." + "\n" + history
        history = "\n>>>Anatomical length calculations." + "\n" + ref_attributes["history"] 
        del ref_attributes['history']

    graph_analysis.output_graph(options.output_name, ref, history, ref_attributes)
    cmd=("\npython /micehome/sghanavati/Desktop/scripts/cerebrovascular_analysis/graph2cylinder.py --clobber --use_label %s" %(options.output_name)) #python /micehome/jgsled/bin/
    os.system(cmd)
        
   