#!/usr/bin/env python

#just open the graph and reduce to intermediaries, calculate length and diameter and save 
#
#  Created July 27, 2011
#  Sahar Ghanavati

    

from vessel_tracking import graph_analysis
from optparse import OptionParser, Option, OptionValueError
from sys import argv
import os, shelve, string ,time
import commands
import copy
#---------------------------------------------------------------------------------
#

program_name = 'simplify_graph.py'


if __name__ == '__main__':
        
    usage = "Usage: "+program_name+" [options] input_filename.db [output_filename.db]\n"+\
            "   or  "+program_name+" --help";
            
    parser = OptionParser(usage)
    
    parser.add_option("--clobber", action="store_true", dest="clobber",
                    default=0, help="overwrite output file")
                    

    options, args = parser.parse_args()
    print '\n\n>>> %s [%s]: %s\n' % (time.ctime(time.time()), os.getcwd(), string.join(argv))
    print "\n\n",options, "\n\n", args , "\n\n"
   

    if len(args) == 1:
        input_file = args[0]	 #eval: convert str to list	#input_file, output_file= args => input_file is list
        output_file= input_file[:-3]+"_simplified.db"

    elif len(args) == 2:	
        input_file, output_file = args
    else:
        parser.error("incorrect number of arguments")
    



    #print("the input file is: %s\n" %(input_file))

    if not options.clobber and os.path.exists(output_file):
        raise SystemExit, \
            "The --clobber option is needed to overwrite an existing file."


    try:
        # open graph data and copy contents to memory
        g, attributes = graph_analysis.input_graph(input_file, ["history", "vertex_offsets"])  	

        print ("Succefully read in the input.db\n")	
    except:
        print("Error reading in the input.db\n")
        
    h=copy.deepcopy(g)
    # simplify h by reducing to vessel segments
    graph_analysis.simplify_graph_retaining_intermediaries(h)
    # estimate lengths and average diameters of each segment
    graph_analysis.estimate_edge_diameters(h)
    graph_analysis.estimate_edge_lengths(h)

    
    history = '\n>>> %s: %s' % (time.ctime(time.time()), string.join(argv))
    if attributes.has_key("history"):
        history = attributes["history"] + "\n" + history
        del attributes['history']
        
    graph_analysis.output_graph(output_file, h, history, attributes)   
    
    print ("Succefully wrote the output_simplified.db\n")
