from vessel_tracking import graph_analysis
from numpy import *
import numpy
import numpy as np
import numpy.matlib as Matlib
from optparse import OptionParser, Option, OptionValueError

from sys import argv
import os, commands, time, string, math
import numpy as np
import scipy
import copy
import operator
from cerebrovascular_analysis import vessel_analysis



program_name = 'keep_connectedcomponents.py'


if __name__ == '__main__':
        
    usage = "Usage: "+program_name+" [options] input_graph.db output_graph.db \n"+\
            "   or  "+program_name+" --help";
            
    parser = OptionParser(usage)
    
    parser.add_option("--clobber", action="store_true", dest="clobber",
                    default=0, help="overwrite output file")
    parser.add_option("--size", dest="size", default=10,
                    help="threshold on the size of the connected components to be deleted",
                    type="int")
    parser.add_option("--lenthresh", dest="lenthresh", default=0,
                    help="threshold on the length of the connected components to be deleted",
                    type="float")
    parser.add_option("--diameterthresh", dest="diameterthresh", default=0.1,
                    help="threshold on the diameter of the connected components to be deleted",
                    type="float")
       
        
                
    options, args = parser.parse_args()
    print '\n\n>>> %s [%s]: %s\n' % (time.ctime(time.time()), os.getcwd(), string.join(argv))
    print "\n\n",options, "\n\n", args , "\n\n"
   
    history = '>>> %s: %s' % (time.ctime(time.time()), string.join(argv))

    if len(args) == 2:
        input_file, output_file = args
    else:
        parser.error("incorrect number of arguments")

    if not options.clobber and os.path.exists(output_file):
        raise SystemExit, \
              "The --clobber option is needed to overwrite an existing file."
    
    g, attributes = graph_analysis.input_graph(input_file, ["history", "vertex_offsets"])
    h = copy.deepcopy(g)
    
    connected_components = vessel_analysis.connected_components(g)
    for c in connected_components:
        print len(c), " ,",
    print ""

    
    for cc in connected_components:
        if len(cc)<options.size:
            for e in cc:
                if g.edge_property(e,'length')<options.lenthresh and g.edge_property(e,'diameter')<options.diameterthresh:
                    g.remove_edge(e) 
                
    if attributes.has_key("history"):
        history = attributes["history"] + "\n" + history
        del attributes['history']

    graph_analysis.output_graph(output_file, g, history, attributes)
    
    cmd=("graph2cylinder.py --use_label %s" %(output_file)) 
    os.system(cmd)

    
