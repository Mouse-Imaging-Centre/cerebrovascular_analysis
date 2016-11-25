
#  Created July 13, 2011
#  Last modified August 7, 2012
#  Sahar Ghanavati


from vessel_tracking import graph_analysis
from morphology import graph
#from numpy import *
from optparse import OptionParser, Option, OptionValueError
from sys import argv
import os, string
import commands
from time import time, ctime
#import copy
#---------------------------------------------------------------------------------
#
program_name = 'radi_correction.py'


if __name__ == '__main__':
        
    usage = "Usage: "+program_name+" [options] input.db outbut.db \n"+\
            "   or  "+program_name+" --help";
            
    parser = OptionParser(usage)

    parser.add_option("--clobber", action="store_true", dest="clobber",
                        default=0, help="overwrite output file")
    parser.add_option("--radi_correction",type="float",dest="radi_correction",default = 1.06,help="The radius correction factor for the vessel_tracking underestimated result, Defult=1.06")
                       

    options, args = parser.parse_args()

    print '>>> %s [%s]: %s\n' % (ctime(time()), os.getcwd(), string.join(argv))
    if not len(args)==2:
        parser.error("incorrect number of arguments")
    else:
        input_file, output_file = args
        
    if not options.clobber and os.path.exists(output_file):
        raise SystemExit, \
            "The --clobber option is needed to overwrite an existing file."
    # reads in the graph for adjust junction and radi corrections
    try:
        g, attributes = graph_analysis.input_graph(input_file, ["history", "vertex_offsets"])
        print ("Succefully read in the %s\n" %input_file) 
    except:
        print("Error reading in the %s\n" %input_file)
        
    # optimize location of branch points
    graph_analysis.adjust_junctions(g, attributes['vertex_offsets'])

    # vertex radii need to be corrected for estimation bias due to blurring of 40% of vessel diameter => radius should be correct by 60%
    for v in g.vertices:
        v.radius = v.radius*options.radi_correction

    # save the corrected graph
    history = '\n>>> %s: corrected branch point locations and radius bias due to blurring by %f' % (ctime(time()), options.radi_correction)     ##argv = ['code.py', 'input.db', 'output']
    if attributes.has_key("history"):
        history = attributes["history"] + "\n" + history
        del attributes['history']
        
    graph_analysis.output_graph(output_file, g, history, attributes)   
    print ("Succefully wrote the %s\n" %output_file)
