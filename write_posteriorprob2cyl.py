from vessel_tracking import graph_analysis
from optparse import OptionParser, Option, OptionValueError
from sys import argv
import os, shelve, string,sys
import commands
from time import time, ctime

program_name = 'write_posteriorprob2cyl.py'


if __name__ == '__main__':
        
    usage = "Usage: "+program_name+" [options] CTimage.mnc postprocess_output_vessel_tracking.db \n"+\
            "   or  "+program_name+" --help";
            
    parser = OptionParser(usage)

    parser.add_option("--clobber", action="store_true", dest="clobber",
                        default=0, help="overwrite output file")

    options, args = parser.parse_args()
    print '>>> %s [%s]: %s\n' % (ctime(time()), os.getcwd(), string.join(argv))



    if len(args) != 2:
        parser.error("incorrect number of arguments")
        
    #print (args)   #list of str = ['inputname', 'outputname']

    input_graph, output_graph = args  #eval: convert str to list #input_file, output_file= args => input_file is list
    output_graph = input_graph[:-3]+"posteriorprob.db"
    
    #### reads in the graph for more corrections
    try:
        g, attributes = graph_analysis.input_graph(input_graph, ["history", "vertex_offsets"])
        print ("Succefully read in the %s\n" %input_graph)  
    except:
        print("Error reading in the %s\n" %input_graph)

    for e in g.edge_list():
        g.set_edge_property(e,'cyl_label', g.edge_property(e,'estimated_label')[1])
    
 
    if attributes.has_key("history"):
        history = attributes["history"] + "\n" + history
        del attributes['history']

    graph_analysis.output_graph(output_file, g, history, attributes)
    
    cmd=("\npython /projects/souris/sghanavati/src/scripts/cerebrovascular_analysis/graph2cylinder.py --use_label %s" %(output_file)) #python /micehome/jgsled/bin/
    os.system(cmd)
   