#  Sahar Ghanavati


from vessel_tracking import graph_analysis
from morphology import graph
#from numpy import *
from optparse import OptionParser, Option, OptionValueError
from sys import argv
import os, string, time
import commands
import shelve, copy

#from time import time, ctime
#import copy
#---------------------------------------------------------------------------------
#
program_name = 'unsimplify_graph.py'


if __name__ == '__main__':
        
    usage = "Usage: "+program_name+" [options] input.db outbut.db \n"+\
            "   or  "+program_name+" --help";
            
    parser = OptionParser(usage)

    parser.add_option("--clobber", action="store_true", dest="clobber",
                        default=0, help="overwrite output file")
                       

    options, args = parser.parse_args()

    print '>>> %s [%s]: %s\n' % (time.ctime(time.time()), os.getcwd(), string.join(argv))
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
        
    h = copy.deepcopy(g)
    for e in h.edge_list():    
        if 'intermediaries' in h.edge_properties(e).keys():
            for i in range(len(h.edge_property(e,'intermediaries'))-1):
                g.add_edge(tuple((h.edge_property(e,'intermediaries')[i],h.edge_property(e,'intermediaries')[i+1])))
                if 'label' in h.edge_properties(e).keys():
                    g.set_edge_property(tuple((h.edge_property(e,'intermediaries')[i],h.edge_property(e,'intermediaries')[i+1])),'label',h.edge_property(e,'label'))
                if 'cyl_label' in h.edge_properties(e).keys():
                    g.set_edge_property(tuple((h.edge_property(e,'intermediaries')[i],h.edge_property(e,'intermediaries')[i+1])),'cyl_label',h.edge_property(e,'cyl_label'))
                if 'estimated_label' in h.edge_properties(e).keys():
                    g.set_edge_property(tuple((h.edge_property(e,'intermediaries')[i],h.edge_property(e,'intermediaries')[i+1])),'estimated_label',h.edge_property(e,'estimated_label'))
                
                g.add_edge(tuple((e[0],h.edge_property(e,'intermediaries')[0])))
                if 'label' in h.edge_properties(e).keys():
                    g.set_edge_property(tuple((e[0],h.edge_property(e,'intermediaries')[0])),'label',h.edge_property(e,'label'))
                if 'cyl_label' in h.edge_properties(e).keys():
                    g.set_edge_property(tuple((e[0],h.edge_property(e,'intermediaries')[0])),'cyl_label',h.edge_property(e,'cyl_label'))
                if 'estimated_label' in h.edge_properties(e).keys():
                    g.set_edge_property(tuple((e[0],h.edge_property(e,'intermediaries')[0])),'estimated_label',h.edge_property(e,'estimated_label'))
                
                g.add_edge(tuple((h.edge_property(e,'intermediaries')[-1],e[1])))
                if 'label' in h.edge_properties(e).keys():
                    g.set_edge_property(tuple((h.edge_property(e,'intermediaries')[-1],e[1])),'label',h.edge_property(e,'label'))
                if 'cyl_label' in h.edge_properties(e).keys():
                    g.set_edge_property(tuple((h.edge_property(e,'intermediaries')[-1],e[1])),'cyl_label',h.edge_property(e,'cyl_label'))
                if 'estimated_label' in h.edge_properties(e).keys():
                    g.set_edge_property(tuple((h.edge_property(e,'intermediaries')[-1],e[1])),'estimated_label',h.edge_property(e,'estimated_label'))
            g.remove_edge(e)
        else:
            #delete the properties
            g.remove_edge(e)
            g.add_edge(e)
            if 'label' in h.edge_properties(e).keys():
                g.set_edge_property(e,'label',h.edge_property(e,'label'))
            if 'cyl_label' in h.edge_properties(e).keys():
                g.set_edge_property(e,'cyl_label',h.edge_property(e,'cyl_label'))
            if 'estimated_label' in h.edge_properties(e).keys():
                g.set_edge_property(e,'estimated_label',h.edge_property(e,'estimated_label'))

    history = '\n>>> %s: %s' % (time.ctime(time.time()), string.join(argv))     ##argv = ['code.py', 'input.db', 'output']
    if attributes.has_key("history"):
        history = attributes["history"] + "\n" + history
        del attributes['history']
        
    graph_analysis.output_graph(output_file, g, history, attributes)   
    
    print ("Succefully wrote the %s\n" %output_file)
            
             