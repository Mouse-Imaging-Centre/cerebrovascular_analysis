

#  Created March 3, 2016
#  Last modified 
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
def distance(x,y):
    d = sqrt((x[0]-y[0])*(x[0]-y[0]) + (x[1]-y[1])*(x[1]-y[1])+(x[2]-y[2])*(x[2]-y[2]))
    return d



program_name = 'delete_cage.py'


if __name__ == '__main__':
        
    usage = "Usage: "+program_name+" [options] input_file.db output_file.db\n"+\
            "   or  "+program_name+" --help";
            
    parser = OptionParser(usage)

    parser.add_option("--clobber", action="store_true", dest="clobber",
                        default=0, help="overwrite output file")
    parser.add_option("--diameter_threshold",type="float",dest="diameter_threshold",default = 100.0,help="The threshold diameter in mm for deleting cage around vessels with diameter larger than threshold. Default value is 100.0mm")
    parser.add_option("--mincfile", type="string", dest="mincfile",
            help="give the name of mincfile to be used as the template for graph2isosurf")
                        

    options, args = parser.parse_args()

    print '>>> %s [%s]: %s\n' % (ctime(time()), os.getcwd(), string.join(argv))


    #if len(args) != 1:
        #parser.error("incorrect number of arguments")
        
    #print (args)   #list of str = ['inputname', 'outputname']

    #if len(args) == 1:
        #input_file = args[0]    #eval: convert str to list #input_file, output_file= args => input_file is list

        #out1= input_file[:-2]+"db"
        #out1=out1.replace('tree','graph')
    #elif len(args) == 2:   
        #input_file, out1 = args
    #else:
        #parser.error("incorrect number of arguments")

    if len(args)==2:
        input_file, output_file = args
    elif len(args)==1:
        input_file = args[0]
        output_file = input_file[:-3]+"_delcage.db"
    else:
        parser.error("incorrect number of arguments")
        

    #print output_file
        
    if not options.clobber and os.path.exists(output_file):
        raise SystemExit, \
            "The --clobber option is needed to overwrite an existing file."
            
    # reads in the graph for adjust junction and radi corrections
    try:
        g, attributes = graph_analysis.input_graph(input_file, ["history", "vertex_offsets"])
        print ("Succefully read in the %s\n" %input_file) 
    except:
        print("Error reading in the %s\n" %input_file)
        
        
    large_edges = []
    for edge in g.edge_list():
        if g.edge_property(edge,'diameter')>options.diameter_threshold:
            large_edges.append(edge)
    
    
    delete_edges = []
    for edge in large_edges:
        for othere in g.edge_list():
            if edge != othere:
                distances = []
                for v1 in [edge[0],edge[1]]+ g.edge_property(edge,'intermediaries'):
                    for v2 in [othere[0],othere[1]]+ g.edge_property(othere,'intermediaries'):
                       distances.append(distance(v1,v2)) 
                if max(distances) <= g.edge_property(edge,'diameter'):
                    delete_edges.append(othere)
                    
    for e in delete_edges:
        g.remove_edge(e)
        
    history = '\n>>> %s: deleted cages around large vessels with diameter larger than %f' % (ctime(time()), options.diameter_threshold)     ##argv = ['code.py', 'input.db', 'output']
    if attributes.has_key("history"):
        history = attributes["history"] + "\n" + history
        del attributes['history']
        
    graph_analysis.output_graph(output_file, g, history, attributes)   
    print ("Succefully wrote the %s\n" %output_file)
    
    cmd=("graph2isosurf.py --contour --save_volume %s %s %s " %(output_file, options.mincfile,output_file[:-3]+"_isosurface.obj")) #python /micehome/jgsled/bin/
    print(cmd)
    os.system(cmd)


            
        
        
        
        
        
