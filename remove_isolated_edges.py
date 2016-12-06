#!/usr/bin/env python

# remove the isolated edges (where each vertex only has one neighbour edge) from the simplified graph
# endpoints edges total number of neighbouring edges are 3 (one vertex has 3 neighbouring edges and one vertex has 1 neighbouring edge)
# for other edges total number of neighbouring edges are 5 (each vertex has 3 neighbouring edges)


#  Created Mar 14, 2012
#  Last modified April 11, 2012
#  Sahar Ghanavati


from vessel_tracking import graph_analysis
from optparse import OptionParser, Option, OptionValueError
from sys import argv
import os, string
import commands
from time import time, ctime
#---------------------------------------------------------------------------------
#

program_name = 'remove_isolated_edges.py'


if __name__ == '__main__':
        
    usage = "Usage: "+program_name+" [options] simplified_graph.db [output_file.db]\n"+\
            "   or  "+program_name+" --help";
            
    parser = OptionParser(usage)
    
    parser.add_option("--clobber", action="store_true", dest="clobber",
                    default=0, help="overwrite output file")
    parser.add_option("--layer",action="store_true",default=0,dest="layer",help="remove isolated edges in the shape of a Y ")
    #parser.add_option("--diameter_threshold",type="float",dest="diameter_threshold",default = 0.2,help="The threshold diameter for masking vessels,Defult=0.2")
    #parser.add_option("--length_threshold",type="float",dest="length_threshold",default = 0.5,help="The threshold length for masking shorter vessels from an endpoint, Defult=0.5")
                    

    options, args = parser.parse_args()
    print '\n\n>>> %s [%s]: %s\n' % (ctime(time()), os.getcwd(), string.join(argv))
    print "\n\n",options, "\n\n", args , "\n\n"
    

    #if len(args) != 1:
        #parser.error("incorrect number of arguments")
        
    #print (args)	#list of str = ['inputname', 'outputname']

    if len(args) == 1:
        input_file = args[0]	 #eval: convert str to list	#input_file, output_file= args => input_file is list
    
        output_file= input_file[:-3]+"_rmisolate.db"

    elif len(args) == 2:	
        input_file, output_file = args
    else:
        parser.error("incorrect number of arguments")
    
    if not options.clobber and os.path.exists(output_file):
        raise SystemExit, \
            "The --clobber option is needed to overwrite an existing file."


    # reads in the graph for more corrections
    try:
        g, attributes = graph_analysis.input_graph(input_file, ["history", "vertex_offsets"])
        print ("Succefully read in the %s\n" %input_file)	
    except:
        print("Error reading in the %s\n" %input_file)
        
    # remove isolate edges
    remove_edge_list =[]
    for e in g.edge_list():
        #n=0     
        #n=len(g.vertices[e[0]].edges)+len(g.vertices[e[1]].edges)
        #if (n!=6 and n!=4):
        if ( (len(g.vertices[e[0]].edges)==1) and (len(g.vertices[e[1]].edges)==1) ):#if n==2:
            #print e
            #print g.vertices[e[0]].edges, g.vertices[e[1]].edges
            remove_edge_list.append(e)
    
    for i in range(len(remove_edge_list)):
        g.remove_edge(remove_edge_list[i])
    print ("%d isolated edges were removed!" %len(remove_edge_list))	

    if options.layer:
        remove_edge_list =[]
        for e in g.edge_list():
            Y_flag=1	#suppose this edge is part of a Y structure unless we find later in for loop that it's not
            if len(g.vertices[e[0]].edges)==1:
                for f in g.vertices[e[1]].edges:	#the centre vertex of Y structure
                    if len(g.vertices[f].edges)!=1:
                        Y_flag=0
                if Y_flag:
                    for f in g.vertices[e[1]].edges:
                        if f>e[1]:
                            remove_edge_list.append(tuple((e[1],f)))
                        else:
                            remove_edge_list.append(tuple((f,e[1])))	
            Y_flag=1	#suppose this edge is part of a Y structure unless we find later in for loop that it's not
            if len(g.vertices[e[1]].edges)==1:
                for f in g.vertices[e[0]].edges:	#the centre vertex of Y structure
                    if len(g.vertices[f].edges)!=1:
                        Y_flag=0
                if Y_flag:
                    for f in g.vertices[e[0]].edges:
                        if f>e[0]:
                            remove_edge_list.append(tuple((e[0],f)))
                        else:
                            remove_edge_list.append(tuple((f,e[0])))
                        
        remove_edge_list=list(set(remove_edge_list))	#unique members				
        for i in range(len(remove_edge_list)):
            g.remove_edge(remove_edge_list[i])
        print ("%d Y-shaped edges were removed!" %len(remove_edge_list))	

                
        
        
    # save the corrected graph
    history = '\n>>> %s: %s' % (ctime(time()), "\n>>>removed the isolated edges (where each vertex only has one neighbour edge)")		##argv = ['code.py', 'input.db', 'output']
    if attributes.has_key("history"):
        history = attributes["history"] + "\n" + history
        del attributes['history']
        
    graph_analysis.output_graph(output_file, g, history, attributes)   
    print ("Succefully wrote the %s\n" %output_file)
















            
            
            
            
            
            
            
            
            
            
            
