## Created on December 8, 2008 by Monique Rennie.  Input file is a .db dataset.  Prunes all terminal vessels (leaves) less than a certain diameter threshold
## Updated on August 25, 2015 by Sahar Ghanavati to add delete by connected components
### Usage: python nibbler.py input.db threshold output.db

from vessel_tracking import path_io, path_analysis, curvature, filters,\
     vessel_tracking, medial_atoms, path_visualize, tubes, graph_analysis
from morphology import cluster_skeleton
import pylab
from Numeric import *
from morphology import object_io
#import py_minc
#from minc_util import minc_view
#from scipy import amax, amin, rand
#from scipy import gplt, linalg
import shelve, copy,time
#from scipy.interpolate import splprep, splev, splrep
#import shelve, copy
#from scipy.optimize import leastsq
from optparse import OptionParser, Option, OptionValueError
import sys, os
import string
from sys import argv
#from time import time, ctime
import vessel_analysis
import numpy as np


program_name = 'nibbler.py'

if __name__ == '__main__':
    
    usage = "Usage: "+program_name+" [options] input.db output.db\n"+\
         "   or  "+program_name+" --help";
    parser = OptionParser(usage)
    parser.add_option("--clobber", action="store_true", dest="clobber", default=0, help="overwrite output file")
    #parser.add_option("--connected_component", action="store_true", dest="connected_component", default=0, help="delete the connected components that have number of edges less than a threshold")
    parser.add_option("--diameter_threshold",type="float", dest="diameter_threshold", default=0, help="the diameter threshold for deleting the leaves in [mm] (the range is usually 0.01 to 0.5mm for mouse brains)")                  
    parser.add_option("--rel_diameter_threshold",type="float", dest="rel_diameter_threshold", default=0, help="the threshold on the relative diameter of leaves w.r.t their parents for deleting the leaves (Default: 0, Ex. 0.1)")                  

    parser.add_option("--cc_threshold",type="int", dest="cc_threshold", default=0, help="the threshold of the number of edges in each connected component that should be deleted")   
    parser.add_option("--cc_length_threshold",type="int", dest="cc_length_threshold", default=0, help="the threshold of the total length of edges in each connected component that should be deleted")   
    parser.add_option("--cc_diameter_threshold",type="int", dest="cc_diameter_threshold", default=0, help="the threshold of the mean diameter of edges in each connected component that should be deleted")   

    options, args = parser.parse_args()
    if len(args) != 2:
         parser.error("incorrect number of arguments")
    
    input, output = args
    	
	
    #correction = string.atof(threshold)
    print "diameter threshold:", options.diameter_threshold
    history = '>>> %s: %s' % (time.ctime(time.time()), string.join(argv))
    
    #cmd=("\npython /hpf/largeprojects/MICe/sghanavati/src/scripts/cerebrovascular_analysis/unsimplify_graph.py --clobber %s %s" %(input, input[:-3]+"_unsimplified.db")) #python /micehome/jgsled/bin/
    #os.system(cmd)
    

    #db = shelve.open(input[:-3]+"_unsimplified.db")
    
    #print "Opening input.db file..."
    
    #g = db['skeletal_graph']
    #vertex_offsets = db['vertex_offsets']
    g, attributes = graph_analysis.input_graph(input, ["history", "vertex_offsets"])    
    vertex_offsets = attributes["vertex_offsets"]
    
    graph_analysis.simplify_graph_retaining_intermediaries(g)
    graph_analysis.estimate_edge_lengths(g)
    graph_analysis.estimate_edge_diameters(g)
    h = copy.deepcopy(g)
    #graph_analysis.adjust_junctions(h, vertex_offsets)
    ##j = copy.deepcopy(h)
    #diameters = [j.edge_property(edge, 'diameter') for edge in j.edge_list()]
    #lengths, diameters = graph_analysis.collect_edge_lengths_and_diameters(j)
    while True:
        cnt = 0 
        del_leaves = []
        for e in g.edge_list():
            if len(g.vertices[e[0]].edges) ==1 and len(g.vertices[e[1]].edges) ==1: #a dangling bit that must be deleted
                del_leaves.append(e)
            elif len(g.vertices[e[0]].edges) ==1 or len(g.vertices[e[1]].edges) ==1:   #a leaf with e[0] the end point or #a leaf with e[1] the end point
                if g.edge_property(e, 'diameter') < options.diameter_threshold:
                    del_leaves.append(e)
        for l in del_leaves:
            g.remove_edge(l)
            cnt += 1
        print "removed ", cnt," leaves"
        if len(del_leaves)==0:
            break
    
    
    
           
    ### determine which vertices are at the leaves (end points)
    #leaves = [index for index in range (len(h.vertices)) if len(h.vertices[index].edges) ==1]
    #print "Popping leaves...\nNumber of initial leaves ", len(leaves)
    #cnt = 0   
    #while len(leaves)>5:
        #print "Number of leaves " , len(leaves)
        #print "Number of disconnected leaves " , cnt
        #leaf=leaves.pop(0)
        #leafotherv = h.vertices[leaf].edges[0]
        #edge = (leaf, leafotherv)
        #if h.edge_property(edge, 'diameter') < options.diameter_threshold and len(h.vertices[leaf].edges) ==1:
            #if options.rel_diameter_threshold:
                #for e in h.vertices[leafotherv].edges:
                    #parent_diameter = max( [h.edge_property(tuple((leafotherv,e)), 'diameter') for e in h.vertices[leafotherv].edges] ) 
                #if (h.edge_property(edge, 'diameter')/parent_diameter) < options.rel_diameter_threshold:
                    #if h.vertices[leaf].edges[0] not in leaves:
                        #leaves.append(h.vertices[leaf].edges[0])
                        #for i in h.edge_property(edge, 'intermediaries'):
                            #g.disconnect_vertex(i)
                        #g.disconnect_vertex(leaf)
                        #cnt += 1
            #else:
                #leaves.append(h.vertices[leaf].edges[0])
                #for i in h.edge_property(edge, 'intermediaries'):
                    #g.disconnect_vertex(i)
                #g.disconnect_vertex(leaf)
                #cnt += 1
    #print "Disconnected ", cnt," leaves"
    

    if options.cc_threshold: 
        ccs = vessel_analysis.connected_components(g)
        print "The number of connected components are " , len(ccs)
        print vessel_analysis.size_of_connected_components(g)
        cnt = 0
        for cc in ccs:
            if len(cc) < options.cc_threshold:
                totlen = sum([h.edge_property(e, 'length') for e in cc])
                meandiameter = np.mean([h.edge_property(e, 'diameter') for e in cc])
                print totlen, meandiameter
                if options.cc_length_threshold==0 or totlen<options.cc_length_threshold:
                    if options.cc_diameter_threshold==0 or meandiameter<options.cc_diameter_threshold:
                        cnt += 1
                        while len(cc)>0:
                            edge = cc.pop(0)
                            for i in h.edge_property(edge, 'intermediaries'):
                                g.disconnect_vertex(i)
                            g.disconnect_vertex(edge[0])
                            g.disconnect_vertex(edge[1])
                
        print "Disconnected ", cnt," connected component of size < ",options.cc_threshold 
    
    #db.close()
    
    #graph_analysis.simplify_graph_retaining_intermediaries(g)
    g = vessel_analysis.modified_simplify_graph(g)
    graph_analysis.estimate_edge_lengths(g)
    graph_analysis.estimate_edge_diameters(g)

    if attributes.has_key("history"):
        history = attributes["history"] + "\n" + history
        del attributes['history']

    graph_analysis.output_graph(output, g, history, attributes)
    
    cmd=("\npython /projects/souris/sghanavati/src/scripts/cerebrovascular_analysis/graph2cylinder.py --use_label %s" %(output)) #python /micehome/jgsled/bin/
    os.system(cmd)

    cmd=("\ngraph2obj.py %s %s" %(output,output[:-3]+".obj")) #python /micehome/jgsled/bin/
    os.system(cmd)
    #s = shelve.open(output)
    #s['skeletal_graph'] = g
    #s['vertex_offsets'] = vertex_offsets
    #s['history'] = "processed with nibbler.py script"
    
    ## graph_analysis.output_graph(output, g, "history", {'vertex_offsets':vertex_offsets})
    
    #s.close()
    
    print "Nibbler applied, ",output," updated"
    