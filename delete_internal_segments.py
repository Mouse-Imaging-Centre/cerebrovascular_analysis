# -*- coding: utf-8 -*-
"""
Created on Wed Nov 17 21:29:09 2015

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
import os, commands, time, string, math
import numpy as np
import scipy
import copy
import operator


def distance(u,v):
    dist = math.sqrt( pow(u[0]-v[0],2) + pow(u[1]-v[1],2) + pow(u[2]-v[2],2))
    return dist
    
    
def detect_internal_segments(g, ls, scale = 1.0, offset = 0.0):
    """
    ## Find all edges whose vertices are inside the radius of another edge
    ## In order to find which vertices should be deleted: 
    ##1. We should determine which edge is the parent. 
    ##2. We should delete the vertices enclosed and leave the rest of the vertices intact, 
    ##3. we should reconnect the intact vertices after deletion retaining the topology   
    """
    
    occlusions = {}
    encompassing_e = []
    
    for e1 in g.edge_list():    #check if e1 is inside e2
        if g.edge_property(e1,'label') in ls:
            occluded_vertices = []
            for e2 in g.edge_list():
                if g.edge_property(e2,'label') in ls:
                    #check if e1 and e2 are not adjacent edges and are not the same edge
                    if not e1[0] in [e2[0],e2[1]] and not e1[1] in [e2[0],e2[1]] :
                        occlusions[e1] = []
                        vertex1 = [e1[0]]+g.edge_property(e1,'intermediaries')+[e1[1]]
                        #r1 = g.edge_property(e1,'diameter')
                        #l1 = g.edge_property(e1,'length')
                        vertex2 = [e2[0]]+g.edge_property(e2,'intermediaries')+[e2[1]]
                        r2 = g.edge_property(e2,'diameter')
                        #l2 = g.edge_property(e2,'length')
                        for v1 in vertex1:
                            for v2 in vertex2:
                                if distance(g.vertices[v1].centre,g.vertices[v2].centre)<=r2:
                                    occluded_vertices.append(v1)
                        if len(occluded_vertices) > 0:
                            occlusions[e1].append([e2,occluded_vertices])
                            if e2 not in encompassing_e:
                                encompassing_e.append(e2)
    
    for k in occlusions.keys():
        if len(occlusions[k]) == 0:
            del occlusions[k]
            
    return occlusions, encompassing_e
    
    
    
#def delete_internal_leaves(g, scale = 1.0, offset = 0.0):
    #"""delete_internal_leaves.py is a script to identify and remove terminal branches
#(leaves) from a vessel geometry in cases where the leaf is fully contained within the remaining structure.

#g is graph which is modified in place
#scale is factor on the radii (default = 1.0)
#offset is added to the radii (default = 0.0)
#"""

    #simplified = copy.deepcopy(g)
    #simplify_graph_retaining_intermediaries(simplified)
    
    ## collect all leaves
    #leaves = []
    #for i in xrange(len(g.vertices)):
        #if len(simplified.vertices[i].edges) == 1:
            #edge = (i, simplified.vertices[i].edges[0])  # connection to branch point
            #leaves.append( [i] + simplified.edge_property(edge, 'intermediaries'))

    ## identify connected vertices
    #active_vertices = [i for i in xrange(len(g.vertices)) \
                       #if len(g.vertices[i].edges) > 0]
    ## create reverse lookup for vertices
    #vertex_map = zeros((len(g.vertices)), int_)
    #for i in range(len(active_vertices)):
        #vertex_map[active_vertices[i]] = i

    #centres = zeros((len(active_vertices), 3), float_)
    #radii2 = zeros((len(active_vertices),), float_)
                    
    #for i in xrange(len(active_vertices)):
        #centres[i, :] = g.vertices[active_vertices[i]].centre
        #radii2[i] = (g.vertices[active_vertices[i]].radius*scale+offset)**2

    #report = progress_report(0, len(leaves), "Testing leaf intersections")
    #removable = []
    #for leaf in leaves:
        #for vertex in leaf:
            ## compute distance from every other active vertex
            #diff = sum((centres - array(g.vertices[vertex].centre)[newaxis, :])**2, 1)
            ## compare to radii
            #test_radius = greater_equal(radii2, diff) 
            ## eliminate current leaf from comparison
            #put(test_radius, take(vertex_map, leaf), 0) 
            #if not any(test_radius):  # if current vertex doesn't overlap with any
                                      ## other then leaf is at least partly outside
                                      ## of remaining structure
                #break
        #else:
            #removable.append(leaf)
    
        #report.update_incr()

    ## delete removable leaves
    #for leaf in removable:
        #for vertex in leaf:
            #g.disconnect_vertex(vertex)




#def adjust_junctions(g, vessel_offsets):
    #"""given a vessel graph as produced by tree2graph adjust positions of
#the branch points so as to minimize the length of the connection while
#reducing the angle between the connection vector and the vessel tangent
#and preserving topology"""

    #def unit_vector(u):
        #return u / norm(u)


    #connections = []

    #report = progress_report(0, len(vessel_offsets), "Moving branch points")
    #for offset in vessel_offsets:
        #for junction in g.vertices[offset].edges:
            #if junction < offset and len(g.vertices[junction].edges) == 3:
                ## create a list of candidates to move connection to
                #candidates = [junction]
                ## consider vertices adjacent to junction
                #for neighbor in g.vertices[junction].edges:
                    ## if not the vertex that starts the branch
                    #if neighbor != offset:
                        ## follow along path from this neighbour
                        #prev_vertex = junction
                        #next_vertex = neighbor
                        #while len(g.vertices[next_vertex].edges) == 2:
                            ## add vertex as a candidate
                            #candidates.append(next_vertex)
                            #index = g.vertices[next_vertex].edges[0] == prev_vertex
                            #prev_vertex = next_vertex
                            #next_vertex = g.vertices[next_vertex].edges[index]

                ## pick best candidate
                #ref_cosine = dot(unit_vector(g.vertices[offset].centre -
                   #g.vertices[junction].centre), g.vertices[offset].tangent)
                #cosines = [dot(unit_vector(g.vertices[offset].centre -
                   #g.vertices[vertex].centre), g.vertices[offset].tangent)
                        #for vertex in candidates]   
                #candidates = [v for v, c in zip(candidates, cosines) if
                              #c >= ref_cosine]

                #distances = [norm(g.vertices[offset].centre -
                        #g.vertices[v].centre) for v in candidates]   

                #selection = candidates[argmin(distances)]
                #p = g.edge_properties((offset, junction))
                #g.remove_edge((offset, junction))
                #g.add_edge((selection, offset), p)

        #report.update_incr()



program_name = 'delete_internal_segments.py'


if __name__ == '__main__':
        
    usage = "Usage: "+program_name+" [options] input_graph.db [output_graph.db] \n"+\
            "   or  "+program_name+" --help";
            
    parser = OptionParser(usage)
    
    parser.add_option("--clobber", action="store_true", dest="clobber",
                    default=0, help="overwrite output file")
    parser.add_option("--scale", type="float", dest="scale", default = 1.0,
                       help="scale vessel radii by given factor (default 1.0)")
    parser.add_option("--offset", type="float", dest="offset", default = 0.0,
                       help="add an offset to the vessel radii (default 0.0)")
    parser.add_option("--labels", dest="labels",
                    help="Label numbers to calculate vascular features for, separated by commas, i.e. --labels 3,11,21.",
                    type="string")
        
        
                
    options, args = parser.parse_args()
    print '\n\n>>> %s [%s]: %s\n' % (time.ctime(time.time()), os.getcwd(), string.join(argv))
    print "\n\n",options, "\n\n", args , "\n\n"
   

        
    if len(args) == 1:
        input_file = args[0]
        output_file = input_file[:-3]+"_delinternalsegs.db"
    elif len(args) == 2:
        input_file, output_file = args
    else:
        parser.error("incorrect number of arguments")

    if not options.clobber and os.path.exists(output_file):
        raise SystemExit, \
              "The --clobber option is needed to overwrite an existing file."
        
    ##1. borrow from delete_internal_leaves to find segments that are inside
    
    ##2. borrow from smooth_graph to decide which one is parent and which is daughter branch that needs to be deleted     
    
    ##3. borrow from adjust_junction to find the nearest connection point with preserving the topology 
    
    
    
    
    history = '>>> %s: %s' % (time.ctime(time.time()), string.join(argv))

    ls = []
    if options.labels:
        ls = options.labels.split(',')
        ls = [int(l) for l in ls]
    
    
    g, attributes = graph_analysis.input_graph(input_file, ["history", "vertex_offsets"])
    
    occlusions, encompassing_e = detect_internal_segments(g, ls, options.scale, options.offset)
    
    
    for e in occlusions.keys():
        if e not in encompassing_e:     #this way only edges that have no other edge inside them will be deleted
            g.remove_edge(e)
    
    #NOTE: what we should actually do is to check if all the vertices of an edge are inside another edge => remove that edge
    # if not, then only disconnect_vertex of those inside, then reconnect the rest of the vertices to the encompassing_e

    if attributes.has_key("history"):
        history = attributes["history"] + "\n" + history
        del attributes['history']

    graph_analysis.output_graph(output_file, g, history, attributes)
    
    cmd=("\npython /projects/souris/sghanavati/src/scripts/cerebrovascular_analysis/graph2cylinder.py --use_label %s" %(output_file)) #python /micehome/jgsled/bin/
    os.system(cmd)

    cmd=("\ngraph2obj.py %s %s %s " %(output_file,output_file[:-3]+"obj", options.clobber)) #python /micehome/jgsled/bin/
    os.system(cmd)


    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
