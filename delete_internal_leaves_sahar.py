#!/usr/bin/python
#
#
#  delete_internal_leaves_sahar.py
#
#
#  Created May 28, 2010
#  John G. Sled
#
#  Modified by Sahar Ghanavati  

from optparse import OptionParser, Option, OptionValueError
import shelve, os, string
from vessel_tracking import graph_analysis
from sys import argv
from time import time, ctime


from numpy import *
import py_minc
import os, string
from scipy import mean
from morphology import cluster_division_opt, voxel_code, graph, \
     cluster_skeleton, ray_trace, histogram, cluster_skeleton_opt, \
     object_io, connected_components, disjoint_sets, shortest_paths,\
     histogram
from minc_util.progress import progress_report
import shelve
import copy, tempfile

def remove_self_loops(g):
    for i in range(len(g.vertices)):
        if g.has_edge((i,i)):
            g.remove_edge((i,i))

def simplify_graph_retaining_intermediaries(g):

    remove_self_loops(g)

    for edge in g.edge_list():
        g.set_edge_property(edge, 'intermediaries', [])

    for i in range(len(g.vertices)):
        if len(g.vertices[i].edges) == 2 and not g.has_edge(tuple(g.vertices[i].edges)):
            n = g.order(tuple(g.vertices[i].edges))
            first = copy.copy(g.edge_property((n[0], i), 'intermediaries'))
            second = copy.copy(g.edge_property((i, n[1]), 'intermediaries'))
            if i < n[0]:
                first.reverse()
            if i > n[1]:
                second.reverse()
            intermediaries = first + [i] + second
            g.add_edge(n, {'intermediaries': intermediaries})
            g.disconnect_vertex(i)


def my_delete_internal_leaves(g, scale = 1.0, offset = 0.0):
    """delete_internal_leaves.py is a script to identify and remove terminal branches
(leaves) from a vessel geometry in cases where the leaf is fully contained within the remaining structure.

g is graph which is modified in place
scale is factor on the radii (default = 1.0)
offset is added to the radii (default = 0.0)

return number of vertices deleted
"""

    simplified = copy.deepcopy(g)
    simplify_graph_retaining_intermediaries(simplified)

    # collect all leaves
    leaves = []
    for i in xrange(len(g.vertices)):
        if len(simplified.vertices[i].edges) == 1:
            edge = (i, simplified.vertices[i].edges[0])  # connection to branch point
            leaves.append( [i] + simplified.edge_property(edge, 'intermediaries'))

    # identify connected vertices
    active_vertices = [i for i in xrange(len(g.vertices)) \
                       if len(g.vertices[i].edges) > 0]
    # create reverse lookup for vertices
    vertex_map = zeros((len(g.vertices)), int_)
    for i in range(len(active_vertices)):
        vertex_map[active_vertices[i]] = i

    centres = zeros((len(active_vertices), 3), float_)
    radii2 = zeros((len(active_vertices),), float_)

    for i in xrange(len(active_vertices)):
        centres[i, :] = g.vertices[active_vertices[i]].centre
        radii2[i] = (g.vertices[active_vertices[i]].radius*scale+offset)**2

    report = progress_report(0, len(leaves), "Testing leaf intersections")
    removable = []
    for leaf in leaves:
        for vertex in leaf:
            # compute distance from every other active vertex
            diff = sum((centres - array(g.vertices[vertex].centre)[newaxis, :])**2, 1)
            # compare to radii
            test_radius = greater_equal(radii2, diff)
            # eliminate current leaf from comparison
            put(test_radius, take(vertex_map, leaf), 0)
            if not any(test_radius):  # if current vertex doesn't overlap with any
                                      # other then leaf is at least partly outside
                                      # of remaining structure
                break
        else:
            removable.append(leaf)

        report.update_incr()

    # delete removable leaves
    for leaf in removable:
        for vertex in leaf:
            g.disconnect_vertex(vertex)

    return len(removable)

#---------------------------------------------------------------------------------
#

"""delete_internal_leaves.py is a script to identify and remove terminal branches
(leaves) from a vessel geometry in cases where the leaf is fully contained within the remaining structure.
"""


program_name = 'delete_internal_leaves_sahar.py'


if __name__ == '__main__':
        
    usage = "Usage: "+program_name+" input.db output.db\n"+\
            "   or  "+program_name+" --help"
    parser = OptionParser(usage)
    parser.add_option("--clobber", action="store_true", dest="clobber",
                       default=0, help="overwrite output file")
    parser.add_option("--scale", type="float", dest="scale", default = 1.0,
                       help="scale vessel radii by given factor (default 1.0)")
    parser.add_option("--offset", type="float", dest="offset", default = 0.0,
                       help="add an offset to the vessel radii (default 0.0)")

    options, args = parser.parse_args()
    if len(args) != 2:
        parser.error("incorrect number of arguments")

    input, output = args

    if not options.clobber and os.path.exists(output):
        raise SystemExit, \
              "The --clobber option is needed to overwrite an existing file."

    history = '>>> %s: %s' % (ctime(time()), string.join(argv))

    g, attributes = graph_analysis.input_graph(input, ["history", "vertex_offsets"])
    
    while True:
        num_del = my_delete_internal_leaves(g, options.scale, options.offset)
	print num_del, " edges has been deleted in this iteration"
	if num_del==0:
	    break

    if attributes.has_key("history"):
        history = attributes["history"] + "\n" + history
        del attributes['history']

    graph_analysis.output_graph(output, g, history, attributes)




