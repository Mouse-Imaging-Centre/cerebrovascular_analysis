#
#  vessel_graph_io.py
#
#  John G. Sled
#  Created: March 1, 2011



""" vessel_graph_io module provides functions for reading and writing
vessel graphs as hdf5 format files.  Vessel graphs are of type
graph.undirected_graph and have vertex annotations and edge annotations
appropriate for representing vessel geometry.  'centre' and 'radius' are
required properties of vertices.  Other properties such as tangent,
contours, and label are optional.


HDF5 vessel_graph file format specification

/

/vessel_graph   : top level group for vessel graph properties

/vessel_graph.history : string with history of processing (optional)

/vessel_graph.description     : string describing data (optional)
/vessel_graph.subject_id      : string identifying the specimen or
                                subject (optional)
/vessel_graph.observation_id  : string identifying the data
                                acquisition (optional)

/vessel_graph.file_format_version   : string in format X.Y.Z indicating
                                vessel_graph file format version

/vessel_graph.distance_unit   : string indicating the unit of all spatial
                                dimensions
                                ["mm" or "um"] (optional)


/vessel_graph.n_vertices : integer -- number of vertices in graph
/vessel_graph.n_edges    : integer -- number of edges in graph

/vessel_graph.vertex_offsets   : array of integer specifying the
                                offsets into the vertex array at which paths
                                start for cases where the graph data
                                structure was originally described as a
                                list of paths


/vessel_graph/edges      : n_edges by 2 unsigned integer array indicating
                            pairs of
                            connected vertices

/vessel_graph/vertex_properties : group containing required and optional
                                properties of each vertex
/vessel_graph/vertex_properties/centre  : n_vertices by 3 float32 array with the
                                coordinates of each vertex
/vessel_graph/vertex_properties/radius  : n_vertices float32 array with the
                                radius at each vertex 
/vessel_graph/vertex_properties/tangent : n_vertices by 3 float32 array with the
                                tangent at each vertex (optional)
/vessel_graph/vertex_properties/contour : n_vertices by 2 by 3 float32 array
                                with the major and minor axes of an elliptic
                                contour around each vertex (optional)
/vessel_graph/vertex_properties/label   : n_vertices unsigned integer array with the
                                label ID of each vertex (optional) 

/vessel_graph/edge_properties           : group containing additional optional
                                properties of each edge
/vessel_graph/edge_properties/length    : n_edges float32 array with length of
                                each edge (optional)
/vessel_graph/edge_properties/diameter  : n_edges float32 array with diameter of
                                each edge (optional)

"""


from numpy import *
import tables
from morphology import graph


# file formal version string
vessel_graph_format_version = "1.0.1"

# top level attributes
known_attributes = ["history", "description", "subject_id", "observation_id",
                    "distance_unit", "vertex_offsets"]

# known properties and associated data types
known_vertex_properties = { "radius" :  (float32, (1,),
                                        "radius at each vertex"),
                            "centre" :  (float32, (3,),
                                        "coordinates of each vertex"),
                            "tangent" : (float32, (3,),
                                        "tangent vector at each vertex"),
                            "contour" : (float32, (2,3),
                                        "major and minor axes of an elliptic contour around each vertex"),
                                        
                            "label" :    (uint32, (1,),
                                        "label ID of each vertex") }
                            
known_edge_properties =   { "length"   : (float32, (1,),
                                        "length of each edge"),
                            "diameter" : (float32, (1,),
                                        "diameter of each edge"),
                            "intermediaries" : (tables.Int32Atom(shape=()), None,
                                                "intemediate vertices comprising a vessel segment"),
                            "cylinder_intermediaries" : (tables.Int32Atom(shape=()), None,
                                                "intemediate vertices comprising a vessel segment"),
                            "cyl_segment_ID" : (tables.Int32Atom(shape=()), None,
                                            "numeric ID for each edge"),
                            "cyl_centreX" : (tables.Float32Atom(shape=()), None,
                                            "X coordinates of centre of simplified segment cylinders"),
                            "cyl_centreY" : (tables.Float32Atom(shape=()), None,
                                            "Y coordinates of centre of simplified segment cylinders"),
                            "cyl_centreZ" : (tables.Float32Atom(shape=()), None,
                                            "Z coordinates of centre of simplified segment cylinders"),
                            "cyl_tangentX" : (tables.Float32Atom(shape=()), None,
                                            "X part of tangent of simplified segment cylinders"),
                            "cyl_tangentY" : (tables.Float32Atom(shape=()), None,
                                            "Y part of tangent of simplified segment cylinders"),
                            "cyl_tangentZ" : (tables.Float32Atom(shape=()), None,
                                            "Z part of tangent of simplified segment cylinders"),
                            "cyl_height" : (tables.Float32Atom(shape=()), None,
                                            "length of simplified segment cylinders"),
                            "cyl_radius" : (tables.Float32Atom(shape=()), None,
                                            "radius of simplified segment cylinders"),
                            "label" : (uint, (1,),
                                            "manually assigned label for vessel segment"),
                            "proximity" : (tables.Float32Atom(shape=(2,)), None,
                                            "pairs of region ID from MR atlas and distances where the distance"+\
                                            "is between edge centre and region centre"),
                            "anglewref" : (tables.Float32Atom(shape=(2,)), None,
                                            "pairs of label in vascalar atlas and the angle made with"+\
                                            " the reference direction"),
                            "estimated_label" : (tables.Float32Atom(shape=(2,)), None,
                                            "pairs of automatically assigned label and posterior probability"),
}

required_vertex_properties = ["radius", "centre"]


def output_vessel_graph(g, filename, history = None,
                        extra_attrs = { "distance_unit" : "mm" },
                        additional_vertex_properties = [],
                        additional_edge_properties = [],
                        guess_unknown_property_format = False):
    """                  
    g                : undirected graph
    filename  : new HDF5 file to create
    history          : string
    extra_attrs      : dictionary with optional keys as follows
        vertex_offsets   : array of integer offsets
        description      : string
        subject_id       : string
        observation_id   : string
        distance_unit    : string (default: "mm")
    additional_vertex_properties : list of vertex property names
    additional_edge_properties   : list of edge property names
    guess_unknown_property_format : if true, attempt guess the correct
                                numerical format of unrecognized
                                edge and vertex properties
    """

    n_vertices = len(g.vertices)
    n_edges = len(g.edge_list())

    d = tables.openFile(filename, mode = "w", title="vessel_graph")

    vg = d.createGroup(d.root, "vessel_graph",
                        "top level group for vessel graph properties")

    # create attributes of top level group
    if history is not None:
        vg._v_attrs.history = history

    # add top level attributes
    for name in extra_attrs.keys():
        if name == "vertex_offsets":
            d.createArray(vg, name, array(extra_attrs[name]),
                        "offsets into vertex list of each path")
        elif name in known_attributes and name != "history":
            setattr(vg._v_attrs, name, extra_attrs[name])
        else:
            print "Warning: unknown attribute %s in output_vessel_graph" % \
                name
            
    vg._v_attrs.file_format_version = vessel_graph_format_version

    vg._v_attrs.n_vertices = n_vertices
    vg._v_attrs.n_edges = n_edges

    # create vertex properties
    vertex_properties = d.createGroup(vg, "vertex_properties",
        "group containing required and optional properties of each vertex")
    edge_properties = d.createGroup(vg, "edge_properties",
        "group containing additional optional properties of each edge")

    d.createArray(vg, "edges", array(g.edge_list()),
        "n_edges by 2 integer array indicating pairs of connected vertices")
    
    # add properties of vertices
    for property in required_vertex_properties + additional_vertex_properties:
        if property in known_vertex_properties:
            # create array of required data type
            ptype, pshape, description = known_vertex_properties[property]
            values = zeros( (n_vertices,) + pshape, ptype )

            for vertex_id in xrange(len(g.vertices)):
                if hasattr(g.vertices[vertex_id], property):
                    values[vertex_id, ...] = getattr(g.vertices[vertex_id],
                                            property)
            # put vertex property array into hdf5 file
            d.createArray(vertex_properties, property, values, description)

        elif guess_unknown_property_format:
            values = []
            pshape = None
            for vertex in g.vertices:
                if hasattr(vertex, property):
                    value = array(getattr(vertex, property))
                    values.append(value)
                    if pshape is None:
                        pshape = value.shape
                    elif pshape != value.shape:
                        raise ValueError, \
                            "vertex property %s is not of uniform dimension" \
                            % property
                # else value is missing
                else:
                    values.append(None)
            # fill in missing values with zeros
            for vertex_id in range(len(values)):
                if values[vertex_id] is None:
                    values[vertex_id] = zeros(pshape, float_)
            #create vertex property array
            d.createArray(vertex_properties, property, array(values))
            
        else:
            raise ValueError, "unknown vertex property: %s" % property


    # add properties of edges
    for property in additional_edge_properties:
        if property in known_edge_properties:
            edges = g.edge_list()

            # create array of required data type
            ptype, pshape, description = known_edge_properties[property]
            if isinstance(ptype, tables.Atom):
                values = d.createVLArray(edge_properties, property, ptype, description)
                for edge_id in xrange(len(edges)):
                    values.append(g.edge_property(edges[edge_id], property))

            else:
                values = zeros( (n_edges,) + pshape, ptype )
                for edge_id in xrange(len(edges)):
                    values[edge_id, ...] = g.edge_property(edges[edge_id], property)
                # put edge property array into hdf5 file
                d.createArray(edge_properties, property, values, description)

        elif guess_unknown_property_format:
            values = []
            pshape = None
            edges = g.edge_list()
            for edge_id in xrange(len(edges)):
                value = array(g.edge_property(edges[edge_id], property))
                values.append(value)
                if pshape is None:
                    pshape = value.shape
                elif pshape != value.shape:
                    raise ValueError, \
                        "edge property %s is not of uniform dimension" \
                        % property

            #create edge property array
            d.createArray(edge_properties, property, array(values))

        else:
            raise ValueError, "unknown edge property: %s" % property


    # finish writing hdf5 file
    d.close()


def input_vessel_graph(filename, extra_attrs = known_attributes):
    """read hdf5 vessel_graph formated and return undirected_graph.  Additional
top level attributes can be requested by name using the extra_attrs options.

    returns:
        g       :   undirected_graph
        extras  :   dictionary of additional attributes
"""
    d = tables.openFile(filename, mode = "r")
    print "input_vessel_graph"

    if not is_known_file_format(d):
        d.close()
        raise ValueError, "unrecognized file format in %s" % filename

    vg = d.root.vessel_graph

    g = graph.undirected_graph(vg._v_attrs.n_vertices)
    g_extras = dict([(attribute,getattr(vg._v_attrs, attribute)) \
                for attribute in extra_attrs if hasattr(vg._v_attrs, attribute)])
    if hasattr(vg, "vertex_offsets"):
        print "vertex_offsets"
        g_extras["vertex_offsets"] = vg.vertex_offsets[:]

    # populate edge list and properties
    edge_property_names = vg.edge_properties._v_children.keys()
    for edge_index in xrange(vg._v_attrs.n_edges):
        properties = dict([(name , getattr(vg.edge_properties, name)[edge_index]) \
                    for name in edge_property_names ])
        g.add_edge(vg.edges[edge_index], properties)
        
    if hasattr(vg, "new_edges"):
        print "new_edges"
        ne = vg.new_edges[:]
        for e in ne:
            g.add_edge(e)
        

    # add properties to vertices
    vertex_property_names = vg.vertex_properties._v_children.keys()
    for vertex_index in xrange(vg._v_attrs.n_vertices):
        for name in vertex_property_names:
            setattr(g.vertices[vertex_index], name,
                    getattr(vg.vertex_properties, name)[vertex_index])

    d.close()

    return g, g_extras

def is_known_file_format(d):
    # note: should consider checking format version number e.g. 
    #   d.root.vessel_graph._v_attrs.file_format_version == \
    #           vessel_graph_format_version
    return hasattr(d.root, "vessel_graph")

def input_history(filename):
    d = tables.openFile(filename, mode = "r")

    if not is_known_file_format(d):
        d.close()
        raise ValueError, "unrecognized file format in %s" % filename

    vg = d.root.vessel_graph

    if hasattr(vg._v_attrs, "history"):
        history = vg._v_attrs.history
    else:
        history = None

    d.close()
    return history

def is_vessel_graph(filename):
    d = tables.openFile(filename, mode = "r")
    known = is_known_file_format(d)
    d.close()

    return known


def find_additional_vertex_properties(g):
    """returns a list of vertex properties found in g that are not required
vertex properties"""

    base_properties = dir(graph.undirected_graph_vertex()) + \
                    required_vertex_properties

    additional_properties = []

    for vertex in g.vertices:
        for p in dir(vertex):
            if (p not in base_properties) and (p not in additional_properties):
                additional_properties.append(p)

    return additional_properties

if __name__ == '__main__':
    "Testing vessel_graph_io.py"

    g = graph.undirected_graph(7)
    g.add_edge((1,3))
    g.add_edge((2,2))
    g.add_edge((3,3))
    g.add_edge((1,2))
    g.add_edge((4,2))
    g.add_edge((0,5))
    g.add_edge((4,5))
    g.add_edge((6,0))
    g.add_edge((6,2))

    for vertex_id in range(len(g.vertices)):
        g.vertices[vertex_id].radius = float(vertex_id)
        g.vertices[vertex_id].centre = float(vertex_id)**2*array([1,1,1])
        g.vertices[vertex_id].tangent = array([1,1,1], float_)
        g.vertices[vertex_id].colour = array([0.5,0.8,0.7], float_)
        
    for edge in g.edge_list():
        g.set_edge_property(edge, "length", 1.0)
        g.set_edge_property(edge, "tortuosity", array([1,2]))

    vertex_offsets = [0,3,4]

    output_vessel_graph(g, "/tmp/test_vessel_graph.h5",
                    history = "testing vessel graph",
                    extra_attrs = { "description" : "arbitrary data",
                                "subject_id": "100",
                                "observation_id" : "100",
                                "vertex_offsets" : vertex_offsets}, 
                    additional_vertex_properties = ["tangent", "colour"],
                    additional_edge_properties = ["length", "tortuosity"],
                    guess_unknown_property_format = True)
                    

    h, h_extras = \
        input_vessel_graph("/tmp/test_vessel_graph.h5",
                        ["history", "description", "subject_id",
                        "observation_id", "vertex_offsets"])

    print h_extras
    print h.vertices
    print h.edge_list()
    print [h.edge_properties(edge) for edge in h.edge_list()]
    
    for vertex_id in range(len(h.vertices)):
        vertex = h.vertices[vertex_id]
        print "id: %u, centre: %s, radius: %s, tangent: %s, colour: %s" % \
            (vertex_id, str(vertex.centre), str(vertex.radius),
                str(vertex.tangent), str(vertex.colour))

    print "additional vertex properties: %s" % \
            str(find_additional_vertex_properties(h))
