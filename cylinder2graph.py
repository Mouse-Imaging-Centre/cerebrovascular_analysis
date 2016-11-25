#!/usr/bin/env python

#add cylinder objects to vessel segments as an attribute to the *.db file hash table and save it in *_cylinder.db
#
#
#  Created May 22, 2012
#  Last modified July 2, 2013
#  Sahar Ghanavati
 
#connect to bianca (system with hardy os 
#export PYTHONPATH=/home/jgsled/lib64_ubuntu_hardy/python:$PYTHONPATH
#python graph2cylinder.py --help
#python graph2cylinder.py ../tree_mask_blur_delleaves.db ../outputname --clobber				

from vessel_tracking import graph_analysis, path_io, path_analysis, curvature, filters,\
     vessel_tracking, medial_atoms, path_visualize, tubes
from morphology import graph, cluster_skeleton, object_io
from numpy import *
from optparse import OptionParser, Option, OptionValueError
from sys import argv
import os, shelve, string ,time,sys
import commands
#from time import time, ctime
import copy
import numpy
##from Multipack import leastsq
#from morphology import cluster_division_opt, voxel_code, \
     #ray_trace, histogram, cluster_skeleton_opt, \
     #connected_components, disjoint_sets, shortest_paths,\
     #histogram
#from priority_queue import PriorityQueue
#import medial_atoms

#---------------------------------------------------------------------------------
#
def this_vessel_path_length(g, edge):
    #"estimate vessel lengths by summing the euclidian distance between intermediate points"
    edge = g.order(edge)
    edge_length = lambda g, edge: \
        linalg.norm(g.vertices[edge[0]].centre - g.vertices[edge[1]].centre)
    points = [edge[0]] + list(g.edge_property(edge, 'intermediaries')) + [edge[1]]
    segments = [(points[i], points[i+1]) for i in range(len(points)-1)]
    return sum([edge_length(g,segment) for segment in segments])


def this_vessel_avg_diameter(g, edge):
    #"estimate vessel diameter by averaging the diameters at intermediate points"
    diameters = [g.vertices[i].radius*2 for i in \
                [edge[0]] + list(g.edge_property(edge, 'intermediaries')) \
                + [edge[1]]]
    return mean(diameters)


def this_estimate_edge_lengths(g, length_fcn=this_vessel_path_length, property='length'):
    for edge in g.edge_list():
        #print edge, " ",
        length = length_fcn(g, edge)
        #print length
        g.set_edge_property(edge, property, numpy.float32(length))


def this_estimate_edge_diameters(g, diameter_fcn=this_vessel_avg_diameter, property='diameter'):
    for edge in g.edge_list():
        #print edge," ",max([edge[0]] + list(g.edge_property(edge, 'intermediaries')) + [edge[1]]), " / ",
        diameter = diameter_fcn(g, edge)
        #print diameter
        g.set_edge_property(edge, property, numpy.float32(diameter))






#simplify the already simplified graph made from cylinder object, based on cyl_labels. 
#The original graph was simplified and converted to cylinder, however, edges are added and deleted in brainview 
#=> needs to do simplified graph using cyl_labels one more time after conversion from cylinder
def modified_simplify_graph(g):

    graph_analysis.remove_self_loops(g)

    for edge in g.edge_list():
        if 'intermediaries' not in g.edge_properties(edge).keys():
            g.set_edge_property(edge, 'intermediaries', [])
        
        
    vertex_list = unique([e[0] for e in g.edge_list()] + [e[1] for e in g.edge_list()])	

    for i in vertex_list:
        if (len(g.vertices[i].edges) == 2):
            if not g.has_edge(tuple(g.vertices[i].edges)):
                new_e = g.order(tuple(g.vertices[i].edges))
                e1 = g.order(tuple((new_e[0], i)))
                e2 = g.order(tuple((i, new_e[1])))
                if len(g.edge_properties(e2).keys())> len(g.edge_properties(e1).keys()):
                    tmp=e2
                    e2=e1
                    e1=tmp
                #print i, " ", new_e, " ", e1, " ", e2
                
                properties={}
                if ('cyl_label' in g.edge_properties(e1).keys()) and ('cyl_label' in g.edge_properties(e2).keys()):
                    g.set_edge_property(e1,'cyl_label',int(g.edge_property(e1,'cyl_label')))	#make sure cyl_label is int
                    g.set_edge_property(e2,'cyl_label',int(g.edge_property(e2,'cyl_label')))
                    #print g.edge_property(e1,'cyl_label') ," ", g.edge_property(e2,'cyl_label')
                    if (g.edge_property(e1,'cyl_label') == g.edge_property(e2,'cyl_label')):
                        properties={'cyl_label': g.edge_property(e1,'cyl_label')}
                    else: 	#if the 2 adjacent edges that should be merged into 1 edge have different cyl_labels then get the label of the longer edge_length
                        if (this_vessel_path_length(g, e1) > this_vessel_path_length(g, e2)):
                            properties={'cyl_label': g.edge_property(e1,'cyl_label')}
                        else:
                            properties={'cyl_label': g.edge_property(e2,'cyl_label')}
                            
                if ('intermediaries' in g.edge_properties(e1).keys()) and ('intermediaries' in g.edge_properties(e2).keys()):
                    first = list(copy.copy(g.edge_property(e1, 'intermediaries')))
                    second = list(copy.copy(g.edge_property(e2, 'intermediaries')))
                    if i < new_e[0]:
                        first.reverse()
                    if i > new_e[1]:
                        second.reverse()
                    intermediaries = first + [i] + second                
                    properties['intermediaries']=intermediaries
                    
                #print g.edge_properties(e1).keys()
                #print g.edge_properties(e2).keys()
                #for n in g.edge_properties(e2).keys():
                    #print n ,":", g.edge_property(e2,n)
                    #print n ,":", g.edge_property(e1,n)
                #print "\n\n"	
                for name in g.edge_properties(e2).keys():
                    if not name=='intermediaries' and not name=='cyl_label':
                        if type(g.edge_property(e1,name))==list or type(g.edge_property(e1,name))==array:
                            print name ,":", g.edge_property(e1,name), " , " , g.edge_property(e2, name)
                        else:
                            properties[name] = []
                            properties[name].append(g.edge_property(e1, name))
                            properties[name].append(g.edge_property(e2, name))	
                g.add_edge(new_e, properties)
                g.disconnect_vertex(i)

    #print len(g.edge_list())		
    return g



        
 
program_name = 'cylinder2graph.py'


if __name__ == '__main__':
        
    usage = "Usage: "+program_name+" [options] labelled.h5 ref_graph.db mr_atlas_file.mnc mr_atlas_centroids.db cba_direction_reference.db [output_file.db]\n"+\
            "   or  "+program_name+" --help";
            
    parser = OptionParser(usage)

    parser.add_option("--clobber", action="store_true", dest="clobber",
                        default=0, help="overwrite output file")
    #parser.add_option("--use_label", action="store_true", dest="use_label",
                        #default=0, help="use labels specified in label edge_property of the input graph. NOTE: input graph needs to have edge_property label, otherwise 0 will be used as the cylinder label!")
    parser.add_option("--tangent_threshold",type="float", dest="tangent_threshold",
                        default=0.174, help="tangent threshold in radian to make cylinders. Defulat: 0.174 radian = 10 degrees")                  
    parser.add_option("--use_estimated_label", action="store_true", dest="use_estimated_label",
                        default=0, help="use labels specified in estimated_label edge_property of the input graph, which is calculated from automatic labeling. ")      
    parser.add_option("--use_cyl_label", action="store_true", dest="use_cyl_label",
                        default=0, help="use labels specified in cyl_label edge_property of the input graph through Brainview.")

    options, args = parser.parse_args()

    print '>>> %s [%s]: %s\n' % (time.ctime(time.time()), os.getcwd(), string.join(argv))
    sys.stdout.flush()

    if len(args) == 6:
        input_file, ref_graph, mr_atlas_file, mr_centroids_file, reference_file, output_file = args
    elif len(args) ==5:
        input_file, ref_graph, mr_atlas_file, mr_centroids_file, reference_file = args
        output_file = (input_file[:-2]+ "db")
    else:	
        parser.error("incorrect number of arguments")


    if not options.clobber and os.path.exists(output_file):
        raise SystemExit, \
            "The --clobber option is needed to overwrite an existing file."

    #if options.use_label and options.use_estimated_label: 
        #raise SystemExit, \
            #"Only one of the options --use_label and --use_estimated_label can be used."
            
    #if options.use_label and options.use_cyl_label: 
        #raise SystemExit, \
            #"Only one of the options --use_label and --use_cyl_label can be used."
            
    if options.use_cyl_label and options.use_estimated_label: 
        raise SystemExit, \
            "Only one of the options --use_cyl_label and --use_estimated_label can be used."


    #cmd=("python /projects/souris/sghanavati/src/scripts/graph2graph_v2.py %s %s --clobber" %(input_file,input_file[:-3]+"graph2graph.db"))	#python /micehome/jgsled/bin/
    cmd=("python /projects/souris/sghanavati/src/scripts/cerebrovascular_analysis/graph2graph_v2.py %s %s " %(input_file,input_file[:-3]+"graph2graph.db")) #python /micehome/jgsled/bin/
    print(cmd)
    sys.stdout.flush()
    os.system(cmd)	

    cmd=("python /projects/souris/sghanavati/src/scripts/cerebrovascular_analysis/modify_db_intermediaries.py %s %s %s" %(input_file[:-3]+"graph2graph.db",ref_graph,output_file[:-3]+"_nofeature.db" ))	#python /micehome/jgsled/bin/
    print(cmd)
    sys.stdout.flush()
    os.system(cmd)	

    #output_file = output_file[:-3]+"_intm.db"

    try:
        g, attributes = graph_analysis.input_graph(output_file[:-3]+"_nofeature.db", ["history", "vertex_offsets"])  	# open graph data and copy contents to memory
        #g = graph_analysis.input_graph("tree_mask_blur_delleaves.db")
        print ("Succefully read in the %s\n" %output_file[:-3]+"_nofeature.db")
        sys.stdout.flush()
    except:
        print("Error reading in the %s\n" %output_file[:-3]+"_nofeature.db")
        sys.stdout.flush()

    #redo the simplified graph and feature extraction for new edges!	
    print len(g.edge_list()), " ", len(g.vertices),"\nsimplifying edges:"
    sys.stdout.flush()
    #for e in g.edge_list():
        #print e, " ",
    #print " "	
    g = modified_simplify_graph(g)
    print len(g.edge_list()), " ", len(g.vertices)
    sys.stdout.flush()
    #for e in g.edge_list():
        #print e, " ",g.edge_property(e, 'intermediaries')
    #print " "	
    this_estimate_edge_diameters(g)
    this_estimate_edge_lengths(g)


    history = '\n>>> %s: %s' % (time.ctime(time.time()), string.join(argv)+"\nmodified_simplify_graph")		##argv = ['code.py', 'input.db', 'output']
    if attributes.has_key("history"):
        history = attributes["history"] + "\n" + history
        del attributes['history']
        
    graph_analysis.output_graph(output_file[:-3]+"_nofeature_simplified.db", g, history, attributes)   

    print ("Succefully wrote the %s\n" %output_file[:-3]+"_nofeature_simplified.db")
    sys.stdout.flush()
    #output_file = output_file[:-3]+"_simp.db"


    #redo feature extraction
    cmnd = ("python /projects/souris/sghanavati/src/scripts/cerebrovascular_analysis/feature_extraction.py %s %s %s %s %s" %(output_file[:-3]+"_nofeature_simplified.db", mr_atlas_file, mr_centroids_file, reference_file,output_file[:-3]+"_featured.db"))		#--new_edges		
    print  "\n", cmnd, "\n"
    sys.stdout.flush()
    os.system(cmnd)
    #output_file=output_file[:-3]+"_features.db"

    #output_file = output_file[:-3]+"_f.db"



    try:
        g, attributes = graph_analysis.input_graph(output_file[:-3]+"_featured.db", ["history", "vertex_offsets"])  	# open graph data and copy contents to memory
        #g = graph_analysis.input_graph("tree_mask_blur_delleaves.db")
        print ("Succefully read in the %s\n" %output_file[:-3]+"_featured.db")
        sys.stdout.flush()
    except:
        print("Error reading in the %s\n" %output_file[:-3]+"_featured.db")
        sys.stdout.flush()

        
    ###borrowed from graph2cylinder code 
    tang_thresh= options.tangent_threshold   
    for e in g.edge_list():
        centreX=[]
        centreY=[]
        centreZ=[]
        tangentX=[]
        tangentY=[]
        tangentZ=[]
        radius=[]
        height=[]
        ##segment=[]
        vertex_list = [e[0]]+ list(g.edge_property(e, 'intermediaries'))+ [e[1]]
        #print e, g.edge_property(e,'intermediaries'), vertex_list
        if len(vertex_list)<2:
            print "ERROR: length of vertex_list is less than 2! ", e, " " , vertex_list, "\nAborted!"
            exit (0)
        for xx in vertex_list:
            if (xx> (len(g.vertices)-1)):
                print "ERROR: ",len(g.vertices), " ", e, " " , vertex_list, "\nAborted!"
                exit (0)
        v0 = vertex_list[0]
        v1 = vertex_list[1]
        r0 = ( g.vertices[vertex_list[0]].radius + g.vertices[vertex_list[1]].radius )/2
        t0= ( g.vertices[vertex_list[1]].centre - g.vertices[vertex_list[0]].centre )
        tang0 = t0/linalg.norm(t0)
        for i in range(1,len(vertex_list)-1):
            r1 = ( g.vertices[vertex_list[i]].radius + g.vertices[vertex_list[i+1]].radius )/2
            t1= ( g.vertices[vertex_list[i+1]].centre - g.vertices[vertex_list[i]].centre )
            tang1 = t1/linalg.norm(t1)
            innerprod= tang0[0]*tang1[0] + tang0[1]*tang1[1] + tang0[2]*tang1[2]
            
            if (innerprod > 1.0):
                print "ERROR: innerprod > 1.0 : v0=", g.vertices[v0].centre, " v1=", g.vertices[v1].centre, " v2=", g.vertices[vertex_list[i+1]].centre, " t0=",t0, " t1=", t1, " tang0=",tang0, " tang1=", tang1, " innerprod=", innerprod
                innerprod=1.0
            if (innerprod < -1.0):
                print "ERROR: innerprod < -1.0 : v0=", g.vertices[v0].centre, " v1=", g.vertices[v1].centre, " v2=", g.vertices[vertex_list[i+1]].centre, " t0=",t0, " t1=", t1, " tang0=",tang0, " tang1=", tang1, " innerprod=", innerprod
                innerprod = -1.0
            
            angle = math.acos (innerprod ) #angle=acos(t1.t2)#in radian #axis= norm(t1 x t2)    
            if ( (abs(r0-r1) < r0/2) and (abs(angle)< tang_thresh) ):   #|r1-r0|<0.5r0=> 0.5r0<r1<1.5r0  &  3*0.17 #10 degrees= 0.174   => merge them into 1 cylinder           
                r0 = (r0+r1)/2
                t0 = t1+t0
                tang0 = t0/linalg.norm(t0)
                v1 = vertex_list[i+1]

            else:   # => make a cylinder obj 
                centre = (g.vertices[v0].centre + g.vertices[v1].centre )/2
                centreX.append(numpy.float32(centre[0]))
                centreY.append(numpy.float32(centre[1]))
                centreZ.append(numpy.float32(centre[2]))
                tangentX.append(numpy.float32(tang0[0]))        #t = g.vertices[v0].centre - g.vertices[v1].centre )    #tang = t/linalg.norm(t)
                tangentY.append(numpy.float32(tang0[1]))
                tangentZ.append(numpy.float32(tang0[2]))
                radius.append(numpy.float32(r0))                #radius = ( g.vertices[v0].radius + g.vertices[v1].radius )/2
                height.append(numpy.float32(linalg.norm(t0)))   #height = linalg.norm(t)
                r0 = r1
                t0 = t1
                tang0 = t0/linalg.norm(t0)
                v0 = vertex_list[i]
                v1 = vertex_list[i+1]
                ##segment.append(segment_ID)
        ##The last part of the edge is between v0 and v1 => need to make a cylinder from it     
        centre = (g.vertices[v0].centre + g.vertices[v1].centre )/2
        centreX.append(numpy.float32(centre[0]))
        centreY.append(numpy.float32(centre[1]))
        centreZ.append(numpy.float32(centre[2]))
        t = (g.vertices[v0].centre - g.vertices[v1].centre )    
        tang = t/linalg.norm(t)
        tangentX.append(numpy.float32(tang[0]))     #
        tangentY.append(numpy.float32(tang[1]))
        tangentZ.append(numpy.float32(tang[2]))
        radius.append( numpy.float32((float(g.vertices[v0].radius) + float(g.vertices[v1].radius) )/2))
        height.append(numpy.float32(linalg.norm(t)))    
        ##segment.append(int(segment_ID))
        
        g.set_edge_property(e,'cyl_centreX', centreX)       #graph.set_edge_property(self,pair,key,value)
        g.set_edge_property(e,'cyl_centreY', centreY)
        g.set_edge_property(e,'cyl_centreZ', centreZ)
        g.set_edge_property(e,'cyl_tangentX', tangentX)  
        g.set_edge_property(e,'cyl_tangentY', tangentY)  
        g.set_edge_property(e,'cyl_tangentZ', tangentZ)  
        g.set_edge_property(e,'cyl_radius', radius)
        g.set_edge_property(e,'cyl_height', height)

        
    #based on use_*_label option => change label property
    for e in g.edge_list():
        if options.use_cyl_label and 'cyl_label' in g.edge_properties(e).keys():	##use the same label as the current edge_property  == making no change
            g.set_edge_property(e,'label', int(g.edge_property(e,'cyl_label')) )
        elif options.use_estimated_label and 'estimated_label' in g.edge_properties(e).keys(): ##use the label in the estimated_label edge_property
            g.set_edge_property(e,'label', int(g.edge_property(e,'estimated_label')[0][0]) )
        elif 'label' not in g.edge_properties(e).keys():
            g.set_edge_property(e,'label', int(0))  #"No label                \0" #add label dataset to edge properties and initialize it to None with null-terminted character '\0'

        
        
    history = '>>> %s: %s' % (time.ctime(time.time()), string.join(argv))		##argv = ['code.py', 'input.db', 'output']
    if attributes.has_key("history"):
        history = attributes["history"] + "\n" + history
        del attributes['history']
        
    graph_analysis.output_graph(output_file, g, history, attributes)   

    print ("Succefully wrote the brainview labelled %s\n" %output_file)
    sys.stdout.flush()

    cmd = ("rm %s %s %s %s" %(output_file[:-3]+"_featured.db",output_file[:-3]+"_nofeature.db ", output_file[:-3]+"_nofeature_simplified.db",input_file[:-3]+"graph2graph.db"))
    print  "\n", cmd, "\n"
    sys.stdout.flush()
    os.system(cmd)

        
        
        

