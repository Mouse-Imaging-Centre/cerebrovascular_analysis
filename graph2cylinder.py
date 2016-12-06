#!/usr/bin/env python

#add cylinder objects to vessel segments as an attribute to the *.db file hash table and save it in *_cylinder.db
#
#
#  Created Mar 24, 2011
#  Last modified April 24, 2012
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
import os, shelve, string,sys
import commands
from time import time, ctime
import copy
import numpy
#---------------------------------------------------------------------------------
#

program_name = 'graph2cylinder.py'


if __name__ == '__main__':
        
    usage = "Usage: "+program_name+" [options] postprocessed_graph(simplified).db [output_file.db]\n"+\
            "   or  "+program_name+" --help";
            
    parser = OptionParser(usage)
    
    parser.add_option("--clobber", action="store_true", dest="clobber",
                    default=0, help="overwrite output file")
    parser.add_option("--tangent_threshold",type="float", dest="tangent_threshold",
                    default=0.174, help="tangent threshold in radian to make cylinders. Defulat: 0.174 radian = 10 degrees")                  
    parser.add_option("--use_label", action="store_true", dest="use_label",
                    default=0, help="use labels specified in label edge_property of the input graph. NOTE: input graph needs to have edge_property label, otherwise 0 will be used as the cylinder label!")
    parser.add_option("--use_estimated_label", action="store_true", dest="use_estimated_label",
                    default=0, help="use labels specified in estimated_label edge_property of the input graph, which is calculated from automatic labeling. NOTE: input graph needs to have edge_property estimated_label, otherwise 0 will be used as the cylinder label!")      

    parser.add_option("--use_error_label", action="store_true", dest="use_error_label",
                    default=0, help="use error of labeling specified [1:correct, -1:incorrect] of the input graph, which is calculated from automatic labeling. NOTE: input graph needs to have edge_property error_label, otherwise 0 will be used as the cylinder label!")      
    
    options, args = parser.parse_args()
    print '\n\n>>> %s [%s]: %s\n' % (ctime(time()), os.getcwd(), string.join(argv))
    print "\n\n",options, "\n\n", args , "\n\n"
    

    if len(args) == 2:
        input_file, output_file = args
    elif len(args) ==1:
        input_file = args[0]
        output_file = (input_file[:-3]+ "_cyl.db")
    else:	
        parser.error("incorrect number of arguments")

    
    tang_thresh= options.tangent_threshold

    #print("the input file is: %s\n" %(input_file))

    if not options.clobber and os.path.exists(output_file):
        raise SystemExit, \
            "The --clobber option is needed to overwrite an existing file."

    if options.use_label and options.use_estimated_label: 
        raise SystemExit, \
            "Only one of the options --use_label and --use_estimated_label can be used."
        
    if options.use_label and options.use_error_label: 
        raise SystemExit, \
            "Only one of the options --use_label and --use_error_label can be used."
    
    if options.use_error_label and options.use_estimated_label: 
        raise SystemExit, \
            "Only one of the options --use_error_label and --use_estimated_label can be used."
    
    try:
        g, attributes = graph_analysis.input_graph(input_file, ["history", "vertex_offsets"])  	# open graph data and copy contents to memory
        print ("Sucessfully read in the %s\n" %input_file)	
    except:
        print("Error reading in the %s\n" %input_file)
        
    
    

    #for edge in g.edge_list():
        ##g.set_edge_property(edge, 'cylinder_intermediaries', [])
        #g.set_edge_property(edge, 'cyl_centreX', [])
        #g.set_edge_property(edge, 'cyl_centreY', [])
        #g.set_edge_property(edge, 'cyl_centreZ', [])
        #g.set_edge_property(edge, 'cyl_tangentX', [])
        #g.set_edge_property(edge, 'cyl_tangentY', [])
        #g.set_edge_property(edge, 'cyl_tangentZ', [])
        #g.set_edge_property(edge, 'cyl_radius', [])
        #g.set_edge_property(edge, 'cyl_height', [])
        ##g.set_edge_property(edge, 'cyl_label', int())
        ##g.set_edge_property(edge, 'cyl_segment_ID', int())
        
    #### for each edge in [e[0],g.edge_property(e, 'intermediaries'),e[1]] do as you did for the vertices down here
    #also create a dimension0 = #of edges
    #dimension1=a list which shows # of cylinders for each edge => to be used in brainview reading h5 file
    #segment_ID = 0
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
            
            angle = math.acos (innerprod ) #angle=acos(t1.t2)#in radian	#axis= norm(t1 x t2)	
            if ( (abs(r0-r1) < r0/2) and (abs(angle)< tang_thresh) ):	#|r1-r0|<0.5r0=> 0.5r0<r1<1.5r0  &  3*0.17 #10 degrees= 0.174	=> merge them into 1 cylinder			
                r0 = (r0+r1)/2
                t0 = t1+t0
                tang0 = t0/linalg.norm(t0)
                v1 = vertex_list[i+1]

            else:	# => make a cylinder obj 
                centre = (g.vertices[v0].centre + g.vertices[v1].centre )/2
                centreX.append(numpy.float32(centre[0]))
                centreY.append(numpy.float32(centre[1]))
                centreZ.append(numpy.float32(centre[2]))
                tangentX.append(numpy.float32(tang0[0]))		#t = g.vertices[v0].centre - g.vertices[v1].centre )	#tang = t/linalg.norm(t)
                tangentY.append(numpy.float32(tang0[1]))
                tangentZ.append(numpy.float32(tang0[2]))
                radius.append(numpy.float32(r0))				#radius = ( g.vertices[v0].radius + g.vertices[v1].radius )/2
                height.append(numpy.float32(linalg.norm(t0)))	#height = linalg.norm(t)
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
        tangentX.append(numpy.float32(tang[0]))		#
        tangentY.append(numpy.float32(tang[1]))
        tangentZ.append(numpy.float32(tang[2]))
        radius.append( numpy.float32((float(g.vertices[v0].radius) + float(g.vertices[v1].radius) )/2))
        height.append(numpy.float32(linalg.norm(t)))	
        ##segment.append(int(segment_ID))
        
        g.set_edge_property(e,'cyl_centreX', centreX)		#graph.set_edge_property(self,pair,key,value)
        g.set_edge_property(e,'cyl_centreY', centreY)
        g.set_edge_property(e,'cyl_centreZ', centreZ)
        g.set_edge_property(e,'cyl_tangentX', tangentX)  
        g.set_edge_property(e,'cyl_tangentY', tangentY)  
        g.set_edge_property(e,'cyl_tangentZ', tangentZ)  
        g.set_edge_property(e,'cyl_radius', radius)
        g.set_edge_property(e,'cyl_height', height)
        #g.set_edge_property(e,'edge_ID', int(segment_ID))

        if options.use_label and 'label' in g.edge_properties(e).keys():	##use the same label as the current edge_property  == making no change
            g.set_edge_property(e,'cyl_label', int(g.edge_property(e,'label')) )
        elif options.use_estimated_label and 'estimated_label' in g.edge_properties(e).keys(): ##use the label in the estimated_label edge_property
            g.set_edge_property(e,'cyl_label', int(g.edge_property(e,'estimated_label')[0][0]) )
        elif options.use_error_label and 'error_label' in g.edge_properties(e).keys():
            g.set_edge_property(e,'cyl_label', int(g.edge_property(e,'error_label')) )		
        else:
            g.set_edge_property(e,'cyl_label', int(0))  #"No label                \0" #add label dataset to edge properties and initialize it to None with null-terminted character '\0'
        ##g.set_edge_property(e,'num_cyl', int(len(segment)))  #"No label                \0" #add label dataset to edge properties and initialize it to None with null-terminted character '\0'

        ##g.set_edge_property(e,'segment_ID', int(segment_ID))
        #segment_ID +=1


    
    history = '>>> %s: %s' % (ctime(time()), string.join(argv))		##argv = ['code.py', 'input.db', 'output']
    if attributes.has_key("history"):
        history = attributes["history"] + "\n" + history
        del attributes['history']
        
    graph_analysis.output_graph(output_file, g, history, attributes)   
    
    print ("Succefully wrote the %s\n" %output_file)
    
    ## save the final cylinder file into h5 file for the Brainview
    output_final=output_file.replace(".db", ".h5")
    cmd=("graph2graph.py %s %s --clobber " %(output_file,output_final))	
    print(cmd)
    os.system(cmd)
