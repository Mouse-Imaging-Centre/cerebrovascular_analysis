#!/usr/bin/env python

# Postprocess the output of the vessel_tracking (tree.h5) before doing any further process on the datasets
# Convert tree.h5 to graph.db, delete_internal_leaves, adjust junction and correct radii due to blurring

#STEPS:
##python /micehome/jgsled/bin/tree2graph.py file_tree.h5 file_graph.db --join_some  => delete_internal_leaves and remember list of disconnected vertices
#python tree2graph.py file_tree.h5 file_graph.db --join_all and disconnect the vertices in join_some delete_internal_leaves
#python /micehome/jgsled/bin/delete_internal_leaves.py input.db output.db
#python /micehome/jgsled/bin/delete_short_leaves.py input.db output.db
#*correct location of vertices : # optimize location of branch points
#graph_analysis.adjust_junctions(g, vertex_offsets)
#*correct radius by 1.06*old_radius : # diameters need to corrected for estimation bias due to blurring
#diameters = [h.edge_property(edge, 'diameter') for edge in h.edge_list()] 
#diameters = array(diameters)*1.06
#smooth_graph by centre w/o aggregate, by radius with aggregate_isolate	
# with --simplify_graph reduce graph to intermediaries, calculate length and diameter and save 
# remove the isolated edges (where each vertex only has one neighbour edge) from the simplified graph, with --layer to get rid of isolated Y and X structures.


#h5ls -r file_tree_graph2graph.h5
#h5ls -rv file_tree_graph2graph.h5


#  Created July 13, 2011
#  Last modified August 7, 2012
#  Sahar Ghanavati


from vessel_tracking import graph_analysis
from morphology import graph
from optparse import OptionParser, Option, OptionValueError
from sys import argv
import os, string
import commands
from subprocess import check_call
from time import time, ctime

#---------------------------------------------------------------------------------
#


def delete_internal_leaves_jsome2jall(jall_name,jsome_name):
    #disconnect vertices from g that are disconnect in h
    try:
        g, attributes = graph_analysis.input_graph(jall_name, ["history", "vertex_offsets"])
        print ("Successfully read in the %s\n" %jall_name)	
    except:
        print("Error reading in the %s\n" %jall_name)

    h = graph_analysis.input_graph(jsome_name)

    if len(g.vertices)!=len(h.vertices):
        print "ERROR: delete_internal_leaves_jsome2jall, the number of vertices of g and h are not equal!!\nAborted!"
        exit(0)
    for i in xrange(len(g.vertices)):
        if len(h.vertices[i].edges)==0:
            g.disconnect_vertex(i)
            
    out_name = jall_name[:-3]+"_delleave.db"
    history = '\n>>> %s: delete_internal_leaves on join_all using the join_some result' % (ctime(time()))		##argv = ['code.py', 'input.db', 'output']
    if attributes.has_key("history"):
        history = attributes["history"] + "\n" + history
        del attributes['history']
        
    graph_analysis.output_graph(out_name, g, history, attributes)   
    print ("Succefully wrote the %s\n" %out_name)
        
    return out_name


program_name = 'postprocess_vessel_track.py'


if __name__ == '__main__':
        
    usage = "Usage: "+program_name+" [options] vessel_tracking_tree1.h5 vessel_tracking_tree2.h5 ... vessel_tracking_treeN.h5 output_file.db\n"+\
            "   or  "+program_name+" --help";
            
    parser = OptionParser(usage)

    parser.add_option("--clobber", action="store_true", dest="clobber",
                        default=0, help="overwrite output file")
    parser.add_option("--mask", type="string", dest="mask", help="mask graph with the given mask")
    parser.add_option("--radi_correction",type="float",dest="radi_correction",default = 1.06,
                      help="The radius correction factor for the vessel_tracking underestimated result, Defult=1.06")
    parser.add_option("--del_leave_factor",type="float",dest="del_leave_factor",default = 1.0,
                      help="scale vessel radii by given factor for delete_internal_leaves.py, Defult=1.0")
    parser.add_option("--offset", type="float", dest="offset", default = 0.06,
                      help="add an offset to the vessel radii for delete_internal_leaves.py distance in mm(default 0.06)")
    parser.add_option("--del_short_leaves",action="store_true",default = 0,dest="del_short_leaves",
                      help="To use delete_short_leaves.py to delete vertices on short leaves (vessels from an endpoint)")
    parser.add_option("--radius_threshold",type="float",dest="radius_threshold",default = 0.1,
                      help="The threshold radius in mm for delete_short_leaves.py. Default value is 0.1mm")
    parser.add_option("--length_threshold",type="float",dest="length_threshold",default = 0.2,
                      help="The threshold length in mm for delete_short_leaves.py,Default is 0.2mm")
    parser.add_option("--intermed_threshold",type="int",dest="intermed_threshold",default = 0,
                      help="The threshold on number of intermediaries for delete_short_leaves.py. Recommended value is 100.")
    parser.add_option("--smooth_vertex",action="store_true",default=0, dest="smooth_vertex",
                      help="smooth graph by vertex property")
    parser.add_option("--smooth_vertex_tol", type="float", dest="smooth_vertex_tol", default = 0.001, 
                      help="specify error tolerance per vertex in smoothed property [default: 0.001]")
    parser.add_option("--smooth_radius",action="store_true",default=0, dest="smooth_radius",
                      help="smooth graph by radius property")
    parser.add_option("--smooth_radius_tol", type="float", dest="smooth_radius_tol", default = 1, 
                      help="specify error tolerance per vertex in smoothed property [default: 1]")
    #parser.add_option("--mask_minor",action="store_true",default=0,dest="mask_minor",help="mask minor vessels smaller than the diameter threshold, otherwise vessels are not masked ")
    #parser.add_option("--diameter_threshold",type="float",dest="diameter_threshold",default = 0.2,help="The threshold diameter for masking vessels, Defult=0.2")
    #parser.add_option("--length_threshold",type="float",dest="length_threshold",default = 0.5,help="The threshold length for masking shorter vessels from an endpoint, Defult=0.5")
    parser.add_option("--remove_interms",action="store_true",default=0, dest="remove_interms",
                      help="remove the intermediate mnc files that were created for calculating seeds in iterations")
    parser.add_option("--simplify_graph",action="store_true",default=0, dest="simplify_graph",
                      help="simplify graph to intermediaries and calculate diameter and length")
    parser.add_option("--remove_isolated",action="store_true",default=0, dest="remove_isolated",
                      help="remove isolated edges from the simplified graph to give a fully-connected graph")
                        

    options, args = parser.parse_args()

    print '>>> %s [%s]: %s\n' % (ctime(time()), os.getcwd(), string.join(argv))

    if len(args)<2:
        parser.error("incorrect number of arguments")
    else:
        input_files = args[:-1]
        output_file = args[-1]

    output_directory = os.path.dirname(output_file)
    print "Writing output files to: %s" % output_directory
    print output_file
    
      
    if not options.clobber and os.path.exists(output_file):
        raise SystemExit, \
            "The --clobber option is needed to overwrite an existing file."
            
    out2s = []
    for input_file in input_files:
        out1_basename = os.path.splitext(os.path.basename(input_file))[0].replace('tree','graph')
        out1=os.path.join(output_directory, out1_basename + ".db")
        print out1
        # convert output of vessel_tracking.h5 to python pickeled object.db file
        out_jsome=os.path.join(output_directory, out1_basename + "_jsome.db")
        cmd=("tree2graph.py %s %s --join_some --clobber" %(input_file,out_jsome))	
        print(cmd)
        check_call(cmd, shell=True)	

        # convert output of vessel_tracking.h5 to python pickeled object.db file
        cmd=("tree2graph.py %s %s --join_all --clobber" %(input_file,out1))
        print(cmd)
        check_call(cmd, shell=True)	

        #mask graph
        if options.mask:
            cmd=("mask_graph.py %s %s %s --negate --clobber" %(out_jsome,options.mask, out_jsome+"_masked.db"))	
            print(cmd)
            check_call(cmd, shell=True)	
            out_jsome = out_jsome[:-3]+"_masked.db"
            
            cmd=("mask_graph.py %s %s %s --negate --clobber" %(out1,options.mask, out1[:-3]+"_masked.db"))	
            print(cmd)
            check_call(cmd, shell=True)	
            out1 = out1[:-3]+"_masked.db"


        # delete internal leaves
        out2=os.path.join(output_directory, out1_basename + "_delleave.db")
        cmd=("delete_internal_leaves_sahar.py %s %s --scale %f --offset %f --clobber" % \
                 (out_jsome,out_jsome[:-3]+"_delleave.db", options.del_leave_factor, options.offset))	
        print(cmd)
        check_call(cmd, shell=True)

        out2 = delete_internal_leaves_jsome2jall(out1,out_jsome[:-3]+"_delleave.db")
        out2s.append(out2)
        
    print type(out2s)
    print len(out2s)
    print out2s

    if len(input_files)>1:
        cmd="merge_graphs.py --clobber "
        for out2 in out2s:
            cmd+=" %s" %out2
        cmd+= " %s" %output_file[:-3]+"_merged.db"
        print(cmd)
        check_call(cmd, shell=True)  
        out2 = output_file[:-3]+"_merged.db"
        
        cmd=("graph2obj.py --clobber %s %s" %(out2,out2[:-2]+"obj")) 
        print(cmd)
        check_call(cmd, shell=True)

    # delete short leaves
    if options.del_short_leaves:
        out2= out2[:-3]+"_delshort.db"
        cmd=("delete_short_leaves.py %s %s --radius_threshold %f --length_threshold %f --intermed_threshold %d --clobber" % \
                 (out1[:-3]+"_delleave.db",out2, options.radius_threshold, options.length_threshold, options.intermed_threshold))	
        print(cmd)
        check_call(cmd, shell=True)



    # reads in the graph for adjust junction and radi corrections
    try:
        g, attributes = graph_analysis.input_graph(out2, ["history", "vertex_offsets"])
        print ("Successfully read in the %s\n" %out2)	
    except:
        print("Error reading in the %s\n" %out2)
        
    # optimize location of branch points
    graph_analysis.adjust_junctions(g, attributes['vertex_offsets'])

    # vertex radii need to be corrected for estimation bias due to blurring of 40% of vessel diameter => radius should be correct by 60%
    for v in g.vertices:
        v.radius = v.radius*options.radi_correction

    # save the corrected graph
    out3=out2[:-3]+("_radi%.2f.db" %options.radi_correction)
    history = '\n>>> %s: corrected branch point locations and radius bias due to blurring by %f' % (ctime(time()), options.radi_correction)		##argv = ['code.py', 'input.db', 'output']
    if attributes.has_key("history"):
        history = attributes["history"] + "\n" + history
        del attributes['history']
        
    graph_analysis.output_graph(out3, g, history, attributes)   
    print ("Succefully wrote the %s\n" %out3)

    if options.smooth_vertex:
        smooth_input=out3
        out3=smooth_input[:-3]+"_smooth.db"
        cmd=("smooth_graph.py --property centre --smooth %f %s %s --clobber" %(options.smooth_vertex_tol, smooth_input,out3)) 
        print(cmd)
        check_call(cmd, shell=True)
        
    if options.smooth_radius:
        cmd=("smooth_graph.py --property radius --smooth %f --aggregate --isolate %s %s --clobber" %(options.smooth_radius_tol, out3, smooth_input[:-3]+"_smooth2.db")) 
        print(cmd)
        check_call(cmd, shell=True)
        out3=smooth_input[:-3]+"_smooth2.db"

    #save the object file of the corrected db file to check if it worked correctly	
    cmd=("graph2obj.py --clobber %s %s" %(out3,out3[:-2]+"obj"))
    print(cmd)
    check_call(cmd, shell=True)

    if options.simplify_graph:
        cmnd = ("simplify_graph.py %s" %out3)
        if options.clobber:
            cmnd += " --clobber"
        print  "\n", cmnd, "\n"
        os.system(cmnd)
        out3 = out3[:-3]+"_simplified.db"
        
    if not options.simplify_graph and options.remove_isolated:
        raise SystemExit, ("--remove_isolated can only be used with --simplify_graph")
    elif options.simplify_graph and options.remove_isolated:
        cmnd = ("remove_isolated_edges.py %s --layer" %(out3))
        print  "\n", cmnd, "\n"
        os.system(cmnd)
        out3 = out3[:-3]+"_rmisolate.db"
        
        cmd=("graph2obj.py --clobber %s %s" %(out3,out3[:-3]+".obj"))	
        print(cmd)
        check_call(cmd, shell=True)
        
    if len(input_files)==1:
        cmd=("mv %s %s" %(out3,output_file))  
        print(cmd)
        check_call(cmd, shell=True)

        
    if options.remove_interms:
        cmnd = "rm "+os.path.dirname(os.path.abspath(out1))+"/*_jsome*"
        print  "\n", cmnd, "\n"
        os.system(cmnd)
        
        cmnd = "rm "+os.path.dirname(os.path.abspath(out1))+"/*_delleave.db"
        print  "\n", cmnd, "\n"
        os.system(cmnd)
        
        #cmnd = "rm "+os.path.dirname(os.path.abspath(out1))+"/*_radi1.06.db"
        #print  "\n", cmnd, "\n"
        #os.system(cmnd)
        
        #cmnd = "rm "+os.path.dirname(os.path.abspath(out1))+"/*_simplified.db"
        #print  "\n", cmnd, "\n"
        #os.system(cmnd)
        
        cmnd = "rm "+os.path.dirname(os.path.abspath(out1))+"/*_smooth.db"
        print  "\n", cmnd, "\n"
        os.system(cmnd)
        
        #cmnd = "rm "+os.path.dirname(os.path.abspath(out1))+"/*_smooth2.db"
        #print  "\n", cmnd, "\n"
        #os.system(cmnd)
        

        
        
        
