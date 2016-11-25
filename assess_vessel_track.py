#!/usr/bin/env python

# To asses how much of the volume of CT vessels are not tracked in mm^3 and % , # of segments and % of tracked
#  Created May 09, 2011
#  Last modified April 11,2012
#  Sahar Ghanavati
 
#PYTHONPATH=/home/jgsled/lib64_ubuntu_hardy/python:$PYTHONPATH


from vessel_tracking import graph_analysis, distribution_analysis
from morphology import graph, cluster_skeleton, object_io
from numpy import *
import py_minc
from optparse import OptionParser, Option, OptionValueError
from sys import argv
import os,commands, shelve, string
from time import time, ctime
import copy

#---------------------------------------------------------------------------------
#

program_name = 'Asses_Vessel_Track.py'


if __name__ == '__main__':
        
    usage = "Usage: "+program_name+" [options] CTimage.mnc postprocess_output_vessel_tracking.db \n"+\
            "   or  "+program_name+" --help";
            
    parser = OptionParser(usage)

    parser.add_option("--clobber", action="store_true", dest="clobber",
                        default=0, help="overwrite output file")
    #parser.add_option("--mask_minor",action="store_true",default=0,dest="mask_minor",help="mask minor vessels smaller than the diameter threshold, otherwise vessels are not masked ")
    parser.add_option("--lower_threshold",type="float",dest="lower_threshold",default = 0.25,help="The lower threshold for binarizing .mnc image,Defult=0.25")
    parser.add_option("--upper_threshold",type="float",dest="upper_threshold",default = 2.0,help="The upper threshold for binarizing .mnc image,Defult=2.0")
    parser.add_option("--scale_factor",type="float",dest="scale_factor",default = 1.2,help="scale factor for isosurface (radius correction factor), Defult=1.2")
    parser.add_option("--remove_interms",action="store_true",default=0, dest="remove_interms",help="remove the intermediate mnc files that were created for calculating seeds in iterations")
    parser.add_option("--graph2dt",action="store_true",default=0, dest="graph2dt",help="To make Diameter based colour rendering")
    parser.add_option("--vertex", action="store_true", dest="vertex", help="graph2dt option: create file with vertex numbers")
    parser.add_option("--property", type="string", dest="property", help="graph2dt option: create file with given vertex property")
    parser.add_option("--edge_property", type="string", dest="edge_property", help="graph2dt option: create file with given edge property")
    parser.add_option("--edge_method", type="string", dest="edge_method", default = "max", help="graph2dt option: method for mapping edge properties to adjacent vertices (choose either: max [default], min, or mean)")
    parser.add_option("--transform", type="string", dest="transform", help="graph2dt option: transform the property or edge_property using the specified function (one of: log, log10, exp, pow10)")
    parser.add_option("--contour", action="store_true", dest="contour", help="graph2dt option: use elliptic contours rather than spheres")

    parser.add_option("--distribution_analysis",action="store_true",default=0, dest="distribution_analysis",help="To generate diameter cumulative distribution")
    parser.add_option("--simplify_graph",action="store_true",default=0, dest="simplify_graph",help="simplify graph to intermediaries and calculate diameter and length, for the distribution_analysis, if the input graph is not simplified already")

    options, args = parser.parse_args()
    print '>>> %s [%s]: %s\n' % (ctime(time()), os.getcwd(), string.join(argv))



    if len(args) != 2:
        parser.error("incorrect number of arguments")
        
    #print (args)	#list of str = ['inputname', 'outputname']

    input_image, input_graph = args	 #eval: convert str to list	#input_file, output_file= args => input_file is list

    l_thresh = options.lower_threshold 
    #if l_thresh <0.2:
        #l_thresh = 0.2 #make sure the minimum entered for lower threshold is 0.2
    u_thresh = options.upper_threshold
            

    #### reads in the graph for more corrections
    try:
        g, attributes = graph_analysis.input_graph(input_graph, ["history", "vertex_offsets"])
        print ("Succefully read in the %s\n" %input_graph)	
    except:
        print("Error reading in the %s\n" %input_graph)
        
    # in postprocessing the vertex radii has been corrected for estimation bias due to blurring of 40% of vessel diameter => radius should be correct by 60% 
    # in order to include all the haze around vessels in the minc file inside the tracked vessels => the calculated not tracked volume then will be only for missed vessels and not the haze around tracked vessels
    # => increase radius by 20% instead of normal 6%
    # save the object file of the corrected db file to check if it worked correctly	
    #path2isosurf gets h5 and return obj, graph2isosurf gets db and returns obj, scale is the radius correction factor
    graph2=input_graph[:-3]+"_isosurface.obj"
    if not options.clobber and os.path.exists(graph2):
        raise SystemExit, \
            "The --clobber option is needed to overwrite an existing file."
            
            
    cmd=("graph2isosurf.py %s %s %s --scale_factor %f --save_volume %s --clobber" %(input_graph,input_image,graph2, options.scale_factor, graph2[:-3]+"mnc"))	#python /micehome/jgsled/bin/
    if options.clobber or (not options.clobber and not os.path.exists(graph2)):
        print(cmd)
        os.system(cmd)

    cmnd = ("mincmath -clobber -lt -const 0.9 %s %s -quiet" %(graph2[:-3]+"mnc",graph2[:-4]+"_binary.mnc") )
    print "\n", cmnd, "\n"
    os.system(cmnd)

    ##### make a mnc out of the corrected.obj
    #cmd=("surface_mask2 -binary_mask %s %s %s " %(input_image,graph3,graph3[:-4]+"_binary.mnc"))	
    #if options.clobber or (not options.clobber and not os.path.exists(graph3[:-4]+"_binary.mnc")):
        #print(cmd)
        #os.system(cmd)

    # make the CT.mnc binary
    #ct=py_minc.ArrayVolume(input_image)	# without this always [z,y,x], with this :dim_names=py_minc.FILE_ORDER_DIMENSION_NAMES order as the image	
    #ct_sizes=ct.get_sizes()

    #### binarize ct image with 0.25 threshold and 2
    cmd=("mincmath -segment -const2 %f %f -clobber -quiet %s %s" %(l_thresh, u_thresh, input_image,input_image[:-4]+"_binary.mnc"))	#python /micehome/jgsled/bin/
    if options.clobber or (not options.clobber and not os.path.exists(input_image[:-4]+"_binary.mnc")):
        print(cmd)
        os.system(cmd)

    #subtract original_image - isosurface
    cmnd = ("minccalc -clobber -expression 'if(A[0]>0.995 && A[1]==0){out=1;} else {out=0}' %s %s %s -quiet" %(input_image[:-4]+"_binary.mnc",graph2[:-4]+"_binary.mnc",graph2[:-4]+"_binary_subtracted.mnc") )
    print "\n", cmnd, "\n"
    os.system(cmnd)
    #os.system("Display %s" %(input_image[:-4]+"_binary_subtracted.mnc"))	#Display subtracted image

    #Open filter on the subtraction
    cmd=("mincmorph -successive OOOO %s %s " %(graph2[:-4]+"_binary_subtracted.mnc", graph2[:-4]+"_binary_subtracted_Open.mnc" ))
    if options.clobber or (not options.clobber and not os.path.exists(graph2[:-4]+"_binary_subtracted_Open.mnc")):
        print(cmd)
        os.system(cmd)	
        
    #calculate the percentage of coverage
    cmd=("mincstats -floor 0.95 %s > %s" %(graph2[:-4]+"_binary_subtracted_Open.mnc",  os.path.dirname(os.path.abspath(input_graph))+"/stats_untracked_"+os.path.basename(input_graph)[:-2]+"txt"))
    #print(cmd)
    if options.clobber or (not options.clobber and not os.path.exists(os.path.dirname(os.path.abspath(input_graph))+"/stats_untracked_"+os.path.basename(input_graph)[:-2]+"txt")):
        print(cmd)
        os.system(cmd)
        
    ##create the polygon of subtract
    #cmd=("marching_cubes %s %s %f" %(graph2[:-4]+"_binary_subtracted.mnc", graph2[:-4]+"_binary_subtracted.obj", 0.95 ))
    #if options.clobber or (not options.clobber and not os.path.exists(graph2[:-4]+"_binary_subtracted.obj")):
        #print(cmd)
        #os.system(cmd)	

        
        
    cmd=("mincstats -floor 0.95 %s > %s" %(input_image[:-4]+"_binary.mnc",  os.path.dirname(os.path.abspath(input_graph))+"/stats_original_mnc_"+os.path.basename(input_image)[:-3]+"txt"))
    #print(cmd)
    if options.clobber or (not options.clobber and not os.path.exists(os.path.dirname(os.path.abspath(input_image))+"/stats_original_mnc_"+os.path.basename(input_image)[:-3]+"txt")):
        os.system(cmd)


    #total number of segments, write in the output "stats_"+input_graph[:-2]+"txt"
    h=copy.deepcopy(g)
    graph_analysis.simplify_graph_retaining_intermediaries(h)
    segment_num =len(h.edge_list())

    filename= os.path.dirname(os.path.abspath(input_graph))+"/stats_segment#_"+os.path.basename(input_graph)[:-2]+"txt"
    f = open(filename, 'a')		#appending mode
    f.write('\nThe number of segments in the vessel_tracking result %s is: %d' %(input_graph,segment_num) )


    if options.remove_interms:
        cmnd = ("rm %s " %(graph2[:-4]+"*mnc"))
        print  "\n", cmnd, "\n"
        os.system(cmnd)

        cmnd = ("rm %s " %(input_image[:-4]+"_binar*.mnc"))
        print  "\n", cmnd, "\n"
        os.system(cmnd)
        
        
    if options.graph2dt:
        cmnd = ("python /home/jgsled/Source/vessel_tracking/graph2dt.py %s %s %s" %(input_graph ,input_image, input_graph[:-3]+"_graph2dt.mnc"))
        if options.vertex:
            cmnd+= " --vertex"
        if options.property:
            cmnd+= (" --property %s" %options.property)
        if options.edge_property:
            cmnd+= (" --edge_property %s" %options.edge_property)
        if options.edge_method:
            cmnd+= (" --edge_method %s" %options.edge_method)
        if options.transform:
            cmnd+= (" --transform %s" %options.transform)
        if options.contour:
            cmnd+= " --contour"

        print  "\n", cmnd, "\n"
        os.system(cmnd)
        
    if options.distribution_analysis:
        cmnd = ("python /projects/souris/sghanavati/src/scripts/distribution_analysis.py %s" %(input_graph))
        if options.simplify_graph:
            cmnd+= " --simplify_graph"
        print  "\n", cmnd, "\n"
        os.system(cmnd)


    print "\n\nyou're done :D\n\n"






















