#!/usr/bin/env python

# -*- coding: utf-8 -*-

#  Run vessel-tracking code in iterations for full coverage of the vessels in the minc image  
#
#  Created Jan 17, 2012
#  modified Nov 15, 2012	#replace os in exec_follow_tree with subprocess and comment import os
#  Last modified Jan 27, 2014  #add run time to each of the iterations
#  Sahar Ghanavati
 
	

from optparse import OptionParser, Option, OptionValueError
from sys import argv
#from numpy import *
import py_minc 
import numpy as np
import os
import subprocess
import operator
import time, string

#---------------------------------------------------------------------------------
#
#variable declared in main will be global whereas variable declared within a func will be local var

def file_exists(in_file):
#Test whether a path exists.  Returns True for broken symbolic links
    try:
        F = open(in_file,'rb')
    except IOError as e:
        print in_file," doesn't exist!"
        return False
    F.close()
    print in_file," exists!"
    return True



def seeding(mnc_bin_file,tag_file,connectivity,all_tags):
    print "****************** seeding : ", tag_file
        
    #create distance_transform (binary input only)
    cmnd = ("mincmorph -clobber -short -successive F %s %s > %s/tmp.log" %(mnc_bin_file,mnc_bin_file[:-4]+"_dt.mnc",os.path.dirname(os.path.abspath(mnc_bin_file))) )
    print cmnd
    os.system(cmnd)
    #p= subprocess.Popen(cmnd, shell=True, stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE, close_fds=True)
    #stdout, stderr = p.communicate()
    ##p.wait()

    #create group file
    cmnd = ("mincmorph -clobber -short -group -%s %s %s > %s/tmp.log" %(connectivity,mnc_bin_file,mnc_bin_file[:-4]+"_group.mnc",os.path.dirname(os.path.abspath(mnc_bin_file))) )
    print cmnd
    os.system(cmnd)
    #p= subprocess.Popen(cmnd, shell=True, stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE, close_fds=True)
    #stdout, stderr = p.communicate()
    ##p.wait()


    dt_img = py_minc.ArrayVolume(mnc_bin_file[:-4]+"_dt.mnc")# without this always [z,y,x], with ,dim_names=py_minc.FILE_ORDER_DIMENSION_NAMES order as the image	
    group_img = py_minc.ArrayVolume(mnc_bin_file[:-4]+"_group.mnc")# without this always [z,y,x], with ,dim_names=py_minc.FILE_ORDER_DIMENSION_NAMES order as the image	

    sizes = dt_img.get_sizes()
    #print "\n\nsize of ", group_filename," is ",sizes, "\n\n" 

    component_maximum = {}
    for z_vox in range(sizes[0]):
        for y_vox in range(sizes[1]):
            for x_vox in range(sizes[2]):
                g_value = int(round(group_img.array[z_vox,y_vox,x_vox]))
                if g_value>0:
                    dt_value = int(round(dt_img.array[z_vox,y_vox,x_vox]))
                    if g_value not in component_maximum.keys():
                        component_maximum[g_value]= [dt_value,z_vox,y_vox,x_vox]
                    else:
                        if dt_value > component_maximum[g_value][0]:
                            component_maximum[g_value]= [dt_value,z_vox,y_vox,x_vox]

    #descending sort tag components based on value (dt_value)
    sorted_component_maximum = sorted(component_maximum.iteritems(), key=operator.itemgetter(1),reverse=True)	#result is not dict, list of tuples=(key,val)
    tags = py_minc.VolumeTags(1)

    for i in range(len(sorted_component_maximum)):
        voxel_location = np.array([sorted_component_maximum[i][1][1],sorted_component_maximum[i][1][2],sorted_component_maximum[i][1][3]])
        world_location = dt_img.convert_voxel_to_world(voxel_location)
        if not list(world_location) in all_tags:
            tags.append(world_location,0.05,1,1,1) 
            all_tags.append(list(world_location))

    if tags.number()>0:
        tags.output(tag_file)
        print "****************** seeding : ", tag_file, " DONE! ",tags.number()," seeds created!\n"
        return True, all_tags
    else:
        print "****************** No seed can be added!\n"
        return False, all_tags



def exec_follow_tree(num,mnc_file, tag_file, h5_output,output_obj,output_avoid,branch_limit,rho,gamma,beta,smooth,contrast,delta,scale,avoidance_radius,extent,blur_ratio,dratio,avoidance_ratio,epsilon,limit,extent_step,store_true_options,output_prev_run=0):
    print"\n\n##################################################################\n ######################## generating tree " , num, " ##################\n##################################################################\n\n"
    localtic = time.time()

    #python /home/jgsled/development/bin/follow_tree.py if export PYTHONPATH=/micehome/jgsled/development/lib/python/:/micehome/jgsled/lib64_ubuntu_lucid/python/:/projects/mice/share/arch/linux-x86_64-eglibc2_11_1/lib/python2.6/site-packages/:/projects/souris/sghanavati/src/bin/bin:/micehome/sghanavati/pydpiperbin/bin
    cmnd = ("follow_tree.py --branch_limit %d %s --save_obj %s %s %s %s --rho %f --epsilon %f --gamma %f --beta %f %f %f %f --smooth %f %f --contrast %f --delta %f --scale %f --avoidance_radius %d %d --extent %d --extent_step %d --limit %d --blur_ratio %f --dratio %f --avoidance_ratio %f --save_avoidance %s --avoidance" %(branch_limit,store_true_options,output_obj,mnc_file, tag_file, h5_output,rho,epsilon,gamma,beta[0],beta[1],beta[2],beta[3],smooth[0],smooth[1],contrast,delta,scale,avoidance_radius[0],avoidance_radius[1],extent,extent_step,limit,blur_ratio,dratio,avoidance_ratio,output_avoid))
    if output_prev_run!=0:
        cmnd += (" --append %s " %output_prev_run)
    print cmnd
    os.system(cmnd)
    #log_file = h5_output[:-2]+"log"
    #p= subprocess.Popen(cmnd, shell=True, stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE, close_fds=True)
    #stdout, stderr = p.communicate()
    #FILE = open(log_file,"w")
    #FILE.write(stdout)
    #FILE.close()
    ##p.wait()

    print"\n\n##################################################################\n ######################## generating tree " , num, " DONE! ##################\n##################################################################\n\n"

    localtoc = time.time()-localtic
    hr = int(localtoc)/3600
    minut = int(localtoc - (hr*3600))/60
    sec = localtoc - (hr*3600) - (minut*60)
    print ("\n\n\nIteration %d of vessel tracking took %d hours and %d minutes and %d seconds\nEND!" %(num,hr,minut,sec)) 		


def binarize_minc_file(mnc_file,bin_threshold,outdir):
    #binarize mnc_file
    mnc_bin_file = outdir+"/"+os.path.basename(mnc_file)[:-4]+"_binarize.mnc"
    cmnd = ("mincmath -clobber -byte -gt -const %f %s %s -quiet" %(bin_threshold,mnc_file,mnc_bin_file) )
    print cmnd
    os.system(cmnd)
    #p= subprocess.Popen(cmnd, shell=True, stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE, close_fds=True)
    #stdout, stderr = p.communicate()
    ##p.wait()
    return mnc_bin_file


def find_untracked(mnc_file,h5_output,scale_factor):
    print "\n****************** find untracked ", h5_output

    #path2isosurf gets h5 and return obj, graph2isosurf gets db and returns obj => --save_volume saves inverse volume where tracked vessels have intensity [0...2.0] and background is 2.0 
    cmnd = ("python /home/jgsled/bin/path2isosurf.py --clobber --scale_factor %f %s %s %s --save_volume %s" %(scale_factor,h5_output,mnc_file,h5_output[:-3]+"_isosurf.obj",h5_output[:-3]+"_isosurf.mnc"))
    print  cmnd
    os.system(cmnd)
    #p= subprocess.Popen(cmnd, shell=True, stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE, close_fds=True)
    #stdout, stderr = p.communicate()
    ##p.wait()

    cmnd = ("mincmath -clobber -byte -lt -const 1.5 %s %s -quiet" %(h5_output[:-3]+"_isosurf.mnc",h5_output[:-3]+"_isosurf_binarize.mnc") )
    print cmnd
    os.system(cmnd)
    #p= subprocess.Popen(cmnd, shell=True, stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE, close_fds=True)
    #stdout, stderr = p.communicate()
    ##p.wait()

    diff_file = h5_output[:-3]+"_diff.mnc"
    cmnd = ("minccalc -clobber -byte -expression 'if(A[0]>%f && A[1]==0){out=1;} else {out=0}' %s %s %s -quiet" %(bin_threshold,mnc_file,h5_output[:-3]+"_isosurf_binarize.mnc",diff_file ) )
    print cmnd
    os.system(cmnd)
    #p= subprocess.Popen(cmnd, shell=True, stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE, close_fds=True)
    #stdout, stderr = p.communicate()
    ##p.wait()

    print "****************** find untracked ", h5_output, " DONE!\n"
    return diff_file



def additional_seeding(mnc_bin_file,tag_file, all_tags):
    print "****************** additional seeding : ", tag_file

    #create distance_transform (binary input only)
    cmnd = ("mincmorph -clobber -short -successive F %s %s > %s/tmp.log" %(mnc_bin_file,mnc_bin_file[:-4]+"_dt.mnc",os.path.dirname(os.path.abspath(mnc_bin_file))) )
    print cmnd
    os.system(cmnd)
    #p= subprocess.Popen(cmnd, shell=True, stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE, close_fds=True)
    #stdout, stderr = p.communicate()
    ##p.wait()

    volume_dt=py_minc.ArrayVolume(mnc_bin_file[:-4]+"_dt.mnc")# without this always [z,y,x], with ,dim_names=py_minc.FILE_ORDER_DIMENSION_NAMES order as the image	

    size=volume_dt.get_sizes()

    tags=py_minc.VolumeTags(1)
                            
    for z in range(size[0]):
        for y in range(size[1]):
            for x in range(size[2]):
                dt_value=volume_dt.array[z,y,x]		

                if dt_value>=0.5: 
                    dt_value_world_location= volume_dt.convert_voxel_to_world(np.array([z,y,x]))		#return world in [x,y,z]
                    if not list(dt_value_world_location) in all_tags:		
                        tags.append(dt_value_world_location,1,1,1,1)

    if tags.number()>0:
        tags.output(tag_file)
        print "****************** additional seeding : ", tag_file, " DONE! ",tags.number()," seeds created!\n"
        return True
    else:
        print "****************** No additional seed can be added!\n"
        return False



    ################################


program_name = 'auto_vessel_tracking_v3.1.py'


if __name__ == '__main__':
        
    usage = "Usage: "+program_name+" [options] input_image.mnc output_tree.h5 [original_seed.tag]\n"+\
            "   or  "+program_name+" --help";
            
    parser = OptionParser(usage)

    #### follow_tree options                      
    parser.add_option("--branch_limit", type="int", dest="branch_limit",
                        metavar="branch_limit", default = 10000, 
                        help="maximum number of branches in tree [default: 10000]")
                        
    parser.add_option("--rho", type="float", dest="rho",
                        metavar="rho", default = 1.5, 
                        help="maximum step size along vessel path [default: 1.5]")                      
    parser.add_option("--limit", type="int", dest="limit",
                        metavar="limit", default = 10000, 
                        help="maximum number of point in path [default: 10000]")

    parser.add_option("--epsilon", type="float", dest="epsilon",
                        metavar="epsilon", default = 0.9, 
                        help="tolerance for eccentric vessels used in finding vessel near seed point [default: 0.9]") 
                        
    parser.add_option("--gamma", type="float", dest="gamma",
                        metavar="gamma", default = 8, 
                        help="tolerance for curvature between steps [default: 8]")
                        
    parser.add_option("--beta", nargs = 4, type="float", dest="beta",
                        metavar="default min max increment",
                        default = (0.2, 0.1, 0.7, 0.05), 
                        help="bounds on step size in voxels for advancing normal plane [default: 0.2 0.01 0.7 0.05]")                      
                        
    parser.add_option("--smooth", type="float", dest="smooth", nargs=2,
                        metavar="smoothness", default = (0.05, 0.01), 
                        help="coefficient on regularization term for path and radius [default: 0.05 0.01]")  
                        
    parser.add_option("--contrast", type="float", dest="contrast",
                        metavar="constrast_limit", default = 0.2, 
                        help="lower limit for edge contrast when continuing a vessel [default: 0.2]")                      
                        
    parser.add_option("--delta", type="float", dest="delta",
                        metavar="delta", default = 0.01, 
                        help="tolerance for non-orthogonal grandient and tangent used in finding vessel near seed point [default: 0.01]")
                            
    parser.add_option("--scale", type="float", dest="scale",
                        metavar="scale", default = 0.1, 
                        help="scale value used in searching near seed point (real space units) [default: 0.1]")
                            
    parser.add_option("--avoidance_radius", type="int", nargs=2,
                        dest="avoidance_radius", default=(0, 10), metavar="min_radius max_radius",
                        help="use a sphere shaped structuring element within given range of radii for avoidance. Radii are in voxel unit (default: [0 10])")
                            
    parser.add_option("--extent", type="int", dest="extent",
                        metavar="nplanes", default = 1, 
                        help="Number of planes ahead and behind normal plane to include in tube detection. [default: 1]")
                            
                        
    parser.add_option("--blur_ratio", type="float", dest="blur_ratio",
                        metavar="ratio", default = 0.4, 
                        help="ratio of blur kernel stddev to vessel radius in multiscale algorithm [default: 0.4]")

    parser.add_option("--dratio", type="float", dest="dratio",
                        metavar="delta", default = 1.0, 
                        help="spoke edge contrast ratio threshold for detecting branches [default: 1.0]")

    parser.add_option("--avoidance_ratio", type="float",
                        dest="avoidance_ratio", default=0.8, metavar="ratio",
                        help="ratio of radius of avoidance sphere to vessel radius (default: 0.8)")

    #parser.add_option("--save_avoidance", type="string", dest="avoidance_file",
                        #metavar="avoidance.mnc", 
                        #help="save avoidance map to given file (Only useful for debugging).")


    parser.add_option("--multiscale", action="store_true", dest="multiscale",default = 0, help="process data in multiscale fashion with variable width blurring kernel")                      
    parser.add_option("--minima_branching", action="store_true",
                        dest="minima_branching",
                        default = 0, help="branch at all minima that are below threshold")
    parser.add_option("--grad_comp", action="store_true", dest="grad_comp",
                        default = 0, help="compensate directional derivatives based on anisotropic psf")
    parser.add_option("--save_blurs", action="store_true", dest="save_blurs_flag",
                        help="record blurring used at every step")
    parser.add_option("--save_contours", action="store_true", dest="save_contours_flag",
                        help="record elliptic contour used at every step")
    #parser.add_option("--avoidance", action="store_true", dest="avoidance",
                        #default=0,
                        #help="use voxel map to avoid intersecting paths")	
    parser.add_option("--voxel_coord", action="store_true", dest="voxel_coord",
                        default=0,
                        help="save result using voxel coordinates")
    parser.add_option("--normalize", action="store_true", dest="normalize",
                        default=0,
                        help="normalize intensity range to [0, 1] as required by the algorithm")
    parser.add_option("--shuffle", action="store_true", dest="shuffle",
                        help="shuffle (reproducibly) the order in which seeds are processed")
    parser.add_option("--clobber", action="store_true", dest="clobber",
                        default=0, help="overwrite output file")

    parser.add_option("--blur_comp", action="store_const", const=1,dest="blur_comp", default = 0, help="compensate contrast by blur width to reduce contrast variation with vessel radius (square_root)")
    parser.add_option("--blur_comp_linear", action="store_const", const=2,dest="blur_comp", help="compensate contrast by blur width to reduce contrast variation with vessel radius (linear)")

    parser.add_option("--extent_step", type="int", dest="extent_step",
                        metavar="nvoxels", default = 1, 
                        help="Number of voxels between planes when extent is greater than 1. [default: 1].  Note, this is at present only implemented in combination with the psf option.")
    parser.add_option("--psf", nargs = 3, type="float", dest="psf",
                        metavar="xspace yspace zspace", 
                        help="Standard deviation of anisotropic points spread function in each of three directions in world space units (mm).")

    parser.add_option("--rescale", dest="rescale", type="float", nargs=2,
                        metavar="a b", help="rescale intensities by mapping specified interval to [0, 1]")
    parser.add_option("--seed_range", type="int", nargs=2, dest="seed_range",
                        metavar="first last", help="specify a subset of seeds to process by seed index (counting from 1), e.g. --seed_range 1 1000 (default: full process all seeds")
                            
    parser.add_option("--noisy", action="store_const", const=2,
                        dest="verbose",
                        help="print hightly verbose progress messages")
    parser.add_option("--verbose", action="store_const", const=1,
                        dest="verbose",  default=1,
                        help="print progress messages (default)")
    parser.add_option("--quiet", action="store_const", const=0,
                        dest="verbose",
                        help="do not print progress messages")
                        
                        
    #### auto_vessel_track options                      
                        
    parser.add_option("--run",type="int", dest="run", default= 3, help="number of times follow tree is executed (<10)") 

    parser.add_option("--add_coverage",action="store_true",default=0, dest="add_coverage",help="running the program through brute force algorithm to get the remaining uncovered vessels thereby adding one additional run to the program" )

    parser.add_option("--append", type="string", dest="append", help="append new branches to specified file")

    parser.add_option("--append_initnum", type="int", dest="append_initnum", metavar="append_initnum", default = 0, help="start iteration number of follow tree, if not equal to 0 then should be used with append option [default: 0]")

    parser.add_option("--bin_threshold",type="float", dest="bin_threshold", default= 0.2, help="intensity threshold for binarizing the image") 

    parser.add_option("--scale_factor",type="float",dest="scale_factor",default = 1.2,help="scale factor for isosurface (radius correction factor)")

    parser.add_option("--connectivity",type="string",dest="connectivity",default="3D06",help="connectivity for mincmorph")

    parser.add_option("--automatic_seed",action="store_true",default=0,dest="automatic_seed",help="input original seed automatically")

    parser.add_option("--remove_interms",action="store_true",default=0, dest="remove_interms",help="remove the intermediate mnc files that were created for calculating seeds in iterations")


    options,args = parser.parse_args() 

    if (len(args) < 2) or (len(args) > 3):
        parser.error("incorrect number of arguments")
        
    if ( (len(args) != 3) and not options.automatic_seed):
        parser.error("either an initial seed file should be given as input or use option automatic_seed")

    if ( (len(args) == 3) and options.automatic_seed):
        parser.error("either an initial seed file should be given as input or use option automatic_seed")

    if options.append and not file_exists(options.append):
        raise SystemExit, ("--append specifed yet file %s does not exist." % options.append)

    if options.append_initnum!=0 and not options.append:
        raise SystemExit, ("--append should be used with a non-zero append_initnum!")

    print '\n\n>>> %s [%s]: %s\n' % (time.ctime(time.time()), os.getcwd(), string.join(argv))
    print "\n\n",options, "\n\n", args , "\n\n"

    tic= time.time()



    branch_limit=options.branch_limit
    rho=options.rho
    limit=options.limit
    epsilon = options.epsilon
    extent_step = options.extent_step
    gamma=options.gamma
    beta=options.beta
    smooth=options.smooth
    contrast=options.contrast
    delta=options.delta
    scale=options.scale
    avoidance_radius=options.avoidance_radius 
    extent=options.extent 
    blur_ratio=options.blur_ratio
    dratio=options.dratio
    avoidance_ratio=options.avoidance_ratio



    store_true_options=""
    if options.multiscale:
        store_true_options+="--multiscale "
    if options.minima_branching:
        store_true_options+="--minima_branching "
    if options.grad_comp:
        store_true_options+="--grad_comp "
    if options.save_blurs_flag:
        store_true_options+="--save_blurs "
    if options.save_contours_flag:
        store_true_options+="--save_contours "
    #if options.avoidance:
        #store_true_options+="--avoidance "
    if options.voxel_coord:
        store_true_options+="--voxel_coord "
    if options.normalize:
        store_true_options+="--normalize "
    if options.shuffle:
        store_true_options+="--shuffle "
    if options.clobber:
        store_true_options+="--clobber "
    if (options.blur_comp==1):
        store_true_options+="--blur_comp "
    if (options.blur_comp==2):
        store_true_options+="--blur_comp_linear "
        
    if options.psf:
        psf=options.psf
        store_true_options+=("--psf %f %f %f " %(psf[0],psf[1],psf[2]))
        
    if options.rescale:	
        rescale=options.rescale
        store_true_options+=("--rescale %f %f " %(rescale[0],rescale[1]))
        
    if options.seed_range is not None:	
        seed_range=options.seed_range
        store_true_options+=("--seed_range %d %d " %(seed_range[0],seed_range[1]))
        
    if (options.verbose==0):
        store_true_options+="--quiet"
    if (options.verbose==1):
        store_true_options+="--verbose"
    if (options.verbose==2):
        store_true_options+="--noisy"
        

    n=options.run
    n = int(n)   	
    scale_factor = options.scale_factor 
    connectivity = options.connectivity 
    bin_threshold = options.bin_threshold

    print "number of runs:",n," ","branch_limit:",branch_limit,"\n","scale_factor:",scale_factor," ","connectivity:",connectivity," ",\
        "blur_ratio:",blur_ratio," ", "contrast:",contrast,"\n"

    print "\n######################## Start of program ########################\n"
    all_tags=[]

    if len(args) == 2:
        print"\n\nautomatic seeding...\n\n"
        minc_file, output_hdf5 = args 
        seed_stat = False
        
    if len(args) == 3:
        print"\n\ninputing seed file manually\n\n"
        minc_file, output_hdf5, seed_file = args 
        seed_stat = True
        manual_seeds = py_minc.VolumeTags(seed_file)		#=> append them to all_tags
        all_tags=[list(i) for i in manual_seeds.locations]

    initnum = options.append_initnum	#default is 0
    initnum = int(initnum) 

    for num in range(initnum,n):
        h5_output = str(output_hdf5[:-3]+"_"+str(num)+".h5")
        output_obj = str(output_hdf5[:-3]+"_"+str(num)+".obj")
        output_avoid = str(output_hdf5[:-3]+"_"+str(num)+"_avoid.mnc")
                
        if num==0:
            if not seed_stat:	#seed_file not inputed => --automatic_seeding: create the seed
                minc_bin_file = binarize_minc_file(minc_file,bin_threshold,os.path.dirname(os.path.abspath(output_hdf5)))
                seed_file = output_hdf5[:-3]+"_0_seeds.tag"
                seed_stat, all_tags = seeding(minc_bin_file,seed_file,connectivity,all_tags)

            if not seed_stat:
                print "No original seed file is inputed or can be created!\nAborted!\n"
                exit(0)
            else:	
                exec_follow_tree(num,minc_file, seed_file, h5_output,output_obj,output_avoid,branch_limit,rho,gamma,beta,smooth,contrast,delta,scale,avoidance_radius,extent,blur_ratio,dratio,avoidance_ratio,epsilon,limit,extent_step,store_true_options)
                seed_stat = False	# => in next iteration creates the seed_file 
                if not file_exists(h5_output):
                    print "\nERROR:",h5_output, " couldn't be created!\nAborted!\n"
                    exit(0)
        else:
            if num == initnum:
                output_prev_run = options.append
            else:
                output_prev_run = str(output_hdf5[:-3]+"_"+str(num-1)+".h5")
                
            if not seed_stat:	#if num==initnum and seed_file inputed => skip if statement and use input seed; if (num>initnum (in the iterations)) or (num==initnum and no seed_file inputed )=> create seed_file
                #create the seed file for this round
                diff_file = find_untracked(minc_file,output_prev_run,scale_factor)
                seed_file = str(output_hdf5[:-3]+"_"+str(num)+"_seeds.tag")
                seed_stat, all_tags = seeding(diff_file,seed_file,connectivity,all_tags)
                
            if seed_stat:	
                exec_follow_tree(num,minc_file, seed_file, h5_output,output_obj,output_avoid,branch_limit,rho,gamma,beta,smooth,contrast,delta,scale,avoidance_radius,extent,blur_ratio,dratio,avoidance_ratio,epsilon,limit,extent_step,store_true_options,output_prev_run)
                seed_stat = False	# => in next iteration creates the seed_file
                if not file_exists(h5_output):
                    print "\nERROR:",h5_output, " couldn't be created!\nIterations aborted!\n"
                    break
            else:
                print "No seed file can be created in this iteration!\nIterations aborted!\n"
                break


    if options.add_coverage:
        output_final = str(output_hdf5[:-3]+"_add_coverage.h5")
        output_obj = str(output_hdf5[:-3]+"_add_coverage.obj")
        output_avoid = str(output_hdf5[:-3]+"_add_coverage_avoid.mnc")
        if n == options.append_initnum:
            output_prev_run = options.append
        else:	
            output_prev_run = str(output_hdf5[:-3]+"_"+str(n-1)+".h5")
        
        diff_file = find_untracked(minc_file,output_prev_run,scale_factor)
        seed_file = str(output_hdf5[:-3]+"_additional_seeding.tag")  

        add_stat = additional_seeding(diff_file,seed_file, all_tags)
        if add_stat:
            exec_follow_tree(n,minc_file, seed_file, output_final,output_obj,output_avoid,branch_limit,rho,gamma,beta,smooth,contrast,delta,scale,avoidance_radius,extent,blur_ratio,dratio,avoidance_ratio,epsilon,limit,extent_step,store_true_options,output_prev_run)
            if not file_exists(output_final):
                print "\nERROR:",output_final, " couldn't be created!\n"
        else:
            print "\nNo additional seed file can be created for add_coverage!\n"
        

    if options.remove_interms:
        cmnd = ("rm %s " %(output_hdf5[:-3]+"*_isosurf*"))	#output_hdf5[:-3]+"*_isosurf.obj output_hdf5[:-3]+"*_isosurf.mnc output_hdf5[:-3]+"*_isosurf_binarize.mnc
        print  "\n", cmnd, "\n"
        os.system(cmnd)
        #p= subprocess.Popen(cmnd, shell=True, stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE, close_fds=True)
        #stdout, stderr = p.communicate()
        ##p.wait()

        cmnd = ("rm %s " %(minc_file[:-4]+"_binarize*.mnc"))
        print  "\n", cmnd, "\n"
        os.system(cmnd)
        #p= subprocess.Popen(cmnd, shell=True, stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE, close_fds=True)
        #stdout, stderr = p.communicate()
        ##p.wait()

        cmnd = ("rm %s/tmp.log" %os.path.dirname(os.path.abspath(output_hdf5)))
        print  "\n", cmnd, "\n"
        os.system(cmnd)
        #p= subprocess.Popen(cmnd, shell=True, stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE, close_fds=True)
        #stdout, stderr = p.communicate()
        ##p.wait()

    print "\n\nyou're done :D\n\n"

    toc = time.time()-tic
    hr = int(toc)/3600
    minut = int(toc - (hr*3600))/60
    sec = toc - (hr*3600) - (minut*60)
    print ("\n\n\nThe vessel tracking took %d hours and %d minutes and %d seconds\nEND!" %(hr,minut,sec)) 		

            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
        	