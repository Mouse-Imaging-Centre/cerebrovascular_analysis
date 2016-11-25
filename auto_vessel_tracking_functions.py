#!/usr/bin/env python
# -*- coding: utf-8 -*-

#  Run follow_tree on subvolumes in iterations for full coverage of the vessels in the minc image and merge subtrees
#
#  Created Jan 17, 2012
#  modified Nov 15, 2012	#replace os in exec_follow_tree with subprocess and comment import os
#  modified Jan 27, 2014  #add run time to each of the iterations
#  Last modified MAy 10, 2014  #add subvolume to automatic seeding and vessel tracking and then in postprocessing we should merge_graphs at the end
#  Sahar Ghanavati
 
    

from optparse import OptionGroup, OptionParser, Option, OptionValueError
from sys import argv

import py_minc 
#import pyminc.volumes.factory as pymincf
#import pyminc.volumes.volumes as pymincv

import numpy as np
import os
import subprocess
import operator
import time, string

#from pydpiper.application import AbstractApplication
#from pydpiper.pipeline import Pipeline, CmdStage, InputFile, OutputFile, LogFile
#import pydpiper.file_handling as fh
from os.path import abspath
#import logging
#import Pyro
import sys
import csv
#import atoms_and_modules.minc_atoms as ma


def addVTOptionGroup(parser):
    group = OptionGroup(parser, "Auto_Vessel_Tracking options", 
                        "Options for performing iterative subvolume vessel tracking with automatic seeding.")
    #### follow_tree options                      
    group.add_option("--branch_limit", type="int", dest="branch_limit",
                        metavar="branch_limit", default = 10000, 
                        help="maximum number of branches in tree [default: 10000]")
                        
    group.add_option("--rho", type="float", dest="rho",
                        metavar="rho", default = 1.5, 
                        help="maximum step size along vessel path [default: 1.5]")                      
    group.add_option("--limit", type="int", dest="limit",
                        metavar="limit", default = 10000, 
                        help="maximum number of point in path [default: 10000]")
    
    group.add_option("--epsilon", type="float", dest="epsilon",
                        metavar="epsilon", default = 0.9, 
                        help="tolerance for eccentric vessels used in finding vessel near seed point [default: 0.9]") 
                        
    group.add_option("--gamma", type="float", dest="gamma",
                        metavar="gamma", default = 2, 
                        help="tolerance for curvature between steps [default: 2] (for /micehome/jgsled/Source/vessel_tracking/follow_tree.py use 8)")
                        
    group.add_option("--beta", nargs = 4, type="float", dest="beta",
                        metavar="default min max increment",
                        default = (0.2, 0.1, 0.7, 0.05), 
                        help="bounds on step size in voxels for advancing normal plane [default: 0.2 0.01 0.7 0.05]")                      
                        
    group.add_option("--smooth", type="float", dest="smooth", nargs=2,
                        metavar="smoothness", default = (0.05, 0.01), 
                        help="coefficient on regularization term for path and radius [default: 0.05 0.01]")  
                        
    group.add_option("--contrast", type="float", dest="contrast",
                        metavar="constrast_limit", default = 0.2, 
                        help="lower limit for edge contrast when continuing a vessel [default: 0.2]")                      
                        
    group.add_option("--delta", type="float", dest="delta",
                        metavar="delta", default = 0.01, 
                        help="tolerance for non-orthogonal grandient and tangent used in finding vessel near seed point [default: 0.01]")
                            
    group.add_option("--scale", type="float", dest="scale",
                        metavar="scale", default = 0.1, 
                        help="scale value used in searching near seed point (real space units) [default: 0.1]")
                            
    group.add_option("--avoidance_radius", type="int", nargs=2,
                        dest="avoidance_radius", default=(0, 10), metavar="min_radius max_radius",
                        help="use a sphere shaped structuring element within given range of radii for avoidance. Radii are in voxel unit (default: [0 10])")
                            
    group.add_option("--extent", type="int", dest="extent",
                        metavar="nplanes", default = 1, 
                        help="Number of planes ahead and behind normal plane to include in tube detection. [default: 1]")
                            
                        
    group.add_option("--blur_ratio", type="float", dest="blur_ratio",
                        metavar="ratio", default = 0.4, 
                        help="ratio of blur kernel stddev to vessel radius in multiscale algorithm [default: 0.4]")
    
    group.add_option("--dratio", type="float", dest="dratio",
                        metavar="delta", default = 1.0, 
                        help="spoke edge contrast ratio threshold for detecting branches [default: 1.0]")
    
    group.add_option("--avoidance_ratio", type="float",
                        dest="avoidance_ratio", default=0.8, metavar="ratio",
                        help="ratio of radius of avoidance sphere to vessel radius (default: 0.8)")
    
    #group.add_option("--save_avoidance", type="string", dest="avoidance_file",
                        #metavar="avoidance.mnc", 
                        #help="save avoidance map to given file (Only useful for debugging).")
    group.add_option("--multiscale", action="store_true", dest="multiscale",default = 0, help="process data in multiscale fashion with variable width blurring kernel")                      
    group.add_option("--minima_branching", action="store_true",
                        dest="minima_branching",
                        default = 0, help="branch at all minima that are below threshold")
    group.add_option("--grad_comp", action="store_true", dest="grad_comp",
                        default = 0, help="compensate directional derivatives based on anisotropic psf")
    group.add_option("--save_blurs", action="store_true", dest="save_blurs_flag",
                        help="record blurring used at every step")
    group.add_option("--save_contours", action="store_true", dest="save_contours_flag",
                        help="record elliptic contour used at every step")
    group.add_option("--avoidance", action="store_true", dest="avoidance",
                        default=0,
                        help="use voxel map to avoid intersecting paths")  
    group.add_option("--voxel_coord", action="store_true", dest="voxel_coord",
                        default=0,
                        help="save result using voxel coordinates")
    group.add_option("--normalize", action="store_true", dest="normalize",
                        default=0,
                        help="normalize intensity range to [0, 1] as required by the algorithm")
    group.add_option("--shuffle", action="store_true", dest="shuffle",
                        help="shuffle (reproducibly) the order in which seeds are processed")
    group.add_option("--clobber", action="store_true", dest="clobber",
                        default=0, help="overwrite output file")
    
    group.add_option("--blur_comp", action="store_const", const=1,dest="blur_comp", default = 0, help="compensate contrast by blur width to reduce contrast variation with vessel radius (square_root)")
    group.add_option("--blur_comp_linear", action="store_const", const=2,dest="blur_comp", help="compensate contrast by blur width to reduce contrast variation with vessel radius (linear)")

    group.add_option("--extent_step", type="int", dest="extent_step",
                        metavar="nvoxels", default = 1, 
                        help="Number of voxels between planes when extent is greater than 1. [default: 1].  Note, this is at present only implemented in combination with the psf option.")

    group.add_option("--extent_ratio", type="float", dest="extent_ratio",
                      metavar="ratio", help="Specify a fixed ratio between vessel radius and extent steps. (optimal value: 0.5)")
    group.add_option("--psf", nargs = 3, type="float", dest="psf",
                      metavar="xspace yspace zspace", 
                      help="Standard deviation of anisotropic points spread function in each of three directions in world space units (mm).")
    group.add_option("--psf_comp", action="store_true", dest="psf_comp",
                      default = 0, help="compensate for directional derivatives based on anisotropic psf (implies --blur_comp_linear)")
    group.add_option("--psf_z_dependence", nargs = 4, type="float", dest="psf_z_dependence",
                      metavar="z1 z2 sigma1 sigma2", 
                      help="Specify psf that varies linearly with z.  Values of sigma1 and sigma2 at depths z1 and z2 all in mm define the linear relation (implies --psf).")
    
    group.add_option("--rescale", dest="rescale", type="float", nargs=2,
                        metavar="a b", help="rescale intensities by mapping specified interval to [0, 1]")
    group.add_option("--seed_range", type="int", nargs=2, dest="seed_range",
                        metavar="first last", help="specify a subset of seeds to process by seed index (counting from 1), e.g. --seed_range 1 1000 (default: full process all seeds")
                            
    group.add_option("--noisy", action="store_const", const=2,
                        dest="verbose",
                        help="print hightly verbose progress messages")
    group.add_option("--verbose", action="store_const", const=1,
                        dest="verbose",  default=1,
                        help="print progress messages (default)")
    group.add_option("--quiet", action="store_const", const=0,
                        dest="verbose",
                        help="do not print progress messages")
    #### auto_vessel_track options                      
                        
    group.add_option("--run",type="int", dest="run", default= 3, help="number of times follow tree is executed (<10)") 
    
    #group.add_option("--add_coverage",action="store_true",default=0, dest="add_coverage",help="running the program through brute force algorithm to get the remaining uncovered vessels thereby adding one additional run to the program" )
    group.add_option("--cmndname", type="string", dest="followtree_cmndname",default="follow_tree.py", help="follow_tree command name (Default:follow_tree in the bin directory [you can use /micehome/jgsled/Source/vessel_tracking/follow_tree.py for old version of follow_tree with gamma 8])")
    
    group.add_option("--append", type="string", dest="append", help="append new branches to specified file")
    
    group.add_option("--append_initnum", type="int", dest="append_initnum", metavar="append_initnum", default = 0, help="start iteration number of follow tree, if not equal to 0 then should be used with append option [default: 0]")
    
    group.add_option("--bin_threshold",type="float", dest="bin_threshold", default= 0.15, help="intensity threshold for binarizing the image (default 0.15)") 
    
    group.add_option("--scale_factor",type="float",dest="scale_factor",default = 1.2,help="scale factor for isosurface (radius correction factor)[default 1.2]")
    
    group.add_option("--connectivity",type="string",dest="connectivity",default="3D06",help="connectivity for mincmorph [default 3D06]")
    
    group.add_option("--automatic_seed",action="store_true",default=0,dest="automatic_seed",help="input original seed automatically")
    
    group.add_option("--remove_interms",action="store_true",default=0, dest="remove_interms",help="remove the intermediate mnc files that were created for calculating seeds in iterations")
    
    group.add_option("--subvol_nums",type="int", dest="subvol_nums", default= 1, help="number of subvolumes to be tracked (default = 1)") 
    group.add_option("--additional_seeding", action="store_true", dest="additional_seeding",
                        default=0, help="Add seeds exhaustively")
    group.add_option("--laplacian", action="store_true", dest="laplacian",
                        default=0, help="Use the laplacian of the distance_transform to force seeds to the centre of the vessels")
    
    parser.add_option_group(group)
    
#---------------------------------------------------------------------------------
#
#variable declared in main will be global whereas variable declared within a func will be local var

def file_exists(in_file):
#Test whether a path exists.  Returns True for broken symbolic links
    try:
        F = open(in_file,'rb')
    except IOError as e:
        #print in_file," doesn't exist!"
        return False
    F.close()
    #print in_file," exists!"
    return True
    
    
def print_timing(tic,toc,n="END"):
    t = toc-tic
    hr = int(t)/3600
    minut = int(t - (hr*3600))/60
    sec = t - (hr*3600) - (minut*60)
    if n == "END":
        print ("\n\n\nThe vessel tracking took %d hours and %d minutes and %d seconds\nEND!" %(hr,minut,sec))
        print ("\n\nyou're done :D\n\n")
    else:
        print ("\n\n\nThe iteration %s of vessel tracking took %d hours and %d minutes and %d seconds\nEND!" %(n,hr,minut,sec))      
        

class SubVolumer:
    """Class that calculate the subvolume dimensions based on the number of subvols,
    followed by resampling the volume to create the subvolumes.
    Currently used in subvolume_vessel_tracking."""
    def __init__(self,
                 minc_file,
                 subvol_nums,clobber=1,
                 overlap=10):
        
        self.minc_file = minc_file
        self.subvol_nums = subvol_nums
        self.clobber = clobber
        self.overlap = overlap
        self.input_subvols = []
        self.spacing = []
        self.nel = []
        self.strts = []
        
        
        self.buildPipeline()
        
        
    def mincResampler(self,outfile):
        cmd = "mincresample -clobber %s %s -step %f %f %f -nelements %d %d %d -start %f %f %f" %(self.minc_file, outfile, self.spacing[2],self.spacing[1],self.spacing[0], self.nel[0],self.nel[1],self.nel[2], self.strts[0],self.strts[1], self.strts[2])   
        print cmd
        os.system(cmd)
        #mincSubvoler = CmdStage(cmd)
        #self.p.addStage(mincSubvoler)      
        #print self.p.cmd
        
    def  buildPipeline(self):   
        #def create_subvolumes (minc_file,subvol_nums, overlap=10):
        """v0min v0max v1min v1max v2min v2max",
        only process a subvolume specified by the bounds in voxel coordinates: v0min <= v0 < v0max, v1min <= v1 < v1max, v2min <= v2 < v2max"   
        create subvolumes 6 parameters to be passed to each parallel process of vesselTracking."""
        v0min = [] #x
        v0max  = []#x
        v1min  = []#y
        v1max = [] #y
        v2min  = []#z
        v2max = []#z
        #img = pymincf.volumeFromFile(self.minc_file)
        #sizes = img.getSizes()      # z_vox(sizes[0]), y_vox(sizes[1]),x_vox(sizes[2])
        #self.spacing = img.getSeparations()# z_vox(sizes[0]), y_vox(sizes[1]),x_vox(sizes[2])
        #origin = img.getStarts()# z_vox(sizes[0]), y_vox(sizes[1]),x_vox(sizes[2])
        img = py_minc.ArrayVolume(self.minc_file)# without this always [z,y,x], with ,dim_names=py_minc.FILE_ORDER_DIMENSION_NAMES order as the image   
        sizes = img.get_sizes()        # z_vox(sizes[0]), y_vox(sizes[1]),x_vox(sizes[2])
        self.spacing = img.get_separations()# z_vox(sizes[0]), y_vox(sizes[1]),x_vox(sizes[2])
        origin = img.get_starts()# z_vox(sizes[0]), y_vox(sizes[1]),x_vox(sizes[2])
        
        if self.subvol_nums > 27:
            self.subvol_nums = 27
        if self.subvol_nums < 1:
            self.subvol_nums = 1
            
        if self.subvol_nums == 1:
            self.input_subvols.append(self.minc_file)
        else:            
            [xdiv, ydiv, zdiv] = self.make_subvol_indx()
            #### crop the mincfile into each subvolume
            subv = 0
            for i in range(xdiv):
                for j in range(ydiv):
                    for k in range(zdiv):
                        subvol_minc_file = self.minc_file[:-4]+"_crop"+str(subv)+".mnc"
                        v0min = max(0,int(i*sizes[2]/xdiv)-self.overlap)
                        v0max = min(sizes[2],int((i+1)*sizes[2]/xdiv)+self.overlap)
                        v1min = max(0,int(j*sizes[1]/ydiv)-self.overlap)
                        v1max = min(sizes[1],int((j+1)*sizes[1]/ydiv)+self.overlap)
                        v2min = max(0,int(k*sizes[0]/zdiv)-self.overlap)
                        v2max = min(sizes[0],int((k+1)*sizes[0]/zdiv)+self.overlap)
                            
                        self.nel = [v0max-v0min, v1max-v1min,v2max-v2min]  #x,y,z
                        self.strts = [origin[2]+(v0min*self.spacing[2]), origin[1]+(v1min*self.spacing[1]), origin[0]+(v2min*self.spacing[0])]  #x,y,z
                        """Resample minc_file into space of subvol_minc_file
                        Here, we add stages of the pipeline that is defined in the main file (subvolume_vessel_tracking)""" 
                        if self.clobber or not file_exists(subvol_minc_file):
                            self.mincResampler(subvol_minc_file)
                        self.input_subvols.append(subvol_minc_file)
                        #resample = ma.mincresample(self.minc_file, subvol_minc_file, likeFile=self.minc_file, argArray=["-clobber", "-step %f %f %f" %(self.spacing[2],self.spacing[1],self.spacing[0]),"-nelements %d %d %d" %(self.nel[0],self.nel[1],self.nel[2]),"-start %f %f %f" %(self.strts[0],self.strts[1], self.strts[2])])
                        #self.p.addStage(resample)
                        #print self.p.cmd
                        #self.input_subvols.append(resample.outputFiles[0])                
                        subv += 1
    
    def make_subvol_indx(self):
        xdiv = 1
        ydiv = 1
        zdiv = 1
        
        if (self.subvol_nums%2 == 0):
            xdiv = 2*xdiv
            self.subvol_nums = self.subvol_nums/2
            if (self.subvol_nums%2 == 0):
                ydiv = 2*ydiv
                self.subvol_nums = self.subvol_nums/2
        if (self.subvol_nums%3 == 0):
            if (xdiv==1):
                xdiv = 3*xdiv
            elif (ydiv==1):
                ydiv = 3*ydiv
            else:   
                zdiv = 3*zdiv
            self.subvol_nums = self.subvol_nums/3
            if (self.subvol_nums%3 == 0):
                zdiv = 3*zdiv
                self.subvol_nums = self.subvol_nums/3
        if (self.subvol_nums%5 == 0):
            if (xdiv==1):
                xdiv = 5*xdiv
            elif (ydiv==1):
                ydiv = 5*ydiv
            else:   
                zdiv = 5*zdiv
            self.subvol_nums = self.subvol_nums/5
        if (xdiv==1):
            xdiv = self.subvol_nums*xdiv
        elif (ydiv==1):
            ydiv = self.subvol_nums*ydiv
        else:   
            zdiv = self.subvol_nums*zdiv
        return xdiv, ydiv, zdiv 

class FollowTreeer:
    """Class that gets all the options for the follow_tree.py,
    and runs it for the current iteration.
    Currently used in subvolume_vessel_tracking."""
    def __init__(self, followtree_cmndname,
                 num,minc_file, tag_file, h5_output,output_obj,output_avoid,
                 branch_limit,rho,gamma,beta,smooth,contrast,delta,scale,
                 avoidance_radius,extent,blur_ratio,dratio,avoidance_ratio,
                 epsilon,limit,extent_step,
                 store_true_options,
                 output_prev_run=0):
        
        self.followtree_cmndname = followtree_cmndname
        self.num = num
        self.minc_file = minc_file
        self.tag_file = tag_file
        self.h5_output = h5_output
        self.output_obj = output_obj
        self.output_avoid = output_avoid
        
        self.branch_limit = branch_limit
        self.rho = rho
        self.gamma = gamma
        self.beta = beta
        self.smooth = smooth
        self.contrast = contrast 
        self.delta = delta
        self.scale = scale
        self.avoidance_radius = avoidance_radius
        self.extent = extent
        self.blur_ratio = blur_ratio
        self.dratio = dratio
        self.avoidance_ratio = avoidance_ratio
        self.epsilon = epsilon
        self.limit = limit
        self.extent_step = extent_step
        self.store_true_options = store_true_options
        self.output_prev_run = output_prev_run
    
        self.buildPipeline()
        
    def  buildPipeline(self):   
    #def exec_follow_tree(num,minc_file, tag_file, h5_output,output_obj,output_avoid,branch_limit,rho,gamma,beta,smooth,contrast,delta,scale,avoidance_radius,extent,blur_ratio,dratio,avoidance_ratio,epsilon,limit,extent_step,store_true_options,output_prev_run=0):
        print"\n\n##################################################################\n ######################## generating tree ", self.num, " subvolume ",self.minc_file," ##################\n##################################################################\n\n"
        log_file = self.h5_output[:-2]+"log"
        cmd = "sge_batch -o %s -e %s -l vf=16G " %(log_file, log_file[:-3]+"e")
        #cmd += " python /micehome/jgsled/Source/vessel_tracking/follow_tree.py " + self.minc_file +" "+ self.tag_file +" "+ self.h5_output
        #cmd += " follow_tree.py " + self.minc_file +" "+ self.tag_file +" "+ self.h5_output
        if not os.path.dirname(self.followtree_cmndname)=='':
            cmd += " python "
        cmd += " " + self.followtree_cmndname + " " + self.minc_file +" "+ self.tag_file +" "+ self.h5_output
        cmd += " --branch_limit %d --rho %f --epsilon %f --gamma %f --beta %f %f %f %f " %(self.branch_limit,self.rho,self.epsilon,self.gamma,self.beta[0],self.beta[1],self.beta[2],self.beta[3])
        cmd += " --smooth %f %f" %(self.smooth[0],self.smooth[1]) + " --contrast %f" %self.contrast + " --delta %f" %self.delta + " --scale %f" %self.scale + " --avoidance_radius %d %d" %(self.avoidance_radius[0],self.avoidance_radius[1]) + " --extent %d" %self.extent
        cmd += " --extent_step %d" %self.extent_step + " --limit %d" %self.limit + " --blur_ratio %f" %self.blur_ratio + " --dratio %f" %self.dratio + " --avoidance_ratio %f" %self.avoidance_ratio       
        cmd += " --save_avoidance %s " %self.output_avoid + self.store_true_options
        cmd += " --save_obj %s " %self.output_obj
        if self.output_prev_run!=0:
            cmd += " --append %s " %self.output_prev_run
        print cmd
        os.system(cmd)   
        #followTree = CmdStage(cmd)
        #followTree.setLogFile(LogFile(log_file))
        #self.p.addStage(followTree)
        #print self.p.cmd
    
        #cmd = ("python /home/jgsled/bin/follow_tree.py --branch_limit %d %s --save_obj %s %s %s %s --rho %f --epsilon %f --gamma %f --beta %f %f %f %f --smooth %f %f --contrast %f --delta %f --scale %f --avoidance_radius %d %d --extent %d --extent_step %d --limit %d --blur_ratio %f --dratio %f --avoidance_ratio %f --save_avoidance %s --avoidance" %(branch_limit,store_true_options,output_obj,minc_file, tag_file, h5_output,rho,epsilon,gamma,beta[0],beta[1],beta[2],beta[3],smooth[0],smooth[1],contrast,delta,scale,avoidance_radius[0],avoidance_radius[1],extent,extent_step,limit,blur_ratio,dratio,avoidance_ratio,output_avoid))
        #if output_prev_run!=0:
            #cmd += (" --append %s " %output_prev_run)


 

class SeedingClass():
    """Class that calculates the uncovered connected components and put seeds in each,
    Currently used in subvolume_vessel_tracking."""
    def __init__(self,
                 minc_file,
                 tag_file,
                 connectivity,
                 bin_threshold,
                 scale_factor,
                 all_tags,
                 remove_interms,additional_seeding=0,laplacian=0,clobber=0,
                 h5_output = 0):
        
        #self.p = Pipeline()
        self.minc_file = minc_file
        self.tag_file = tag_file
        self.h5_output = h5_output
        self.connectivity = connectivity
        self.bin_threshold = bin_threshold
        self.scale_factor = scale_factor
        self.all_tags = all_tags
        self.remove_interms = remove_interms
        self.additional_seeding = additional_seeding
        self.laplacian = laplacian
        self.clobber = clobber
        
        
        self.mnc_bin_file = self.minc_file[:-4]+"_binarize.mnc"
        try: # the attempt to access the binary minc volume will fail if it doesn't yet exist at pipeline creation
            self.fileH = pymincf.volumeFromFile(self.mnc_bin_file)
        except:
            # if it indeed failed, create it
            self.binarize_file(self.minc_file,self.mnc_bin_file,self.bin_threshold,"gt")
             
        #intermediary files to be deleted after this iteration of vesselTracking is done 
        if self.h5_output == 0:
            self.mnc_dt_file = self.mnc_bin_file[:-4]+"_dt.mnc"
            self.mnc_group_file = self.mnc_bin_file[:-4]+"_group.mnc"            
        else:    
            self.graph_iso_mnc = self.h5_output[:-3]+"_iso.mnc"
            self.graph_iso_obj = self.h5_output[:-3]+"_iso.obj"
            self.graph_iso_bin_mnc = self.graph_iso_mnc[:-4]+"_bin.mnc"
            self.mnc_diff_file = self.graph_iso_bin_mnc[:-4]+"_diff.mnc"
            self.mnc_dt_file = self.mnc_diff_file[:-4]+"_dt.mnc"
            self.mnc_group_file = self.mnc_diff_file[:-4]+"_group.mnc"
        
        self.buildPipeline()
 
    def binarize_file(self,binin,binout,binthresh,lg,dtype="byte"):
        #binarize minc_file "" + 
        cmd = "mincmath" + " -clobber -quiet" + " -%s" %dtype + " -%s -const %f %s %s" %(lg,binthresh,binin,binout)
        #print cmd
        os.system(cmd)
        #binarizer = CmdStage(cmd)
        #self.p.addStage(binarizer)       
        #print self.p.cmd        
        
    def path2IsoSurfacer(self):
        #path2isosurf gets h5 and return obj, graph2isosurf gets db and returns obj => --save_volume saves inverse volume where tracked vessels have intensity [0...2.0] and background is 2.0 
        cmd = "path2isosurf.py" + " --clobber" + " --scale_factor %f %s %s %s" %(self.scale_factor,self.h5_output,self.minc_file,self.graph_iso_obj) + " --save_volume %s" %(self.graph_iso_mnc)
        #print cmd
        os.system(cmd)
        #isoSurfer = CmdStage(cmd)
        #self.p.addStage(isoSurfer)
        #print self.p.cmd
        
    def mincDifferencer(self):
        cmd = "minccalc" + " -clobber" + " -byte" + " -quiet" + " -expression 'if(A[0]>%f && A[1]==0){out=1;} else {out=0}' %s %s %s" %(self.bin_threshold,self.minc_file,self.graph_iso_bin_mnc,self.mnc_diff_file ) 
        #print cmd
        os.system(cmd)
        #mincdiff = CmdStage(cmd)
        #self.p.addStage(mincdiff)      
        #print self.p.cmd
        
    def find_untracked(self):
        print "\n****************** find untracked ", self.h5_output
        if self.clobber or not file_exists(self.graph_iso_mnc):
            self.path2IsoSurfacer()
        if self.clobber or not file_exists(self.graph_iso_bin_mnc):
            self.binarize_file(self.graph_iso_mnc,self.graph_iso_bin_mnc,1.5,"lt")
        self.mincDifferencer()
        print "****************** find untracked ", self.h5_output, " DONE!\n"
        
    def distanceTransfor(self,infile,outfile):
        #create distance transform (binary input only)
        cmd = "mincmorph" + " -clobber" + " -short" + " -successive F %s %s" %(infile,outfile)
        #print cmd
        os.system(cmd)
        #distanceTr = CmdStage(cmd)
        #self.p.addStage(distanceTr)
        #print self.p.cmd
    
    def connectedCompontentFinder(self,infile, outfile,connectiv):
        #create group file
        cmd = "mincmorph" + " -clobber" + " -short" + " -group -%s %s %s" %(connectiv,infile,outfile)
        #print cmd
        os.system(cmd)
        #createGroup = CmdStage(cmd)
        #self.p.addStage(createGroup)     
        #print self.p.cmd
        
    def seeder(self, dtfile, groupfile, seedfile, alltags):
        #create seed file
        cmd = "auto_vessel_tag.py --clobber --bin_threshold %f" %self.bin_threshold + " %s %s %s " %(dtfile,groupfile,seedfile) #+alltags  
        if self.additional_seeding:
            cmd += " --additional_seeding " 
        if self.laplacian:
            cmd += " --laplacian "
        #print cmd
        os.system(cmd)
        #createSeeds = CmdStage(cmd)
        #self.p.addStage(createSeeds)     
        
    def fileRemover(self):
        if self.remove_interms and self.h5_output!=0:
            cmd = "rm" + " tmp.log %s %s %s" %(self.graph_iso_obj, self.h5_output[:-3]+"*_iso*.mnc", self.minc_file[:-4]+"_binarize*.mnc")
            #print cmd
            os.system(cmd)
            #rmf = CmdStage(cmd)
            #self.p.addStage(rmp)
            #print self.p.cmd
            
        
    #def seeding(self,mnc_bin_file,tag_file,connectivity,all_tags):
    def  buildPipeline(self):   
        print "****************** seeding : ", self.tag_file
        if self.h5_output==0:
            if self.clobber or not file_exists(self.mnc_dt_file):
                self.distanceTransfor(self.mnc_bin_file,self.mnc_dt_file)
            if self.clobber or not file_exists(self.mnc_group_file):
                self.connectedCompontentFinder(self.mnc_bin_file,self.mnc_group_file,self.connectivity)            
        else:
            if self.clobber or not file_exists(self.mnc_diff_file):
                self.find_untracked()
            if self.clobber or not file_exists(self.mnc_dt_file):
                self.distanceTransfor(self.mnc_diff_file,self.mnc_dt_file)
            if self.clobber or not file_exists(self.mnc_group_file):
                self.connectedCompontentFinder(self.mnc_diff_file,self.mnc_group_file,self.connectivity)
        if self.clobber or not file_exists(self.tag_file):
            self.seeder(self.mnc_dt_file, self.mnc_group_file, self.tag_file, self.all_tags)
        """
        tags = py_minc.VolumeTags(1)
        #dt_img = py_minc.ArrayVolume(mnc_bin_file[:-4]+"_dt.mnc")# without this always [z,y,x], with ,dim_names=py_minc.FILE_ORDER_DIMENSION_NAMES order as the image  
        #group_img = py_minc.ArrayVolume(mnc_bin_file[:-4]+"_group.mnc")# without this always [z,y,x], with ,dim_names=py_minc.FILE_ORDER_DIMENSION_NAMES order as the image    
        #dt_img = pymincf.volumeFromFile(self.mnc_dt_file)# without this always [z,y,x], with ,dim_names=py_minc.FILE_ORDER_DIMENSION_NAMES order as the image  
        #group_img = pymincf.volumeFromFile(self.mnc_group_file)# without this always [z,y,x], with ,dim_names=py_minc.FILE_ORDER_DIMENSION_NAMES order as the image    
        dt_img = pymincf.volumeFromFile(self.minc_file)# without this always [z,y,x], with ,dim_names=py_minc.FILE_ORDER_DIMENSION_NAMES order as the image   
        group_img = pymincf.volumeFromFile(self.minc_file)# without this always [z,y,x], with ,dim_names=py_minc.FILE_ORDER_DIMENSION_NAMES order as the image 
        dt_img.openFile()
        group_img.openFile()
        
        #sizes = dt_img.get_sizes()
        sizes = dt_img.getSizes()
        
        component_maximum = {}
        for z_vox in range(sizes[0]):
            for y_vox in range(sizes[1]):
                for x_vox in range(sizes[2]):
                    #g_value = int(round(group_img.array[z_vox,y_vox,x_vox]))
                    g_value = int(round(group_img[z_vox,y_vox,x_vox]))
                    if g_value>0:
                        #dt_value = int(round(dt_img.array[z_vox,y_vox,x_vox]))
                        dt_value = int(round(dt_img[z_vox,y_vox,x_vox]))
                        if g_value not in component_maximum.keys():
                            component_maximum[g_value]= [dt_value,z_vox,y_vox,x_vox]
                        else:
                            if dt_value > component_maximum[g_value][0]:
                                component_maximum[g_value]= [dt_value,z_vox,y_vox,x_vox]
                        if dt_value>0.5 and self.n == "additional_seeding": 
                            #dt_value_world_location= dt_img.convert_voxel_to_world(np.array([z,y,x]))       #return world in [x,y,z]
                            dt_value_world_location= dt_img.convertVoxelToWorld(np.array([z_vox,y_vox,x_vox]))       #return world in [x,y,z]
                            if not list(dt_value_world_location) in self.all_tags:      
                                tags.append(dt_value_world_location,1,1,1,1)
                                self.all_tags.append(list(dt_value_world_location))
        

        #descending sort tag components based on value (dt_value)
        sorted_component_maximum = sorted(component_maximum.iteritems(), key=operator.itemgetter(1),reverse=True)   #result is not dict, list of tuples=(key,val)

        for i in range(len(sorted_component_maximum)):
            voxel_location = np.array([sorted_component_maximum[i][1][1],sorted_component_maximum[i][1][2],sorted_component_maximum[i][1][3]])
            #world_location = dt_img.convert_voxel_to_world(voxel_location)
            world_location = dt_img.convertVoxelToWorld(voxel_location)
            if not list(world_location) in self.all_tags:
                tags.append(world_location,0.05,1,1,1) 
                self.all_tags.append(list(world_location))

        dt_img.closeVolume()
        group_img.closeVolume()
        if tags.number()>0:
            tags.output(self.tag_file)
            print "****************** seeding : ", self.tag_file, " DONE! ",tags.number()," seeds created!\n"
        else:
            print "****************** No seed can be added!\n"
        """            
        #self.fileRemover()
        

#class SeedingClass_v2():
    #"""Class that calculates the uncovered connected components and put seeds in each,
    #Currently used in subvolume_vessel_tracking."""
    #def __init__(self,
                 #minc_file,
                 #tag_file,
                 #connectivity,
                 #bin_threshold,
                 #additional_seeding,
                 #clobber=1):
        
        ##self.p = Pipeline()
        #self.minc_file = minc_file
        #self.tag_file = tag_file
        #self.connectivity = connectivity
        #self.bin_threshold = bin_threshold
        #self.additional_seeding = additional_seeding
        #self.clobber = clobber
        
        
        #self.mnc_bin_file = self.minc_file[:-4]+"_binarize.mnc"
        #try: # the attempt to access the binary minc volume will fail if it doesn't yet exist at pipeline creation
            #self.fileH = pymincf.volumeFromFile(self.mnc_bin_file)
        #except:
            ## if it indeed failed, create it
            #self.binarize_file(self.minc_file,self.mnc_bin_file,self.bin_threshold,"gt")
             
        #self.mnc_dt_file = self.mnc_bin_file[:-4]+"_dt.mnc"
        #self.mnc_group_file = self.mnc_bin_file[:-4]+"_group.mnc"            
        
        #self.buildPipeline()
 
    #def binarize_file(self,binin,binout,binthresh,lg,dtype="byte"):
        ##binarize minc_file "" + 
        #cmd = "mincmath" + " -clobber -quiet" + " -%s" %dtype + " -%s -const %f %s %s" %(lg,binthresh,binin,binout)
        ##print cmd
        #os.system(cmd)
        ##binarizer = CmdStage(cmd)
        ##self.p.addStage(binarizer)       
        ##print self.p.cmd        
        
        
    #def distanceTransfor(self,infile,outfile):
        ##create distance transform (binary input only)
        #cmd = "mincmorph" + " -clobber" + " -short" + " -successive F %s %s" %(infile,outfile)
        ##print cmd
        #os.system(cmd)
        ##distanceTr = CmdStage(cmd)
        ##self.p.addStage(distanceTr)
        ##print self.p.cmd
    
    #def connectedCompontentFinder(self,infile, outfile,connectiv):
        ##create group file
        #cmd = "mincmorph" + " -clobber" + " -short" + " -group -%s %s %s" %(connectiv,infile,outfile)
        ##print cmd
        #os.system(cmd)
        ##createGroup = CmdStage(cmd)
        ##self.p.addStage(createGroup)     
        ##print self.p.cmd
        
    #def seeder(self, dtfile, groupfile, seedfile):
        ##create seed file
        #cmd = "auto_vessel_tag_v2.py --clobber --bin_threshold %f --connectivity %s" %(self.bin_threshold,self.connectivity) + " %s %s %s " %(dtfile,groupfile,seedfile) #+alltags  
        #if self.additional_seeding:
            #cmd += " --self.additional_seeding"
        ##print cmd
        #os.system(cmd)
        ##createSeeds = CmdStage(cmd)
        ##self.p.addStage(createSeeds)     
        
    #def fileRemover(self):
        #if self.remove_interms and self.h5_output!=0:
            #cmd = "rm" + " tmp.log %s %s %s" %(self.graph_iso_obj, self.h5_output[:-3]+"*_iso*.mnc", self.minc_file[:-4]+"_binarize*.mnc")
            ##print cmd
            #os.system(cmd)
            ##rmf = CmdStage(cmd)
            ##self.p.addStage(rmp)
            ##print self.p.cmd
            
        
    ##def seeding(self,mnc_bin_file,tag_file,connectivity,all_tags):
    #def  buildPipeline(self):   
        #print "****************** seeding : ", self.tag_file
        #if self.clobber or not file_exists(self.mnc_dt_file):
            #self.distanceTransfor(self.mnc_bin_file,self.mnc_dt_file)
        #if self.clobber or not file_exists(self.mnc_group_file):
            #self.connectedCompontentFinder(self.mnc_bin_file,self.mnc_group_file,self.connectivity)            
        #if self.clobber or not file_exists(self.tag_file):
            #self.seeder(self.mnc_dt_file, self.mnc_group_file, self.tag_file)  
        
          
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
                    	
