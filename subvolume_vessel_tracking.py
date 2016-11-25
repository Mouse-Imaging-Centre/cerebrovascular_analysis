#!/usr/bin/env python

#  Run follow_tree on subvolumes in iterations for full coverage of the vessels in the minc image and merge subtrees
#
#  Created Jan 17, 2012
#  modified Nov 15, 2012	#replace os in exec_follow_tree with subprocess and comment import os
#  modified Jan 27, 2014  #add run time to each of the iterations
#  Last modified May 10, 2014  #add subvolume to automatic seeding and vessel tracking and then in postprocessing we should merge_graphs at the end
#  Sahar Ghanavati


from optparse import OptionGroup, OptionParser, Option, OptionValueError
from sys import argv
#from pydpiper.application import AbstractApplication
#import pydpiper.file_handling as fh
#import Pyro
from os.path import abspath, join, isdir
#import logging
import glob
import fnmatch
import sys
import os,time, string
from sys import argv

#from pydpiper.pipeline import CmdStage, InputFile, OutputFile, LogFile
import auto_vessel_tracking_functions as avt
#import atoms_and_modules.minc_modules as mm
#import atoms_and_modules.minc_atoms as ma
#import atoms_and_modules.stats_tools as st
#import atoms_and_modules.option_groups as og
from optparse import OptionGroup
import csv

program_name = 'subvolume_auto_vessel_tracking.py'

class VesselTrackingApplication():
    def setup_options(self):        
        usage = "Usage: "+program_name+" [options] input_image.mnc output_tree.h5 [original_seed.tag]\n"+\
            "   or  "+program_name+" --help"                
        self.parser = OptionParser(usage)
        #group = OptionGroup(self.parser, "Automatic seeding options", 
                        #"Options for the one round or iterative automatic seeding.")
        #group.add_option("--additional_seeding", action="store_true", dest="additional_seeding",
                        #default=0, help="Add seeds exhaustively")
        #group.add_option("--gradient_seeding", action="store_const", const=1,
                            #dest="seed_type",  default=1,
                            #help="extensive seeding using the second gradient of distance transform for 1 iteration follow_tree (default)")
        #group.add_option("--simple_seeding", action="store_const", const=0,
                            #dest="seed_type",
                            #help="simple seeding for iterative follow_tree")
        #self.parser.add_option_group(group)
        """Add option groups from specific modules"""
        avt.addVTOptionGroup(self.parser)
        
        self.parser.set_usage("%prog [options] input_image.mnc output_tree.h5 [original_seed.tag]\n   or %prog  --help") 
        

        self.options,self.args = self.parser.parse_args() 
        print '\n\n>>> %s [%s]: %s\n' % (time.ctime(time.time()), os.getcwd(), string.join(argv))
        
        if (len(self.args) < 2) or (len(self.args) > 3):
            self.parser.error("incorrect number of arguments")
            
        if ( (len(self.args) != 3) and self.options.automatic_seed== None):
            self.parser.error("either an initial seed file should be given as input or use option automatic_seed")
        
        if ( (len(self.args) == 3) and self.options.automatic_seed):
            self.parser.error("either an initial seed file should be given as input or use option automatic_seed")
        
        if self.options.append and not avt.file_exists(self.options.append):
            raise SystemExit, ("--append specifed yet file %s does not exist." % self.options.append)

        if self.options.append_initnum!=0 and not self.options.append:
            raise SystemExit, ("--append should be used with a non-zero append_initnum!")
        if int(self.options.run) <1:
            raise SystemExit, ("--run should be >=1!")
        
        # check option consistency
        if self.options.normalize and self.options.rescale != None:
            raise SystemExit,("--normalize and --rescale options are incompatible")
        if self.options.extent_step != 1 and self.options.extent_ratio is not None:
            raise SystemExit,("--extent_step and --extent_ratio options are incompatible")
        # set additional flags
        if self.options.psf_z_dependence is not None:
            if self.options.psf_z_dependence[0] == self.options.psf_z_dependence[1]:
                raise SystemExit, \
                "Depth values for --psf_z_dependence must not be identical."
            if self.options.psf is None:
                self.options.psf = array([0.0,0.0,0.0])
                
        self.store_true_options=""
        if self.options.multiscale:
            self.store_true_options+="--multiscale "
        if self.options.minima_branching:
            self.store_true_options+="--minima_branching "
        if self.options.grad_comp:
            self.store_true_options+="--grad_comp "
        if self.options.save_blurs_flag:
            self.store_true_options+="--save_blurs "
        if self.options.save_contours_flag:
            self.store_true_options+="--save_contours "
        if self.options.avoidance:
            self.store_true_options+="--avoidance "
        if self.options.voxel_coord:
            self.store_true_options+="--voxel_coord "
        if self.options.normalize:
            self.store_true_options+="--normalize "
        if self.options.shuffle:
            self.store_true_options+="--shuffle "
        if self.options.clobber:
            self.store_true_options+="--clobber "
        if (self.options.blur_comp==1):
            self.store_true_options+="--blur_comp "
        if (self.options.blur_comp==2):
            self.store_true_options+="--blur_comp_linear "
        
        if self.options.extent_ratio is not None:
            self.store_true_options+=("--extent_ratio %f " %self.options.extent_ratio)
            
        if self.options.psf is not None:
            self.store_true_options+=("--psf %f %f %f " %(self.options.psf[0],self.options.psf[1],self.options.psf[2]))
        if self.options.psf_z_dependence is not None:
            self.store_true_options+=("--psf_z_dependence %f %f %f %f " %(self.options.psf_z_dependence[0],self.options.psf_z_dependence[1],self.options.psf_z_dependence[2],self.options.psf_z_dependence[3]))
        if self.options.psf_comp:
            self.store_true_options+="--psf_comp "
            
        if self.options.rescale:    
            self.store_true_options+=("--rescale %f %f " %(self.options.rescale[0],self.options.rescale[1]))
            
        if self.options.seed_range is not None: 
            self.store_true_options+=("--seed_range %d %d " %(self.options.seed_range[0],self.options.seed_range[1]))
            
        if (self.options.verbose==0):
            self.store_true_options+="--quiet"
        if (self.options.verbose==1):
            self.store_true_options+="--verbose"
        if (self.options.verbose==2):
            self.store_true_options+="--noisy"
                
            

    def setup_appName(self):
        appName = "VesselTracking"
        return appName

    def run(self):
        self.setup_options()
        print "\n######################## Start of program ########################\n"
        print '\n\n>>> %s: %s\n' % (time.ctime(time.time()), string.join(argv))
        #print "\n\n",self.options, "\n\n", self.args , "\n\n"
        
        begintime = time.time()
        all_tags=[]
        
        if len(self.args) == 2:
            print"\n\nautomatic seeding...\n\n"
            minc_file, output_hdf5 = self.args 
            seed_stat = False
            
        elif len(self.args) == 3:
            print"\n\ninputing seed file manually\n\n"
            minc_file, output_hdf5, seed_file = self.args 
            seed_stat = True
            manual_seeds = py_minc.VolumeTags(seed_file)        #=> append them to all_tags
            all_tags=[list(i) for i in manual_seeds.locations]
        
        #### run vesseltracking for each subvolume
        # Create SubVolume MODULE
        # Add the pipeline for creating subvolumes
        subvolsModule = avt.SubVolumer(minc_file,self.options.subvol_nums,self.options.clobber)
        #self.pipeline.addPipeline(subvolsModule.p)
        #subvolume_inputs = avt.create_subvolumes (minc_file,self.options.subvol_nums)
        subvolume_inputs = subvolsModule.input_subvols
        
        for num in range(int(self.options.append_initnum),int(self.options.run)):
            iterationtic = time.time()
            for subv in range(len(subvolume_inputs)):
                subvol_minc_file = subvolume_inputs[subv]
                h5_output = str(output_hdf5[:-3]+"_crop"+str(subv)+"_"+str(num)+".h5")
                output_obj = str(output_hdf5[:-3]+"_crop"+str(subv)+"_"+str(num)+".obj")
                output_avoid = str(output_hdf5[:-3]+"_crop"+str(subv)+"_"+str(num)+"_avoid.mnc")
                
                output_prev_run = 0
                if num > 0:
                    if num == self.options.append_initnum:
                        output_prev_run = self.options.append
                    else:
                        output_prev_run = str(output_hdf5[:-3]+"_crop"+str(subv)+"_"+str(num-1)+".h5")  
                    if not avt.file_exists(output_prev_run):
                        output_prev_run = 0    
                    
                if not seed_stat:   #seed_file not inputed => --automatic_seeding: create the seed file for this round
                    while output_prev_run!=0:     #### wait for the previous vessel tracking to complete execution
                        if  avt.file_exists(output_prev_run):
                            break  
                    seed_file = str(output_hdf5[:-3]+"_crop"+str(subv)+"_"+str(num)+"_seeds.tag")
                    seedingModule = avt.SeedingClass(subvol_minc_file,seed_file,
                                                    self.options.connectivity,self.options.bin_threshold,
                                                    self.options.scale_factor,all_tags,
                                                    self.options.remove_interms,self.options.additional_seeding,
                                                    self.options.laplacian,self.options.clobber,
                                                    output_prev_run)
                    #self.pipeline.addPipeline(seedingModule.p)
                    seed_stat = avt.file_exists(seed_file)
                    all_tags = seedingModule.all_tags
                    
                if seed_stat:
                    # Create followtree MODULE
                    # Add the pipeline for running follow_tree.py
                    #avt.exec_follow_tree(num,subvol_minc_file, seed_file, h5_output,output_obj,output_avoid,self.options.branch_limit,self.options.rho,self.options.gamma,self.options.beta,self.options.smooth,self.options.contrast,self.options.delta,self.options.scale,self.options.avoidance_radius,self.options.extent,self.options.blur_ratio,self.options.dratio,self.options.avoidance_ratio,self.options.epsilon,self.options.limit,self.options.extent_step,self.store_true_options,output_prev_run)
                    if self.options.clobber or not avt.file_exists(h5_output):
                        followTreeModule = avt.FollowTreeer(self.options.followtree_cmndname,num,subvol_minc_file, seed_file, h5_output,output_obj,output_avoid,
                                                        self.options.branch_limit,self.options.rho,self.options.gamma,
                                                        self.options.beta,self.options.smooth,
                                                        self.options.contrast,self.options.delta,self.options.scale,
                                                        self.options.avoidance_radius,self.options.extent,
                                                        self.options.blur_ratio,self.options.dratio,self.options.avoidance_ratio,
                                                        self.options.epsilon, self.options.limit,self.options.extent_step,
                                                        self.store_true_options,
                                                        output_prev_run)
                    #self.pipeline.addPipeline(followTreeModule.p)
                    seed_stat = False   # => in next iteration creates the seed_file 
                    #if not avt.file_exists(h5_output):
                        #print "\nERROR:",h5_output, " couldn't be created!\nAborted!\n"
                        #exit(0)
                #else:
                    #print "No original seed file is inputed or can be created!\nAborted!\n"
                    #exit(0)
            avt.print_timing(iterationtic,time.time(),str(num))   
            
    
        avt.print_timing(begintime,time.time())



#################################


if __name__ == "__main__":

    application = VesselTrackingApplication()
    application.run()
                

        #### merge all at the end
            #in the post_processing convert each crop to db 
            #merge_graphs.py input1.db input2.db ... output.db
