#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#
#
#
#  Created Aug 1, 2013
#  Last modified April 27, 2015 
#  Sahar Ghanavati


from sys import argv
import math, os, shelve, string ,time, sys
import commands
#import matplotlib.pyplot as plt
#from pylab import *		#needed for command find

from scipy import special
from scipy.optimize import fmin,fmin_ncg
import scipy.linalg

from vessel_tracking import graph_analysis
from numpy import *
import numpy
import numpy.matlib as Matlib
from optparse import OptionParser, Option, OptionValueError
from minc_util.progress import progress_report
#import matplotlib.pyplot as plt
from pylab import *		#needed for command find
from operator import itemgetter
import random
import copy

def cyl_volume_calculation(r,h):
    vol = h*(pi*r*r)
    return vol


############################################################################################################################################0
program_name = 'compare_w_groundtruth.py'


if __name__ == '__main__':
        
    #usage = "Usage: "+program_name+" [options] groundtruth.h5 autolabeled.h5 ref_graph.db mr_atlas_file.mnc mr_atlas_centroids.db cba_direction_reference.db output_errormap.db\n"+\
    usage = "Usage: "+program_name+" [options] groundtruth.db autolabeled.db output_errormap.db\n"+\
       "   or  "+program_name+" --help";
            
    parser = OptionParser(usage)
    
    
    parser.add_option("--verbose", action="store_true", dest="verbose", default=0, help="spit out detailed execution report to shell")	

    options, args = parser.parse_args()
    

    if len(args) < 3:
        parser.error("incorrect number of arguments")
        
    print '\n\n>>> %s [%s]: %s\n' % (time.ctime(time.time()), os.getcwd(), string.join(argv))
    print "\n\n",options, "\n\n", args , "\n\n"
    sys.stdout.flush()
    
    tic= time.time()
    
    #groundturht_file, test_file, ref_graph, mr_atlas_file, mr_centroids_file, reference_file, output_file = args
    groundturht_file, test_file, output_file = args


    #cmnd = ("python /projects/souris/sghanavati/src/scripts/cerebrovascular_analysis/cylinder2graph.py %s %s %s %s %s %s --use_cyl_label --clobber" %(groundturht_file,ref_graph, mr_atlas_file, mr_centroids_file, reference_file,groundturht_file[:-2]+"db"))					
    #print "\n", cmnd
    #os.system(cmnd)


    #cmnd = ("python /projects/souris/sghanavati/src/scripts/cerebrovascular_analysis/cylinder2graph.py %s %s %s %s %s %s --use_cyl_label --clobber" %(test_file,ref_graph, mr_atlas_file, mr_centroids_file, reference_file,test_file[:-2]+"db"))					
    #print "\n", cmnd
    #os.system(cmnd)

    ref = graph_analysis.input_graph(groundturht_file[:-2]+"db")
    g = graph_analysis.input_graph(test_file[:-2]+"db")
    h = copy.deepcopy(g)
    
    labeled_num = 0
    error_num = 0
    total_volume = 0
    error_volume = 0
    availabel_labels = []
    
    for e in g.edge_list():
        if e in ref.edge_list():
            if int(ref.edge_property(e,'cyl_label')) not in availabel_labels:
                availabel_labels.append(int(ref.edge_property(e,'cyl_label')))
            labeled_num += 1
            e_radius = ref.edge_property(e,'cyl_radius')
            e_height = ref.edge_property(e,'cyl_height')
            e_volume = 0
            for j in range(len(e_radius)):
                e_volume += cyl_volume_calculation(e_radius[j],e_height[j])
            total_volume += e_volume
            #if int(ref.edge_property(e,'label')) != int(g.edge_property(e,'estimated_label')):
            if int(ref.edge_property(e,'cyl_label')) != int(g.edge_property(e,'cyl_label')):                
                h.set_edge_property(e,'error_label',2)	#error labelling	
                error_volume += e_volume
                error_num += 1
            else:
                h.set_edge_property(e,'error_label',1)	#correct labelling	
        else:
            h.set_edge_property(e,'error_label',2)  #error labelling (new edge in g that does not exist in ref
    graph_analysis.output_graph(output_file, h)   
    
    
    cmd=("\npython /projects/souris/sghanavati/src/scripts/cerebrovascular_analysis/graph2cylinder.py %s %s --use_error_label --clobber " %(output_file,output_file[:-3]+"_error_cyl.h5"))	#python /micehome/jgsled/bin/
    print(cmd)
    os.system(cmd)
    
    if (labeled_num > 0) and (total_volume > 0): #otherwise test set was not labeled
        print ("Out of %d edges that were labeled, %d were labelled correctly with second round of manual labelling. The error is %f%% (volumetric error %f%%). Recognition Rate is %f%% (volumetric RR %f%%)." %(labeled_num,labeled_num-error_num, 100*float(error_num)/float(labeled_num), 100*float(error_volume)/float(total_volume), 100-100*float(error_num)/float(labeled_num), 100-100*float(error_volume)/float(total_volume) )) 

                 
   #### error by label
    if 0 in availabel_labels:
        availabel_labels.remove(0)
    each_label_err={}
    each_label_num={}
    each_label_total_volume={}
    each_label_error_volume={}
    for l in availabel_labels:
        each_label_err[l] = 0
        each_label_num[l] = 0
        each_label_total_volume[l] = 0
        each_label_error_volume[l] = 0

    for e in g.edge_list():
        if e in ref.edge_list():
            l = int(ref.edge_property(e,'cyl_label'))   
            each_label_num[l] += 1
            e_radius = ref.edge_property(e,'cyl_radius')
            e_height = ref.edge_property(e,'cyl_height')
            e_volume = 0
            for j in range(len(e_radius)):
                e_volume += cyl_volume_calculation(e_radius[j],e_height[j])
            each_label_total_volume[l] += e_volume
            if int(ref.edge_property(e,'cyl_label')) != int(g.edge_property(e,'cyl_label')):
                each_label_err[l] += 1
                each_label_error_volume[l] += 1
            
    for l in availabel_labels:
        print  ("\nOut of %d edges (volume %f) with label %d, %d were labelled correctly (correct volume %f). Error rate is %f%% (volumetric error %f%%). Recognition Rate is %f%% (volumetric RR %f%%)." %(each_label_num[l],each_label_total_volume[l],l,each_label_num[l]-each_label_err[l],each_label_total_volume[l]-each_label_error_volume[l],100*float(each_label_err[l])/float(each_label_num[l]),100*float(each_label_error_volume[l])/float(each_label_total_volume[l]),100-100*float(each_label_err[l])/float(each_label_num[l]),100-100*float(each_label_error_volume[l])/float(each_label_total_volume[l]) ))                                                 
        

            



























