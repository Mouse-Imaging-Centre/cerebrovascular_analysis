#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#
#
#
#  Created Dec 15, 2015
#  Last modified  
#  Sahar Ghanavati


from sys import argv
import math, os, shelve, string ,time, sys
import commands
#import matplotlib.pyplot as plt
#from pylab import *        #needed for command find

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
from pylab import *     #needed for command find
from operator import itemgetter
import random
import copy

def cyl_volume_calculation(r,h):
    vol = h*(pi*r*r)
    return vol


############################################################################################################################################0
program_name = 'compare_w_groundtruth_tocsv.py'


if __name__ == '__main__':
        
    usage = "Usage: "+program_name+" [options] groundtruth.db autolabeled_initial.db autolabeled_GD.db autolabeled_S.db output.csv\n"+\
       "   or  "+program_name+" --help";
            
    parser = OptionParser(usage)
    
    
    parser.add_option("--combine_labels",type="string", dest="combine_labels",
                        help="combine mother and level1 labels. This should be used without labelLU option, otherwise it's ignored")
    parser.add_option("--labels", dest="labels",
                    help="Label numbers to calculate vascular features for, separated by commas, i.e. --labels 3,11,21.",
                    type="string")
    parser.add_option("--labelLU", type="string", dest="labelLU",
            help="give the name of labelLU.config file")
    parser.add_option("--append", action="store_true", dest="append",
                        default=0, help="add to the output file")

    options, args = parser.parse_args()
    

    if len(args) < 5:
        parser.error("incorrect number of arguments")
        
    print '\n\n>>> %s [%s]: %s\n' % (time.ctime(time.time()), os.getcwd(), string.join(argv))
    print "\n\n",options, "\n\n", args , "\n\n"
    sys.stdout.flush()
    
    tic= time.time()
    
    #groundturht_file, test_file, ref_graph, mr_atlas_file, mr_centroids_file, reference_file, output_file = args
    groundturht_file, test_file,  test_file_gd, test_file_s, output_file = args

    options, args = parser.parse_args()
    
    print '>>> %s [%s]: %s\n' % (time.ctime(time.time()), os.getcwd(), string.join(argv))
    sys.stdout.flush()
    
    if options.append:
        f = open(output_file,'a')
    else:
        f = open(output_file,'w')
        f.write( "strain,sample, label, label_Name, initial_RR,GD_RR,SA_RR, initialvol_RR,GDvol_RR,SAvol_RR") 
     
    if not options.labelLU:
        labelNumeric2Name = {0:"None", 2:"Internal Carotid A. R.", 43:"Internal Carotid A. L.", 102:"Internal Carotid A. R.", \
            143:"Internal Carotid A. L.", 122:"Thalamoperforating A.", 8:"PCA R.", 5:"PCA L.", 108:"PCA R.", \
            105:"PCA L.", 35:"ACA", 13:"Azygos Anterior C. A.", 236:"Azygos Anterior C. A.", 191:"MCA R.", 190:"MCA L.", \
            91:"MCA R.", 90:"MCA L.", 200:"Posterior Comm. A. R.", 9:"Posterior Comm. A. L.", 7:"Vertebral A. R.", \
            10:"Vertebral A. L.", 196:"Basilar A.", 96:"Basilar A.", 198:"Ventral spinal A.", 68:"SCA R.", 227:"SCA L.", \
            168:"SCA R.", 169:"SCA L.", 46:"AICA R.", 12:"AICA L.", 49:"Internal Auditory A. R.", \
            45:"Internal Auditory A. L.", 14:"Anterior Spinal A.", 15:"Pontine Arteries", 17:"Medial Orbitofrontal A. R.", \
            18:"Medial Orbitofrontal A. L.", 3:"Paraolivary A. R.", 4:"Paraolivary A. L.", 11:"Superior Saggital Sinus", \
            111:"Superior Saggital Sinus", 6:"Great Cerebral V. of Galen", 206:"Great Cerebral V. of Galen", \
            30:"Transverse Sinus R.", 246:"Transverse Sinus L.", 230:"Transverse Sinus R.", 231:"Transverse Sinus L.", \
            192:"Caudal Rhinal V. R.", 34:"Caudal Rhinal V. L.", 20:"Rostral Rhinal V. R.", 21:"Rostral Rhinal V. L.", \
            120:"Rostral Rhinal V. R.", 121:"Rostral Rhinal V. L.", 22:"Superior Olfactory Sinus", \
            101:"Sigmoid Sinus R.", 24:"Sigmoid Sinus L.", 58:"Longitudinal Hippocampal V. R.", 57:"Longitudinal Hippocampal V. L.", \
            158:"Longitudinal Hippocampal V. R.", 157:"Longitudinal Hippocampal V. L.", 56:"Thalamostriate V. R.", \
            54:"Thalamostriate V. L.", 1:"Medial Collicular V. R.", 16:"Medial Collicular V. L.", 170:"Medial Cerebellar  sinus", \
            171:"Lateral collicular V. R.", 172:"Lateral collicular V. L.", 250:"Lateral Venrtal Cerebellar sinus L.", 249:"Lateral Venrtal Cerebellar sinus R.", \
            135:"Not Labeled" }
    else:
        f = open(options.labelLU, 'r')
        lines=[]
        for line in f:
            if ((not line[0]=='#') and (not line=='')):
                lines.append(line)
        labelNumeric2Name ={}
        for l in lines:
            i0=l.index(';')
            i1=l[i0+1:].index(';')+i0+1
            labelNumeric2Name[int(l[0:i0])]=l[i0+1:i1]
        f.close()   
        
    if options.labels:
        ls = options.labels.split(',')
        ls = [int(l) for l in ls]
    else:
        ls = labelNumeric2Name.keys()
        
    if options.combine_labels:
        lpairs = options.combine_labels.split(',')
        labelpairs = {}
        for pair in lpairs:
            plabels = pair.split(':')
            plabels = [int(l) for l in plabels]
            labelpairs[plabels[0]]=plabels[1]
        
    ref = graph_analysis.input_graph(groundturht_file)
    g = graph_analysis.input_graph(test_file)
    g_gd = graph_analysis.input_graph(test_file_gd)
    g_s = graph_analysis.input_graph(test_file_s)
   
    for e in g.edge_list():
        if int(g.edge_property(e,'cyl_label')) in  labelpairs.keys():
            g.set_edge_property(e,'cyl_label',labelpairs[int(g.edge_property(e,'cyl_label'))])
            
    for e in g_gd.edge_list():
        if int(g_gd.edge_property(e,'cyl_label')) in  labelpairs.keys():            
            g_gd.set_edge_property(e,'cyl_label',labelpairs[int(g_gd.edge_property(e,'cyl_label'))])
            
    for e in g_s.edge_list():
        if int(g_s.edge_property(e,'cyl_label')) in  labelpairs.keys():
            g_s.set_edge_property(e,'cyl_label',labelpairs[int(g_s.edge_property(e,'cyl_label'))])
            
    for e in ref.edge_list():
        if int(ref.edge_property(e,'cyl_label')) in  labelpairs.keys():
            ref.set_edge_property(e,'cyl_label',labelpairs[int(ref.edge_property(e,'cyl_label'))])
            
    labeled_num = {}
    labeled_vol = {}
    rr_labeled_num = {}
    rr_labeled_vol = {}
    rr_labeled_num_gd = {}
    rr_labeled_vol_gd = {}
    rr_labeled_num_s = {}
    rr_labeled_vol_s = {}
    
    for l in ls:
        labeled_num[l] = 0
        labeled_vol[l] = 0
        rr_labeled_num[l] = 0
        rr_labeled_vol[l] = 0
        rr_labeled_num_gd[l] = 0
        rr_labeled_vol_gd[l] = 0
        rr_labeled_num_s[l] = 0
        rr_labeled_vol_s[l] = 0
        
    
    for e in g.edge_list():
        if e in ref.edge_list():
            if int(ref.edge_property(e,'cyl_label')) in ls:
                labeled_num[int(ref.edge_property(e,'cyl_label'))] += 1
                e_radius = ref.edge_property(e,'cyl_radius')
                e_height = ref.edge_property(e,'cyl_height')
                e_volume = 0
                for j in range(len(e_radius)):
                    e_volume += cyl_volume_calculation(e_radius[j],e_height[j])
                labeled_vol[int(ref.edge_property(e,'cyl_label'))]  += e_volume
                if int(ref.edge_property(e,'cyl_label')) == int(g.edge_property(e,'cyl_label')):                
                    rr_labeled_num[int(ref.edge_property(e,'cyl_label'))] += 1
                    rr_labeled_vol[int(ref.edge_property(e,'cyl_label'))]  += e_volume
                if int(ref.edge_property(e,'cyl_label')) == int(g_gd.edge_property(e,'cyl_label')):                
                    rr_labeled_num_gd[int(ref.edge_property(e,'cyl_label'))] += 1
                    rr_labeled_vol_gd[int(ref.edge_property(e,'cyl_label'))]  += e_volume
                if int(ref.edge_property(e,'cyl_label')) == int(g_s.edge_property(e,'cyl_label')):                
                    rr_labeled_num_s[int(ref.edge_property(e,'cyl_label'))] += 1
                    rr_labeled_vol_s[int(ref.edge_property(e,'cyl_label'))]  += e_volume
    
    for l in ls:
        if labeled_num[l] > 0:
            f.write( "\n"+os.path.basename(groundturht_file)[:3]+ ","+os.path.basename(groundturht_file)[:14]+ ","+ str(l)+ ","+ labelNumeric2Name[l]+ ","+ str(100.0*float(rr_labeled_num[l])/float(labeled_num[l]))+","+str(100.0*float(rr_labeled_vol[l])/float(labeled_num[l])) +","+ str(100.0*float(rr_labeled_num_gd[l])/float(labeled_num[l])) + ","+ str(100.0*float(rr_labeled_vol_gd[l])/float(labeled_vol[l]))+","+str(100.0*float(rr_labeled_num_s[l])/float(labeled_vol[l])) +","+ str(100.0*float(rr_labeled_vol_s[l])/float(labeled_vol[l])) ) 
     
    f.write( "\n"+os.path.basename(groundturht_file)[:3]+ ","+os.path.basename(groundturht_file)[:14]+ ",0,total,"+ str(100.0*float(sum(rr_labeled_num.values()))/float(sum(labeled_num.values())))+","+str(100.0*float(sum(rr_labeled_vol.values()))/float(sum(labeled_num.values()))) +","+ str(100.0*float(sum(rr_labeled_num_gd.values()))/float(sum(labeled_num.values()))) + ","+ str(100.0*float(sum(rr_labeled_vol_gd.values()))/float(sum(labeled_vol.values())))+","+str(100.0*float(sum(rr_labeled_num_s.values()))/float(sum(labeled_vol.values()))) +","+ str(100.0*float(sum(rr_labeled_vol_s.values()))/float(sum(labeled_vol.values()))) ) 
    f.close()  
    
           