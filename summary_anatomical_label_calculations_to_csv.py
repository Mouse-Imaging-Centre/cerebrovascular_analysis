# -*- coding: utf-8 -*-
"""
Created on Wed Nov 17 21:29:09 2015

@author: sghanavati
"""

#!/usr/bin/env python

# Get manuallabel_graph.db  and corresponding CT.mnc
# initialize the CT.mnc to 0 everywhere
# for each vertex in the graph give the label to the mnc voxels in a sphere around the graph vertex of vertex radius 
# it will save the labeled CT.mnc
#  




from vessel_tracking import graph_analysis, distribution_analysis
#from graph_labeling import *
from optparse import OptionParser, Option, OptionValueError
from sys import argv
import os, commands, time, string, math
import numpy as np
import scipy
import operator
from optparse import OptionParser, Option, OptionValueError
#from time import time, ctime
from sys import argv
import shelve, os, string
from copy import deepcopy
from numpy import *
from scipy.linalg import norm
from scipy.interpolate import splprep, splev, splrep
from vessel_tracking import path_io
from vessel_tracking import graph_analysis
import vessel_analysis
import copy

program_name = 'summary_anatomical_label_calculations_to_csv.py'


if __name__ == '__main__':
        
    usage = "Usage: "+program_name+" [options]  input1.db input2.db ... inputN.db [output_graph.csv] \n"+\
            "   or  "+program_name+" --help";
            
    parser = OptionParser(usage)
    
    parser.add_option("--clobber", action="store_true", dest="clobber",
                    default=0, help="overwrite output file")
    parser.add_option("--labels", dest="labels",
                    help="Label numbers to keep in the output_file, separated by commas, i.e. --labels 3,11,21.",
                    type="string")
    parser.add_option("--labelLU", type="string", dest="labelLU",
            help="give the name of labelLU.config file")
    parser.add_option("--output_name", type="string", dest="output_name",
            help="give the name for the output csv file")
        
        
                
    options, args = parser.parse_args()
    print '\n\n>>> %s [%s]: %s\n' % (time.ctime(time.time()), os.getcwd(), string.join(argv))
    print "\n\n",options, "\n\n", args , "\n\n"
   
    input_files = args
    #output_file = input_files[0][:-3]+"_anatomicalcalcs.csv"

        
    #if len(args) == 1:
        #input_file = args[0]
    #elif len(args) == 2:
        #input_file, output_file = args
    #else:
        #parser.error("incorrect number of arguments")

    if not options.clobber and os.path.exists(options.output_name):
        raise SystemExit, \
              "The --clobber option is needed to overwrite an existing file."
          
    
    
    history = '\n>>> %s: %s' % (time.ctime(time.time()), string.join(argv))

    ls = []
    if options.labels:
        ls = options.labels.split(',')
        ls = [int(l) for l in ls]
    
    if options.labelLU and not os.path.exists(options.labelLU):
        raise SystemExit, ("--labelLU specifed yet file %s does not exist." % options.labelLU)
    if not options.labelLU:
        #labelNumeric2Name ={0:"No label",35:"Anterior Cerebral Artery", 191:"R. Middle Cerebral Artry", 190:"L. Middle Cerebral Artry", 2:"R. Intern Carotid Artery", 43:"L. Intern Carotid Artery", 200:"R. Posterior Comm. Artry", 9:"L. Posterior Comm. Artry", 8:"R. Posterior Cereb Artry", 5:"L. Posterior Cereb Artry", 68:"R. Superior Cereb Artery", 227:"L. Superior Cereb Artery", 46:"R. Ant. Inf. Cereb Artry", 12:"L. Ant. Inf. Cereb Artry", 196:"Basilar Artery", 7:"Vertebral Artery", 49:"R. Internal Audit Artery", 45:"L. Internal Audit Artery", 7: "Vertebral Artery" , 3: "R. Paraolivary Artry", 4: "L. Paraolivary Artry" , 11:"Superior Saggital Sinus", 6:"Great Cerbral Vein Galen", 30:"R. Transverse Sinus", 246:"L. Transverse Sinus", 192:"R. Caudal Rhinal Vein", 34:"L. Caudal Rhinal Vein", 20:"R. Rostral Rhinal Vein", 21:"L. Rostral Rhinal Vein", 101:"R. Sigmoid Sinus", 24:"L. Sigmoid Sinus", 58:"R. Longitud. Hippo. Vein", 57:"L. Longitud. Hippo. Vein", 56:"R. Thalamostriate Vein",\
        #54:"L. Thalamostriate Vein", 1:"R. Medial Colicular Vein", 16:"L. Medial Colicular Vein", 170:"Unknown Sinus/Vein #01", 171:"R. Lateral collicular V.", 172:"L. Lateral collicular V.", 250:"L. Unknown Sinus/Vein #2", 251:"R. Unknown Sinus/Vein #2"}                                            
        labelNumeric2Name = {0:"None", 2:"Internal Carotid A. R.", 43:"Internal Carotid A. L.", 102:"Internal Carotid A. R. level1", \
            143:"Internal Carotid A. L. level1", 122:"Thalamoperforating A.", 8:"PCA R.", 5:"PCA L.", 108:"PCA R. level1", \
            105:"PCA L. level1", 35:"ACA", 13:"Olfactory A.", 236:"Azygos Anterior C. A.", 191:"MCA R.", 190:"MCA L.", \
            91:"MCA R. level1", 90:"MCA L. level1", 200:"Posterior Comm. A. R.", 9:"Posterior Comm. A. L.", 7:"Vertebral A. R.", \
            10:"Vertebral A. L.", 196:"Basilar A.", 96:"Basilar A. level1", 198:"Ventral spinal A.", 68:"SCA R.", 227:"SCA L.", \
            168:"SCA R. level1", 169:"SCA L. level1", 46:"AICA R.", 12:"AICA L.", 49:"Internal Auditory A. R.", \
            45:"Internal Auditory A. L.", 14:"Anterior Spinal A.", 15:"Pontine Arteries", 17:"Medial Orbitofrontal A. R.", \
            18:"Medial Orbitofrontal A. L.", 3:"Paraolivary A. R.", 4:"Paraolivary A. L.", 11:"Superior Saggital Sinus", \
            111:"Superior Saggital Sinus level1", 6:"Great Cerebral V. of Galen", 206:"Great Cerebral V. of Galen level1", \
            30:"Transverse Sinus R.", 246:"Transverse Sinus L.", 230:"Transverse Sinus R. level1", 231:"Transverse Sinus L. level1", \
            192:"Caudal Rhinal V. R.", 34:"Caudal Rhinal V. L.", 20:"Rostral Rhinal V. R.", 21:"Rostral Rhinal V. L.", \
            120:"Rostral Rhinal V. R. level1", 121:"Rostral Rhinal V. L. level1", 22:"Superior Olfactory Sinus", \
            101:"Sigmoid Sinus R.", 24:"Sigmoid Sinus L.", 58:"Longitudinal Hippocampal V. R.", 57:"Longitudinal Hippocampal V. L.", \
            158:"Longitudinal Hippocampal V. R. level1", 157:"Longitudinal Hippocampal V. L. level1", 56:"Thalamostriate V. R.", \
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
    
    f = open(options.output_name,'w')
    f.write( "strain,sample, label, label_Name, number_of_segments,total_distlength_mm,mean_diameter_nooutlier_um, total_correctedvol_mm3, total_volumemm3, total_length_mm,min_diamter_um, max_diameter_um,mean_diameter_um, sd_diameter_um, mean_curvature ,mean_tortuosity") 
    
    for input_file in input_files:
        g, attributes = graph_analysis.input_graph(input_file, ["history", "vertex_offsets"])    
        for l in ls:
            #f.write( "\n"+os.path.basename(input_file)[:3]+ ","+os.path.basename(input_file)[:14]+ ","+ str(l)+ ","+ labelNumeric2Name[l] )
            f.write( "\n"+os.path.basename(input_file)[:3]+ ","+os.path.basename(input_file)[:14]+ ","+ str(l)+ ","+ labelNumeric2Name[l]+ ","+ str(vessel_analysis.num_segment_by_label(g,l,'label'))+","+str(vessel_analysis.total_feature_by_label(g,l,'dist_length', 'label')) +","+ str(vessel_analysis.mean_feature_by_label(g,l, 'nooutlier_diameter', 'label')) +","+ str(vessel_analysis.corrected_volume_by_label(g,l,'label'))+","+ str(vessel_analysis.volume_by_label(g,l,'label'))+","+str(vessel_analysis.total_feature_by_label(g,l,'length', 'label')) +","+str(vessel_analysis.min_feature_by_label(g,l, 'diameter', 'label'))+","+str(vessel_analysis.max_feature_by_label(g,l, 'diameter', 'label'))+","+ str(vessel_analysis.mean_feature_by_label(g,l, 'diameter', 'label'))+"," +str(vessel_analysis.sd_feature_by_label(g,l, 'diameter', 'label')) ) #+","+ str(vessel_analysis.mean_feature_by_label(g,l, 'curvature', 'label'))+","+ str(vessel_analysis.mean_feature_by_label(g,l, 'tortuosity', 'label')) ) 
        
    f.close()     
