#!/usr/bin/env python
from vessel_tracking import graph_analysis
from cerebrovascular_analysis import vessel_analysis
from optparse import OptionParser, Option, OptionValueError
from sys import argv
import os, shelve, string ,time,sys
import commands
import py_minc
import numpy as np
#---------------------------------------------------------------------------------
#


program_name = 'selectlabelsofgraphs.py'


if __name__ == '__main__':
        
    usage = "Usage: "+program_name+" [options] input1.db input2.db ... inputN.db\n"+\
            "   or  "+program_name+" --help";
            
    parser = OptionParser(usage)

    parser.add_option("--clobber", action="store_true", dest="clobber",
                        default=0, help="overwrite output file")
    
    #parser.add_option("--mask", type="string", dest="mask",
                         #help="the mask minc file")

    parser.add_option("--labels", dest="labels",
                    help="Label numbers to calculate vascular features for, separated by commas, i.e. --labels 3,11,21.",
                    type="string")

    parser.add_option("--output_name", type="string", dest="output_name",
            help="give the name to be added to the end of the input_file name (e.g. label11_30_246 ) ")
                         
    #parser.add_option("--combine_labels", action="store_true", dest="combine_labels",
    #                    default=0, help="combine mother and level1 labels. This should be used without labelLU option, otherwise it's ignored")

    options, args = parser.parse_args()
    
    print '>>> %s [%s]: %s\n' % (time.ctime(time.time()), os.getcwd(), string.join(argv))
    sys.stdout.flush()
    
    input_files = args
    
    #if not options.combine_labels:
        #labelNumeric2Name = {0:"None", 2:"Internal Carotid A. R.", 43:"Internal Carotid A. L.", 102:"Internal Carotid A. R. level1", \
            #143:"Internal Carotid A. L. level1", 122:"Thalamoperforating A.", 8:"PCA R.", 5:"PCA L.", 108:"PCA R. level1", \
            #105:"PCA L. level1", 35:"ACA", 13:"Olfactory A.", 236:"Azygos Anterior C. A.", 191:"MCA R.", 190:"MCA L.", \
            #91:"MCA R. level1", 90:"MCA L. level1", 200:"Posterior Comm. A. R.", 9:"Posterior Comm. A. L.", 7:"Vertebral A. R.", \
            #10:"Vertebral A. L.", 196:"Basilar A.", 96:"Basilar A. level1", 198:"Ventral spinal A.", 68:"SCA R.", 227:"SCA L.", \
            #168:"SCA R. level1", 169:"SCA L. level1", 46:"AICA R.", 12:"AICA L.", 49:"Internal Auditory A. R.", \
            #45:"Internal Auditory A. L.", 14:"Anterior Spinal A.", 15:"Pontine Arteries", 17:"Medial Orbitofrontal A. R.", \
            #18:"Medial Orbitofrontal A. L.", 3:"Paraolivary A. R.", 4:"Paraolivary A. L.", 11:"Superior Saggital Sinus", \
            #111:"Superior Saggital Sinus level1", 6:"Great Cerebral V. of Galen", 206:"Great Cerebral V. of Galen level1", \
            #30:"Transverse Sinus R.", 246:"Transverse Sinus L.", 230:"Transverse Sinus R. level1", 231:"Transverse Sinus L. level1", \
            #192:"Caudal Rhinal V. R.", 34:"Caudal Rhinal V. L.", 20:"Rostral Rhinal V. R.", 21:"Rostral Rhinal V. L.", \
            #120:"Rostral Rhinal V. R. level1", 121:"Rostral Rhinal V. L. level1", 22:"Superior Olfactory Sinus", \
            #101:"Sigmoid Sinus R.", 24:"Sigmoid Sinus L.", 58:"Longitudinal Hippocampal V. R.", 57:"Longitudinal Hippocampal V. L.", \
            #158:"Longitudinal Hippocampal V. R. level1", 157:"Longitudinal Hippocampal V. L. level1", 56:"Thalamostriate V. R.", \
            #54:"Thalamostriate V. L.", 1:"Medial Collicular V. R.", 16:"Medial Collicular V. L.", 170:"Medial Cerebellar  sinus", \
            #171:"Lateral collicular V. R.", 172:"Lateral collicular V. L.", 250:"Lateral Venrtal Cerebellar sinus L.", 249:"Lateral Venrtal Cerebellar sinus R.", \
            #135:"Not Labeled" }
    #else:
        #labelNumeric2Name = {0:"None", 2:"Internal Carotid A. R.", 43:"Internal Carotid A. L.", 102:"Internal Carotid A. R.", \
            #143:"Internal Carotid A. L.", 122:"Thalamoperforating A.", 8:"PCA R.", 5:"PCA L.", 108:"PCA R.", \
            #105:"PCA L.", 35:"ACA", 13:"Azygos Anterior C. A.", 236:"Azygos Anterior C. A.", 191:"MCA R.", 190:"MCA L.", \
            #91:"MCA R.", 90:"MCA L.", 200:"Posterior Comm. A. R.", 9:"Posterior Comm. A. L.", 7:"Vertebral A. R.", \
            #10:"Vertebral A. L.", 196:"Basilar A.", 96:"Basilar A.", 198:"Ventral spinal A.", 68:"SCA R.", 227:"SCA L.", \
            #168:"SCA R.", 169:"SCA L.", 46:"AICA R.", 12:"AICA L.", 49:"Internal Auditory A. R.", \
            #45:"Internal Auditory A. L.", 14:"Anterior Spinal A.", 15:"Pontine Arteries", 17:"Medial Orbitofrontal A. R.", \
            #18:"Medial Orbitofrontal A. L.", 3:"Paraolivary A. R.", 4:"Paraolivary A. L.", 11:"Superior Saggital Sinus", \
            #111:"Superior Saggital Sinus", 6:"Great Cerebral V. of Galen", 206:"Great Cerebral V. of Galen", \
            #30:"Transverse Sinus R.", 246:"Transverse Sinus L.", 230:"Transverse Sinus R.", 231:"Transverse Sinus L.", \
            #192:"Caudal Rhinal V. R.", 34:"Caudal Rhinal V. L.", 20:"Rostral Rhinal V. R.", 21:"Rostral Rhinal V. L.", \
            #120:"Rostral Rhinal V. R.", 121:"Rostral Rhinal V. L.", 22:"Superior Olfactory Sinus", \
            #101:"Sigmoid Sinus R.", 24:"Sigmoid Sinus L.", 58:"Longitudinal Hippocampal V. R.", 57:"Longitudinal Hippocampal V. L.", \
            #158:"Longitudinal Hippocampal V. R.", 157:"Longitudinal Hippocampal V. L.", 56:"Thalamostriate V. R.", \
            #54:"Thalamostriate V. L.", 1:"Medial Collicular V. R.", 16:"Medial Collicular V. L.", 170:"Medial Cerebellar  sinus", \
            #171:"Lateral collicular V. R.", 172:"Lateral collicular V. L.", 250:"Lateral Venrtal Cerebellar sinus L.", 249:"Lateral Venrtal Cerebellar sinus R.", \
            #135:"Not Labeled" }

    ls = options.labels.split(',')
    ls = [int(l) for l in ls]
    history = '>>> %s: %s' % (time.ctime(time.time()), string.join(argv))

    f = open(options.output_name+".csv",'w')
    f.write( "strain,sample, label, seg_num, total_length_mm, mean_diameter_um, sd_diameter_um") 

    for input_file in input_files:
        output_file = input_file[:-3]+"_"+options.output_name+".db"
        g, attributes = graph_analysis.input_graph(input_file, ["history", "vertex_offsets"])
        edge2rm = []
        for e in g.edge_list():
            if g.edge_property(e,'label') not in ls:
                edge2rm.append(e)
        for e in edge2rm:
            g.remove_edge(e)
        
        if attributes.has_key("history"):
            history = attributes["history"] + "\n" + history
            del attributes['history']

        graph_analysis.output_graph(output_file, g, history, attributes)
        
        cmd=("graph2obj.py %s %s %s " %(output_file,output_file[:-3]+"obj", options.clobber)) 
        os.system(cmd)

        cmd=("graph2cylinder.py --use_label %s %s %s " %(output_file,output_file[:-3]+".h5", options.clobber)) 
        os.system(cmd)

        num_seg = {}
        diameter = {}
        length = {}
        for l in ls:
            num_seg[l]=0
            diameter[l]=[]
            length[l]=0
        
        for e in g.edge_list():
            l = g.edge_property(e,'label')
            num_seg[l] +=1
            diameter[l].append(g.edge_property(e,'diameter'))
            length[l] += g.edge_property(e,'length')
            
        for l in ls:
            f.write( "\n"+os.path.basename(input_file)[:3]+ ","+os.path.basename(input_file)[:14]+ ","+ str(l) +","+ str(num_seg[l])+","+str(length[l]) +","+str(np.mean(diameter[l]))+","+str(np.std(diameter[l])))
 
    f.close()
            
        
 
 


 
 
 
 
 
