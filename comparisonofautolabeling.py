#!/usr/bin/env python
from vessel_tracking import graph_analysis
import vessel_analysis
from optparse import OptionParser, Option, OptionValueError
from sys import argv
import os, shelve, string ,time,sys
import commands
import py_minc
from numpy import *

#---------------------------------------------------------------------------------
#
def cyl_volume_calculation(r,h):
    vol = h*(pi*r*r)
    return vol


def totallabelcomparison(g,h,labelNumeric2Name,property='cyl_label'):
    labeled_num = 0
    error_num = 0
    total_volume = 0
    error_volume = 0
    
    for e in g.edge_list():
        if e in h.edge_list():
            labeled_num += 1
            e_radius = h.edge_property(e,'cyl_radius')
            e_height = h.edge_property(e,'cyl_height')
            e_volume = 0
            for j in range(len(e_radius)):
                e_volume += cyl_volume_calculation(e_radius[j],e_height[j])
            total_volume += e_volume
            if labelNumeric2Name[int(h.edge_property(e,property))] != labelNumeric2Name[int(g.edge_property(e,property))]:        #for combining the labels        
                error_volume += e_volume
                error_num += 1
    if labeled_num>0 and total_volume> 0:
        return [100.0*(1.0-(float(error_num)/float(labeled_num))), 100.0*(1.0-(float(error_volume)/float(total_volume)))]
    else:
        return [0,0]

def labelcomparison(g,h,l,labelNumeric2Name,property='cyl_label'):
    labeled_num = 0
    error_num = 0
    total_volume = 0
    error_volume = 0
    
    for e in g.edge_list():
        if e in h.edge_list():
            if int(g.edge_property(e,property))==l:
                labeled_num += 1
                e_radius = h.edge_property(e,'cyl_radius')
                e_height = h.edge_property(e,'cyl_height')
                e_volume = 0
                for j in range(len(e_radius)):
                    e_volume += cyl_volume_calculation(e_radius[j],e_height[j])
                total_volume += e_volume
                if labelNumeric2Name[int(h.edge_property(e,property))] != labelNumeric2Name[int(g.edge_property(e,property))]:        #for combining the labels                 
                    error_volume += e_volume
                    error_num += 1
    if labeled_num>0 and total_volume> 0:
        return [100.0*(1.0-(float(error_num)/float(labeled_num))), 100.0*(1.0-(float(error_volume)/float(total_volume)))]
    else:
        return [0,0]
    
    
    

program_name = 'comparisonofautolabeling.py'


if __name__ == '__main__':
        
    usage = "Usage: "+program_name+" [options] input1.db input2.db ... inputN.db\n"+\
            "   or  "+program_name+" --help";
            
    parser = OptionParser(usage)

    #parser.add_option("--clobber", action="store_true", dest="clobber",
                        #default=0, help="overwrite output file")
    
    parser.add_option("--labels", dest="labels",
                    help="Label numbers to calculate vascular features for, separated by commas, i.e. --labels 3,11,21.",
                    type="string")
    parser.add_option("--labelLU", type="string", dest="labelLU",
            help="give the name of labelLU.config file")
    parser.add_option("--output_name", type="string", dest="output_name",
            help="give the name for the output csv file")
                         
    parser.add_option("--combine_labels", action="store_true", dest="combine_labels",
                        default=0, help="combine mother and level1 labels. This should be used without labelLU option, otherwise it's ignored")

    options, args = parser.parse_args()
    
    print '>>> %s [%s]: %s\n' % (time.ctime(time.time()), os.getcwd(), string.join(argv))
    sys.stdout.flush()
    
    input_files = args


    if not options.labelLU:
        #labelNumeric2Name ={0:"No label",35:"Anterior Cerebral Artery", 191:"R. Middle Cerebral Artry", 190:"L. Middle Cerebral Artry", 2:"R. Intern Carotid Artery", 43:"L. Intern Carotid Artery", 200:"R. Posterior Comm. Artry", 9:"L. Posterior Comm. Artry", 8:"R. Posterior Cereb Artry", 5:"L. Posterior Cereb Artry", 68:"R. Superior Cereb Artery", 227:"L. Superior Cereb Artery", 46:"R. Ant. Inf. Cereb Artry", 12:"L. Ant. Inf. Cereb Artry", 196:"Basilar Artery", 7:"Vertebral Artery", 49:"R. Internal Audit Artery", 45:"L. Internal Audit Artery", 7: "Vertebral Artery" , 3: "R. Paraolivary Artry", 4: "L. Paraolivary Artry" , 11:"Superior Saggital Sinus", 6:"Great Cerbral Vein Galen", 30:"R. Transverse Sinus", 246:"L. Transverse Sinus", 192:"R. Caudal Rhinal Vein", 34:"L. Caudal Rhinal Vein", 20:"R. Rostral Rhinal Vein", 21:"L. Rostral Rhinal Vein", 101:"R. Sigmoid Sinus", 24:"L. Sigmoid Sinus", 58:"R. Longitud. Hippo. Vein", 57:"L. Longitud. Hippo. Vein", 56:"R. Thalamostriate Vein",\
        #54:"L. Thalamostriate Vein", 1:"R. Medial Colicular Vein", 16:"L. Medial Colicular Vein", 170:"Unknown Sinus/Vein #01", 171:"R. Lateral collicular V.", 172:"L. Lateral collicular V.", 250:"L. Unknown Sinus/Vein #2", 251:"R. Unknown Sinus/Vein #2"}                                            
        if not options.combine_labels:
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
        ls = []
        ls_name = []
        for l in labelNumeric2Name.keys():
            if labelNumeric2Name[l] not in ls_name:
                ls.append(l)
                ls_name.append(labelNumeric2Name[l])
    #ls = labelNumeric2Name.keys()
        
    f = open(options.output_name,'w')
    f.write( "strain,sample, label, label_Name, initial_RR,GD_RR,SA_RR, initialvol_RR,GDvol_RR,SAvol_RR") 
    ls = [int(l) for l in ls]
    for input_file in input_files:
        if os.path.exists(input_file) and os.path.exists(input_file[:-3]+"_autolabel_initial_cyl.db") and os.path.exists(input_file[:-3]+"_autolabel_GD_cyl.db") and os.path.exists(input_file[:-3]+"_autolabel_S_cyl.db"):
            g = graph_analysis.input_graph(input_file)
            ginit = graph_analysis.input_graph(input_file[:-3]+"_autolabel_initial_cyl.db")
            ggd = graph_analysis.input_graph(input_file[:-3]+"_autolabel_GD_cyl.db")
            gsa = graph_analysis.input_graph(input_file[:-3]+"_autolabel_S_cyl.db")
            f.write( "\n"+os.path.basename(input_file)[:3]+ ","+os.path.basename(input_file)[:14]+ ",0,total,"+ str(totallabelcomparison(g,ginit,labelNumeric2Name,'cyl_label')[0])+","+str(totallabelcomparison(g,ggd,labelNumeric2Name,'cyl_label')[0]) +","+ str(totallabelcomparison(g,gsa,labelNumeric2Name,'cyl_label')[0]) +"," + str(totallabelcomparison(g,ginit,labelNumeric2Name,'cyl_label')[1])+","+str(totallabelcomparison(g,ggd,labelNumeric2Name,'cyl_label')[1]) +","+ str(totallabelcomparison(g,gsa,labelNumeric2Name,'cyl_label')[1])) 
            for l in ls:
                f.write( "\n"+os.path.basename(input_file)[:3]+ ","+os.path.basename(input_file)[:14]+ ","+ str(l)+ ","+ labelNumeric2Name[l]+ ","+ str(labelcomparison(g,ginit,l,labelNumeric2Name,'cyl_label')[0])+","+str(labelcomparison(g,ggd,l,labelNumeric2Name,'cyl_label')[0]) +","+ str(labelcomparison(g,gsa,l,labelNumeric2Name,'cyl_label')[0]) + ","+ str(labelcomparison(g,ginit,l,labelNumeric2Name,'cyl_label')[1])+","+str(labelcomparison(g,ggd,l,labelNumeric2Name,'cyl_label')[1]) +","+ str(labelcomparison(g,gsa,l,labelNumeric2Name,'cyl_label')[1]) ) 
    f.close()  
    
