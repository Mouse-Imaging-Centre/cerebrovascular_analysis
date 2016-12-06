#!/usr/bin/env python
from vessel_tracking import graph_analysis
from cerebrovascular_analysis import vessel_analysis
from optparse import OptionParser, Option, OptionValueError
from sys import argv
import os, shelve, string ,time,sys
import commands
#---------------------------------------------------------------------------------
#

#122;Thalamoperforating A.;50.0;50.0;50.0;   //unknown
#8;PCA R.;255.0; 170.0; 255.0;    //pink 
#5;PCA L.;255.0; 85.0; 127.0;    //pink 
#108;PCA R. level1; 200.0; 170.0; 255.0;    //pink
#105;PCA L. level1;200.0; 85.0; 127.0;    //pink
#35;ACA;255.0; 0.0; 0.0;    //red
#135;ACA level1;255.0; 0.0; 50.0;    //red //(this should branch into all in the paper of adrienne!)
#13;Olfactory A.;;0.0; 85.0; 0.0;    //dark green
#236;Azygos Anterior C. A.;200.0;100.0;100.0;    //unknown
#191;MCA R.;255.0; 255.0; 150.0;    //yellow
#190;MCA L.;255.0; 255.0; 0.0;    //yellow 
#91;MCA R. level1;200.0; 255.0; 150.0;    //yellow
#90;MCA L. level1;200.0; 255.0; 0.0;    //yellow
#68;SCA R.;255.0; 170.0; 127.0;    //orange 
#227;SCA L.;255.0; 100.0; 0.0;     //orange
#168;SCA R. level1;220.0; 170.0; 127.0;    //orange
#169;SCA L. level1;220.0; 100.0; 0.0;     //orange
#46;AICA R.;0.0; 80.0; 0.0;    //dark green 
#12;AICA L.;0.0; 85.0; 80.0;     //dark green
#49;Internal Auditory A. R.;100.0; 100.0; 150.0;     //light green
#45;Internal Auditory A. L.;120.0; 120.0; 180.0;     //light green
#14;Anterior Spinal A.;;0.0; 0.0; 255.0;     //blue
#15;Pontine Arteries;185.0; 255.0; 255.0;    //light blue
#17;Medial Orbitofrontal A. R.;200.0; 0.0; 0.0;    //red
#18;Medial Orbitofrontal A. L.;200.0; 30.0; 0.0;    //red
#3;Paraolivary A. R.;150.0;150.0;150.0;      //uknown
#4;Paraolivary A. L.;100.0;100.0;100.0;  //unknown
#111;Superior Saggital Sinus level1;170.; 85.0; 0.0;    //red brown

#remap the labels:
#8->108
#5->105
#135->35
#111->35
#236->35
#191->91
#190->90
#68->168
#227->169

program_name = 'perfusionterritory.py'


if __name__ == '__main__':
        
    usage = "Usage: "+program_name+" [options] input_graph mr_atlas_file [output_mnc]\n"+\
            "   or  "+program_name+" --help";
            
    parser = OptionParser(usage)

    parser.add_option("--clobber", action="store_true", dest="clobber",
                        default=0, help="overwrite output file")
    
    parser.add_option("--veins", action="store_true", dest="veins",
                        default=0, help="Create drainage maps for veins")
    
    parser.add_option("--mask", type="string", dest="mask",
                         help="the mask minc file")

    parser.add_option("--remap", dest="remap",
                    help="Label numbers to replace with new label numbers. For example, to turn label 3 into label 11, it would be --remap 3:11. More than one can be specified at once, separated by commas, i.e. --remap 3:11,21:55.",
                    type="string")
    parser.add_option("--labels", dest="labels",
                    help="Label numbers to keep in the output_file, separated by commas, i.e. --labels 3,11,21.",
                    type="string")
                         
    options, args = parser.parse_args()

    print '>>> %s [%s]: %s\n' % (time.ctime(time.time()), os.getcwd(), string.join(argv))
    sys.stdout.flush()

    if len(args) == 3:
        input_file, mr_atlas_file, output_file = args
    elif len(args) ==2:
        input_file, mr_atlas_file = args
        if options.veins:
            output_file = (input_file[:-3]+ "_venousdrainagemap.mnc")
        else:
            output_file = (input_file[:-3]+ "_arteriolperfusionmap.mnc")
        output_file = output_file.replace("_add_coverage","")
        output_file = output_file.replace("_delleave","")
        output_file = output_file.replace("_radi1.06","")
        output_file = output_file.replace("_smooth2","")
        output_file = output_file.replace("_simplified","")
        output_file = output_file.replace("_rmisolate","")
        output_file = output_file.replace("_features","")
        output_file = output_file.replace("_manuallabel","")
        output_file = output_file.replace("_groundtruthlabel","")
        output_file = output_file.replace("_cyl","")
    else:  
        print len(args), ": ", args
        parser.error("incorrect number of arguments")

    if options.remap:
        output_file = output_file[:-4]+"_remap.mnc"
    if not options.mask:
       output_file = output_file[:-4]+"_notmasked.mnc" 
    print   output_file  
    if not options.clobber and os.path.exists(output_file):
        raise SystemExit, \
            "The --clobber option is needed to overwrite an existing file."

    g= graph_analysis.input_graph(input_file)
    #g= graph_analysis.input_graph("/micehome/sghanavati/Documents/data/journal1data_c57_feb13/LOO_final_revised/autolabel_c57_1_4_S_cyl.db")
    #g= graph_analysis.input_graph("./c57_1_9Nov12/c57_1_graph_add_coverage_delleave_radi1.06_smooth2_simplified_rmisolate_features_groundtruthlabel_cyl_revised.db")

    mr_atlas_file = "/projects/souris/sghanavati/data/c57_male_female_mri_atlas/c57_brain_segmentation_mnc2_merged_reg2vascular.mnc"
    ##mr_atlas_file = "/projects/souris/sghanavati/data/c57_male_female_mri_atlas/c57_brain_segmentation_mnc2_merged_reg2vascular_downsample_merge.mnc"

    ##labels = [2,43,102,143,122,8,5,108,105,35,135,13,236,111,191,190,91,90,200,9,7,10,196,96,198,68,227,168,169,46,12,49,45,14,15,17,18,3,4]
    #arteries :    8:108,5:105,135:35,111:35,236:35,191:91,190:90,68:168,227:169 
    if not options.labels:
        labels = [122, 8, 5, 108, 105, 35, 135, 13, 236, 111, 191, 190, 91, 90, 68, 227, 168, 169, 46, 12, 49, 45, 14, 15, 17, 18, 3, 4]
        if options.veins:
            labels = [11,111,30,246,230,231, 192,34, 20,21,120,121,22, 101,24, 57,58,157,158, 56,54, 1,16, 170, 171,172,250,249] #6,206, 
    else:
        labels = []
        labels = options.labels.split(',')
        labels = [int(l) for l in labels]
        
    if not options.mask:
        vessel_analysis.perfusion_region_nomask(g, mr_atlas_file, labels,options.remap, output_file)
    else:
        vessel_analysis.perfusion_region(g, mr_atlas_file, labels, options.mask, options.remap, output_file)
    
    
        
    
