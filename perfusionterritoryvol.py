#!/usr/bin/env python
#from vessel_tracking import graph_analysis
import vessel_analysis
from optparse import OptionParser, Option, OptionValueError
from sys import argv
import os, shelve, string ,time,sys
import commands
#import py_minc
import pyminc.volumes.factory as pymincf
import pyminc.volumes.volumes as pymincv
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

label2name = {0:"None", 2:"Internal Carotid A. R.", 43:"Internal Carotid A. L.", 102:"Internal Carotid A. R. level1", \
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

program_name = 'perfusionterritoryvol.py'


if __name__ == '__main__':
        
    usage = "Usage: "+program_name+" [options] input1_mnc input2_mnc input3_mnc ...  inputN_mnc  output.csv\n"+\
            "   or  "+program_name+" --help";
            
    parser = OptionParser(usage)

    #parser.add_option("--clobber", action="store_true", dest="clobber",
                        #default=0, help="overwrite output file")
    
    #parser.add_option("--mask", type="string", dest="mask",
                         #help="the mask minc file")

    #parser.add_option("--remap", dest="remap",
                    #help="Label numbers to replace with new label numbers. For example, to turn label 3 into label 11, it would be --remap 3:11. More than one can be specified at once, separated by commas, i.e. --remap 3:11,21:55.",
                    #type="string")
                         
    options, args = parser.parse_args()

    print '>>> %s [%s]: %s\n' % (time.ctime(time.time()), os.getcwd(), string.join(argv))
    sys.stdout.flush()
    
    perfusion_files = args[0:-1]
    outputfile = args[-1]
    
    f = open(outputfile,'w')
    f.write( "sample, strain, label, label_Name ,number_of_voxels, perfusion_volume[mm3], brain_volume_voxels, brain_volume_mm3, perfusion_percent")
    print( "sample, strain, label, label_Name ,number_of_voxels, perfusion_volume[mm3], brain_volume_voxels, brain_volume_mm3, perfusion_percent")
    
    for perfusion_file in perfusion_files:
        print perfusion_file
        
        #if len(args) == 3:
            #input_file, mr_atlas_file, output_file = args
        #elif len(args) ==2:
            #input_file, mr_atlas_file = args
            #output_file = (input_file[:-3]+ "_arteriolperfusionmap.mnc")
            #output_file = output_file.replace("_add_coverage","")
            #output_file = output_file.replace("_delleave","")
            #output_file = output_file.replace("_radi1.06","")
            #output_file = output_file.replace("_smooth2","")
            #output_file = output_file.replace("_simplified","")
            #output_file = output_file.replace("_rmisolate","")
            #output_file = output_file.replace("_features","")
            #output_file = output_file.replace("_manuallabel","")
            #output_file = output_file.replace("_groundtruthlabel","")
            #output_file = output_file.replace("_cyl","")
        #else:   
            #parser.error("incorrect number of arguments")

        #if options.remap:
            #output_file = output_file[:-4]+"_remap.mnc"
            
        #if not options.clobber and os.path.exists(output_file):
            #raise SystemExit, \
                #"The --clobber option is needed to overwrite an existing file."

        #g= graph_analysis.input_graph(input_file)
        ##g= graph_analysis.input_graph("/micehome/sghanavati/Documents/data/journal1data_c57_feb13/LOO_final_revised/autolabel_c57_1_4_S_cyl.db")
        ##g= graph_analysis.input_graph("./c57_1_9Nov12/c57_1_graph_add_coverage_delleave_radi1.06_smooth2_simplified_rmisolate_features_groundtruthlabel_cyl_revised.db")

        #mr_atlas_file = "/projects/souris/sghanavati/data/c57_male_female_mri_atlas/c57_brain_segmentation_mnc2_merged_reg2vascular.mnc"
        ###mr_atlas_file = "/projects/souris/sghanavati/data/c57_male_female_mri_atlas/c57_brain_segmentation_mnc2_merged_reg2vascular_downsample_merge.mnc"

        ###labels = [2,43,102,143,122,8,5,108,105,35,135,13,236,111,191,190,91,90,200,9,7,10,196,96,198,68,227,168,169,46,12,49,45,14,15,17,18,3,4]
        #labels = [122, 8, 5, 108, 105, 35, 135, 13, 236, 111, 191, 190, 91, 90, 68, 227, 168, 169, 46, 12, 49, 45, 14, 15, 17, 18, 3, 4]

        labels = {}
        #perfusion_pattern = py_minc.ArrayVolume(perfusion_file, nc_data_type=py_minc.NC_BYTE)# without this always [z,y,x], with ,dim_names=py_minc.FILE_ORDER_DIMENSION_NAMES order as the image   
        #sizes = perfusion_pattern.get_sizes()
        #spacing = perfusion_pattern.get_separations()
        perfusion_pattern = pymincf.volumeFromFile(perfusion_file, dtype='ubyte')  #ushort
        perfusion_pattern.openFile()
        sizes = perfusion_pattern.getSizes()   
        spacing = perfusion_pattern.getSeparations()
        #print sizes[0], sizes[1], sizes[2],":",spacing[0],spacing[1],spacing[2]
        for z_vox in range(sizes[0]):
            for y_vox in range(sizes[1]):
                for x_vox in range(sizes[2]):
                    value = perfusion_pattern[z_vox,y_vox,x_vox]
                    #value = int(round(perfusion_pattern.array[z_vox,y_vox,x_vox]))
                    if value>0:
                        if value not in labels.keys():
                            labels[value] = 1
                        else:
                            labels[value] += 1
        print  labels.keys()          
        for l in labels.keys():
            print l, labels[l]
            print( "\n",os.path.basename(perfusion_file)[:-4], " , ",os.path.basename(perfusion_file)[:3], " , ",str(l), " , ", label2name[l] , " , " , str(labels[l]),",", str(labels[l]*spacing[0]*spacing[1]*spacing[2]) )
            f.write( "\n"+os.path.basename(perfusion_file)[:-4]+ " , "+os.path.basename(perfusion_file)[:3]+ " , "+str(l)+ " , "+ label2name[l] + " , " + str(labels[l])+","+ str(labels[l]*spacing[0]*spacing[1]*spacing[2]) )
        
        
        #print  os.path.basename(perfusion_file)[:-4], " , ",os.path.basename(perfusion_file)[:3],", total_voxels , ", sum(labels.values()),", total_volume , ", sum(labels.values())*spacing[0]*spacing[1]*spacing[2]   
        f.write( "\n"+os.path.basename(perfusion_file)[:-4]+ " , "+os.path.basename(perfusion_file)[:3]+ " , total , total , , ," +str(sum(labels.values()))+","+ str(sum(labels.values())*spacing[0]*spacing[1]*spacing[2]) )
    f.close()
    
    
    
    