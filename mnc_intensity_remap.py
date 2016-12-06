#!/usr/bin/env python
from vessel_tracking import graph_analysis
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


program_name = 'mnc_intensity_remap.py'


if __name__ == '__main__':
        
    usage = "Usage: "+program_name+" [options] mnc_file [output_mnc]\n"+\
            "   or  "+program_name+" --help";
            
    parser = OptionParser(usage)

    parser.add_option("--remap", dest="remap",
                    help="Label numbers to replace with new label numbers. For example, to turn label 3 into label 11, it would be --remap 3:11. More than one can be specified at once, separated by commas, i.e. --remap 5:105,8:108,135:35,111:35,236:35,191:91,190:90,68:168,227:169.",
                    type="string")
    parser.add_option("--clobber", action="store_true", dest="clobber",
                        default=0, help="overwrite output file")
    
    options, args = parser.parse_args()

    print '>>> %s [%s]: %s\n' % (time.ctime(time.time()), os.getcwd(), string.join(argv))
    sys.stdout.flush()

    if len(args) == 2:
        input_file, output_file = args
    elif len(args) ==1:
        input_file = args[0]
        output_file = input_file[:-4]+ "_remap.mnc"
    else:   
        parser.error("incorrect number of arguments")


    if not options.clobber and os.path.exists(output_file):
        raise SystemExit, \
            "The --clobber option is needed to overwrite an existing file."

    if options.remap:
        lpairs = options.remap.split(',')
        labelpairs = []
        tomap_l = {}
        for pair in lpairs:
            plabels = pair.split(':')
            plabels = [int(l) for l in plabels]
            labelpairs.append(plabels)
        for l in labelpairs:    
            tomap_l[int(l[0])] = int(l[1]) 
        print tomap_l

        
    mincfile = pymincf.volumeFromFile(input_file, dtype='ubyte') #ushort
    outfile = pymincf.volumeFromInstance(mincfile, output_file, dtype="ubyte")
    mincfile.openFile()
    size = mincfile.getSizes()
    #spacing = mincfile.getSeparations()
    for i in range(size[0]):
        for j in range(size[1]):
            for k in range(size[2]):
                if mincfile[ i,  j,  k] in tomap_l.keys():
                    outfile.data[ i,  j,  k] = tomap_l[mincfile[ i,  j,  k]]
                else:
                    outfile.data[ i,  j,  k] = mincfile[ i,  j,  k]
                    
    #volume = counter*spacing[0]*spacing[1]*spacing[2]
    mincfile.closeVolume()
    outfile.writeFile()
    outfile.closeVolume()
        
    '''    
    #mincfile = py_minc.ArrayVolume(input_file, nc_data_type=py_minc.NC_BYTE)# without this always [z,y,x], with ,dim_names=py_minc.FILE_ORDER_DIMENSION_NAMES order as the image    
    #sizes = mincfile.get_sizes()    
         
    #for z_vox in range(sizes[0]):
        #for y_vox in range(sizes[1]):
            #for x_vox in range(sizes[2]):
                #if round(mincfile.array[z_vox,y_vox,x_vox]) == 8:
                    #mincfile.array[z_vox,y_vox,x_vox] = 108
                #if round(mincfile.array[z_vox,y_vox,x_vox]) == 5:
                    #mincfile.array[z_vox,y_vox,x_vox] = 105
                #if round(mincfile.array[z_vox,y_vox,x_vox] )== 135:
                    #mincfile.array[z_vox,y_vox,x_vox] = 35
                #if round(mincfile.array[z_vox,y_vox,x_vox] )== 111:
                    #mincfile.array[z_vox,y_vox,x_vox] = 35
                #if round(mincfile.array[z_vox,y_vox,x_vox] )== 236:
                    #mincfile.array[z_vox,y_vox,x_vox] = 35
                #if round(mincfile.array[z_vox,y_vox,x_vox]) == 191:
                    #mincfile.array[z_vox,y_vox,x_vox] = 91
                #if round(mincfile.array[z_vox,y_vox,x_vox]) == 190:
                    #mincfile.array[z_vox,y_vox,x_vox] = 90
                #if round(mincfile.array[z_vox,y_vox,x_vox] )== 68:
                    #mincfile.array[z_vox,y_vox,x_vox] = 168
                #if round(mincfile.array[z_vox,y_vox,x_vox]) == 227:
                    #mincfile.array[z_vox,y_vox,x_vox] = 169
                      
    #mincfile.output(output_file)
    '''
    
    
    
