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
program_name = 'perfusionterritoryvariation.py'


if __name__ == '__main__':
        
    usage = "Usage: "+program_name+" [options] input1_mnc input2_mnc input3_mnc ...  inputN_mnc  average.mnc\n"+\
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
    averagefile = args[-1]
    
    averageminc = pymincf.volumeFromFile(averagefile, dtype='ubyte')
    averageminc.openFile()
    
    for perfusion_file in perfusion_files:
        print perfusion_file
        
        perfusion_pattern = pymincf.volumeFromFile(perfusion_file, dtype='ubyte')  #ushort
        perfusion_pattern.openFile()
        sizes = perfusion_pattern.getSizes()   
        spacing = perfusion_pattern.getSeparations()
        outminc = pymincf.volumeFromInstance(perfusion_pattern,perfusion_file[:-4]+"_variation.mnc" , dtype="ubyte")
        #print sizes[0], sizes[1], sizes[2],":",spacing[0],spacing[1],spacing[2]
        for z_vox in range(sizes[0]):
            for y_vox in range(sizes[1]):
                for x_vox in range(sizes[2]):
                    value = perfusion_pattern[z_vox,y_vox,x_vox]
                    #value = int(round(perfusion_pattern.array[z_vox,y_vox,x_vox]))
                    if int(round(value))==int(round(averageminc[z_vox,y_vox,x_vox])):
                        outminc.data[z_vox,y_vox,x_vox] = 0
                    else:
                        outminc.data[z_vox,y_vox,x_vox] = 1
        outminc.writeFile()
        outminc.closeVolume()
        perfusion_pattern.closeVolume()
    averageminc.closeVolume()