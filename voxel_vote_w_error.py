#!/usr/bin/python

import pyminc.volumes.factory as pymincf
import pyminc.volumes.volumes as pymincv
from pyminc.volumes.factory import *
from numpy import *
from scipy.stats import *
from optparse import OptionParser
from sys import argv
import os, shelve, string ,time,sys
import commands

program_name = 'voxel_w_error_v3.py'


if __name__ == "__main__":

    usage = "Usage text"
    description = "Description text"
    
    parser = OptionParser(usage=usage, description=description)
    parser.add_option("--clobber", dest="clobber",
                    help="clobber output file",
                    type="string")
    parser.add_option("--reference", type="string", dest="reference_name",
            help="give the name for the reference minc to be used for output")

    (options, args) = parser.parse_args()
    print '\n\n>>> %s [%s]: %s\n' % (time.ctime(time.time()), os.getcwd(), string.join(argv))
    print "\n\n",options, "\n\n", args , "\n\n"

    if len(args) < 3:
        parser.error("Incorrect number of arguments")

    outfilename = args[-1]
    
    referencevol = pymincf.volumeFromFile(options.reference_name, dtype='ubyte')  #ushort
    #referencevol.openFile()
    # clobber check should go here
    
    volhandles = []

    # dtype='ubyte' or 'float' is writing a byte (floating point) but we want the voxel intensities to be int labels
    # dtype='short' or 'ushort' is writing a short integer data and we want the voxel intensities to be int labels
    nfiles = len(args)-1
    for i in range( nfiles ):
        volhandles.append(volumeFromFile(args[i], dtype='ubyte'))  #

    outfile = volumeFromInstance(referencevol, outfilename, dtype="ubyte")
    erroutfile = volumeFromInstance(referencevol, outfilename[:-4]+"_Error.mnc", dtype="ubyte")

    '''
    sliceArray = zeros( (nfiles,
                        volhandles[0].sizes[1],
                        volhandles[0].sizes[2]))
                        
    for i in range(volhandles[0].sizes[0]):
        for j in range(nfiles):
            t = volhandles[j].getHyperslab((i,0,0),
                                        (1,volhandles[0].sizes[1],
                                            volhandles[0].sizes[2]))
            t.shape = (volhandles[0].sizes[1], volhandles[0].sizes[2])
            sliceArray[j::] = t
        
        outfile.data[i::] = mode(sliceArray)[0]
        erroutfile.data[i::] = 100.0*mode(sliceArray)[1]/float(nfiles)
    '''
    sizes = volhandles[0].getSizes()   
    for z_vox in range(sizes[0]):
        for y_vox in range(sizes[1]):
            for x_vox in range(sizes[2]):
                #perfuselist = [int(round(vol[z_vox,y_vox,x_vox])) for vol in volhandles]
                perfuselist = [vol[z_vox,y_vox,x_vox] for vol in volhandles]
                outfile.data[z_vox,y_vox,x_vox] = mode(perfuselist)[0]
                erroutfile.data[z_vox,y_vox,x_vox] = 100.0*mode(perfuselist)[1]/float(nfiles)

    
    outfile.writeFile()
    outfile.closeVolume()
    erroutfile.writeFile()
    erroutfile.closeVolume()
    #referencevol.closeVolume()

                                                        
