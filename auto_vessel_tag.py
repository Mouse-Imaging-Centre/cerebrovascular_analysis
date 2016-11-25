#!/usr/bin/env python

from optparse import OptionParser, Option, OptionValueError
import py_minc 
import pyminc.volumes.factory as pymincf
import pyminc.volumes.volumes as pymincv
import numpy as np
from numpy import *
import operator
import os        
from morphology import voxel_code
from vessel_tracking import medial_atoms
import time, string
from sys import argv

################################
object_dtype = np.int16   # integer type for binary and distance transform images

'''
def set_periphery(data, value):
    "an 'in place' operation to set voxels at edge of data volume to zero"

    data[0,:,:] = value
    data[-1,:,:] = value
    data[:,0,:] = value
    data[:,-1,:] = value
    data[:,:,0] = value
    data[:,:,-1] = value
    
def flat_array_strides(array):
    "compute the number of steps for a flat array corresponding to a step in each index"
    sizes = shape(array)

    strides = [1]
    for i in range(1,len(sizes)):
        strides.insert(0,multiply.reduce(sizes[-i:]))

    return asarray(strides)

    
def fev_neighborhood_offsets(dim, connectivity=26):
    "compute linear offset corresponding to F, E, and V neighbors in an array of shape dim"
    if connectivity not in [6, 18, 26]:
        raise ValueError, "Connectivity must be 6, 18, or 26"

    neighborhoods = [];
    steps = [1, dim[2], dim[1]*dim[2]]
    # face neighbors
    neighborhoods.append(list(ravel(map(None, negative(steps), steps))))

    if connectivity >= 18:
        # edge neighbors
        neighborhoods.append( list(ravel(map(
            lambda x, y: [-x-y, -x+y, x-y, x+y],
            steps, [steps[1], steps[2], steps[0]]))))

    if connectivity == 26:
        # vertice neighbors
        v = []
        for i in [-steps[0], steps[0]]:
            for j in [-steps[1], steps[1]]:
                for k in [-steps[2], steps[2]]:
                    v.append(i+j+k)
        neighborhoods.append(v)

    return neighborhoods


def voxel_code(data, seeds, metrics, final=[], prep_flag=0):
    """Use voxel coding to return a labelling for the given data starting from the
given seeds and using the given distance metric.  Optionally stop before
completion if the final point(s) is specified.  Note that this implementation
is a mix of python and C++.

label = voxel_code(data, seeds, metrics, final=[], prep_flag=0)

Notes:

If prep_flag is true then the data array is assumed to have come from the function
voxel_code_prep.
"""

    max_value = 2**15-1

    # create array for voxel codes
    if prep_flag:
        label = copy(data)
    else:
        label = where(greater(data, 0), asarray(max_value,object_dtype),
                      asarray(0,object_dtype))

    # set values at edge to be background
    set_periphery(label, asarray(0, object_dtype))

    # assign neighborhoods to metrics
    neighborhoods = fev_neighborhood_offsets(shape(label))
    neighbors = []
    distances = []
    for i in range(len(neighborhoods)):
        for n in neighborhoods[i]:
            neighbors.append(n)
            distances.append(metrics[i])
            
    # set values of initial seeds
    s = [ x[0] for x in seeds ];
    v = [ x[1] for x in seeds ];

    label, index =  voxel_code_opt.voxel_code_inner_loops(label.flat,
              asarray(neighbors, Int),
              asarray(distances, object_dtype), asarray(s, Int),
              asarray(v, object_dtype), asarray(final, Int))

    return reshape(label, shape(data)), index

    
def find_boundary_seeds(data, metrics):
    "find boundary seed points used for distance transform"
    
    if metrics[2] < metrics[1] or metrics[1] < metrics[0]:
        raise ValueError, \
          "metrics argument must be a tuple of three integers in non-decreasing order"

    # assign neighborhoods to metrics
    neighborhoods = map(None, metrics, fev_neighborhood_offsets(shape(data)))

    # create a list of boundary points to use as seeds
    seeds = []

    # iterate over all points except volume boundaries
    sizes = shape(data)

    strides = flat_array_strides(data)
    for i in range(strides[0], (sizes[0]-1)*strides[0], strides[0]):
        for j in range(strides[1], (sizes[1]-1)*strides[1], strides[1]):
            for k in range(strides[2], (sizes[2]-1)*strides[2], strides[2]):
                index = i+j+k

                # only consider non-zero points
                if data.flat[index] > 0:
                    distance = 0
                    # for each of F, E, V neighborhoods
                    # (note: nearest neighbors are first)
                    for neighborhood in neighborhoods:
                        (metric, neighbors) = neighborhood
                        for neighbor in neighbors:
                            # if neighbor is a background voxel
                            if data.flat[index+neighbor] == 0:
                                distance = metric
                                break
                            # if a neighbor has already been found
                        if distance > 0:
                            # add this point as a boundary point with
                            #  given distance value
                            seeds.append((index, distance))
                            break

    return seeds

def distance_transform(data, metrics):
    "compute distance transform for given 3D binary array"

    # get boundary points
    seeds = find_boundary_seeds(data, metrics)

    # using boundary points as seeds, voxel code the entire volume
    label, index = voxel_code(data, seeds, metrics)

    return label
'''    
    
program_name = 'auto_vessel_tag.py'


if __name__ == '__main__':
        
    usage = "Usage: "+program_name+" [options] input_dt_image.mnc input_group_image.mnc seed_file.tag \n"+\
            "   or  "+program_name+" [options] input_image.mnc seed_file.tag \n"+\
            "   or  "+program_name+" --help";
            
    parser = OptionParser(usage)
    parser.add_option("--clobber", action="store_true", dest="clobber",
                        default=0, help="overwrite output file")
    parser.add_option("--bin_threshold",type="float", dest="bin_threshold", default= 0.2, help="intensity threshold for binarizing the image (default 0.15)") 
    parser.add_option("--additional_seeding", action="store_true", dest="additional_seeding",
                        default=0, help="Add seeds exhaustively")
    parser.add_option("--laplacian", action="store_true", dest="laplacian",
                        default=0, help="Use the laplacian of the distance_transform to force seeds to the centre of the vessels")
    parser.add_option("--weight", action="store_true", dest="weight",
                        default=0, help="Weighting option for the seeds for new follow_tree.py version. Otherwise, the weight is set to 1.")
    parser.add_option("--connectivity",type="string",dest="connectivity",default="3D06",help="connectivity for mincmorph")

    options,args = parser.parse_args() 
    start_time = time.time()
    print '\n\n>>> %s [%s]: %s\n' % (time.ctime(time.time()), os.getcwd(), string.join(argv))
    print "\n\n",options, "\n\n", args , "\n\n"

    if len(args)==2:
        mnc_file, tag_file = args 
        mnc_bin_file = mnc_file[:-4]+"_binarize.mnc"
        mnc_group_file = mnc_bin_file[:-4]+"_group.mnc"
        mnc_dt_file = mnc_bin_file[:-4]+"_dt.mnc"
        #create binarize file
        cmd = "mincmath" + " -clobber -quiet" + " -byte" + " -gt -const %f %s %s" %(options.bin_threshold,mnc_file,mnc_bin_file)
        os.system(cmd)
         #create dt file
        cmd = "mincmorph" + " -clobber" + " -short" + " -successive F %s %s" %(mnc_bin_file,mnc_dt_file)
        os.system(cmd)
        #create group file
        cmd = "mincmorph" + " -clobber" + " -short" + " -group -%s %s %s" %(options.connectivity,mnc_bin_file,mnc_group_file)
        os.system(cmd)
    elif len(args)==3:
        mnc_dt_file, mnc_group_file, tag_file = args 
    else:	
        parser.error("incorrect number of arguments")
        

    if not options.clobber and os.path.exists(tag_file):
        raise SystemExit, \
            "The --clobber option is needed to overwrite an existing file."

    threshold=options.bin_threshold
    if not options.laplacian:
        print "No Laplacian options"
    if options.additional_seeding:
        print "additional_seeding options"

    tags = py_minc.VolumeTags(1)
    #dt_img = py_minc.ArrayVolume(mnc_dt_file)# without this always [z,y,x], with ,dim_names=py_minc.FILE_ORDER_DIMENSION_NAMES order as the image	
    #group_img = py_minc.ArrayVolume(mnc_group_file)# without this always [z,y,x], with ,dim_names=py_minc.FILE_ORDER_DIMENSION_NAMES order as the image	
    dt_img = pymincf.volumeFromFile(mnc_dt_file)# without this always [z,y,x], with ,dim_names=py_minc.FILE_ORDER_DIMENSION_NAMES order as the image	
    group_img = pymincf.volumeFromFile(mnc_group_file)# without this always [z,y,x], with ,dim_names=py_minc.FILE_ORDER_DIMENSION_NAMES order as the image	
    dt_img.openFile()
    group_img.openFile()

    
   
    #sizes = dt_img.get_sizes()
    sizes = dt_img.getSizes()
    dt_scale_factor = min(dt_img.getSeparations())/3.0  # used to convert distance transform values world space units 
    #dt_scale_factor = min(dt_img.get_separations())/3.0  # used to convert distance transform values world space units

    usable_bounds = medial_atoms.bounding_box_from_shape(sizes)
    #usable_bounds.deflate(1)  # periphery is unusable
    
    # extract hyperslab to process
    dt_hyperslab = dt_img.getHyperslab(usable_bounds.start(), usable_bounds.shape(), dtype = dt_img.dtype)
    
    
    if options.weight:
        #object_map = np.greater_equal(dt_img, threshold).astype(object_dtype)
        object_map = greater_equal(dt_hyperslab, threshold).astype(object_dtype)
        if amax(object_map.flat) <=0:
            print "object_map is all 0. No distance_map nor seed can be created."
        else:
            distance_map = voxel_code.distance_transform(object_map, (3,4,5))
            #distance_map = distance_transform(object_map, (3,4,5))
            del object_map

    if options.laplacian:
        cmd = "make_gradient_volume %s %s 2" %(mnc_dt_file, mnc_dt_file[:-4]+"_grad2.mnc") 
        os.system(cmd)
        grad2_img = pymincf.volumeFromFile(mnc_dt_file[:-4]+"_grad2.mnc")# without this always [z,y,x], with ,dim_names=py_minc.FILE_ORDER_DIMENSION_NAMES order as the image  
        grad2_img.openFile()
        grad_hyperslab = grad2_img.getHyperslab(usable_bounds.start(), usable_bounds.shape(), dtype = grad2_img.dtype)      
        grad2_threshold = amax(grad_hyperslab.flat)/2.0                
        component_maximum_gard2 = {}
            
    component_maximum = {}
    for z_vox in range(sizes[0]):
        for y_vox in range(sizes[1]):
            for x_vox in range(sizes[2]):
                g_value = int(round(group_img[z_vox,y_vox,x_vox]))
                if g_value>0:
                    dt_value = int(round(dt_img[z_vox,y_vox,x_vox]))
                    if dt_value>0:
                        weight = -1
                        if options.weight:
                            weight = distance_map[z_vox,y_vox,x_vox] * dt_scale_factor   # compute distance from edge in world space units
                        if g_value not in component_maximum.keys():
                            component_maximum[g_value]= [dt_value,z_vox,y_vox,x_vox,weight]
                            if options.laplacian:
                                grad_value = int(round(grad2_img[z_vox,y_vox,x_vox]))
                                component_maximum_gard2[g_value]= [grad_value,z_vox,y_vox,x_vox,weight]
                        else:
                            if dt_value > component_maximum[g_value][0]:
                                component_maximum[g_value]= [dt_value,z_vox,y_vox,x_vox,weight]
                            if options.laplacian and grad_value > component_maximum_gard2[g_value][0]:
                                    component_maximum_gard2[g_value]= [grad_value,z_vox,y_vox,x_vox,weight]
                        if dt_value>0.5 and options.additional_seeding: 
                            if not options.laplacian:
                                dt_value_world_location= dt_img.convertVoxelToWorld(np.array([z_vox,y_vox,x_vox]))       #return world in [x,y,z]
                                tags.append(dt_value_world_location,weight, -1, -1, "")
                            else:
                                if grad2_img[z_vox,y_vox,x_vox]>=grad2_threshold:
                                    dt_value_world_location= dt_img.convertVoxelToWorld(np.array([z_vox,y_vox,x_vox]))       #return world in [x,y,z]
                                    tags.append(dt_value_world_location,weight, -1, -1, "")
 

        #for group in component_maximum.keys():
            #voxel_location = np.array([component_maximum[group][1],component_maximum[group][2],component_maximum[group][3]])
            #world_location = dt_img.convertVoxelToWorld(voxel_location)
            #tags.append(world_location,component_maximum[group][4], -1, -1, "") 
        
        #if options.laplacian:
            #for group in component_maximum_gard2.keys():
                #voxel_location = np.array([component_maximum_gard2[group][1],component_maximum_gard2[group][2],component_maximum_gard2[group][3]])
                #world_location = dt_img.convertVoxelToWorld(voxel_location)
                #tags.append(world_location,component_maximum_gard2[group][4], -1, -1, "") 

        #descending sort tag components based on value (dt_value)
        sorted_component_maximum = sorted(component_maximum.iteritems(), key=operator.itemgetter(1),reverse=True)  #result is not dict, list of tuples=(key,val)
        
        for i in range(len(sorted_component_maximum)):
            voxel_location = np.array([sorted_component_maximum[i][1][1],sorted_component_maximum[i][1][2],sorted_component_maximum[i][1][3]])
            #world_location = dt_img.convert_voxel_to_world(voxel_location)
            world_location = dt_img.convertVoxelToWorld(voxel_location)
            tags.append(world_location,sorted_component_maximum[i][1][4], -1, -1, "") 
               
                
        if options.laplacian:
            #descending sort tag components based on value (dt_value)
            sorted_component_maximum_gard2 = sorted(component_maximum_gard2.iteritems(), key=operator.itemgetter(1),reverse=True)  #result is not dict, list of tuples=(key,val)
            for i in range(len(sorted_component_maximum_gard2)):
                voxel_location = np.array([sorted_component_maximum_gard2[i][1][1],sorted_component_maximum_gard2[i][1][2],sorted_component_maximum_gard2[i][1][3]])
                #world_location = dt_img.convert_voxel_to_world(voxel_location)
                world_location = dt_img.convertVoxelToWorld(voxel_location)
                tags.append(world_location,sorted_component_maximum_gard2[i][1][4], -1, -1, "") 
                

    dt_img.closeVolume()
    group_img.closeVolume()
    if tags.number()>0:
        tags.output(tag_file)
        print "****************** seeding : ", tag_file, " DONE! ",tags.number()," seeds created!\n"
    else:
        print "****************** No seed can be added!\n"
        
    # report completion
    end_time = time.time()
    elapsed_time = int(round(end_time - start_time))
    print '>>> Processing complete: %s\nElapsed time: %u:%02u:%02u' % (time.ctime(end_time), \
                                  elapsed_time / 3600, (elapsed_time % 3600) / 60, elapsed_time % 60)
      