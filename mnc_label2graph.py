#!/usr/bin/env python

# Get corrected_graph.db (before simplified) and labeled_CT.mnc
# for each vertex it finds the corresponding label from the mnc file (with labels initialzied to 0 and final labels are [0...255])
# it will then simplify the graph
# then it will run edge_labeling which includes label_propage to replace label=0 with a correct one
# and add the label as an edge_property and save the graph
# it will save the labeled simplified graph and graph-to-graph h5 for Brainview
#  
#  Created Feb 6, 2012
#  Last modified Apr 25, 2012
#  Sahar Ghanavati




from vessel_tracking import graph_analysis, distribution_analysis
#from graph_labeling import *
from optparse import OptionParser, Option, OptionValueError
from sys import argv
import os, commands, time, string
import py_minc 
from numpy import *
import scipy
import copy
import operator

#---------------------------------------------------------------------------------
#
#def mean_list(numberList):
    ##floatNums = [float(x) for x in numberList]
    #return sum(numberList) / float(len(numberList))	
    
def create_dist_voxels(N):
    #if x,y,z goes from [-N..N], create a list of all possible coordinates and sort them by distance:
    indxs =[]
    for i in range(-N,N+1):
        for j in range(-N,N+1):
            for k in range(-N,N+1):
                dist=calc_dist(i,j,k)
                indxs.append([i,j,k,dist])
    indx_l = sorted(indxs,key=lambda indxs:indxs[3])
    return indx_l
    
def calc_dist(x,y,z):
    d=sqrt(x*x+y*y+z*z)
    return d


def edge_labeling(h,labeled_ct):
    sizes=labeled_ct.get_sizes()
    neighb=1
    dist_indx=create_dist_voxels(neighb)
    ###for each edge, we find the class intermediaries belong to in mnc file, and edge label will be the one repeated the most
    for e in h.edge_list():
        vertex_list = [e[0]]+ g.edge_property(e, 'intermediaries')+ [e[1]]
        tmp_label={0:0}	#{label:number of intermediaries with this label}
        for v in vertex_list:
            v_world= g.vertices[v].centre	##get world in [x,y,z] 
            v_voxel= labeled_ct.convert_world_to_voxel(v_world) #convert it to image voxels according to dimension order of image which is read to be [z,y,x].			
            found=0
            l=0
            for i in dist_indx:		#search in radius 1 of the vertex for labels
                if found==0:
                    if (((v_voxel[0]+i[0])>0) and ((v_voxel[0]+i[0])<sizes[0]) and ((v_voxel[1]+i[1])>0) and ((v_voxel[1]+i[1])<sizes[1]) and ((v_voxel[2]+i[2])>0) and ((v_voxel[2]+i[2])<sizes[2]) ):
                        if int(labeled_ct.array[v_voxel[0]+i[0],v_voxel[1]+i[1],v_voxel[2]+i[2]])!=0:
                            l = int(labeled_ct.array[v_voxel[0]+i[0],v_voxel[1]+i[1],v_voxel[2]+i[2]])   #convert to list =list(intermed_voxel) or map(None,intermed_voxel)
                            found=1

            if not tmp_label.has_key(l):
                tmp_label[l]=1
            if tmp_label.has_key(l):
                tmp_label[l]+=1
        if len(tmp_label.keys())>1: 
            del tmp_label[0] #choose the non-zero label that most intermediaries has and assign it to the edge label
            label=max(tmp_label.iteritems(), key=operator.itemgetter(1))[0]

        elif len(tmp_label.keys())==1:
            label=tmp_label.keys()[0]
        else:
            print ("Error in reading CT labels for edge (%d,%d)!\nAborted!" %(e[0], e[1]) )
            exit(0)
        #else: 	#if edge has no intermediaries
            ##print ("edge (%d,%d) has no intermediaries " %(edge[0], edge[1]))
            #if h.vertices[e[0]].radius > h.vertices[e[1]].radius:		#assign the label of mother branch to the segment => which should come from the larger endpoint of segment edge   
                #label=h.vertices[e[0]].label
            #else:
                #label=h.vertices[e[1]].label
                
        h.set_edge_property(e,"label", label)	
        
    return h
    

def label_propagate(h,label):
    #gets the graph h and label that needs to be replaced by a label from higher archial vessel
    num0=1	#number of edges with label, initialize num0 to 1 so while will be executed the first time
    counter=0
    while (num0!=0 and counter<20):	#len(h.edge_list()) ):
        counter=counter+1
        num0=0		
        for edge in h.edge_list():
            if h.edge_property(edge,"label") == label:
                num0=num0+1
                #create list of all connecting edges to this edge
                list_edge=[[(min(edge[0],e),max(edge[0],e)),h.edge_property((min(edge[0],e),max(edge[0],e)),'diameter')] for e in h.vertices[edge[0]].edges] \
                + [ [(min(edge[1],e),max(edge[1],e)),h.edge_property((min(edge[1],e),max(edge[1],e)),'diameter')] for e in h.vertices[edge[1]].edges]
                
                #list_e = sorted(list_edge,key=lambda list_edge:list_edge[1])	#sort ascending by edge diameters
                list_e = sorted(list_edge,key=lambda list_edge:list_edge[1], reverse=True)	#sort descending by edge diameters
                
                for e in list_e:
                    l=h.edge_property(e[0],'label')
                    if l!=label:
                        h.set_edge_property(edge,'label',l)		#will replace edge_label==0 with label of largest diameter higher hierarchial vessel
        print ("number of segments with label=%d: %d" %(label, num0))
        
    if counter==len(h.edge_list()):	
        print ("Error: couldn't replace some %d labels!\nAborted!" %label )
        #exit(0)
        
    return h	




program_name = 'mnc_label2graph.py'


if __name__ == '__main__':
        
    usage = "Usage: "+program_name+" [options] postprocess_vessel_tracking(simplified).db LabeledCT.mnc [output_labeled.db] \n"+\
            "   or  "+program_name+" --help";
            
    parser = OptionParser(usage)
    
    parser.add_option("--clobber", action="store_true", dest="clobber",
                    default=0, help="overwrite output file")
    parser.add_option("--remove_interms",action="store_true",default=0, dest="remove_interms",help="remove the intermediate mnc files that were created for calculating seeds in iterations")
        
        
    tic = time.time()
    #ct_labeling = ct_labeling()		#initialize the class ct_labeling, make an instance called ct_labeling
                
    options, args = parser.parse_args()
    print '\n\n>>> %s [%s]: %s\n' % (time.ctime(time.time()), os.getcwd(), string.join(argv))
    print "\n\n",options, "\n\n", args , "\n\n"
   

        
    if len(args) == 2:
        input_file, labeled_CT_file = args
        output_file = input_file[:-3]+"_mnclabel.db"
    elif len(args) == 3:
        input_file, labeled_CT_file, output_file = args
    else:
        parser.error("incorrect number of arguments")
        
    
    if not options.clobber and os.path.exists(output_file):
        raise SystemExit, \
            "The --clobber option is needed to overwrite an existing file."
            
        
    # reads in the graph to be labeled
    try:
        g, attributes = graph_analysis.input_graph(input_file, ["history", "vertex_offsets"])
        print ("Succefully read in the %s\n" %input_file)	
    except:
        print("Error reading in the %s\n" %input_file)
        
    

    # read in the labeled CT : used to make labling training dataset
    #labeled_ct=py_minc.ArrayVolume(labeled_CT_file,dim_names=py_minc.FILE_ORDER_DIMENSION_NAMES) #which is in [x,y,z]
    labeled_ct=py_minc.ArrayVolume(labeled_CT_file)	# without this always [z,y,x], with this :dim_names=py_minc.FILE_ORDER_DIMENSION_NAMES order as the image	
    #
    
    
    g = edge_labeling(g,labeled_ct)
    #history = history + '\n>>> %s: %s' % (time.ctime(time.time()), "edge_labeling")		##argv = ['code.py', 'input.db', 'output']
    history = '\n>>> %s: %s' % (time.ctime(time.time()), string.join(argv))		##argv = ['code.py', 'input.db', 'output']
    if attributes.has_key("history"):
        history = attributes["history"] + "\n" + history
        del attributes['history']
    graph_analysis.output_graph(output_file, g, history, attributes)   
    print ("Succefully wrote the %s\n" %output_file)
    
    g = label_propagate(g,0)
    #history = history + '\n>>> %s: %s' % (time.ctime(time.time()), "edge_labeling")		##argv = ['code.py', 'input.db', 'output']
    history += '\n>>> %s: %s' % (time.ctime(time.time()), "label_propagate")		##argv = ['code.py', 'input.db', 'output']
    if attributes.has_key("history"):
        #history = attributes["history"] + "\n" + history
        del attributes['history']
    graph_analysis.output_graph(output_file[:-3]+"_propagate.db", g, history, attributes)   
    print ("Succefully wrote the %s\n" %output_file[:-3]+"_propagate.db")


#####
    ###if not os.path.exists(input_file[:-3]+"_dt.mnc"):
    ##ct_labeling.graph2dt_exec (input_file,labeled_CT_file)
    ##graph_dt_minc=py_minc.ArrayVolume(input_file[:-3]+"_dt.mnc")
    ##sizes_dt=graph_dt_minc.get_sizes()

    ####for each vertex in the graph find its corresponding CTlabel and save it as an attribute in g.vertices[].CTlabel
    #print "begin calculating the vertices labels "
    #h = ct_labeling.vertex_labeling(g,labeled_ct)
    #print "successfully calculated vertices labels"
    
    #num=0
    #for v in range(len(h.vertices)):
        #if h.vertices[v].label==0:
            #num=num+1
            ##print ("vertex %d has label 0" %v)
    ##print("%d vertices out of %d have label 0" %(num,len(h.vertices)))

    ####for each edge in the simplified graph find its corresponding CTlabel and save it as an edge_property label
    #print "begin calculating the edge labels "
    #h = ct_labeling.edge_labeling(h)
    #print "successfully calculated edge labels"
#####
#####Label propagation: for those with label 0 we should find label by upper branches labels
#h=self.label_propagate(h,0)


##check if any edge has label 0 by mistake:
#for edge in h.edge_list():
    #if h.edge_property(edge,"label") <1:	#label 0 and -1
        #print("edge (%d,%d) has label %d" %(edge[0], edge[1],h.edge_property(edge,"label")) )
    
    
#for edge in h.edge_list():
##correct the label of all intermediaries to the label of edge
    #for i in h.edge_property(edge,'intermediaries'):
        #if h.vertices[i].label!=h.edge_property(edge,'label'):
            #h.vertices[i].label=h.edge_property(edge,'label')
#####


            

    
    ## save the corrected final h5 file for the Brainview
    #cmd=("graph2graph.py %s %s --clobber" %(output_file,output_file[:-2]+"h5"))	#python /micehome/jgsled/bin/
    #print(cmd)
    #os.system(cmd)
    
    if options.remove_interms:
        cmnd = "rm "+output_file
        print  "\n", cmnd, "\n"
        os.system(cmnd)
        
    
    toc = time.time()
    print ("total elapsed time: %f" %(toc - tic))