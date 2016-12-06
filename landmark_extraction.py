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

#landmarks at labels intersections:


def find_pair_label_vertices(g, labelpairs, unique_l ): 
    landmarks_g = []
    for lpair in labelpairs:
        found_v = []
        for e in g.edge_list():
            if g.edge_property(e, 'label') in unique_l:
                e1_neighbs = vessel_analysis.find_neighbor_edges_asymmetrical(g,e,0)
                e2_neighbs = vessel_analysis.find_neighbor_edges_asymmetrical(g,e,1)
                for e1n in e1_neighbs:
                    if set([g.edge_property(e, 'label'), g.edge_property(e1n, 'label')])==set(lpair):  
                        if list(g.vertices[e[0]].centre) not in found_v:
                            found_v.append(list(g.vertices[e[0]].centre))   
                            break
                for e2n in e2_neighbs:
                    if set([g.edge_property(e, 'label'), g.edge_property(e2n, 'label')])==set(lpair): 
                        if list(g.vertices[e[0]].centre) not in found_v:
                            found_v.append(list(g.vertices[e[1]].centre)) 
                            break
                if len(found_v)> 0:
                    break
        if len(found_v)> 0:
            #found_vertex
            landmarks_g.append(found_v[0]) 
        else:
            #None,none,none
            landmarks_g.append(['None','None','None']) 
    return landmarks_g   

   
    
    
def find_triple_label_vertices(g, labeltriples, unique_triple_l ): 
    landmarks_g = []
    for ltriple in labeltriples:        #ltriple is a set
        found_v = []
        for e in g.edge_list():
            if g.edge_property(e, 'label') in unique_triple_l:
                if len(g.vertices[e[0]].edges) >= 3:
                    e1_neighbs = vessel_analysis.find_neighbor_edges_asymmetrical(g,e,0)
                    for e1n in e1_neighbs:
                        for e2n in e1_neighbs:
                            if set([g.edge_property(e, 'label'), g.edge_property(e1n, 'label'), g.edge_property(e2n, 'label')])==ltriple:            
                                if list(g.vertices[e[0]].centre) not in found_v:
                                    found_v.append(list(g.vertices[e[0]].centre)) 
                                    break
                            if len(found_v)> 0:
                                break
                        if len(found_v)> 0:
                            break              
                elif len(g.vertices[e[1]].edges) >= 3:
                    e2_neighbs = vessel_analysis.find_neighbor_edges_asymmetrical(g,e,1)
                    for e1n in e2_neighbs:
                        for e2n in e2_neighbs:
                            if set([g.edge_property(e, 'label'), g.edge_property(e1n, 'label'), g.edge_property(e2n, 'label')])==ltriple:            
                                if list(g.vertices[e[0]].centre) not in found_v:
                                    found_v.append(list(g.vertices[e[0]].centre)) 
                                    break
                            if len(found_v)> 0:
                                break
                        if len(found_v)> 0:
                            break              
                if len(found_v)> 0:
                    break
        if len(found_v)> 0:
            #found_vertex
            landmarks_g.append(found_v[0]) 
        else:
            #None,none,none
            landmarks_g.append(['None','None','None']) 
    return landmarks_g   
    
    
    
program_name = 'landmark_extraction.py'


if __name__ == '__main__':
        
    usage = "Usage: "+program_name+" [options] input_graph1 [input_graph2] output.tag\n"+\
            "   or  "+program_name+" --help";
            
    parser = OptionParser(usage)

    parser.add_option("--clobber", action="store_true", dest="clobber",
                        default=0, help="overwrite output file")
    
    parser.add_option("--landmark", dest="landmark",
                    help="Label pair numbers where intersection landmarks should have. For example, to create a landmark where label 3 and label 11 intersect, it would be --landmark 3:11. More than one can be specified at once, separated by commas, i.e. --landmark 3:11,21:55.",
                    type="string")
                         
    parser.add_option("--landmark3", dest="landmark3",
                    help="Label triple numbers where intersection landmarks should have. For example, to create a landmark where label 3, label 11 and 15 intersect, it would be --landmark3 3:11:15. More than one can be specified at once, separated by commas, i.e. --landmarks3 3:11:12,21:55:246",
                    type="string")

    parser.add_option("--lsq6", dest="lsq6", help="6 parameter (scale=1.0) least-squares linear transformation.",
                    type="string")
    parser.add_option("--lsq7", dest="lsq7", help="7 parameter (one scale) least-squares linear transformation.",
                    type="string")
    parser.add_option("--lsq9", dest="lsq9", help="9 parameter least-squares linear transformation.",
                    type="string")
    parser.add_option("--lsq10", dest="lsq10", help="10 parameter least-squares linear transformation.",
                    type="string")
    parser.add_option("--lsq12", dest="lsq12", help="12 parameter least-squares linear transformation.",
                    type="string")
 
    parser.add_option("--inverse", dest="inverse", help="Swap tags, then compute transform (default=FALSE).",
                    type="string")

    options, args = parser.parse_args()

    print '>>> %s [%s]: %s\n' % (time.ctime(time.time()), os.getcwd(), string.join(argv))
    sys.stdout.flush()

    if len(args) == 2:
        input_file, output_file = args

    #elif len(args) ==1:
        #input_file = args[0]
        #output_file = (input_file[:-3]+ "_landmarks.tag")
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

    elif len(args) == 3:
        input_file1, input_file2, output_file = args   
    else:   
        parser.error("incorrect number of arguments")

    if not options.clobber and os.path.exists(output_file):
        raise SystemExit, \
            "The --clobber option is needed to overwrite an existing file."

            
    if options.landmark:
        lpairs = options.landmark.split(',')
        labelpairs = []
        for pair in lpairs:
            labels = pair.split(':')
            labels = [int(l) for l in labels]
            labelpairs.append(labels)
        all_l = [l[0] for l in labelpairs] + [l[1] for l in labelpairs]
        unique_l = []
        for l in all_l:
            if l not in unique_l:
                unique_l.append(l)
                
    if options.landmark3:
        ltriples = options.landmark3.split(',')
        labeltriples = []
        for pair in ltriples:
            labels = pair.split(':')
            labels = [int(l) for l in labels]
            labeltriples.append(set(labels))
        all_l = [list(l)[0] for l in labeltriples] + [list(l)[1] for l in labeltriples] + [list(l)[2] for l in labeltriples]
        unique_triple_l = []
        for l in all_l:
            if l not in unique_triple_l:
                unique_triple_l.append(l)

                
    if len(args) == 2:
        g = graph_analysis.input_graph(input_file)        
        landmarks = []
        according_l = []
            
        #for e in g.edge_list():
            #if len(vessel_analysis.find_neighbor_edges_asymmetrical(g,e,0)) not in [0,2]:
                #print e,":", vessel_analysis.find_neighbor_edges_asymmetrical(g,e,0)
            #if len(vessel_analysis.find_neighbor_edges_asymmetrical(g,e,1)) not in [0,2]:
                #print e,":", vessel_analysis.find_neighbor_edges_asymmetrical(g,e,1)
                
                
        if options.landmark:        
            for e in g.edge_list():
                if g.edge_property(e, 'label') in unique_l:
                    e1_neighbs = vessel_analysis.find_neighbor_edges_asymmetrical(g,e,0)
                    e2_neighbs = vessel_analysis.find_neighbor_edges_asymmetrical(g,e,1)
                    for e1n in e1_neighbs:
                        if [g.edge_property(e, 'label'), g.edge_property(e1n, 'label')] in labelpairs or [g.edge_property(e1n, 'label'), g.edge_property(e, 'label')] in labelpairs:            
                            if list(g.vertices[e[0]].centre) not in landmarks:
                                landmarks.append(list(g.vertices[e[0]].centre)) 
                                according_l.append([g.edge_property(e, 'label'), g.edge_property(e1n, 'label')])
                    for e2n in e2_neighbs:
                        if [g.edge_property(e, 'label'), g.edge_property(e2n, 'label')] in labelpairs or [g.edge_property(e2n, 'label'), g.edge_property(e, 'label')] in labelpairs:            
                            if list(g.vertices[e[1]].centre) not in landmarks:
                                landmarks.append(list(g.vertices[e[1]].centre))  
                                according_l.append([g.edge_property(e, 'label'), g.edge_property(e2n, 'label')])
            #for i in range(len(landmarks)):
                #print according_l[i],":",landmarks[i]
        
        if options.landmark3:
            for e in g.edge_list():
                if g.edge_property(e, 'label') in unique_triple_l:
                    if len(g.vertices[e[0]].edges) >= 3:
                        e1_neighbs = vessel_analysis.find_neighbor_edges_asymmetrical(g,e,0)
                        for e1n in e1_neighbs:
                            for e2n in e1_neighbs:
                                if set([g.edge_property(e, 'label'), g.edge_property(e1n, 'label'), g.edge_property(e2n, 'label')]) in labeltriples:            
                                    if list(g.vertices[e[0]].centre) not in landmarks:
                                        landmarks.append(list(g.vertices[e[0]].centre))
                                        according_l.append([g.edge_property(e, 'label'), g.edge_property(e1n, 'label'), g.edge_property(e2n, 'label')])
                    if len(g.vertices[e[1]].edges) >= 3:
                        e2_neighbs = vessel_analysis.find_neighbor_edges_asymmetrical(g,e,1)
                        for e1n in e2_neighbs:
                            for e2n in e2_neighbs:
                                if set([g.edge_property(e, 'label'), g.edge_property(e1n, 'label'), g.edge_property(e2n, 'label')]) in labeltriples:            
                                    if list(g.vertices[e[1]].centre) not in landmarks:
                                        landmarks.append(list(g.vertices[e[1]].centre))                                          
                                        according_l.append([g.edge_property(e, 'label'), g.edge_property(e1n, 'label'), g.edge_property(e2n, 'label')])
        for i in range(len(landmarks)):
            print according_l[i],":",landmarks[i]
            
       
    if len(args) ==3:
        g = graph_analysis.input_graph(input_file1)        
        h = graph_analysis.input_graph(input_file2)        
        landmarks_g = []
        landmarks_h = []
        
        #### write in .tags
        f = open(output_file, 'w')
        f.write( "MNI Tag Point File\nVolumes = 2;\n")
        special_letter = "%"
        f.write("%s" %special_letter)
        f.write("Volume1: %s\n" %input_file1)
        f.write("%s" %special_letter)
        f.write("Volume2: %s\n" %input_file2)
        f.write("%s" %special_letter)
        f.write( "Tag file created by landmark_extraction.py\n\nPoints =\n")

        
        if options.landmark:  
            landmarks_g = find_pair_label_vertices(g, labelpairs, unique_l ) 
            landmarks_h = find_pair_label_vertices(h, labelpairs, unique_l )
            if len(landmarks_g)!= len(landmarks_h):
                print "ERROR! the length of landmarks are different in the 2 graphs: ", len(landmarks_g), " and ", len(landmarks_h)
            else:
                for i in range(len(landmarks_g)):
                    print labelpairs[i],":",landmarks_g[i], ",", landmarks_h[i]
                    if landmarks_g[i]!=['None','None','None'] and landmarks_h[i]!=['None','None','None']:
                        f.write( "%f %f %f %f %f %f\n"%(landmarks_g[i][0],landmarks_g[i][1],landmarks_g[i][2],landmarks_h[i][0],landmarks_h[i][1],landmarks_h[i][2]))
            
            
        if options.landmark3:
            landmarks_g = find_triple_label_vertices(g, labeltriples, unique_triple_l ) 
            landmarks_h = find_triple_label_vertices(h, labeltriples, unique_triple_l )
            if len(landmarks_g)!= len(landmarks_h):
                print "ERROR! the length of landmarks are different in the 2 graphs: ", len(landmarks_g), " and ", len(landmarks_h)
            else:
                for i in range(len(landmarks_g)):
                    print list(labeltriples[i]),":",landmarks_g[i], ",", landmarks_h[i]
                    if landmarks_g[i]!=['None','None','None'] and landmarks_h[i]!=['None','None','None']:
                        f.write( "%f %f %f %f %f %f\n"%(landmarks_g[i][0],landmarks_g[i][1],landmarks_g[i][2],landmarks_h[i][0],landmarks_h[i][1],landmarks_h[i][2]))
            
        f.write( ";")
        f.close()
        
        #### convert with correct lsq and inverse to xfm
        cmd = "tagtoxfm  " + output_file +" "+ output_file[:-3]+".xfm "
        if options.inverse: 
            cmd += "-inverse "       
        if options.lsq6:
            cmd += "-lsq6 "
        if options.lsq7:
            cmd += "-lsq7 "
        if options.lsq9:
            cmd += "-lsq9 "
        if options.lsq10:
            cmd += "-lsq10 "
        if options.lsq12:
            cmd += "-lsq12 "
        if options.clobber:
            cmd += "-clobber "

        os.system(cmd)
        
        
        
