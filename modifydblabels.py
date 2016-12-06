# -*- coding: utf-8 -*-
"""
Created on Thu Sep 10 11:10:07 2015

@author: sghanavati
"""

from vessel_tracking import graph_analysis, distribution_analysis
from sys import argv
import os, shelve, string ,time,sys
import commands
from optparse import OptionParser, Option, OptionValueError

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
# basilar branch 1 -> basilar
# internal auditory -> pontine
#ICA branch 1 -> ICA
#ACA branch 1 -> ACA
#Galen branch 1 -> Galen

#96->196
#45->15
#49->15
#102->2
#143->43
#135->35
#206->6


program_name = 'modifydblabels.py'


if __name__ == '__main__':
        
    usage = "Usage: "+program_name+" [options] input_graph  [output_graph]\n"+\
            "   or  "+program_name+" --help";
            
    parser = OptionParser(usage)

    parser.add_option("--clobber", action="store_true", dest="clobber",
                        default=0, help="overwrite output file")
    
    parser.add_option("--remap", dest="remap",
                    help="Label numbers to replace with new label numbers. For example, to turn label 3 into label 11, it would be --remap 3:11. More than one can be specified at once, separated by commas, i.e. --remap 3:11,21:55.",
                    type="string")
    
    parser.add_option("--remove_nonlabel", action="store_true", dest="remove_nonlabel",
                        default=0, help="remove edges with no label of label None")
                         
    options, args = parser.parse_args()

    print '>>> %s [%s]: %s\n' % (time.ctime(time.time()), os.getcwd(), string.join(argv))
    sys.stdout.flush()

    if len(args) == 2:
        input_file, output_file = args
    elif len(args) ==1:
        input_file = args[0]
        output_file = input_file[:-3]+ "_remapped.db"
    else:   
        parser.error("incorrect number of arguments")
      
    if not options.clobber and os.path.exists(output_file):
        raise SystemExit, \
            "The --clobber option is needed to overwrite an existing file."

    try:
        g, attributes = graph_analysis.input_graph(input_file, ["history", "vertex_offsets"])  	# open graph data and copy contents to memory
        print ("Succefully read in the %s\n" %input_file)	
        sys.stdout.flush()
    except:
        print("Error reading in the %s \n" %input_file)
        sys.stdout.flush()


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

        for e in g.edge_list():
            if g.edge_property(e,'label') in tomap_l.keys():
                oldl = g.edge_property(e,'label')
                g.set_edge_property(e,'label',tomap_l[oldl])
                g.set_edge_property(e,'cyl_label',tomap_l[oldl])   #,int(labelpairs[tomap_l.index(g.edge_property(e,'label'))][1]))

    if options.remove_nonlabel:
        for e in g.edge_list():
            if (g.edge_property(e,'label')==0):
                g.remove_edge(e)

    history = '>>> %s: %s' % (time.ctime(time.time()), string.join(argv))           ##argv = ['code.py', 'input.db', 'output']
    if attributes.has_key("history"):
            history = attributes["history"] + "\n" + history
            del attributes['history']
    
    graph_analysis.output_graph(output_file, g, history, attributes)   
    
    print ("Succefully wrote the %s\n" %output_file)
    sys.stdout.flush()
            
    cmd=("graph2cylinder.py %s --use_label --clobber " %output_file)	
    os.system(cmd)


    
