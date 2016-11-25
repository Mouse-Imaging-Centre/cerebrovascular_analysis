#!/usr/bin/env/python

from distutils.core import setup

setup(name='cerebrovascular_analysis',
      version='1.0',
      license='Modified BSD',
      description='Python code for cerebrovascular analysis',
      long_description = 'Python code for cerebrovascular analysis', 
      author='Sahar Ghanavati, John Sled',
      maintainer_email='sahar.ghanavati@mouseimaging.ca',
      url='https://github.com/sghanavati/',
      platforms="any",
      #packages=['pydpiper', 'applications', 'atoms_and_modules'], 
      scripts=['assess_vessel_track.py', 'subvolume_vessel_tracking.py', 'auto_vessel_tag.py', 'auto_vessel_tracking_functions.py', 
               'auto_vessel_tracking.py', 'postprocess_vessel_track.py', 'feature_extraction.py', 
               'bayesClassifier.py', 'vessel_analysis.py', '/micehome/jgsled/Source/vessel_tracking/lib/vessel_tracking/graph_analysis.py', 
               'vessel_graph_io.py'])
