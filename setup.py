#!/usr/bin/env python
#
#  Setup script for cerebrovascular_analysis package
#
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import os, sys

setup(name="cerebrovascular_analysis",
      version="0.9",
      description="automatic labelling of mouse cerebral vasculature",
      author="Sahar Ghanavati",
      author_email="john.sled@utoronto.ca",
      packages=['cerebrovascular_analysis'],
      package_dir = {'': 'lib'},
      scripts = ('anatomical_label_calculations_combineall.py', 'assess_vessel_track.py', 'auto_vessel_tag.py', 'auto_vessel_tracking_functions.py', 'auto_vessel_tracking.py', 'bayesClassifier.py', 'bayesianLearner.py', 'bayesTrainer.py', 'compare_w_groundtruth.py', 'compare_w_groundtruth_tocsv.py', 'comparisonofautolabeling.py', 'cylinder2graph.py', 'delete_cage.py', 'delete_internal_leaves_sahar.py', 'delete_internal_segments.py', 'delete_short_leaves.py', 'delete_small_connected_components.py', 'feature_extraction.py', 'feature_MRI_labels.py', 'graph2cylinder.py', 'graph_label2mnc.py', 'graphs_label_average.py', 'keep_connectedcomponents.py', 'keep_labels.py', 'landmark_extraction.py', 'mnc_intensity_remap.py', 'mnc_label2graph.py', 'modify_db_intermediaries.py', 'modifydblabels.py', 'nibbler.py', 'PCA.py', 'perfusionterritory.py', 'perfusionterritoryvariation.py', 'perfusionterritoryvol.py', 'postprocess_vessel_track.py', 'radi_correction.py', 'remove_isolated_edges.py', 'selectlabelsofgraphs.py', 'setup.py', 'simplify_graph.py', 'subvolume_vessel_tracking.py', 'summary_anatomical_label_calculations_to_csv.py', 'unsimplify_graph.py', 'vessel_graph_io.py', 'vessel_num_by_MR_region.py', 'voxel_vote_w_error.py', 'write_posteriorprob2cyl.py')
     )

