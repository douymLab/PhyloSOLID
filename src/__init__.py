#!/usr/bin/env python3

"""
PhyloSOLID core modules
"""

from .__version__ import __version__

# Export core functions to maintain original import style
from .data_loader import load_all
from .germline_filter import identify_germline_variants
from .mutation_integrator import attach_mutations_to_current_tree
from .scaffold_builder import build_scaffold_tree
from .scrna_classifier import real_time_classifier_predict as scrna_predict
from .scdna_classifier import real_time_classifier_predict as scdna_predict

__all__ = [
    '__version__',
    'load_all',
    'identify_germline_variants',
    'attach_mutations_to_current_tree',
    'build_scaffold_tree',
    'scrna_predict',
    'scdna_predict',
]