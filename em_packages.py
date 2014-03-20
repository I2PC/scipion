#!/usr/bin/env python

import os

# EM-Software root
os.environ['EM_SOFTWARE'] = os.path.join(os.environ['PW_HOME'], 'software', 'em')

# Xmipp 3
XMIPP_HOME = os.path.join(os.environ['EM_SOFTWARE'], 'xmipp')
# Spider 
SPIDER_DIR = os.path.join(os.environ['EM_SOFTWARE'], 'spider')
# Eman2
EMAN2DIR = os.path.join(os.environ['EM_SOFTWARE'], 'EMAN2-1_alpha1')
# frealign 8.11
FREALIGN_HOME = os.path.join(os.environ['EM_SOFTWARE'], 'brandeis', 'frealign_v8.11')
# Ctffind 3
CTFFIND_HOME = os.path.join(os.environ['EM_SOFTWARE'], 'brandeis', 'ctffind')
# Relion
RELION_HOME = os.path.join(os.environ['EM_SOFTWARE'], 'relion-1.2')
# xmipp_optical_alignment
OPT_ALIGN_HOME = os.path.join(os.environ['EM_SOFTWARE'], 'opt_flow_alignment')
# BSOFT
BSOFT_HOME = os.path.join(os.environ['EM_SOFTWARE'], 'bsoft')

