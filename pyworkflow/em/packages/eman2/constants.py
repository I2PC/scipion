# **************************************************************************
# *
# *  Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk)
# *
# * MRC Laboratory of Molecular Biology (MRC-LMB)
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
"""
This modules contains constants related to EMAN2 protocols
"""

#------------------ Constants values --------------------------------------

# ctf processing type
HIRES = 0
MIDRES = 1
LORES = 2

# centering algorithms
XFORM_NOCENTER = 0
XFORM_CENTER = 1
XFORM_CENTERACF = 2
XFORM_CENTEROFMASS = 3
XFORM_CENTER_NONE = 4

# comparators
CMP_CCC = 0
CMP_DOT = 1
CMP_FRC = 2
CMP_LOD = 3
CMP_OPTSUB = 4
CMP_OPTVARIANCE = 5
CMP_PHASE = 6
CMP_QUADMINDOT = 7
CMP_SQEUCLIDEAN = 8
CMP_VERTICAL = 9
CMP_NONE = 10

# aligners
ALN_FRM2D = 0
ALN_ROTATE_FLIP = 1
ALN_ROTATE_FLIP_ITERATIVE = 2
ALN_ROTATE_PRECENTER = 3
ALN_ROTATE_TRANS_FLIP_SCALE = 4
ALN_ROTATE_TRANS_FLIP_SCALE_ITER = 5
ALN_ROTATE_TRANS_SCALE_ITER = 6
ALN_ROTATE_TRANSLATE = 7
ALN_ROTATE_TRANSLATE_FLIP = 8
ALN_ROTATE_TRANSLATE_FLIP_ITERATIVE = 9
ALN_ROTATE_TRANSLATE_FLIP_RESAMPLE = 10
ALN_ROTATE_TRANSLATE_ITERATIVE = 11
ALN_ROTATE_TRANSLATE_RESAMPLE = 12
ALN_ROTATE_TRANSLATE_SCALE = 13
ALN_ROTATE_TRANSLATE_TREE  =14
ALN_ROTATIONAL = 15
ALN_ROTATIONAL_ITERATIVE = 16
ALN_RTF_EXHAUSTIVE = 17
ALN_RTF_SLOW_EXHAUSTIVE = 18
ALN_SCALE = 19
ALN_SYMALIGN = 20
ALN_SYMALIGNQUAT = 21
ALN_TRANSLATIONAL = 22
ALN_NONE = 23

RALN_NONE = 0
RALN_REFINE = 1
RALN_REFINE_3D = 2
RALN_REFINE_3D_GRID = 3

# averagers
AVG_ABSMAXMIN = 0
AVG_CTF_AUTO = 1
AVG_CTF_WEIGHT = 2
AVG_CTF_WEIGHT_AUTOFILT = 3
AVG_CTFW_AUTO = 4
AVG_ITERATION = 5
AVG_LOCALWEIGHT = 6
AVG_MEAN = 7
AVG_MEAN_TOMO = 8
AVG_MINMAX = 9
AVG_SIGMA = 10
AVG_WEIGHTEDFOURIER = 11

#processors normalize
PROC_NORMALIZE = 0
PROC_NORMALIZE_BYMASS = 1
PROC_NORMALIZE_CIRCLEMEAN = 2
PROC_NORMALIZE_EDGEMEAN = 3
PROC_NORMALIZE_LOCAL = 4
PROC_NORMALIZE_LREDGE = 5
PROC_NORMALIZE_MASK = 6
PROC_NORMALIZE_MAXMIN = 7
PROC_NORMALIZE_RAMP_NORMVAR = 8
PROC_NORMALIZE_ROWS = 9
PROC_NORMALIZE_TOIMAGE = 10
PROC_NORMALIZE_UNITLEN = 11
PROC_NORMALIZE_UNITSUM = 12
PROC_NONE = 13

# Reconstruction methods
RECON_BACKPROJ = 0
RECON_FOURIER = 1
RECON_FOURIER_PLANE = 2
RECON_FOURIER_SIMPLE = 3
RECON_NN4 = 4
RECON_NN4_CTF = 5
RECON_NN4_CTF_RECT = 6
RECON_NN4_CTFW = 7
RECON_NN4_CTFWS = 8
RECON_NN4_RECT = 9
RECON_NNSSNR = 10
RECON_NNSSNR_CTF = 11
RECON_WIENER_FOURIER = 12

# modes to reconstruct with fourier method
FOURIER_NEIGHBOR = 0
FOURIER_GAUSS2 = 1
FOURIER_GAUSS3 = 2
FOURIER_GAUSS5 = 3
FOURIER_GAUSS5_SLOW = 4
FOURIER_GYPERGEOM5 = 5
FOURIER_EXPERIMENTAL = 6

# speed
SPEED_1 = 0
SPEED_2 = 1
SPEED_3 = 2
SPEED_4 = 3
SPEED_5 = 4
SPEED_6 = 5
SPEED_7 = 6

# Keep parameter for e2make3d.py
KEEP_PERCENTAGE = 0
KEEP_STDDEV = 1
KEEP_ABSQUAL = 2

# Amplitude correction type for e2refine_easy
AMP_AUTO = 0
AMP_STRUCFAC = 1
AMP_FLATTEN = 2
AMP_NONE = 3

# tophat filter for e2refine_easy
TOPHAT_NONE = 0
TOPHAT_LOCAL = 1
TOPHAT_GLOBAL = 2

# e2boxer autopick modes
AUTO_LOCAL = 0
AUTO_REF = 1
AUTO_GAUSS = 2
AUTO_CONVNET = 3

WIKI_URL = "[[http://blake.bcm.edu/emanwiki/EMAN2][Wiki]]"
