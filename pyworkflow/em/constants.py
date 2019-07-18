# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
This modules contains constants related to EM
"""
# ------------------ Constants values --------------------------------------

NO_INDEX = 0  # This is the index value for single images

# Sampling rate input modes
SAMPLING_FROM_IMAGE = 0
SAMPLING_FROM_SCANNER = 1

# This is the name for track which data is the source of others
RELATION_SOURCE = 'relation_datasource'
RELATION_TRANSFORM = 'relation_transform'
RELATION_CTF = 'relation_ctf'

UNIT_PIXEL = 'px'
UNIT_PIXEL_FOURIER = '1/px'
UNIT_ANGSTROM = 'A'
UNIT_ANGSTROM_FOURIER = '1/A'

# Fourier Filter options
FILTER_LOW_PASS = 0
FILTER_HIGH_PASS = 1
FILTER_BAND_PASS = 2
FILTER_GAUSSIAN = 3
FILTER_LOW_PASS_NO_DECAY = 4
FILTER_NO_DECAY = 5

# Transform
ALIGN_NONE = 'None'
ALIGN_2D = '2D'            # 2D image alignment
ALIGN_3D = '3D'            # 3D map alignment
ALIGN_PROJ = 'Projection'  # relate projections with 3d map

ALIGNMENTS = [ALIGN_NONE, ALIGN_2D, ALIGN_3D, ALIGN_PROJ]

#Constants related with colormaps for viewers
# Color maps
COLOR_JET = 0
COLOR_TERRAIN = 1
COLOR_GIST_EARTH = 2
COLOR_GIST_NCAR = 3
COLOR_GNU_PLOT = 4
COLOR_GNU_PLOT2 = 5
COLOR_OTHER = 6

COLOR_CHOICES = ['jet', 'terrain',  'gist_earth', 'gist_ncar', 'gnuplot', 'gnuplot2', 'other']

# Axis code
AX_X = 0
AX_Y = 1
AX_Z = 2


# SYMMETRY, conventions described atr: https://scipion-em.github.io/docs/docs/developer/symmetries
# Some notes:
# Icosahedral
#   xmipp and relion define I1, I2, I3 and I4 that corerspond to
#   I222, I222r, In25 and In25z
#   cryosparc has I1 and I2 
#   EMAN always puts the highest symmetry axis on Z. 
#   and a 2 fold axis in X (so it uses I2n5 or I2n5r, not sure)
# DIEDRAL: first axis Z, second X (DX) except cryosparc that uses y (DY)
# Tetraedral, most of the programs use T222, EMAN uses Tz3
SYM_CYCLIC = 0  # CN cyclic symmetry Cn around z axis
SYM_DIHEDRAL = 1  # DN cyclic symmetry plus and extra 2fold symmetry axis around X
SYM_DIHEDRAL_X = SYM_DIHEDRAL
SYM_DIHEDRAL_Y = 2 # DN cyclic symmetry plus and extra 2fold symmetry axis around Y (cryoSparc)
SYM_TETRAHEDRAL = 3  # T 222 Tetahedral symmetry with two-fold symmetry axes along the X, Y, 
                     # and Z axes, a three-fold along axis (1,1,1)
SYM_TETRAHEDRAL_Z3 = 4  # T  three-fold symmetry axis along Z, another three-fold axis 
                        # in the YZ plane such that rotation about the X axis by ~110° 
                        # is a symmetry operation (EMAN convention)

SYM_OCTAHEDRAL = 5  # O

# icosahedric IXXX
# (no crowther 222 and standard in heyman et al 2005 article).
# 2-fold axes on x,y,z axes. With the positive z-axis pointing at the viewer,
# the front-most 5-fold vertices are in yz plane, and the front-most 3-fold
# axes are in the xz plane.
SYM_I222 = 6

# (crowther) 2-fold axes on x,y,z axes. With the positive z-axis pointing at
# the viewer, the front-most 5-fold vertices are in xz plane,
# and the front-most 3-fold axes are in the yz plane.
# (222 rotated 90 degrees around Z)
SYM_I222r = 7 

# '2-fold symmetry along y and 5-fold along z 
SYM_In25 = 8

# 'n25' with 180 degree rotation about x
SYM_In25r = 9

SYM_I2n3 = 10 # Two-fold symmetry along X and 3-fold along Z
SYM_I2n3r = 11 # Idem but rotated 180 degree about Y
SYM_I2n5 = 12 # Two-fold symmetry along Y and 5-fold along Z
SYM_I2n5r = 13 # Idem but rotated 180 degree about X

# Symmetry dictionary
SCIPION_SYM_NAME = {}
SCIPION_SYM_NAME[SYM_CYCLIC] = 'Cn'
SCIPION_SYM_NAME[SYM_DIHEDRAL_X] = 'Dxn'
SCIPION_SYM_NAME[SYM_DIHEDRAL_Y] = 'Dyn'
SCIPION_SYM_NAME[SYM_TETRAHEDRAL] = 'T222'
SCIPION_SYM_NAME[SYM_TETRAHEDRAL_Z3] = 'Tz3'
SCIPION_SYM_NAME[SYM_OCTAHEDRAL] = 'O'
SCIPION_SYM_NAME[SYM_I222] = 'I222'
SCIPION_SYM_NAME[SYM_I222r] = 'I222r'
SCIPION_SYM_NAME[SYM_In25] = 'In25'
SCIPION_SYM_NAME[SYM_In25r] = 'In25r'
SCIPION_SYM_NAME[SYM_I2n3] = 'I2n3'
SCIPION_SYM_NAME[SYM_I2n3r] = 'I2n3r'
SCIPION_SYM_NAME[SYM_I2n5] = 'I2n5'
SCIPION_SYM_NAME[SYM_I2n5r] = 'I2n5r'
