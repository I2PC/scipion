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
#------------------ Constants values --------------------------------------

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
ALIGN_2D   = '2D'         # 2D image alignment
ALIGN_3D   = '3D'         # 3D map alignment
ALIGN_PROJ = 'Projection' # relate projections with 3d map

ALIGNMENTS = [ALIGN_NONE, ALIGN_2D, ALIGN_3D, ALIGN_PROJ]

#SYMMETRY, follows chimera conventions https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/midas/sym.html
counter = 0
#cyclic CN
SYM_CYCLIC  = counter # symmetry Cn around axis
counter += 1
#icosahedric IXXX
SYM_I222  = counter # (no crowther 222 and standard in heyman et al 2005 article).
              # 2-fold axes on x,y,z axes. With the positive z-axis pointing at
              # the viewer, the front-most 5-fold vertices are in yz plane,
              #  and the front-most 3-fold axes are in the xz plane.
counter += 1
SYM_I222r = counter # (crowther) 2-fold axes on x,y,z axes. With the positive z-axis
              # pointing at the viewer, the front-most 5-fold vertices are in xz plane,
              # and the front-most 3-fold axes are in the yz plane.
counter += 1
SYM_In25   = counter # '2-fold symmetry along y and 5-fold along z
counter += 1
SYM_In25r   = counter # 'n25' with 180 degree rotation about x

#Diccionario sym
SCIPION_SYM_NAME={}
SCIPION_SYM_NAME[SYM_CYCLIC]  = 'Cn'
SCIPION_SYM_NAME[SYM_I222]  = 'I222'
SCIPION_SYM_NAME[SYM_I222r] = 'I222r'
SCIPION_SYM_NAME[SYM_In25]   = 'In25'
SCIPION_SYM_NAME[SYM_In25r]   = 'In25r'