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

