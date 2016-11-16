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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************
"""
This modules contains constants related to Xmipp3 protocols
"""

#------------------ Constants values --------------------------------------

MASK_FILL_VALUE = 0
MASK_FILL_MIN = 1
MASK_FILL_MAX = 2
MASK_FILL_AVG = 3

PROJECT_FOURIER = 0
PROJECT_REALSPACE = 1

KERNEL_NEAREST = 0
KERNEL_LINEAR = 1
KERNEL_BSPLINE = 2

SELECT_NONE = 0
SELECT_MAXCC = 1
SELECT_PERCENTAGE = 2
SELECT_CLASSPERCENTAGE = 3

RECONSTRUCT_FOURIER = 0
RECONSTRUCT_ART = 1
RECONSTRUCT_WBP = 2

FILTER_SPACE_FOURIER = 0
FILTER_SPACE_REAL = 1
FILTER_SPACE_WAVELET = 2

# Rotational spectra find center mode
ROTSPECTRA_CENTER_MIDDLE = 0
ROTSPECTRA_CENTER_FIRST_HARMONIC = 1

# 3D Geometrical Mask
MASK3D_SPHERE = 0
MASK3D_BOX = 1
MASK3D_CROWN = 2
MASK3D_CYLINDER = 3
MASK3D_GAUSSIAN = 4
MASK3D_RAISED_COSINE = 5
MASK3D_RAISED_CROWN = 6

# Mask Choice
SOURCE_GEOMETRY=0
SOURCE_MASK=1

# 2D Geometrical Mask
MASK2D_CIRCULAR = 0
MASK2D_BOX = 1
MASK2D_CROWN = 2
MASK2D_GAUSSIAN = 3
MASK2D_RAISED_COSINE = 4
MASK2D_RAISED_CROWN = 5

# Threshold value substitute
FILL_VALUE = 0
FILL_BINARIZE = 1
FILL_AVG = 2

# Reconstruction methods
RECONSTRUCT_FOURIER = 0
RECONSTRUCT_WSLART = 1

# Micrograph type constants for particle extraction
SAME_AS_PICKING = 0
OTHER = 1
