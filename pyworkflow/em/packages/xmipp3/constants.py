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

MASK_NONE = 0
MASK_RAISED_COSINE = 0
MASK_CIRCULAR = 1
MASK_FILE = 2

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

FILTER_LOW_PASS = 0
FILTER_HIGH_PASS = 1
FILTER_BAND_PASS = 2

# Rotational spectra find center mode
ROTSPECTRA_CENTER_MIDDLE = 0
ROTSPECTRA_CENTER_FIRST_HARMONIC = 1