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

import os
from os.path import join, exists

CTFFIND3 = 'ctffind3.exe'
CTFFIND4 = 'ctffind'
FREALIGN = 'frealign_v9.exe'
FREALIGNMP = 'frealign_v9_mp.exe'
MAGDISTEST = 'mag_distortion_estimate_openmp_1_12_16.exe'
MAGDISTCORR = 'mag_distortion_correct_openmp_8_18_15.exe'
CALC_OCC = 'calc_occ.exe'
RSAMPLE = 'rsample.exe'
UNBLUR = 'unblur'
SUMMOVIE = 'sum_movie_openmp_7_17_15.exe'

def _getCtffind4():
    ctffind4 = join(os.environ['CTFFIND4_HOME'], 'bin', CTFFIND4)
    if exists(ctffind4):
        return ctffind4
    else:
        return join(os.environ['CTFFIND4_HOME'], CTFFIND4)
    
def _getHome(key, default):
    """ Get the required home path, if not present..
    the default value will be used from EM_ROOT.
    """
    return os.environ.get(key, join(os.environ['EM_ROOT'], default))

CTFFIND_PATH = join(os.environ['CTFFIND_HOME'], CTFFIND3)
CTFFIND4_PATH = _getCtffind4()

FREALING_HOME = _getHome('FREALIGN_HOME', 'frealign')
FREALIGN_PATH = join(FREALING_HOME, 'bin', FREALIGN)
FREALIGNMP_PATH = join(FREALING_HOME, 'bin', FREALIGNMP)
CALC_OCC_PATH = join(FREALING_HOME, 'bin', CALC_OCC)
RSAMPLE_PATH = join(FREALING_HOME, 'bin', RSAMPLE)

MAGDIST_HOME = _getHome('MAGDIST_HOME', 'mag_distortion')
MAGDISTEST_PATH = join(MAGDIST_HOME, 'bin', MAGDISTEST)
MAGDISTCORR_PATH = join(MAGDIST_HOME, 'bin', MAGDISTCORR)

UNBLUR_PATH  = join(_getHome('UNBLUR_HOME', 'unblur'), 'bin', UNBLUR)
SUMMOVIE_PATH  = join(_getHome('SUMMOVIE_HOME', 'summovie'), 'bin', SUMMOVIE)
