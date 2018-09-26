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

# Program common names
FREALIGN_APP = 'FREALIGN'
CTFFIND4_APP = 'CTFFIND4'
V4_0_15 = '4.0.15'
V4_1_10 = '4.1.10'

MAGDIST_HOME_VAR = 'MAGDIST_HOME'

SUMMOVIE_HOME = 'SUMMOVIE_HOME'

UNBLUR_HOME = 'UNBLUR_HOME'

FREALIGN_HOME_VAR = 'FREALIGN_HOME'

CTFFIND4_HOME = 'CTFFIND4_HOME'
CTFFIND_HOME = 'CTFFIND_HOME'

CTFFIND3 = 'ctffind3.exe'
CTFFIND3MP = 'ctffind3_mp.exe'
CTFFIND4 = 'ctffind'
CTFTILT = 'ctftilt.exe'
CTFTILTMP = 'ctftilt_mp.exe'
FREALIGN = 'frealign_v9.exe'
FREALIGNMP = 'frealign_v9_mp.exe'
MAGDISTEST = 'mag_distortion_estimate_openmp.exe'
MAGDISTCORR = 'mag_distortion_correct_openmp.exe'
CALC_OCC = 'calc_occ.exe'
RSAMPLE = 'rsample.exe'
UNBLUR = 'unblur_openmp.exe'
SUMMOVIE = 'sum_movie_openmp.exe'


def getVersion(var=FREALIGN_APP):
    varHome = var + '_HOME'
    path = os.environ[varHome]
    for v in getSupportedVersions(var):
        if v in path or v in os.path.realpath(path):
            return v
    return ''


def getSupportedVersions(var=FREALIGN_APP):
    if var == 'UNBLUR':
        return ['1.0_150529', '1.0.2']
    elif var == CTFFIND4_APP:
        return ['4.0.15', '4.1.5', '4.1.8', V4_1_10]
    else:  # FREALIGN
        return ['9.07']


def _getCtffind4():
    ctffind4 = join(os.environ[CTFFIND4_HOME], 'bin', CTFFIND4)
    if exists(ctffind4):
        return ctffind4
    else:
        return join(os.environ[CTFFIND4_HOME], CTFFIND4)


def _getHome(key, default):
    """ Get the required home path, if not present..
    the default value will be used from EM_ROOT.
    """
    return os.environ.get(key, join(os.environ['EM_ROOT'], default))

CTFFIND_PATH = join(os.environ[CTFFIND_HOME], CTFFIND3)
CTFFINDMP_PATH = join(os.environ[CTFFIND_HOME], CTFFIND3MP)
CTFFIND4_PATH = _getCtffind4()

CTFTILT_PATH = join(os.environ[CTFFIND_HOME], CTFTILT)
CTFTILTMP_PATH = join(os.environ[CTFFIND_HOME], CTFTILTMP)

FREALIGN_HOME = _getHome(FREALIGN_HOME_VAR, 'frealign')
FREALIGN_PATH = join(FREALIGN_HOME, 'bin', FREALIGN)
FREALIGNMP_PATH = join(FREALIGN_HOME, 'bin', FREALIGNMP)
CALC_OCC_PATH = join(FREALIGN_HOME, 'bin', CALC_OCC)
RSAMPLE_PATH = join(FREALIGN_HOME, 'bin', RSAMPLE)

MAGDIST_HOME = _getHome(MAGDIST_HOME_VAR, 'mag_distortion')
MAGDISTEST_PATH = join(MAGDIST_HOME, 'bin', MAGDISTEST)
MAGDISTCORR_PATH = join(MAGDIST_HOME, 'bin', MAGDISTCORR)

UNBLUR_PATH = join(_getHome(UNBLUR_HOME, 'unblur'), 'bin', UNBLUR)
SUMMOVIE_PATH = join(_getHome(SUMMOVIE_HOME, 'summovie'), 'bin', SUMMOVIE)


def validateMagDistorsionInstallation():
    """ Check if the installation of this protocol is correct.
    Can't rely on package function since this is a "multi package" package
    Returning an empty list means that the installation is correct
    and there are not errors. If some errors are found, a list with
    the error messages will be returned.
    """
    missingPaths = []

    if not os.path.exists(MAGDIST_HOME):
        missingPaths.append("%s : %s" % (MAGDIST_HOME_VAR, MAGDIST_HOME))
    return missingPaths
