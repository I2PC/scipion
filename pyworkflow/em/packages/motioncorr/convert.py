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

from pyworkflow.utils import Environ

def _getHome(key, default):
    """ Get the required home path, if not present..
    the default value will be used from EM_ROOT.
    """
    return os.environ.get(key, join(os.environ['EM_ROOT'], default))


MOTIONCORR = 'dosefgpu_driftcorr'
MOTIONCOR2 = 'motioncor2'
MOTIONCORR_PATH = join(_getHome('MOTIONCORR_HOME', 'motioncorr'),
                       'bin', MOTIONCORR)
MOTIONCOR2_PATH = join(_getHome('MOTIONCOR2_HOME', 'motioncor2'),
                       'bin', MOTIONCOR2)


def getEnviron():
    """ Return the environ settings to run motioncorr programs. """
    environ = Environ(os.environ)

    if exists(MOTIONCORR_PATH):
        environ.update({'PATH': join(os.environ['MOTIONCORR_HOME'], 'bin')},
                       position=Environ.BEGIN)

    if exists(MOTIONCOR2_PATH):
        environ.update({'PATH': join(os.environ['MOTIONCOR2_HOME'], 'bin')},
                       position=Environ.BEGIN)

    #FIXME: do we need separate libs for motioncor2?
    environ.update({'LD_LIBRARY_PATH': join(os.environ.get('MOTIONCORR_CUDA_LIB', ''))},
                   position=Environ.BEGIN)

    return environ


def parseMovieAlignment(logFile):
    """ Get the first and last frames together with the shifts
    between frames aligned. Motioncorr old version
    """
    f = open(logFile, 'a+')
    first = None
    xshifts = []
    yshifts = []

    for line in f:
        l = line.strip()
        if 'Add Frame #' in l:
            parts = l.split()
            if first is None: # read the first frame number
                first = int(parts[2][1:]) # take the id from #000 format
            # take the id from the last two colums of the line
            xshifts.append(float(parts[-2]))
            yshifts.append(float(parts[-1]))
    f.close()
    return xshifts, yshifts


def parseMovieAlignment2(logFile):
    """ Get the first and last frames together with the shifts
    between frames aligned. Motioncor2 new version
    """
    f = open(logFile, 'a+')
    first = None
    xshifts = []
    yshifts = []

    for line in f:
        l = line.strip()
        if not '#' in l and len(l) > 0:
            parts = l.split()
            if first is None: # read the first frame number
                first = int(parts[0]) # take the id from first column
            # take the shifts from the last two columns of the line
            xshifts.append(float(parts[1]))
            yshifts.append(float(parts[2]))
    f.close()
    return xshifts, yshifts
