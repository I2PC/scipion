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
from os.path import join

from pyworkflow.utils import Environ
#from pyworkflow.em.data import MovieAlignment


def getEnviron():
    """ Return the envirion settings to run motioncorr programs. """
    """ Setup the environment variables needed to launch Motioncorr. """
    MOTIONCORR_HOME = os.environ.get('MOTIONCORR_HOME', 
                                     join(os.environ['EM_ROOT'], 'motioncorr'))
    environ = Environ(os.environ)
    environ.update({
            'PATH': join(MOTIONCORR_HOME, 'bin'),
            'LD_LIBRARY_PATH': join(os.environ.get('MOTIONCORR_CUDA_LIB', ''))                                    
            }, position=Environ.BEGIN)
    return environ


def getVersion():
    path = os.environ['MOTIONCORR_HOME']
    for v in getSupportedVersions():
        if v in path:
            return v
    return ''


def getSupportedVersions():
    return ['2.1', '2(new)']


def validateVersion(protocol, errors):
    """ Validate if motioncorr version is set properly according
     to installed version and the one set in the config file.
     Params:
        protocol: the input protocol calling to validate
        errors: a list that will be used to add the error message.
    """
    protocol.validatePackageVersion('MOTIONCORR_HOME', errors)


def parseMovieAlignment(logFile):
    """ Get the first and last frames together with the shifts
    between frames aligned. Motioncorr old version
    """
    f = open(logFile)
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
    f = open(logFile)
    first = None
    xshifts = []
    yshifts = []

    for line in f:
        l = line.strip()
        if not '# full frame alignment' in l:
            parts = l.split()
            if first is None: # read the first frame number
                first = int(parts[0]) # take the id from first column
            # take the shifts from the last two columns of the line
            xshifts.append(float(parts[1]))
            yshifts.append(float(parts[2]))
    f.close()
    return xshifts, yshifts
