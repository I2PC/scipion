# **************************************************************************
# *
# * Authors:     Grigory Sharov (sharov@igbmc.fr)
# *
# * L'Institut de genetique et de biologie moleculaire et cellulaire (IGBMC)
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
This module contains converter functions that will serve to:
1. Write from base classes to Gctf specific files
2. Read from Gctf files to base classes
"""

import os
import re

from pyworkflow.em.packages.gctf import GCTF_HOME
from pyworkflow.object import Float
import pyworkflow.utils as pwutils


def getVersion():
    path = os.environ[GCTF_HOME]
    for v in getSupportedVersions():
        if v in path or v in os.path.realpath(path):
            return v
    return ''


def getSupportedVersions():
    return ['0.50', '1.06']


def parseGctfOutput(filename):
    """ Retrieve defocus U, V, angle, crossCorrelation
    and ctfResolution from the output file of the Gctf execution.
    """

    if os.path.exists(filename):
        # Create an empty list with: defU, defV, angle, CC and resolution
        result = [0.] * 6
        f = open(filename)
        for line in f:
            if 'Final Values' in line:
                parts = line.strip().split()
                # line = DefocusU, DefocusV, Angle, crossCorrelation, Final, Values
                # OR
                # line = DefocusU, DefocusV, Angle, ctfPhaseShift, crossCorrelation, Final, Values
                length = len(line.split())
                # Always map defU, defV and angle
                result[0:3] = map(float, parts[0:3])

                if parts[4] == 'Final': # no ctfPhaseShift
                    result[3] = float(parts[3])
                else:
                    result[3] = float(parts[4]) # CC is now in position 4
                    result[4] = float(parts[3]) # get ctfPhaseShift
            if 'Resolution limit estimated by EPA' in line:
                # Take ctfResolution as a tuple
                # that is the last value in the line
                # but remove escape characters first
                result[5] = float(line.strip().split()[-2])
                break
        f.close()
    else:
        result = None
        print "Warning: Missing file: ", filename

    return result


def setWrongDefocus(ctfModel):
    ctfModel.setDefocusU(-999)
    ctfModel.setDefocusV(-1)
    ctfModel.setDefocusAngle(-999)
    

def readCtfModel(ctfModel, filename, ctf4=False):        
    result = parseGctfOutput(filename)
    if result is None:
        setWrongDefocus(ctfModel)
        ctfFit, ctfResolution, ctfPhaseShift = -999, -999, -999
    else:
        defocusU, defocusV, defocusAngle, ctfFit, ctfPhaseShift, ctfResolution = result
        ctfModel.setStandardDefocus(defocusU, defocusV, defocusAngle)
    ctfModel.setFitQuality(ctfFit)
    ctfModel.setResolution(ctfResolution)
    ctfModel.setPhaseShift(ctfPhaseShift)


def getEnviron():
    """ Return the environ settings to run gautomatch programs. """
    environ = pwutils.Environ(os.environ)

    # Take Scipion CUDA library path
    cudaLib = environ.getFirst(('GCTF_CUDA_LIB', 'CUDA_LIB'))
    environ.addLibrary(cudaLib)

    return environ

