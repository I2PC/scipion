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
from pyworkflow.object import Float
import pyworkflow.utils as pwutils


def getVersion():
    path = os.environ['GCTF_HOME']
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
    result = None
    result1 = ()
    ansi_escape = re.compile(r'\x1b[^m]*m')
    if os.path.exists(filename):
        f = open(filename)
        for line in f:
            if 'Final Values' in line:
                # Take DefocusU, DefocusV, Angle and crossCorrelation as a tuple
                # that are the first four values in the line
                result1 = tuple(map(float, line.split()[:4]))
            if 'Resolution limit estimated by EPA' in line:
                # Take ctfResolution as a tuple
                # that is the last value in the line
                # but remove escape characters first
                result2 = ansi_escape.sub('', line)
                result3 = tuple(map(float,result2.split()[-1:]))
                result = result1 + result3
                break
        f.close()
    return result


def setWrongDefocus(ctfModel):
    ctfModel.setDefocusU(-999)
    ctfModel.setDefocusV(-1)
    ctfModel.setDefocusAngle(-999)
    

def readCtfModel(ctfModel, filename, ctf4=False):        
    result = parseGctfOutput(filename)
    if result is None:
        setWrongDefocus(ctfModel)
        ctfFit, ctfResolution = -999, -999 
    else:
        defocusU, defocusV, defocusAngle, ctfFit, ctfResolution = result
        ctfModel.setStandardDefocus(defocusU, defocusV, defocusAngle)
    ctfModel._gctf_crossCorrelation = Float(ctfFit)
    ctfModel._gctf_ctfResolution = Float(ctfResolution)


def getEnviron():
    """ Return the environ settings to run gautomatch programs. """
    environ = pwutils.Environ(os.environ)

    # Take Scipion CUDA library path
    cudaLib = environ.getFirst(('GCTF_CUDA_LIB', 'CUDA_LIB'))
    environ.addLibrary(cudaLib)

    return environ

