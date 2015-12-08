# **************************************************************************
# *
# * Authors:     Grigory Sharov (sharov@igbmc.fr)
# *
# * L'Institut de génétique et de biologie moléculaire et cellulaire (IGBMC)
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
This module contains converter functions that will serve to:
1. Write from base classes to Gctf specific files
2. Read from Gctf files to base classes
"""

import os
import re
from pyworkflow.object import Float


def parseGctfOutput(filename):
    """ Retrieve defocus U, V, angle, crossCorrelation
    and ctfResolution from the output file of the Gctf execution.
    """
    result = None
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
                # that is a 5-th value in the line
                # but remove escape characters first
                result2 = ansi_escape.sub('', line)
                result3 = tuple(map(float,result2.split()[5:6]))
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
    ctfModel._ctffind4_crossCorrelation = Float(ctfFit)
    ctfModel._ctffind4_ctfResolution = Float(ctfResolution)

