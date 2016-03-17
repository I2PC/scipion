# **************************************************************************
# *
# * Authors:     Josue Gomez Blanco (jgomez@cnb.csic.es)
# *
# *
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
This module contains converter functions that will serve to:
1. Write from base classes to Grigorieff packages specific files
2. Read from Grigo packs files to base classes
"""

import os

from itertools import izip
import numpy
from collections import OrderedDict
from numpy import rad2deg
from numpy.linalg import inv

from pyworkflow.object import Float
import pyworkflow.em as em


HEADER_COLUMNS = ['INDEX', 'PSI', 'THETA', 'PHI', 'SHX', 'SHY', 'MAG',
                  'FILM', 'DF1', 'DF2', 'ANGAST', 'OCC',
                  '-LogP', 'SIGMA', 'SCORE', 'CHANGE']


class FrealignParFile(object):
    """ Handler class to read/write frealign metadata."""
    def __init__(self, filename, mode='r'):
        self._file = open(filename, mode)
        self._count = 0

    def __iter__(self):
        """PSI   THETA     PHI       SHX       SHY     MAG  FILM      DF1      DF2  ANGAST     OCC     -LogP      SIGMA   SCORE  CHANGE
        """
        for line in self._file:
            line = line.strip()
            if not line.startswith('C'):
                row = OrderedDict(zip(HEADER_COLUMNS, line.split()))
                yield row

    def close(self):
        self._file.close()


def readSetOfParticles(inputSet, outputSet, parFileName):
    """
     Iterate through the inputSet and the parFile lines
     and populate the outputSet with the same particles
     of inputSet, but with the angles and shift (3d alignment)
     updated from the parFile info.
     It is assumed that the order of iteration of the particles
     and the lines match and have the same number.
    """
    #create dictionary that matches input particles with param file
    samplingRate = inputSet.getSamplingRate()
    parFile = FrealignParFile(parFileName)
    partIter = iter(inputSet.iterItems(orderBy=['_micId', 'id'], direction='ASC'))
    
    for particle, row in izip(partIter, parFile):        
        particle.setTransform(rowToAlignment(row, samplingRate))
        # We assume that each particle have ctfModel
        # in order to be processed in Frealign
        # JMRT: Since the CTF will be set, we can setup
        # an empty CTFModel object
        if not particle.hasCTF():
            particle.setCTF(em.CTFModel())
        rowToCtfModel(row, particle.getCTF())
        outputSet.append(particle)
    outputSet.setAlignment(em.ALIGN_PROJ)


def rowToAlignment(alignmentRow, samplingRate):
    """
    Return an Transform object representing the Alignment
    from a given parFile row.
    """
    angles = numpy.zeros(3)
    shifts = numpy.zeros(3)
    alignment = em.Transform()
    # PSI   THETA     PHI       SHX       SHY
    angles[0] = float(alignmentRow.get('PSI'))
    angles[1] = float(alignmentRow.get('THETA'))
    angles[2] = float(alignmentRow.get('PHI'))
    shifts[0] = float(alignmentRow.get('SHX'))/samplingRate
    shifts[1] = float(alignmentRow.get('SHY'))/samplingRate

    M = matrixFromGeometry(shifts, angles)
    alignment.setMatrix(M)

    return alignment


def matrixFromGeometry(shifts, angles):
    """ Create the transformation matrix from a given
    2D shifts in X and Y...and the 3 euler angles.
    """
    inverseTransform = True
    from pyworkflow.em.transformations import euler_matrix
    from numpy import deg2rad
    radAngles = -deg2rad(angles)

    M = euler_matrix(radAngles[0], radAngles[1], radAngles[2], 'szyz')
    if inverseTransform:
        M[:3, 3] = -shifts[:3]
        M = inv(M)
    else:
        M[:3, 3] = shifts[:3]

    return M


def rowToCtfModel(ctfRow, ctfModel):
    defocusU = float(ctfRow.get('DF1'))
    defocusV = float(ctfRow.get('DF2'))
    defocusAngle = float(ctfRow.get('ANGAST'))
    ctfModel.setStandardDefocus(defocusU, defocusV, defocusAngle)


#-------------- Old fuctions (before using EMX matrix for alignment) ------
def parseCtffindOutput(filename):
    """ Retrieve defocus U, V and angle from the
    output file of the ctffind3 execution.
    """
    result = None
    if os.path.exists(filename):
        f = open(filename)
        for line in f:
            if 'Final Values' in line:
                # Take DefocusU, DefocusV and Angle as a tuple
                # that are the first three values in the line
                result = tuple(map(float, line.split()[:3]))
                break
        f.close()
    return result


def parseCtffind4Output(filename):
    """ Retrieve defocus U, V and angle from the
    output file of the ctffind3 execution.
    """
    result = None
    if os.path.exists(filename):
        f = open(filename)
        for line in f:
            if not line.startswith("#"):
                result = tuple(map(float, line.split()[1:]))
        f.close()
    return result


def ctffindOutputVersion(filename):
    """ Detect the ctffind version (3 or 4) that produced
    the given filename.
    """
    f = open(filename)
    for line in f:
        if 'Output from CTFFind version 4.' in line:
            return 4
    return 3


def setWrongDefocus(ctfModel):
    ctfModel.setDefocusU(-999)
    ctfModel.setDefocusV(-1)
    ctfModel.setDefocusAngle(-999)
    
    
def readCtfModel(ctfModel, filename, ctf4=False):        
    if not ctf4:
        result = parseCtffindOutput(filename)
        if result is None:
            setWrongDefocus(ctfModel)
        else:
            defocusU, defocusV, defocusAngle = result
            ctfModel.setStandardDefocus(defocusU, defocusV, defocusAngle)
    else:
        result =  parseCtffind4Output(filename)
        if result is None:
            setWrongDefocus(ctfModel)
            ctfFit, ctfResolution = -999, -999 
        else:
            defocusU, defocusV, defocusAngle, _, ctfFit, ctfResolution = result
            ctfModel.setStandardDefocus(defocusU, defocusV, defocusAngle)
        ctfModel._ctffind4_crossCorrelation = Float(ctfFit)
        ctfModel._ctffind4_ctfResolution = Float(ctfResolution)


def geometryFromMatrix(matrix):
    from pyworkflow.em.transformations import translation_from_matrix, euler_from_matrix
    inverseTransform = True

    if inverseTransform:
        matrix = inv(matrix)
        shifts = -translation_from_matrix(matrix)
    else:
        shifts = translation_from_matrix(matrix)
    angles = -rad2deg(euler_from_matrix(matrix, axes='szyz'))
    return shifts, angles


def geometryFromAligment(alignment):
    shifts, angles = geometryFromMatrix(alignment.getMatrix(),True)#####

    return shifts, angles


def _createErrorCtf4Fn(self, filename):
            f = open(filename, 'w+')
            lines = """# Error report file 
  -999       -999       -999       -999       -999      -999       -999     
"""
            f.write(lines)
            f.close()


def _createErrorCtf3Fn(self, filename):
    f = open(filename, 'w+')
    lines = """-999    -999       -999     -999  Final Values"""
    f.write(lines)
    f.close()


def writeShiftsMovieAlignment(movie, shiftsFn, s0, sN):
    movieAlignment=movie.getAlignment()
    shiftListX, shiftListY = movieAlignment.getShifts()
    
    # Generating metadata for global shifts
    a0, aN = movieAlignment.getRange()
    alFrame = a0
    
    if s0 < a0:
        diff = a0 - s0
        initShifts = "0.0000 " * diff
    else:
        initShifts = ""
    
    if sN > aN:
        diff = sN - aN
        finalShifts = "0.0000 " * diff
    else:
        finalShifts = ""
    
    shiftsX = ""
    shiftsY = ""
    for shiftX, shiftY in izip(shiftListX, shiftListY):
        if alFrame >= s0 and alFrame <= sN:
            shiftsX = shiftsX + "%0.4f " % shiftX
            shiftsY = shiftsY + "%0.4f " % shiftY
        alFrame += 1
    
    f=open(shiftsFn,'w')
    shifts = (initShifts + shiftsX + " " + finalShifts + "\n" 
              + initShifts + shiftsY + " " + finalShifts)
    f.write(shifts)
    f.close()
