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

import numpy
from collections import OrderedDict
from numpy import rad2deg
from numpy.linalg import inv
#from itertools import izip

from pyworkflow.em.packages.grigoriefflab.grigoriefflab import *

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
    partIter = iter(inputSet.get().iterItems(orderBy=['_micId', 'id'],
                                         direction='ASC'))
    for particle, row in izip(partIter, parFile):
        import sys
        print >> sys.stderr, particle.getObjId(), particle.get_CTF(), row
        particle.setTransform(rowToAlignment(row, samplingRate))
        # We assume that each particle have ctfModel
        # in order to be processed in Frealign
        rowToCtfModel(row, particle.getCTF())
        outputSet.append(particle)
    outputSet.setAlignment(ALIGN_PROJ)

def rowToAlignment(alignmentRow, samplingRate):
    """
    Return an Transform object representing the Alignment
    from a given parFile row.
    """
    angles = numpy.zeros(3)
    shifts = numpy.zeros(3)
    alignment = Transform()
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
        from numpy.linalg import inv
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

def readSetOfClasses3D(classes3DSet, fileparList, volumeList):
    """read from frealign .par.
    """
    imgSet = classes3DSet.getImages()

    for ref, volFn in enumerate(volumeList):
        class3D = Class3D()
        class3D.setObjId(ref+1)
        vol = Volume()
        vol.copyObjId(class3D)
        vol.setLocation(volFn)

        class3D.setRepresentative(vol)
        classes3DSet.append(class3D)

        file1 = fileparList[ref]
        f1 = open(file1)
        for l in f1:
            if not l.startswith('C'):
                values = l.split()
                prob = float(values[11])
                if prob > 0:
                    objId = int(values[7])
                    img = imgSet[objId]
                    class3D.append(img)
        f1.close()

        # Check if write function is necessary
        class3D.write()


def parseCtffindOutput(filename):
    """ Retrieve defocus U, V and angle from the
    output file of the ctffind3 execution.
    """
    f = open(filename)
    result = None
    for line in f:
        if 'Final Values' in line:
            # Take DefocusU, DefocusV and Angle as a tuple
            # that are the first three values in the line
            result = tuple(map(float, line.split()[:3]))
            break
    f.close()
    return result

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
    matrix = alignment.getMatrix()
    shifts, angles = geometryFromMatrix(alignment.getMatrix(),True)#####

    return shifts, angles


    alignmentRow.setValue(xmipp.MDL_SHIFT_X, shifts[0])
    alignmentRow.setValue(xmipp.MDL_SHIFT_Y, shifts[1])

    alignmentRow.setValue(xmipp.MDL_SHIFT_Z, shifts[2])
    alignmentRow.setValue(xmipp.MDL_ANGLE_ROT,  angles[0])
    alignmentRow.setValue(xmipp.MDL_ANGLE_TILT, angles[1])
    alignmentRow.setValue(xmipp.MDL_ANGLE_PSI,  angles[2])
    alignmentRow.setValue(xmipp.MDL_FLIP, flip)
