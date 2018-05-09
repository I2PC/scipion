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
import numpy
from collections import OrderedDict

from pyworkflow.em.packages.gctf import GCTF_HOME
from pyworkflow.em.constants import ALIGN_2D, ALIGN_3D, ALIGN_PROJ, ALIGN_NONE
from pyworkflow.object import ObjectWrap
import pyworkflow.em.metadata as md
import pyworkflow.utils as pwutils


CTF_DICT = OrderedDict([
       ("_defocusU", md.RLN_CTF_DEFOCUSU),
       ("_defocusV", md.RLN_CTF_DEFOCUSV),
       ("_defocusAngle", md.RLN_CTF_DEFOCUS_ANGLE)
       ])


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
        ansi_escape = re.compile(r'\x1b[^m]*m')
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
                resol = ansi_escape.sub('', line)
                result[5] = float(resol.strip().split()[-1])
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

    # Avoid creation of phaseShift
    if ctfPhaseShift != 0:
        ctfModel.setPhaseShift(ctfPhaseShift)


def getEnviron():
    """ Return the environ settings to run gautomatch programs. """
    environ = pwutils.Environ(os.environ)

    # Take Scipion CUDA library path
    cudaLib = environ.getFirst(('GCTF_CUDA_LIB', 'CUDA_LIB'))
    environ.addLibrary(cudaLib)

    return environ


def writeSetOfCoordinates(coordDir, coordSet, micsSet):
    """ Write a star file on metadata format for each micrograph
    on the coordSet.
    Params:
        coordDir: the directory where the .star files will be written.
        coordSet: the SetOfCoordinates that will be read.
        micsSet: the SetOfMicrographs that will be read.
    """
    header = """
data_

loop_
_rlnCoordinateX #1
_rlnCoordinateY #2
"""

    # Create a dictionary with the pos filenames for each micrograph
    posDict = {}
    for mic in micsSet:
        micBase = pwutils.removeBaseExt(mic.getFileName())
        posDict[mic.getObjId()] = pwutils.join(coordDir, micBase,
                                               micBase + '_coords.star')

    f = None
    lastMicId = None

    # Iterate only once over the whole SetOfCoordinates, but ordering by
    # micrograph Id, so we can detect when there are coordinates from a
    # new micrographs to write the new star file
    for coord in coordSet.iterItems(orderBy='_micId'):
        micId = coord.getMicId()

        if micId != lastMicId:  # Detect there is a new micrograph
            if f:  # we need to close previous opened file
                f.close()
            f = open(posDict[micId], 'w')
            f.write(header)
            lastMicId = micId

        f.write("%d %d\n" % coord.getPosition())

    if f:
        f.close()


def rowToCtfModel(ctfRow, ctfModel):
    """ Create a CTFModel from a row of a meta """
    if ctfRow.containsAll(CTF_DICT):
        for attr, label in CTF_DICT.iteritems():
            value = ctfRow.getValue(label)
            if not hasattr(ctfModel, attr):
                setattr(ctfModel, attr, ObjectWrap(value))
            else:
                getattr(ctfModel, attr).set(value)

        ctfModel.standardize()
    else:
        ctfModel = None

    return ctfModel


def getShifts(transform, alignType):
    """
    is2D == True-> matrix is 2D (2D images alignment)
            otherwise matrix is 3D (3D volume alignment or projection)
    invTransform == True  -> for xmipp implies projection
                          -> for xmipp implies alignment
    """
    if alignType == ALIGN_NONE:
        return None

    inverseTransform = alignType == ALIGN_PROJ
    # only flip is meaningful if 2D case
    # in that case the 2x2 determinant is negative
    flip = False
    matrix = transform.getMatrix()
    if alignType == ALIGN_2D:
        # get 2x2 matrix and check if negative
        flip = bool(numpy.linalg.det(matrix[0:2, 0:2]) < 0)
        if flip:
            matrix[0, :2] *= -1.  # invert only the first two columns keep x
            matrix[2, 2] = 1.  # set 3D rot
        else:
            pass

    elif alignType == ALIGN_3D:
        flip = bool(numpy.linalg.det(matrix[0:3, 0:3]) < 0)
        if flip:
            matrix[0, :4] *= -1.  # now, invert first line including x
            matrix[3, 3] = 1.  # set 3D rot
        else:
            pass

    else:
        pass
        # flip = bool(numpy.linalg.det(matrix[0:3,0:3]) < 0)
        # if flip:
        #    matrix[0,:4] *= -1.#now, invert first line including x
    shifts = geometryFromMatrix(matrix, inverseTransform)

    return shifts


def geometryFromMatrix(matrix, inverseTransform):
    from pyworkflow.em.transformations import translation_from_matrix
    if inverseTransform:
        matrix = numpy.linalg.inv(matrix)
        shifts = -translation_from_matrix(matrix)
    else:
        shifts = translation_from_matrix(matrix)
    return shifts
