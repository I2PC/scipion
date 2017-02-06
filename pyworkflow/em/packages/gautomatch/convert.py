# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# * Authors:     Grigory Sharov (sharov@igbmc.fr)
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
"""
This module contains converter functions that will serve to:
1. Write from base classes to Gautomatch specific files
2. Read from Gautomatch files to base classes
"""

import os
from collections import OrderedDict

import pyworkflow.em as em
import pyworkflow.em.metadata as md
from pyworkflow.object import ObjectWrap
import pyworkflow.utils as pwutils



COOR_DICT = OrderedDict([
             ("_x", md.RLN_IMAGE_COORD_X),
             ("_y", md.RLN_IMAGE_COORD_Y)
             ])

COOR_EXTRA_LABELS = [
    # Additional autopicking-related metadata
    md.RLN_PARTICLE_AUTOPICK_FOM,
    md.RLN_PARTICLE_CLASS,
    md.RLN_ORIENT_PSI
    ]


def rowToObject(row, obj, attrDict, extraLabels=[]):
    """ This function will convert from a XmippMdRow to an EMObject.
    Params:
        row: the XmippMdRow instance (input)
        obj: the EMObject instance (output)
        attrDict: dictionary with the map between obj attributes(keys) and
            row MDLabels in Xmipp (values).
        extraLabels: a list with extra labels that could be included
            as properties with the label name such as: _rlnSomeThing
    """
    obj.setEnabled(row.getValue(md.RLN_IMAGE_ENABLED, 1) > 0)

    for attr, label in attrDict.iteritems():
        value = row.getValue(label)
        if not hasattr(obj, attr):
            setattr(obj, attr, ObjectWrap(value))
        else:
            getattr(obj, attr).set(value)

    attrLabels = attrDict.values()

    for label in extraLabels:
        if label not in attrLabels and row.hasLabel(label):
            labelStr = md.label2Str(label)
            setattr(obj, '_' + labelStr, row.getValueAsObject(label))


def rowToCoordinate(coordRow):
    """ Create a Coordinate from a row of a meta """
    # Check that all required labels are present in the row
    if coordRow.containsAll(COOR_DICT):
        coord = em.Coordinate()
        rowToObject(coordRow, coord, COOR_DICT, extraLabels=COOR_EXTRA_LABELS)

        # FIXME: THE FOLLOWING IS NOT CLEAN
        try:
            coord.setMicId(int(coordRow.getValue(md.RLN_MICROGRAPH_NAME)))
        except Exception:
            pass
    else:
        coord = None

    return coord


def readSetOfCoordinates(workDir, micSet, coordSet, suffix=None):
    """ Read from coordinates from Gautomatch .star files.
    For a micrograph: mic1.mrc, the expected coordinate file is:
    mic1_automatch.star
    Params:
        workDir: where the Gautomatch output files are located.
        micSet: the SetOfMicrographs.
        coordSet: the SetOfCoordinates that will be populated.
        suffix: input coord file suffix
    """
    if suffix == None:
        suffix = '_automatch.star'

    for mic in micSet:
        micBase = pwutils.removeBaseExt(mic.getFileName())
        fnCoords = os.path.join(workDir, micBase + suffix)
        readCoordinates(mic, fnCoords, coordSet)


def readCoordinates(mic, fileName, coordsSet):
    if os.path.exists(fileName):
        for row in md.iterRows(fileName):
            coord = rowToCoordinate(row)
            coord.setX(coord.getX())
            coord.setY(coord.getY())
            coord.setMicrograph(mic)
            coordsSet.append(coord)


def writeSetOfCoordinates(workDir, coordSet, isGlobal=False):
    """ Write set of coordinates from md to star file.
    Used only for exclusive picking. Creates .star files with
    bad coordinates (for each mic) and/or a single .star file with
    global detector defects.
    """
    # Gautomatch cannot read default star header (with # XMIPP_STAR_1 *),
    # so we write directly to file

    header = """
data_

loop_
_rlnCoordinateX #1
_rlnCoordinateY #2
_rlnAnglePsi #3
_rlnClassNumber #4
_rlnAutopickFigureOfMerit #5
"""

    if isGlobal:
        fnOut = os.path.join(workDir, 'defects.star')
        f = open(fnOut, 'w')
        f.write(header)
        for coord in coordSet:
            f.write("%0.6f %0.6f %0.6f %d %0.6f\n" %
                    (coord.getX(), coord.getY(), 0, 0, 0))
        f.close()

    else:
        for mic in coordSet.iterMicrographs():
            micId = mic.getObjId()
            micBase = pwutils.removeBaseExt(mic.getFileName())
            fnCoords = os.path.join(workDir, micBase + '_rubbish.star')
            f = open(fnCoords, 'w')
            f.write(header)
            for coord in coordSet:
                if coord.getMicId() == micId:
                    f.write("%0.6f %0.6f %0.6f %d %0.6f\n"
                            % (coord.getX(), coord.getY(),
                               coord.getAttributeValue('_rlnAnglePsi', 0.0),
                               coord.getAttributeValue('_rlnClassNumber', 0),
                               coord.getAttributeValue('_rlnAutopickFigureOfMerit', 0.0)))
            f.close()


def getEnviron():
    """ Return the environ settings to run gautomatch programs. """
    environ = pwutils.Environ(os.environ)

    # Take Scipion CUDA 7.5 library path
    cudaLib = environ.getFirst(('GAUTOMATCH_CUDA_LIB', 'CUDA_LIB'))
    environ.addLibrary(cudaLib)

    return environ


def getProgram():
    """ Return the program binary that will be used. """
    if (not 'GAUTOMATCH' in os.environ or
        not 'GAUTOMATCH_HOME' in os.environ):
        return None

    return os.path.join(os.environ['GAUTOMATCH_HOME'], 'bin',
                           os.path.basename(os.environ['GAUTOMATCH']))


def runGautomatch(micName, refStack, workDir, args, env=None):
    # We convert the input micrograph on demand if not in .mrc
    outMic = os.path.join(workDir, pwutils.replaceBaseExt(micName, 'mrc'))
    if micName.endswith('.mrc'):
        pwutils.createLink(micName, outMic)
    else:
        em.ImageHandler().convert(micName, outMic)
    if refStack is not None:
        args = ' %s --T %s %s' % (outMic, refStack, args)
    else:
        args = ' %s %s' % (outMic, args)

    pwutils.runJob(None, getProgram(), args, env=env)
    # After picking we can remove the temporary file.
    pwutils.cleanPath(outMic)
