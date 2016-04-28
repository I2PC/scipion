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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************
"""
This module contains converter functions that will serve to:
1. Write from base classes to Gempicker specific files
2. Read from Gempicker files to base classes
"""

import os
from os.path import join

import pyworkflow.em as em
import pyworkflow.em.metadata as md
import pyworkflow.utils as pwutils



def readSetOfCoordinates(workingDir, micSet, coordSet):
    """ Read from coordinates from Gempicker .box files.
    For a micrograph: mic1.mrc, the expected coordinate file is: mic1.box
    Params:
        workingDir: where the Gempicker .box output files are located.
        micSet: the input SetOfMicrographs
        coordSet: the SetOfCoordinates that will be populated.
    """
    for mic in micSet:
        micBase = pwutils.removeBaseExt(mic.getFileName())
        fnCoords = join(workingDir, 'pik_box', micBase + '.box')
        readCoordinates(mic, fnCoords, coordSet)


def readCoordinates(mic, fileName, coordsSet):
    if os.path.exists(fileName):
        mdCoords = md.MetaData()
        mdCoords.readPlain(fileName, "xcoor ycoor xSize ySize")
        for row in md.iterRows(mdCoords):
            coord = em.Coordinate(x = row.getValue("xcoor") + row.getValue("xSize")/2,
                                  y = row.getValue("ycoor") + row.getValue("ySize")/2)
            coord.setMicrograph(mic)
            coordsSet.append(coord)


def getProgram(useGPU):
    """ Return the program binary that will be used. """
    if useGPU:
        binary = os.environ['GEMPICKER_CUDA']
    else:
        binary = os.environ['GEMPICKER']

    program = join(os.environ['GEMPICKER_HOME'], os.path.basename(binary))
    return program


def runGempicker(micName, workingDir, useGPU, args):
    # We convert the input micrograph on demand if not in .mrc
    outMic = os.path.join(workingDir, pwutils.replaceBaseExt(micName, 'mrc'))
    if micName.endswith('.mrc'):
        pwutils.createLink(micName, outMic)
    else:
        em.ImageHandler().convert(micName, outMic)

    refDir = join(workingDir, 'templates')
    maskSchDir = join(workingDir, 'maskSch')
    args += ' --dirTgt=%s --dirSch=%s --dirMskSch=%s ' % (workingDir,
                                                          refDir, maskSchDir)
    # Run Gempicker:
    for mode in [0, 1]:
        pwutils.runJob(None, getProgram(useGPU), args + ' --mode=%d' % mode)
    # After picking we can remove the temporary file.
    pwutils.cleanPath(outMic)
