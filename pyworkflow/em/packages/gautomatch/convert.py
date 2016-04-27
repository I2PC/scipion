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
1. Write from base classes to Eman specific files
2. Read from Eman files to base classes
"""

import os

import pyworkflow.em as em
import pyworkflow.utils as pwutils



def readSetOfCoordinates(workDir, micSet, coordSet):
    """ Read from coordinates from Gautomatch .box files.
    For a micrograph: mic1.mrc, the expected coordinate file is: mic1.box
    It is expected a file named: base.json under the workDir.
    Params:
        workDir: where the Gautomatch boxer output files are located.
        micSet: the SetOfMicrographs to associate the .json, which 
            name should be the same of the micrographs.
        coordSet: the SetOfCoordinates that will be populated.
    """
    for mic in micSet:
        micBase = pwutils.removeBaseExt(mic.getFileName())
        fnCoords = os.path.join(workDir, micBase + '_automatch.box')
        readCoordinates(mic, fnCoords, coordSet)


def readCoordinates(mic, fileName, coordsSet):
    with open(fileName, "r") as source:
        for line in source:
            tokens = map(int, line.split())
            coord = em.Coordinate(x=tokens[0] + tokens[2]/2,
                                  y=tokens[1] + tokens[3]/2)
            coord.setMicrograph(mic)
            coordsSet.append(coord)
