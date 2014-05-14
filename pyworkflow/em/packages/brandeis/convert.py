# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              Laura del Cano (ldelcano@cnb.csic.es)
# *              Josue Gomez Blanco (ldelcano@cnb.csic.es)
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

import subprocess
from data import *
from brandeis import *
import glob
from pyworkflow.em.constants import NO_INDEX
from os.path import abspath



def objectToRow(obj, row, attrDict):
    pass


def _rowToObject(row, obj, attrDict):
    pass

    
def rowToObject(md, objId, obj, attrDict):
    pass

    
def readCTFModel(filename):
    pass


def writeCTFModel(ctfObj, filename):
    """ Write a CTFModel object as Xmipp .ctfparam"""
    pass


def readMicrograph(md, objId):
    """ Create a Micrograph object from a row of Xmipp metadata. """
    pass


def locationToEman(index, filename):
    pass


def micrographToRow(mic, micRow):
    pass


def rowToCoordinate(md, objId):
    """ Create a Coordinate from a json. """
    pass


def readSetOfMicrographs(filename):
    pass


def writeSetOfMicrographs(micSet, filename, rowFunc=None):
    pass

    
def readPosCoordinates(posFile):
    pass

            
def readSetOfCoordinates(workDir, micSet, coordSet):
    pass


def writeSetOfCoordinates():
    pass


def readSetOfParticles(filename):
    pass


def writeSetOfParticles(partSet, filename, inputdir):
    pass


def writeSetOfClasses3D(classes3DSet, filename, classesBlock='classes'):    
    pass


def readSetOfClasses3D(classes3DSet, fileparList, volumeList):
    """read from frealign .par.
    """
    import re
    imgSet = classes3DSet.getImages()
    samplingRate = imgSet.getSamplingRate()
    averages = classes3DSet.createRepresentatives()
    
    for ref, volFn in enumerate(volumeList):
        class3D = Class3D()
        class3D.setObjId(ref+1)
        vol = volume()
        vol.copyObjId(class3D)
        vol.setLocation(volFn)
        
        class3D.setRepresentative(vol)
        averages.append(vol)
        classes3DSet.append(class3D)
        
        file1 = fileparList[ref]
        f1 = open(file1)
        for l in f1:
            if not l.startswith('C'):
                values = re.split(" +", l)
                prob = float(values[11])
                if prob > 0:
                    img = Particle()
                    particle = imgSet[int(values[7])]
                    img.copy(particle)
                    class3D.append(img)
        f1.close()
        
        # Check if write function is necessary
        class3D.write()

