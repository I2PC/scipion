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
1. Write from base classes to Eman specific files
2. Read from Eman files to base classes
"""

import subprocess
from data import *
from eman2 import *

from pyworkflow.em.constants import NO_INDEX
from os.path import abspath

# LABEL_TYPES = { 
#                xmipp.LABEL_SIZET: long,
#                xmipp.LABEL_DOUBLE: float,
#                xmipp.LABEL_INT: int,
#                xmipp.LABEL_BOOL: bool              
#                }

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
    """ Read from Eman .json files.
    It is expected a file named: base.json under the workDir.
    Params:
        workDir: where the Eman boxer output files are located.
        micSet: the SetOfMicrographs to associate the .json, which 
            name should be the same of the micrographs.
        coordSet: the SetOfCoordinates that will be populated.
    """
    # Read the boxSize from the e2boxercache/base.json
    jsonFnbase = join(workDir, 'e2boxercache', 'base.json')
    jsonBoxDict = loadJson(jsonFnbase)
    size = int(jsonBoxDict["box_size"])
    
    for mic in micSet:
        micFnroot = removeBaseExt(mic.getFileName()) + '_info.json'
        micPosRelFn = join("info", micFnroot)
        micPosFn = join(workDir, micPosRelFn)
        
        if exists(micPosFn):
            jsonPosDict = loadJson(micPosFn)
            boxes = jsonPosDict["boxes"]

            for box in boxes:
                x, y = box[:2] 
                coord = Coordinate()
                coord.setPosition(x, y)
                coord.setMicrograph(mic)
                coordSet.append(coord)
    coordSet.setBoxSize(size)

def writeSetOfCoordinates():
    pass

def readSetOfParticles(filename):
    pass

def writeSetOfParticles(partSet, filename):
    """ Convert the imgSet particles to a single .hdf file as expected by Eman. 
    This function should be called from a current dir where
    the images in the set are available.
    """
    cwd = os.getcwd()
    loadEnvironment()
    program = pw.join('em', 'packages', 'eman2', 'e2converter.py')        
    cmd = getEmanCommand(program, filename)
    
#    gcmd = greenStr(cmd)
    print "** Running: '%s'" % cmd
    proc = subprocess.Popen(cmd, shell=True, 
                            stdin=subprocess.PIPE,
                            stdout=subprocess.PIPE)
    
    for part in partSet:
        index, fn = part.getLocation()
        # Write the e2converter.py process from where to read the image
        print "sending: ", part.getId(), index, fn
        print >> proc.stdin, part.getId(), index, join(cwd, fn)
        proc.stdin.flush()
        response = proc.stdout.readline()
        print "response: ", response
    #proc.wait()
