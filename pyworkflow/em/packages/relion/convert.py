# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              Laura del Cano (ldelcano@cnb.csic.es)
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
1. Write from base classes to Relion specific files
2. Read from Relion files to base classes
"""

import os
from pyworkflow.object import String
from pyworkflow.utils.path import createLink
from pyworkflow.em import ImageHandler
from pyworkflow.em.data import Class2D, SetOfClasses2D
from constants import *
# Since Relion share several conventions with Xmipp, we will reuse 
# the xmipp3 tools implemented for Scipion here
import xmipp



def locationToRelion(index, filename):
    """ Convert an index and filename location
    to a string with @ as expected in Relion.
    """
    from pyworkflow.em.packages.xmipp3.convert import locationToXmipp
    return locationToXmipp(index, filename)


def relionToLocation(filename):
    from pyworkflow.em.packages.xmipp3.convert import xmippToLocation
    return xmippToLocation(filename)


_xmippLabelsDict = {} # Dictionary to store mappings replaced

            
def restoreXmippLabels():
    global _xmippLabelsDict
    for k, v in _xmippLabelsDict.iteritems():
        xmipp.addLabelAlias(k, v, True)
    _xmippLabelsDict = {}

            
class ParticleAdaptor():
    """ Class used to convert a set of particles for Relion.
    It will write an stack in Spider format and also
    modify the output star file to point to the new stack.
    """
    def __init__(self, stackFile):
        self._rowCount = 1
        self._ih = ImageHandler()
        self._stackFile = stackFile
        
        import pyworkflow.em.packages.xmipp3 as xmipp3
        self._particleToRow = xmipp3.particleToRow
        
    def setupRow(self, img, imgRow):
        """ Convert image and modify the row. """
        newLoc = (self._rowCount, self._stackFile)
        self._ih.convert(img.getLocation(), newLoc)
        img.setLocation(newLoc)
        # Re-write the row with the new location
        self._particleToRow(img, imgRow)
        self._rowCount += 1
    

def addRelionLabels(replace=False, extended=False):
    """ Add relion labels as aliases for Xmipp metadata. """
    global _xmippLabelsDict
    _xmippLabelsDict = {}
    for k, v in XMIPP_RELION_LABELS.iteritems():
        _xmippLabelsDict[k] = xmipp.label2Str(k) # store original label string
        xmipp.addLabelAlias(k, v, replace)
    if extended:
        for k, v in XMIPP_RELION_LABELS_EXTRA.iteritems():    
            _xmippLabelsDict[k] = xmipp.label2Str(k) # store original label string
            xmipp.addLabelAlias(k, v, replace)

            
def addRelionLabelsToEnviron(env):
    """ create an string that can be used for XMIPP_EXTRA_ALIASES
    for adding the labels of Relion.
    """
    pairs = []
    for k, v in XMIPP_RELION_LABELS.iteritems():
        pairs.append('%s=%s' % (xmipp.label2Str(k), v))
    for k, v in XMIPP_RELION_LABELS_EXTRA.iteritems():
        pairs.append('%s=%s' % (xmipp.label2Str(k), v))        
    varStr = ';'.join(pairs)
    env['XMIPP_EXTRA_ALIASES'] = varStr
    
    
def writeSetOfParticles(imgSet, starFile, stackFile):
    """ This function will write a SetOfImages as Relion metadata.
    Params:
        imgSet: the SetOfImages instance.
        filename: the filename where to write the metadata.
    """
    import pyworkflow.em.packages.xmipp3 as xmipp3
    addRelionLabels(replace=True)
    pa = ParticleAdaptor(stackFile)
    xmipp3.writeSetOfParticles(imgSet, starFile, rowFunc=pa.setupRow)
    imgSet._relionStar = String(starFile)
    restoreXmippLabels()
    
    
def createRelionInputParticles(imgSet, starFile, stackFile): 
    """ Ensure that in 'filename' it is a valid STAR files with particles.
    If the imgSet comes from Relion, just create a link.
    If not, then write the proper file.
    """
    imgsStar = getattr(imgSet, '_relionStar', None)
    if imgsStar is None:
        writeSetOfParticles(imgSet, starFile, stackFile)
    else:
        imgsFn = imgsStar.get()
        createLink(imgsFn, imgsStar.get())


def createClassesFromImages(inputImages, inputStar, classesFn, classLabel, classFnTemplate, iter):
    """ From an intermediate dataXXX.star file produced by relion, create
    the set of classes in which those images are classified.
    Params:
        inputImages: the SetOfImages that were initially classified by relion. 
        inputStar: the filename with relion star format.
        classesFn: filename where to write the classes.
        classLabel: label that is the class reference.
        classFnTemplate: the template to get the classes averages filenames
        iter: the iteration number, just used in Class template
    """
    import pyworkflow.em.packages.xmipp3 as xmipp3
    addRelionLabels(replace=True, extended=True)
    # We asume here that the volumes (classes3d) are in the same folder than imgsFn
    # rootDir here is defined to be used expanding locals()
    rootDir = os.path.dirname(inputStar)
    
    md = xmipp.MetaData(inputStar)
    clsDict = {} # Dictionary to store the (classId, classSet) pairs
    clsSet = SetOfClasses2D(filename=classesFn)
    clsSet.setImages(inputImages)
    
    for objId in md:
        ref = md.getValue(classLabel, objId)
        if not ref in clsDict: # Register a new class set if the ref was not found.
            cls = Class2D(objId=ref)
            refLocation = relionToLocation(classFnTemplate % locals())
            rep = cls.ITEM_TYPE()
            rep.setLocation(refLocation)
            cls.setRepresentative(rep)
            
            clsDict[ref] = cls
            clsSet.append(cls)
        class2D = clsDict[ref] # Try to get the class set given its ref number
        # Set images attributes from the md row values
        img = xmipp3.rowToParticle(md, objId, hasCtf=False)
        class2D.append(img)
        
    for class2D in clsDict.values():
        clsSet.update(class2D)
        
    clsSet.write()
    restoreXmippLabels()    
    
#def createClassesFromImages(inputImages, inputStar, classesFn, ):
    

def readSetOfClasses2D(classes2DSet, filename, classesBlock='classes', **args):
    """ Fill a set of classes 2d from a Relion star file.
    Params:
        classesStar: the set of classes that will be populated.
        filename: the path to the input Relion star file.
        classesBlock: the name of the block in the star file.
    """
    clsSet = SetOfClasses2D(filename=filename)
    for cls in clsSet:
        newCls = Class2D()
        newCls.copyInfo(cls)
        classes2DSet.append(newCls)
        for img in cls:
            newCls.append(img)
        classes2DSet.update(newCls)

def readSetOfClasses3D(classes3DSet, filename, classesBlock='classes', **args):
    """ Fill a set of classes 3d from a Relion star file.
    Params:
        classesStar: the set of classes that will be populated.
        filename: the path to the input Relion star file.
        classesBlock: the name of the block in the star file.
    """
    import pyworkflow.em.packages.xmipp3 as xmipp3
    addRelionLabels(replace=True, extended=True)
    xmipp3.readSetOfClasses3D(classes3DSet, filename, classesBlock)
    restoreXmippLabels()


def sortImagesByLL(imgStar, imgSortedStar):
    """ Given a Relion images star file, sort by LogLikelihood 
    and save a new file. 
    """
    addRelionLabels(replace=True, extended=True)
    print "Sorting particles by likelihood, input: ", imgStar
    fn = 'images@'+ imgStar
    md = xmipp.MetaData(fn)
    md.sort(xmipp.MDL_LL, False)
    md.write('images_sorted@' + imgSortedStar)   
    restoreXmippLabels()
    
