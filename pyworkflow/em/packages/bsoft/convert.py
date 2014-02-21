'''
Created on Feb 21, 2014

@author: airen
'''
import os
import xmipp
from pyworkflow.utils.path import join, dirname, replaceBaseExt
from protlib_import import addBsoftLabelAliases
from pyworkflow.em import *
from pyworkflow.em.packages.xmipp3.convert import rowToCoordinate


def readSetOfCoordinates(outputDir, micSet, coordSet):
    """ Read from Bsoft .star files.
    Params:
        outputDir: the directory where the .star files are.
           
        micSet: the SetOfMicrographs to associate the .star, which 
            name should be the same of the micrographs.
        coordSet: the SetOfCoordinates that will be populated.
    """

    addBsoftLabelAliases()
    for mic in micSet:
        outputFile = join(outputDir, replaceBaseExt(mic.getFileName(), 'star'))
        scipionPosFile = join(outputDir, "scipion_" + replaceBaseExt(mic.getFileName(), 'pos'))
        if exists(outputFile):
            posMd = xmipp.MetaData(outputFile)
        else:
            posMd = xmipp.MetaData()
        
        for objId in posMd:
            coord = rowToCoordinate(posMd, objId)
            coord.setMicrograph(mic)
            coordSet.append(coord)      
            # Add an unique ID that will be propagated to particles
            posMd.setValue(xmipp.MDL_ITEM_ID, long(coord.getObjId()), objId)
        if not posMd.isEmpty():
            boxSize = posMd.getValue(xmipp.MDL_PICKING_PARTICLE_SIZE, posMd.firstObject())
            posMd.write("particles@%s"  % scipionPosFile)
    
    coordSet.setBoxSize(boxSize)