
'''
Created on Jan 27, 2014

@author: airen
'''
import os, sys
from pyworkflow.manager import Manager

from pyworkflow.em import *
import pyworkflow.em.packages.xmipp3 as xmipp3
from xmipp import *
from pyworkflow.utils.path import moveTree
from pyworkflow.em.packages.xmipp3 import readSetOfCoordinates



if __name__ == '__main__':

    protlabel = sys.argv[1]
    outputdir = sys.argv[2]
    type = 'Coordinates'
    projectid = sys.argv[3]
    inputid = sys.argv[4]
    protstate = sys.argv[5]
    

    
    
    project = Manager().loadProject(projectid)
    
    prot = ProtUserSubSet(label=protlabel, inputType=type, outputType=type)
    
    inputset = project.mapper.selectById(int(inputid))
    prot.createInputPointer(inputset)
    
    project._setupProtocol(prot)
    prot.makePathsAndClean()
    
    extradir = prot._getExtraPath()
    moveTree(outputdir, extradir)
    outputset = prot._createSetOfCoordinates(inputset.getMicrographs())
    readSetOfCoordinates(extradir, outputset.getMicrographs(), outputset)
    

    prot.createOutputSet(outputset)
    prot.setStatus(STATUS_FINISHED)
    project._storeProtocol(prot)
    