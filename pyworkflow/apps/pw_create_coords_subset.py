#!/usr/bin/env python
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

    
    outputdir = sys.argv[1]
    type = 'Coordinates'
    projectid = sys.argv[2]
    inputid = sys.argv[3]

    
    
    project = Manager().loadProject(projectid)
    
    prot = ProtUserSubSet(setType=type)
    
    inputset = project.mapper.selectById(int(inputid))
    prot.createInputSet(inputset)
    
    project._setupProtocol(prot)
    prot.makePathsAndClean()
    
    extradir = prot._getExtraPath()
    moveTree(outputdir, extradir)
    outputset = prot._createSetOfCoordinates()
    outputset.setMicrographs(inputset.getMicrographs())
    readSetOfCoordinates(extradir, outputset.getMicrographs(), outputset)
    

    prot.createOutputSet(outputset)
    prot.setStatus(STATUS_FINISHED)
    project._storeProtocol(prot)
    