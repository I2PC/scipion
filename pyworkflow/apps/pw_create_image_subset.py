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



if __name__ == '__main__':

    mdfile = sys.argv[1]
    setType = sys.argv[2]
    projectid = sys.argv[3]
    inputid = sys.argv[4]
    inputimagesid = sys.argv[5]
    
    
    project = Manager().loadProject(projectid)
    
    prot = ProtUserSubSet(setType=setType)
    inputset = project.mapper.selectById(int(inputid))
    inputimages = project.mapper.selectById(int(inputimagesid))
    
    prot.createInputSet(inputset)
    project._setupProtocol(prot)
    prot.makePathsAndClean()
    moveFile(mdfile, prot._getExtraPath())
    mdfile = join(prot._getExtraPath(), "selection.xmd")
    createSetFun = getattr(prot, '_createSetOf' + setType)
    if setType != 'Classes2D':
        outputset = createSetFun()
        readSetFun = getattr(xmipp3, 'readSetOf' + setType )
        outputset.copyInfo(inputimages)    
        readSetFun(mdfile, outputset, outputset.hasCTF())
    else:
        outputset = createSetFun(inputimages)
        readSetFun = getattr(xmipp3, 'readSetOf' + setType )       
        readSetFun(outputset, mdfile, inputImages=inputimages)

        
    prot.createOutputSet(outputset)
   
    prot.setStatus(STATUS_FINISHED)
    project._storeProtocol(prot)
    
    
  
    
    
    