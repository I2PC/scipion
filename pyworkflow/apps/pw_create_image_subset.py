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
    inputType = sys.argv[2]
    setType = sys.argv[3]
    projectid = sys.argv[4]
    inputid = sys.argv[5]
    inputimagesid = sys.argv[6]
    mdname = mdfile[mdfile.rfind(os.sep) + 1:]
    
    
    project = Manager().loadProject(projectid)
    
    prot = ProtUserSubSet(inputType=inputType, outputType=setType)
    input = project.mapper.selectById(int(inputid))
    inputimages = project.mapper.selectById(int(inputimagesid))
    
    prot.createInputPointer(input)
    project._setupProtocol(prot)
    prot.makePathsAndClean()
    moveFile(mdfile, prot._getExtraPath())
    mdfile = join(prot._getExtraPath(), mdname)


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
    
    
  
    
    
    