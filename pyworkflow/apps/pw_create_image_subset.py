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

    protlabel = sys.argv[1]
    mdfile = sys.argv[2]
    inputType = sys.argv[3]
    setType = sys.argv[4]
    projectid = sys.argv[5]
    inputid = sys.argv[6]
    inputimagesid = sys.argv[7]
    mdname = mdfile[mdfile.rfind(os.sep) + 1:]
    
    
    project = Manager().loadProject(projectid)
    
    prot = ProtUserSubSet(label=protlabel, inputType=inputType, outputType=setType)
    input = project.mapper.selectById(int(inputid))
    inputimages = project.mapper.selectById(int(inputimagesid))
    
    prot.createInputPointer(input)
    project._setupProtocol(prot)
    prot.makePathsAndClean()
    moveFile(mdfile, prot._getExtraPath())
    mdfile = join(prot._getExtraPath(), mdname)


    createSetFun = getattr(prot, '_createSetOf' + setType)
    if not 'Classes' in setType:
        outputset = createSetFun()
        readSetFun = getattr(xmipp3, 'readSetOf' + setType )
        outputset.copyInfo(inputimages)    
        readSetFun(mdfile, outputset, outputset.hasCTF())
    else:
        outputset = createSetFun(inputimages)
        readSetFun = getattr(xmipp3, 'readSetOf' + setType )       
        readSetFun(outputset, mdfile)

        
    prot.createOutputSet(outputset)
   
    prot.setStatus(STATUS_FINISHED)
    project._storeProtocol(prot)
    
    
  
    
    
    