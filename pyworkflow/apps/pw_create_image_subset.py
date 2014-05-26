#!/usr/bin/env python
'''
Created on Jan 27, 2014

@author: airen
'''
import os, sys
from pyworkflow.manager import Manager

from pyworkflow.em import *
from pyworkflow.utils import moveFile

import pyworkflow.em.packages.xmipp3 as xmipp3
from xmipp import *

    
def createSetFromMd(selfile, setType, inputimages):
    createSetFun = getattr(prot, '_createSetOf' + setType)
        
    if not 'Classes' in setType:
        outputset = createSetFun()
        readSetFun = getattr(xmipp3, 'readSetOf' + setType )
        outputset.copyInfo(inputimages)    
        readSetFun(selfile, outputset, outputset.hasCTF())
    else:
        outputset = createSetFun(inputimages)
        readSetFun = getattr(xmipp3, 'readSetOf' + setType )       
        readSetFun(outputset, selfile)
    return outputset


def createSetFromSqlite(selfile, setType, inputimages):
    cls = getattr(my_import('pyworkflow.em.data'), 'SetOf' + setType)
    outputset = cls(filename=selfile)
    if not 'Classes' in setType:
        outputset.copyInfo(inputimages)    
    else:
        outputset.setImages(inputimages)
    return outputset


def my_import(name):
    mod = __import__(name)
    components = name.split('.')
    for comp in components[1:]:
        mod = getattr(mod, comp)
    return mod

if __name__ == '__main__':

    protlabel = sys.argv[1]
    selfile = sys.argv[2]
    inputType = sys.argv[3]
    setType = sys.argv[4]
    projectid = sys.argv[5]
    inputid = sys.argv[6]
    inputimagesid = sys.argv[7]
    mdname = selfile[selfile.rfind(os.sep) + 1:]
    
    project = Manager().loadProject(projectid)
    prot = ProtUserSubSet(label=protlabel, inputType=inputType, outputType=setType)
    input = project.mapper.selectById(int(inputid))
    inputimages = project.mapper.selectById(int(inputimagesid))
    
    prot.createInputPointer(input)
    project._setupProtocol(prot)
    prot.makePathsAndClean()
    moveFile(selfile, prot._getExtraPath())
    selfile = join(prot._getExtraPath(), mdname)
    outputset = createSetFromSqlite(selfile, setType, inputimages) if selfile.endswith('.sqlite') else createSetFromMd(selfile, setType, inputimages)
    prot.createOutputSet(outputset)
   
    prot.setStatus(STATUS_FINISHED)
    project._storeProtocol(prot)
    
    
  
    

    