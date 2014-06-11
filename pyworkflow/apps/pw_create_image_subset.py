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




#def createSetFromSqlite(selfile, setType, inputimages):
#    cls = getattr(my_import('pyworkflow.em.data'), 'SetOf' + setType)
#    outputset = cls(filename=selfile)
#    if not 'Classes' in setType:
#        outputset.copyInfo(inputimages)    
#    else:
#        outputset.setImages(inputimages)
#    return outputset

def my_import(name):
    mod = __import__(name)
    components = name.split('.')
    for comp in components[1:]:
        mod = getattr(mod, comp)
    return mod

if __name__ == '__main__':

    projectid = sys.argv[1]
    inputid = sys.argv[2]
    sqlitefile = sys.argv[3] # Sqlite file with disabled objects
    outputType = sys.argv[4] # Output set type
    protlabel = sys.argv[5] #Protocol label 

    inputType = None#to be readed
    
    selfilename = sqlitefile[sqlitefile.rfind(os.sep) + 1:]
    
    project = Manager().loadProject(projectid)
    prot = ProtUserSubSet(label=protlabel, inputType=inputType, outputType=outputType)
    input = project.mapper.selectById(int(inputid))

    
    prot.createInputPointer(input)
    project._setupProtocol(prot)
    prot.makePathsAndClean()
    moveFile(sqlitefile, prot._getExtraPath())
    selfile = join(prot._getExtraPath(), selfilename)
    
    readInputSetFun = getattr(prot, '_createSetOf' + inputType)
    
    inputset = readInputSetFun(filename=selfile)
    outputcls = getattr(my_import('pyworkflow.em.data'), 'SetOf' + outputType)
    outputset = outputcls()
    for obj in inputset:
        if(obj.isEnabled()):
            outputset.add(obj)
            
    if isinstance(input, SetOfImages):
        inputImages = input
    elif isinstance(input, SetOfClasses):
        inputImages = input.getImages()
    elif isinstance(input, EMProtocol):
        inputImages = input.getAttribute('inputParticles', None)
        
    if not 'Classes' in outputType:
        outputset.copyInfo(inputImages)    
    else:
        outputset.setImages(inputImages)

    prot.defineOutputSet(outputset)
   
    prot.setStatus(STATUS_FINISHED)
    project._storeProtocol(prot)
    
    
  
    

    