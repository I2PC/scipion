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
    type = sys.argv[2]
    projectid = sys.argv[3]
    inputid = sys.argv[4]

    
    
    project = Manager().loadProject(projectid)
    
    prot = ProtUserSubSet(setType=type)
    inputset = project.mapper.selectById(int(inputid))
    prot.createInputSet(inputset)
    
    project._setupProtocol(prot)
    prot.makePathsAndClean()
    moveFile(mdfile, prot._getExtraPath())
    mdfile = join(prot._getExtraPath(), "selection.xmd")
    createSetFun = getattr(prot, '_createSetOf' + type)
    if isinstance(inputset, SetOfClasses2D):
            inputset = inputset.getImages()
    if type != 'Classes2D':
         # For inputset of classes, we need to copy properties from
         # the images used in the classification.
        
        outputset = createSetFun()
        readSetFun = getattr(xmipp3, 'readSetOf' + type )
        outputset.copyInfo(inputset)    
        readSetFun(mdfile, outputset, outputset.hasCTF())
    else:
        outputset = createSetFun(inputset)
        readSetFun = getattr(xmipp3, 'readSetOf' + type )       
        readSetFun(outputset, mdfile)

        
    prot.createOutputSet(outputset)
   
    prot.setStatus(STATUS_FINISHED)
    project._storeProtocol(prot)
    
    print "%s Subset created...\n(You may need to refresh the main window to visualize output)"%(type)
  
    
    
    