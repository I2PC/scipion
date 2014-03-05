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
    createSetFun = getattr(prot, '_createSetOf' + type)
    outputset = createSetFun()
    readSetFun = getattr(xmipp3, 'readSetOf' + type )
    
    # For inputset of classes, we need to copy properties from
    # the images used in the classification.
    if isinstance(inputset, SetOfClasses2D):
        inputset = inputset.getImages()
    outputset.copyInfo(inputset)    
    readSetFun(mdfile, outputset, outputset.hasCTF())
    prot.createOutputSet(outputset)
   
    prot.setStatus(STATUS_FINISHED)
    project._storeProtocol(prot)
    
    
    