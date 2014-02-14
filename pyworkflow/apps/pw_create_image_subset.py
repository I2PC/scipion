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
    
    prot = ProtUserSubSet(type)
    inputset = project.mapper.selectById(int(inputid))
    prot.createInputSet(inputset)
    
    
    createSetFun = getattr(prot, '_createSetOf' + type)
    outputset = createSetFun()
    readSetFun = getattr(xmipp3, 'readSetOf' + type )
    
    readSetFun(mdfile, outputset)
    createSetFun = getattr(prot, '_createSetOf' + type)
    outputset.setSamplingRate(inputset.getSamplingRate())
    project._setupProtocol(prot)
    prot.makePathsAndClean()
    prot.createOutputSet(outputset)
   
    prot.setStatus(STATUS_FINISHED)
    project._storeProtocol(prot)
    
    
    