#!/usr/bin/env python
'''
Created on Jan 27, 2014

@author: airen
'''
import os, sys
from pyworkflow.project import Project

from pyworkflow.em import *
import pyworkflow.em.packages.xmipp3 as xmipp3
from xmipp import *



if __name__ == '__main__':

    mdfile = sys.argv[1]
    set = sys.argv[2]
    inputid = sys.argv[3]

    
    
    project = Project(".")
    project.load()
    prot = ProtUserSubSet()
    
    inputset = project.mapper.selectById(int(inputid))
    prot.inputset = Pointer()
    prot.inputset.set(inputset)
    
    createFun = getattr(prot, '_create' + set)
    outputset = createFun() #prot._createSetOfParticles() 
    readFun = getattr(xmipp3, 'read' + set)
    readFun(mdfile, outputset)
    #readSetOfParticles(mdfile, outputParticles)
    #this value is not always read it on readSetOfParticles
    outputset.setSamplingRate(inputset.getSamplingRate())
    project._setupProtocol(prot)
    prot.makePathsAndClean()

    prot._defineOutputs(outputset=outputset)
    prot.setStatus(STATUS_FINISHED)
    project._storeProtocol(prot)
    
    
    print 'done'