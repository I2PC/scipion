#!/usr/bin/env python
'''
Created on Jan 27, 2014

@author: airen
'''
import os, sys
from pyworkflow.manager import Manager
from pyworkflow.em import *
from pyworkflow.em.packages.xmipp3 import *
from xmipp import *



if __name__ == '__main__':

    mdfile = sys.argv[1]
    projectid = sys.argv[2]
    protid = sys.argv[3]
    dbpath = sys.argv[4]
    
    manager = Manager()
    project = manager.loadProject(projectid)
    
    
    prot = ProtUserSelection()
    project._setupProtocol(prot)
    imgSet = prot._createSetOfParticles()
    
    readSetOfParticles(mdfile, imgSet)
#    prot._defineOutputs()

    args = {}
    outputSet = prot._getOutputSet()
    args[outputSet] = imgSet
    prot._defineOutputs(**args)
    project._storeProtocol(prot)
    print 'done'