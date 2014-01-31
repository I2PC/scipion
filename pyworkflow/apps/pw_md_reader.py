#!/usr/bin/env python
'''
Created on Jan 27, 2014

@author: airen
'''
import os, sys
from pyworkflow.project import Project

from pyworkflow.em import *
from pyworkflow.em.packages.xmipp3 import *
from xmipp import *



if __name__ == '__main__':

    mdfile = sys.argv[1]
    inputid = sys.argv[2]

    
    
    project = Project(".")
    project.load()
    prot = ProtParticlesSubset()
    
    particles = project.mapper.selectById(int(inputid))
    prot.inputParticles = Pointer()
    prot.inputParticles.set(particles)
    
    prot.outputParticles = prot._createSetOfParticles() 
    readSetOfParticles(mdfile, prot.outputParticles)
    
    
    project._setupProtocol(prot)
    prot.makePathsAndClean()
    prot.setStatus(STATUS_FINISHED)
    project._storeProtocol(prot)
    
    
    print 'done'