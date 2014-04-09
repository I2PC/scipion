
'''
Created on Jan 27, 2014

@author: airen
'''
import os, sys
from pyworkflow.manager import Manager

from pyworkflow.em import *
import pyworkflow.em.packages.xmipp3 as xmipp3
from xmipp import *
from pyworkflow.utils.path import moveTree
from pyworkflow.em.packages.xmipp3 import readSetOfCoordinates



if __name__ == '__main__':

    
    outputdir = sys.argv[1]
    type = 'Coordinates'
    projectid = sys.argv[2]
    inputid = sys.argv[3]
    protid = sys.argv[4]
    
    
    project = Manager().loadProject(projectid)
    prot = project.mapper.selectById(protid)
    inputset = project.mapper.selectById(int(inputid))
    extradir = prot._getExtraPath()
    
    if(prot is None):
        print 'prot is finished'
        protlabel = protid
        prot = ProtUserSubSet(label=protlabel, inputType=type, outputType=type)        
        prot.createInputPointer(inputset)
        project._setupProtocol(prot)
        prot.makePathsAndClean()
        extradir = prot._getExtraPath()
        moveTree(outputdir, extradir)
        prot.setStatus(STATUS_FINISHED)
        outputName = 'output' + type
        outputset = prot._createSetOfCoordinates(inputset.getMicrographs())
    else:
        print 'prot output will be updated'
        count = 0;
        for key, output in prot.iterOutputAttributes(EMObject):
            print key
            count += 1
        suffix = str(count) if count > 1 else ''
        outputName = 'output' + type + suffix
        outputset = prot._createSetOfCoordinates(inputset)#micrographs are the input set if protocol is not finished
    print outputset.printAll()
    readSetOfCoordinates(extradir, outputset.getMicrographs(), outputset)
    outputs = {outputName: outputset}
    prot._defineOutputs(**outputs)
    prot._defineSourceRelation(inputset, outputset)
    project._storeProtocol(prot)
    