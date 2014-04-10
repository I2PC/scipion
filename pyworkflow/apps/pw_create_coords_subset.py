
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
    updateprot = sys.argv[2].endswith('.db')
    if updateprot:
        dbpath = sys.argv[2]
        protid = sys.argv[3]
        prot = getProtocolFromDb(dbpath, protid)
        inputset = prot.inputMicrographs.get()
        extradir = prot._getExtraPath()
        count = 0;
        for key, output in prot.iterOutputAttributes(EMObject):
            print key
            count += 1
        print 'count %s'%(count)
        suffix = str(count + 1) if count > 0 else ''
        outputName = 'output' + type + suffix
        outputset = prot._createSetOfCoordinates(inputset, suffix=suffix)#micrographs are the input set if protocol is not finished
    else:
        projectid = sys.argv[2]
        inputid = sys.argv[3]
        protlabel = sys.argv[4]
        project = Manager().loadProject(projectid)
        inputset = project.mapper.selectById(int(inputid))
        
        prot = ProtUserSubSet(label=protlabel, inputType=type, outputType=type)        
        prot.createInputPointer(inputset)
        project._setupProtocol(prot)
        prot.makePathsAndClean()
        extradir = prot._getExtraPath()
        moveTree(outputdir, extradir)
        prot.setStatus(STATUS_FINISHED)
        outputName = 'output' + type
        outputset = prot._createSetOfCoordinates(inputset.getMicrographs())
    
    
    readSetOfCoordinates(extradir, outputset.getMicrographs(), outputset)
    outputs = {outputName: outputset}
    prot._defineOutputs(**outputs)
    prot._defineSourceRelation(inputset, outputset)
    if not updateprot:
        project._storeProtocol(prot)
    else:
        prot._store()
        
    