
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
    dbpath = sys.argv[2]
    protid = sys.argv[3]
    print outputdir, dbpath, protid
    prot = getProtocolFromDb(dbpath, protid)
    prot.printAll()
    inputset = prot.inputMicrographs.get()
    extradir = prot._getExtraPath()
    count = 0;
    for key, output in prot.iterOutputAttributes(EMObject):
        count += 1
    
    suffix = str(count + 1) if count > 0 else ''
    outputName = 'outputCoordinates' + suffix
    outputset = prot._createSetOfCoordinates(inputset, suffix=suffix)#micrographs are the input set if protocol is not finished
    readSetOfCoordinates(extradir, outputset.getMicrographs(), outputset)
    outputs = {outputName: outputset}
    prot._defineOutputs(**outputs)
    prot._defineSourceRelation(inputset, outputset)
    prot._store()
        
    