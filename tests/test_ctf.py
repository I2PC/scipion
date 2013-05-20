import sys

from pyworkflow.manager import Manager
from pyworkflow.project import *
from pyworkflow.em import *
from pyworkflow.em.packages.xmipp3 import *


projName = sys.argv[1]
pattern = sys.argv[2]

manager = Manager()
proj = manager.createProject(projName) # Now it will be loaded if exists

from tests.tester import *

protImport = ProtImportMicrographs(pattern=pattern, samplingRate=1, voltage=200, tiltPairs=True)
proj.launchProtocol(protImport, wait=True)

l = proj.mapper.selectByClass('SetOfMicrographs')

for p in l:
    p.printAll()
    

if len(l):
    print "Continue after import..."   
    protCTF = XmippProtCTFMicrographs()
        
    #prot4.inputMicrographs.set(l[0])
    protCTF.inputMicrographs.set(protImport.outputMicrographs)

    proj.launchProtocol(protCTF, wait=True)

    l = proj.mapper.selectByClass('SetOfMicrographs')

    for p in l:
        if p.hasCTF():
            for mic in p:
                if mic.hasCTF():
                    mic.ctfModel.printAll() 
    
else:
    print "Not micrographs found"
