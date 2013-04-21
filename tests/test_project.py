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

#prot = MyProtocol(name="test3", n=2, workingDir=proj.getPath('Runs', 'T3'))
#proj.launchProtocol(prot)


prot2 = ProtImportMicrographs(pattern=pattern, samplingRate=1, voltage=200)
proj.launchProtocol(prot2, wait=True)

l = proj.mapper.selectByClass('SetOfMicrographs')

sys.exit(0)

for p in l:
    p.printAll()
    

if len(l):
    print "Continue after import..."
    prot3 = XmippProtPreprocessMicrographs(workingDir=proj.getPath('Runs', 'PreprocessDown'), 
                                       doDownsample=True, downFactor=3, doCrop=True)
    prot3.inputMicrographs.set(prot2.outputMicrographs)

    proj.launchProtocol(prot3)
    
    prot4 = XmippProtCTFMicrographs(workingDir=proj.getPath('Runs', 'Ctf'))
        
    #prot4.inputMicrographs.set(l[0])
    prot4.inputMicrographs.set(prot3.outputMicrographs)

    proj.launchProtocol(prot4)

    l = proj.mapper.selectByClass('SetOfMicrographs')

    for p in l:
        if p.hasCTF():
            for mic in p:
                if mic.hasCTF():
                    mic.ctfModel.printAll() 
    
else:
    print "Not micrographs found"
