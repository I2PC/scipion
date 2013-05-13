import sys

from pyworkflow.manager import Manager
from pyworkflow.project import *
from pyworkflow.em import *
from pyworkflow.em.packages.xmipp3 import *


#projName = sys.argv[1]
#pattern = sys.argv[2]
projName = "myproject"
pattern = "/home/laura/Scipion_Projects/InputData/*.mrc"


manager = Manager()
proj = manager.createProject(projName) # Now it will be loaded if exists

from tests.tester import *

#prot = MyProtocol(name="test3", n=2, workingDir=proj.getPath('Runs', 'T3'))
#proj.launchProtocol(prot)


prot2 = ProtImportMicrographs(pattern=pattern, samplingRate=1, voltage=300, tiltPairs=True)
proj.launchProtocol(prot2, wait=True)

#l = proj.mapper.selectByClass('SetOfMicrographs')
   

if prot2.outputMicrographs is not None:
    for (mT,mU) in prot2.outputMicrographs.iterTiltPairs():
        print 'Tilted micrograph:', mT.getFileName()
        print 'Untilted micrograph:', mU.getFileName()
        
    print "Downsampling..."
    prot3 = XmippProtPreprocessMicrographs(doDownsample=True, downFactor=3)
    prot3.inputMicrographs.set(prot2.outputMicrographs)

    proj.launchProtocol(prot3, wait=True)
    
else:
    print "Not micrographs found"
