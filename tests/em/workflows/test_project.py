import sys

from pyworkflow.manager import Manager
from pyworkflow.project import *
from pyworkflow.em import *
from pyworkflow.em.packages.xmipp3 import *
from pyworkflow.em.packages.brandeis.protocol_ctffind3 import ProtCTFFind


projName = sys.argv[1]
pattern = sys.argv[2]

manager = Manager()
proj = manager.createProject(projName) # Now it will be loaded if exists

from tests.tester import *

#prot = MyProtocol(name="test3", n=2, workingDir=proj.getPath('Runs', 'T3'))
#proj.launchProtocol(prot)

#pattern = "/home/josem/work/data/InputMicrographs/*6.mrc"

prot2 = ProtImportMicrographs(pattern=pattern, samplingRate=1, voltage=300)
proj.launchProtocol(prot2, wait=True)

l = proj.mapper.selectByClass('SetOfMicrographs')
#
#for p in l:
#    p.printAll()
    

if len(l):
    print "Continue after import..."
    prot3 = XmippProtPreprocessMicrographs(doDownsample=True, downFactor=3, doCrop=False)
    prot3.inputMicrographs.set(prot2.outputMicrographs)
    proj.launchProtocol(prot3, wait=True)
    
    prot4 = XmippProtCTFMicrographs()
    prot4.inputMicrographs.set(prot3.outputMicrographs)
    proj.launchProtocol(prot4, wait=True)
    
    prot5 = ProtCTFFind()
    prot5.inputMicrographs.set(prot3.outputMicrographs)
    proj.launchProtocol(prot5, wait=True)
      
    
else:
    print "Not micrographs found"
