import sys

from pyworkflow.manager import Manager
from pyworkflow.project import *
from pyworkflow.em import *
from pyworkflow.em.packages.xmipp3 import *


projName = sys.argv[1]

manager = Manager()
proj = manager.createProject(projName) # Now it will be loaded if exists

from tests.tester import *

#prot2 = ProtImportMicrographs(pattern=pattern, samplingRate=1, voltage=200)
#proj.launchProtocol(prot2, wait=True)

coords = proj.mapper.selectByClass('XmippSetOfCoordinates')

mics = proj.mapper.selectByClass('SetOfMicrographs')

#for p in mics:
#    p.printAll()
#
#sys.exit(0)   

if len(coords) and len(mics):
    print "We have a set of coordinates..."
    protExtract = XmippProtExtractParticles(boxSize=258, downsampleType=1, downFactor=3)
    protExtract.inputCoordinates.set(coords[-1])
    protExtract.inputMicrographs.set(mics[-1])

    proj.launchProtocol(protExtract, wait=True)
    
else:
    print "Not coordinates found"
