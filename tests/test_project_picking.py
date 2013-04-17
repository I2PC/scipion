import sys

from pyworkflow.manager import Manager
from pyworkflow.project import *
from pyworkflow.em import *
from pyworkflow.em.packages.xmipp3.protocol_particle_pick import XmippProtParticlePicking


projName = sys.argv[1]
pattern = sys.argv[2]

manager = Manager()
proj = manager.createProject(projName)  # Now it will be loaded if exists

from tests.tester import *

# prot = MyProtocol(name="test3", n=2, workingDir=proj.getPath('Runs', 'T3'))
# proj.launchProtocol(prot)

prot2 = ProtImportMicrographs(workingDir=proj.getPath('Runs', 'Import1'),
                              pattern=pattern, samplingRate=1, voltage=200)
proj.launchProtocol(prot2)

l = proj.mapper.selectByClass('SetOfMicrographs')
# l = proj.mapper.selectAll()

#for p in l:
#    p.printAll()
l[0].printAll()
    
if len(l):    
    # Particle Picking-------------
    print "testing Particle Picking XMIPP"
    protPartPick = XmippProtParticlePicking(workingDir=proj.getPath('Runs', 'ParticlePicking'))
    protPartPick.inputMicrographs.set(l[0])

    proj.launchProtocol(protPartPick)
    sys.exit(0)

else:
    print "Not micrographs found"
