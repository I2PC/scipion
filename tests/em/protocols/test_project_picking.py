import sys

from pyworkflow.manager import Manager
from pyworkflow.project import *
from pyworkflow.em import *
from pyworkflow.em.packages.xmipp3.protocol_particle_pick import XmippProtParticlePicking


projName = sys.argv[1]
pattern = sys.argv[2]
importFolder = sys.argv[3]

manager = Manager()
proj = manager.createProject(projName)  # Now it will be loaded if exists

from tests.tester import *

prot2 = ProtImportMicrographs(workingDir=proj.getPath('Runs', 'Import1'),
                              pattern=pattern, samplingRate=1.237, voltage=200, tiltPairs=True)
proj.launchProtocol(prot2, wait=True)

l = proj.mapper.selectByClass('SetOfMicrographs')
# l = proj.mapper.selectAll()

for p in l:
    if p.getObjId() == prot2.outputMicrographs.getObjId():
        break

if len(l):    
    # Particle Picking-------------
    print "Launching particle picking..."   
    protPP = XmippProtParticlePicking(importFolder=importFolder)
            
    protPP.inputMicrographs.set(prot2.outputMicrographs)
    
    proj.launchProtocol(protPP, wait=True)

else:
    print "Not micrographs found"
