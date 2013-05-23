import sys

from pyworkflow.manager import Manager
from pyworkflow.project import *
from pyworkflow.em import *
from pyworkflow.em.packages.xmipp3.protocol_particle_pick import XmippProtParticlePicking


projName = 'myproject'#sys.argv[1]
pattern = '/home/laura/Scipion_Projects/InputData/*.mrc' #sys.argv[2]

manager = Manager()
proj = manager.createProject(projName)  # Now it will be loaded if exists

#l = proj.mapper.selectByClass('SetOfCoordinates')
#
#for c in l[1].iterCoordinates():
#    print "%d, %d" % (c.x, c.y) 
#    print " from mic: ", c.getMicrograph().getFileName()
#    
#sys.exit(0)



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
    print "testing Particle Picking XMIPP"
    protPartPick = XmippProtParticlePicking(workingDir=proj.getPath('Runs', 'ParticlePicking'))
    #protPartPick.inputMicrographs.set(prot2.outputMicrographs)
    protPartPick.inputMicrographs.set(p)
    proj.launchProtocol(protPartPick, wait=True)
    sys.exit(0)

else:
    print "Not micrographs found"
