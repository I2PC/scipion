import sys

from pyworkflow.manager import Manager
from pyworkflow.project import *
from pyworkflow.em import *
from pyworkflow.em.packages.xmipp3 import *

projName = sys.argv[1]

manager = Manager()
proj = manager.createProject(projName) # Now it will be loaded if exists

from tests.tester import *

#prot = MyProtocol(name="test3", n=2, workingDir=proj.getPath('Runs', 'T3'))
#proj.launchProtocol(prot)

prot2 = ProtImportMicrographs(workingDir=proj.getPath('Runs', 'Import1'), 
                              pattern="/home/josem/pyworkflow/InputData/*.mrc", sampling=1, voltage=200)
proj.launchProtocol(prot2)

l = proj.mapper.select(classname='SetOfMicrographs')
#l = proj.mapper.getAll()

for p in l:
    p.printAll()

if len(l):
    print "Running downsampling..."
    prot3 = XmippProtDownsampleMicrographs(workingDir=proj.getPath('Runs', 'Downsample1'),
                                       inputMicrographs=l[0], downFactor=2)

    proj.launchProtocol(prot3)

    l = proj.mapper.select(classname='SetOfMicrographs')
#l = proj.mapper.getAll()

    for p in l:
      p.printAll()
else:
    print "Not micrographs found"
