import sys

from pyworkflow.manager import Manager
from pyworkflow.project import *
from pyworkflow.em import *

projName = sys.argv[1]

manager = Manager()
proj = manager.createProject(projName) # Now it will be loaded if exists

from tests.tester import *

#prot = MyProtocol(name="test3", n=2, workingDir=proj.getPath('Runs', 'T3'))
#proj.launchProtocol(prot)

prot = ProtImportMicrographs(workingDir=proj.getPath('Runs', 'Import3'), 
                              pattern="InputData/*.mrc", sampling=1, voltage=200)
proj.launchProtocol(prot)

l = proj.mapper.select(classname='SetOfMicrographs')
#l = proj.mapper.getAll()

for p in l:
    p.printAll()

