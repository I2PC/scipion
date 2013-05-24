import sys

from pyworkflow.manager import Manager
from pyworkflow.project import *
from pyworkflow.em import *
from pyworkflow.em.packages.brandeis import *


projName = sys.argv[1]

manager = Manager()
proj = manager.createProject(projName) # Now it will be loaded if exists

from tests.tester import *

l = proj.mapper.selectByClass('SetOfMicrographs')

for p in l:
    p.printAll()
    

if len(l):
    protCTF = ProtCTFFind()
        
    protCTF.inputMicrographs.set(l[0])

    proj.launchProtocol(protCTF, wait=True)

    l = proj.mapper.selectByClass('SetOfMicrographs')

    for p in l:
        if p.hasCTF():
            for mic in p:
                if mic.hasCTF():
                    mic.ctfModel.printAll() 
    
else:
    print "Not micrographs found"
