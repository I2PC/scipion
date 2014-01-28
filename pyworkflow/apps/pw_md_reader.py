#!/usr/bin/env python
'''
Created on Jan 27, 2014

@author: airen
'''
import os, sys
from pyworkflow.em import *
from pyworkflow.em.packages.xmipp3 import *
from xmipp import *


if __name__ == '__main__':
    mdfile = sys.argv[1]
    prot = ProtUserSelection()
    imgSet = prot._createSetOfParticles()
    md = MetaData(mdfile)
    md.save("test.xmd")
    print md
    readSetOfParticles(mdfile, imgSet)