#!/usr/bin/env xmipp_python

import os, glob

xmippHome=os.getenv('XMIPP_HOME')
testPrograms = glob.glob(os.path.join(xmippHome,'bin/xmipp_test_*'))
for program in testPrograms:
    print('Testing '+program+' --------------------------')
    os.system(program)
