#!/usr/bin/python

import os, sys, shutil
from shutil import copy

#First we look for the address of the xmipp protocols and we store this path in scriptdir
scriptdir = os.popen('which xmipp_protocols | sed "s|/bin/xmipp_protocols|/protocols|"').read()
scriptdir = scriptdir.rstrip('\n')


if (not os.path.exists("xmipp_protocol_setup.py")):
    print "Xmipp_protocols was never executed before in this directory. \n",\
          "Make sure you run xmipp_protocols only in your project directory. \n",\
          "You are in directory: ", str(os.environ.get('PWD')), "\n"
          
    answer = raw_input("Do you want to setup Xmipp protocols here? [y/n]:")
    
    if answer == "y":
        print "Making a local copy of xmipp_protocol_setup.py ... "
        copy(scriptdir + "/xmipp_protocol_setup.py",'.')
    else:
        exit(1)
else:
    print "Reading existing copy of xmipp_protocol_setup.py ..."

print "Launching GUI ... "
    


if len(sys.argv) == 1:
    os.system(scriptdir + "/xmipp_protocol_gui.py xmipp_protocol_setup.py")
else:
    os.system(scriptdir + "/xmipp_protocol_gui.py "+ sys.argv[1])
