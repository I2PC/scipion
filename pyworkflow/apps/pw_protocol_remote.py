#!/usr/bin/env python
# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
"""
This script is intended to be invoked using ssh from a master.
So, it will just do the same job of launch._launchLocal and print
the jobId to be tracked from the machine that was invoked.
"""

import os
import sys
import subprocess


def usage(msg=''):
    print "Usage: pw_protocol_remote.py [run|stop] project protDbPath protID\n%s" % msg
    sys.exit(1)
    
    
if __name__ == '__main__':
    n = len(sys.argv)
    if n < 5:
        usage("Received only %d arguments" % n)
        
    mode = sys.argv[1]
    
    if mode not in ['run', 'stop']:
        usage("Mode should be 'run' or 'stop'. Received: '%s'" % mode)

    projectPath = os.path.join(os.environ['SCIPION_USER_DATA'], 'projects', sys.argv[2])
    print "projectPath: ", projectPath
    protDbPath = sys.argv[3]
    protId = int(sys.argv[4])

    from pyworkflow.protocol import getProtocolFromDb
    from pyworkflow.protocol.launch import launch, stop

    protocol = getProtocolFromDb(projectPath, protDbPath, protId, chdir=False)
    # We need to change the hostname to localhost, since it will 
    # be considered 'local' from now on to either run or stop
    protocol.setHostName('localhost')
    
    if mode == 'run':        
        FNULL = open(os.devnull, 'w')
        jobId = launch(protocol, stdin=None, stdout=FNULL, stderr=subprocess.STDOUT)        
        print "Scipion remote jobid: %d" % jobId
    elif mode == 'stop':
        stop(protocol)
    else:
        usage()
        
        
        

