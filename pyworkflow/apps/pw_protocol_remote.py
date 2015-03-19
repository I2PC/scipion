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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************
"""
This script is intended to be invoked using ssh from a master.
So, it will just do the same job of launch._launchLocal and print
the jobId to be tracked from the machine that was invoked.
"""

import sys


if __name__ == '__main__':
    if len(sys.argv) > 2:
        projectPath = sys.argv[1]
        protDbPath = sys.argv[2]
        protId = int(sys.argv[3])
        
        from pyworkflow.protocol import getProtocolFromDb
        from pyworkflow.protocol.launch import launch 
        
        protocol = getProtocolFromDb(projectPath, protDbPath, protId, chdir=False)
        # We need to change the hostname to localhost, since it will run here. 
        protocol.setHostName('localhost')
        jobId = launch(protocol)
        
        print "Scipion remote jobid: ", jobId
    else:
        from os.path import basename
        print "usage: %s dbPath protocolID" % basename(sys.argv[0])
