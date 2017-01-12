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
This module is responsible for launching protocol executions.
"""
import sys
from pyworkflow.em import *
from pyworkflow.config import *

# Add callback for remote debugging if available.
try:
    from rpdb2 import start_embedded_debugger
    from signal import signal, SIGUSR2
    signal(SIGUSR2, lambda sig, frame: start_embedded_debugger('a'))
except ImportError:
    pass


if __name__ == '__main__':
    if len(sys.argv) > 2:
        projPath = sys.argv[1]
        dbPath = sys.argv[2]
        protId = int(sys.argv[3])
        from pyworkflow.protocol import runProtocolMain
        runProtocolMain(projPath, dbPath, protId)
        
    else:
        from os.path import basename
        print "usage: %s dbPath protocolID" % basename(sys.argv[0])
