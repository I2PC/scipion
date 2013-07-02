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
This module will be executed in remote host in order
to launch a protocol remotely. 
It will receive as params:
    - path to mapper to read protocol
    - id of the protocol object in the mapper
    - wait True force to wait until protocol has finished 
"""
import sys
from os.path import basename
from pyworkflow.protocol import getProtocolFromDb
from pyworkflow.em import *
from pyworkflow.apps.config import *
from pyworkflow.protocol.launch import _launchLocalProtocol


if __name__ == '__main__':
    if len(sys.argv) > 3:
        dbPath = sys.argv[1]
        protId = int(sys.argv[2])
        wait = bool(sys.argv[3])
        protocol = getProtocolFromDb(dbPath, protId, globals())
        jobId = _launchLocalProtocol(protocol, wait)
        # Print the status to be read 
        print "OUTPUT jobId %s" % jobId
    else:
        print "usage: %s dbPath protocolID wait" % basename(sys.argv[0])
