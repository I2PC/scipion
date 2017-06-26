#!/usr/bin/env python
# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *
# * [1] SciLifeLab, Stockholm University
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

import os
import time
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
        from pyworkflow.protocol import runProtocolMain, getProtocolFromDb

        # Enter to the project directory and load protocol from db
        protocol = getProtocolFromDb(projPath, dbPath, protId, chdir=True)
        mapper = protocol.getMapper()

        log = open(protocol._getLogsPath('schedule.log'), 'w')
        pid = os.getpid()
        protocol.setPid(pid)

        def _log(msg):
            print >> log, "%s: %s" % (pwutils.prettyTimestamp(), msg)
            log.flush()

        _log("Scheduling protocol %s, pid: %s" % (protId, pid))

        mapper.store(protocol)
        mapper.commit()
        mapper.close()

        while True:
            protocol = getProtocolFromDb(projPath, dbPath, protId, chdir=True)
            project = protocol.getProject()

            # Check if there are missing inputs
            missing = False
            # Keep track of the last time the protocol was checked and
            # its modification date to avoid unnecesary db opening
            lastCheckedDict = {}

            for key, attr in protocol.iterInputAttributes():
                if attr.hasValue() and attr.get() is None:
                    missing = True
                    inputProt = attr.getObjValue()
                    protDb = inputProt.getDbPath()

                    # One of the required input protocols has not even started
                    if not os.path.exists(protDb):
                        continue


                    inputProtId = inputProt.getObjId()
                    lastChecked = lastCheckedDict.get(inputProtId, None)
                    lastModified = pwutils.getFileLastModificationDate(protDb)

                    if lastChecked is None or (lastModified > lastChecked):
                        project._updateProtocol(inputProt,
                                                skipUpdatedProtocols=False)
                        _log("Updated protocol %s" % inputProtId)

            if not missing:
                break

            project.mapper.commit()
            project.mapper.close()

            _log("Still missing input, sleeping...")
            time.sleep(15)

        _log("Launching the protocol >>>>")
        log.close()
        project.launchProtocol(protocol, scheduled=True)
    else:
        from os.path import basename
        print "usage: %s dbPath protocolID" % basename(sys.argv[0])
