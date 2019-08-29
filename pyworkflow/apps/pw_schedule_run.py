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
import json
import argparse

from pyworkflow.em import *
from pyworkflow.protocol import (getProtocolFromDb,
                                 STATUS_FINISHED, STATUS_ABORTED, STATUS_FAILED)

# Add callback for remote debugging if available.
try:
    from rpdb2 import start_embedded_debugger
    from signal import signal, SIGUSR2

    signal(SIGUSR2, lambda sig, frame: start_embedded_debugger('a'))
except ImportError:
    pass


class RunScheduler():
    """ Check that all dependencies are met before launching a run. """

    def _parseArgs(self):
        parser = argparse.ArgumentParser()
        _addArg = parser.add_argument  # short notation

        _addArg("projPath", metavar='PROJECT_NAME',
                help="Project database path.")

        _addArg("dbPath", metavar='DATABASE_PATH',
                help="Protocol database path.")

        _addArg("protId", type=int, metavar='PROTOCOL_ID',
                help="Protocol ID.")

        _addArg("--sleep_time", type=int, default=15,
                dest='sleepTime', metavar='SECONDS',
                help="Sleeping time (in seconds) between updates.")

        _addArg("--wait_for", nargs='*', type=int, default=[],
                dest='waitProtIds', metavar='PROTOCOL_ID',
                help="List of protocol ids that should be not running "
                     "(i.e, finished, aborted or failed) before this "
                     "run will be executed.")

        self._args = parser.parse_args()

    def _loadProtocol(self):
        return getProtocolFromDb(self._args.projPath,
                                 self._args.dbPath,
                                 self._args.protId, chdir=True)

    def main(self):
        self._parseArgs()

        stopStatuses = [STATUS_FINISHED, STATUS_ABORTED, STATUS_FAILED]

        # Enter to the project directory and load protocol from db
        protocol = self._loadProtocol()
        mapper = protocol.getMapper()

        log = open(protocol._getLogsPath('schedule.log'), 'w')
        pid = os.getpid()
        protocol.setPid(pid)

        prerequisites = map(int, protocol.getPrerequisites())

        def _log(msg):
            print >> log, "%s: %s" % (pwutils.prettyTimestamp(), msg)
            log.flush()

        _log("Scheduling protocol %s, pid: %s, prerequisites: %s"
             % (protocol.getObjId(), pid, prerequisites))

        mapper.store(protocol)
        mapper.commit()
        mapper.close()

        # Keep track of the last time the protocol was checked and
        # its modification date to avoid unnecessary db opening
        lastCheckedDict = {}

        def _updateProtocol(protocol, project):
            protDb = protocol.getDbPath()

            if os.path.exists(protDb):
                protId = protocol.getObjId()
                lastChecked = lastCheckedDict.get(protId, None)
                lastModified = pwutils.getFileLastModificationDate(protDb)

                if lastChecked is None or (lastModified > lastChecked):
                    project._updateProtocol(protocol,
                                            skipUpdatedProtocols=False)
                    _log("Updated protocol %s" % protId)

        def _getProtocolFromPointer(pointer):
            """
            The function return a protocol from an attribute

               A) When the pointer points to a protocol

               B) When the pointer points to another object (INDIRECTLY).
                  - The pointer has an _extended value (new parameters
                    configuration in the protocol)

               C) When the pointer points to another object (DIRECTLY).
                  - The pointer has not an _extended value (old parameters
                    configuration in the protocol)
            """
            output = pointer.get()
            if isinstance(output, Protocol):  # case A
                protocol = output
            else:
                if pointer.hasExtended():  # case B
                    protocol = pointer.getObjValue()
                else:  # case C
                    protocol = self.getProject().getProtocol(
                        output.getObjParentId())
            return protocol

        while True:
            protocol = self._loadProtocol()
            project = protocol.getProject()

            # Check if there are missing inputs
            missing = False
            # Check if there are input protocols failed or aborted
            failedInputProtocols = False

            _log("Checking input data...")
            # FIXME: This does not cover all the cases:
            # When the user registers new coordinates after clicking the
            # "Analyze result" button, this action is registered in the project.sqlite
            # and not in it's own run.db and never gets updated. It is not critical and
            # will only affect a combination of json import with extended that will
            # appear after clicking on the "Analyze result" button.
            for key, attr in protocol.iterInputAttributes():
                inputProt = _getProtocolFromPointer(attr)
                _updateProtocol(inputProt, project)
                if (inputProt.getStatus() == STATUS_ABORTED or
                        inputProt.getStatus() == STATUS_FAILED):
                    failedInputProtocols = True

            if not failedInputProtocols:
                _updateProtocol(protocol, project)
                if len(protocol.validate()) > 0:
                        missing = True
                elif not protocol.worksInStreaming():
                    for key, attr in protocol.iterInputAttributes():
                        pointer = attr.get()
                        if isinstance(pointer, Set) and pointer.isStreamOpen():
                            missing = True
                            break

                if not missing:
                    inputProtocolDict = protocol.inputProtocolDict()
                    for prot in inputProtocolDict.values():
                        _updateProtocol(prot, project)

            _log("Checking prerequisited...")
            wait = False  # Check if we need to wait for required protocols
            for protId in prerequisites:
                prot = project.getProtocol(protId)
                if prot is not None:
                    _updateProtocol(prot, project)
                    if prot.getStatus() not in stopStatuses:
                        wait = True

            if not missing and not wait:
                break

            project.mapper.commit()
            project.mapper.close()

            _log("Still missing input, sleeping %s seconds..."
                 % self._args.sleepTime)
            time.sleep(self._args.sleepTime)

        _log("Launching the protocol >>>>")
        log.close()
        project.launchProtocol(protocol, scheduled=True, force=True)


if __name__ == '__main__':
    scheduler = RunScheduler()
    scheduler.main()
