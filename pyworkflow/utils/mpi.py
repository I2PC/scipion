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
This module contains some MPI utilities
"""

import os
from time import sleep
from process import buildRunCommand, runCommand


TAG_RUN_JOB = 1000

def runJobMPI(log, programname, params, mpiComm, mpiDest,
              numberOfMpi=1, numberOfThreads=1,
              runInBackground=False, hostConfig=None, cwd=None):
    """ Send the command to the MPI node in which it will be executed. """

    print "runJobMPI: hostConfig: ", hostConfig
    command = buildRunCommand(log, programname, params,
                              numberOfMpi, numberOfThreads, 
                              runInBackground, hostConfig)

    # Send directory (cwd) in a non-blocking way (with isend())
    if cwd is not None:  # we only mention it if we actually change directory
        print "Asking to change to directory: %s" % cwd
    req_cd = mpiComm.isend(str(cwd), dest=mpiDest, tag=TAG_RUN_JOB+mpiDest)
    while not req_cd.test()[0]:
        sleep(1)

    # Send command in a non-blocking way (with isend())
    print "Sending command: %s to %d" % (command, mpiDest)
    req_send = mpiComm.isend(command, dest=mpiDest, tag=TAG_RUN_JOB+mpiDest)
    while not req_send.test()[0]:
        sleep(1)

    # Receive the exit code in a non-blocking way (with irecv())
    req_recv = mpiComm.irecv(dest=mpiDest, tag=TAG_RUN_JOB+mpiDest)
    while True:
        done, result = req_recv.test()
        if done:
            break
        sleep(1)

    # Our convention: if we get a string, an error happened.
    # Else, it is the return code of our command, which we return too.
    if isinstance(result, str):
        raise Exception(result)
    else:
        return result


def runJobMPISlave(mpiComm):
    """ This slave will be receiving commands to execute
    until 'None' is received. 
    """
    rank = mpiComm.Get_rank()
    print "Running runJobMPISlave: ", rank

    # Listen for commands until we get 'None'
    while True:
        # Receive dir (cwd) and command in a non-blocking way
        req_cd = mpiComm.irecv(dest=0, tag=TAG_RUN_JOB+rank)
        while True:
            done, cwd = req_cd.test()
            if done:
                break
            sleep(1)

        req_recv = mpiComm.irecv(dest=0, tag=TAG_RUN_JOB+rank)
        while True:
            done, command = req_recv.test()
            if done:
                break
            sleep(1)

        print "Slave %d, received command: %s" % (rank, command)
        if command == 'None':
            break

        # Run the command and get the result (exit code or exception)
        try:
            if cwd != 'None':
                os.chdir(cwd)
            result = runCommand(command)
        except Exception, e:
            result = str(e)

        # Send result in a non-blocking way
        req_send = mpiComm.isend(result, dest=0, tag=TAG_RUN_JOB+rank)
        while not req_send.test()[0]:
            sleep(1)

    print "finishing slave...", rank
