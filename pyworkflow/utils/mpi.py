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
from time import time, sleep
from process import buildRunCommand, runCommand


TIMEOUT = 60  # seconds trying to send/receive data thru a socket

TAG_RUN_JOB = 1000


def send(command, comm, dest, tag):
    """ Send command in a non-blocking way and return the exit code. """

    print "Sending command to %d: %s" % (dest, command)

    # Send command with isend()
    req_send = comm.isend(command, dest=dest, tag=tag)
    t0 = time()
    while not req_send.test()[0]:
        sleep(1)
        if time() - t0 > TIMEOUT:
            raise Exception("Timeout, cannot send command to slave.")

    # Receive the exit code in a non-blocking way too (with irecv())
    req_recv = comm.irecv(dest=dest, tag=tag)
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


def runJobMPI(log, programname, params, mpiComm, mpiDest,
              numberOfMpi=1, numberOfThreads=1,
              runInBackground=False, hostConfig=None, 
              env=None, cwd=None):
    """ Send the command to the MPI node in which it will be executed. """

    command = buildRunCommand(log, programname, params,
                              numberOfMpi, numberOfThreads,
                              runInBackground, hostConfig)
    if cwd is not None:
        send("cwd=%s" % cwd, mpiComm, mpiDest, TAG_RUN_JOB+mpiDest)

    return send(command, mpiComm, mpiDest, TAG_RUN_JOB+mpiDest)


def runJobMPISlave(mpiComm):
    """ This slave will be receiving commands to execute
    until 'None' is received. 
    """
    rank = mpiComm.Get_rank()
    print "Running runJobMPISlave: ", rank

    # Listen for commands until we get 'None'
    cwd = None  # We run without changing directory by default
    env = None
    
    while True:
        # Receive command in a non-blocking way
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
            if command.startswith("cwd="):
                cwd = command.split("=", 1)[-1]
                result = 0
            else:
                result = runCommand(command, cwd=cwd)
                cwd = None  # unset directory and environment
                env = None
        except Exception, e:
            result = str(e)

        # Send result in a non-blocking way
        req_send = mpiComm.isend(result, dest=0, tag=TAG_RUN_JOB+rank)
        t0 = time()
        while not req_send.test()[0]:
            sleep(1)
            if time() - t0 > TIMEOUT:
                print "Timeout, cannot send result to master."
                break

    print "finishing slave...", rank
