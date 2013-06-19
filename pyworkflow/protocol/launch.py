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
This module is responsible for launching protocol executions.
There are two main scenarios: local execution and remote execution.

A. Local execution: 
This will depend on the 'localhost' configuration
1- Check if the protocol will be launched with MPI or not (using MPI template from config)
2- Check if the protocol will be submitted to a queue (using Queue template from config)
3- Build the command that will be launched.

B. Remote execution:
1- Establish a connection with remote host for protocol execution
2- Copy necesary files to remote host.
3- Run a local proccess (for local execution, see case A) in the remote host
4- Get the result back after launching remotely
"""
import re
from subprocess import Popen, PIPE
from pyworkflow.utils import buildRunCommand, redStr, greenStr, makeFilePath
from pyworkflow.protocol import STEPS_PARALLEL

UNKNOWN_JOBID = -1


def launchProtocol(protocol, wait=False):
    """ This is the entry point to launch a protocol
    This function will decide wich case, A or B will be used.
    """
    if protocol.getHostName() == 'localhost':
        return _launchLocalProtocol(protocol, wait)
    else:
        return _launchRemoteProtocol(protocol, wait)
    
    
def _launchLocalProtocol(protocol, wait):
    bg = not wait
    # Check first if we need to launch with MPI or not
    if (protocol.stepsExecutionMode == STEPS_PARALLEL and
        protocol.numberOfMpi > 1):
        program = 'pw_protocol_mpirun.py'
        mpi = protocol.numberOfMpi.get() + 1
    else:
        program = 'pw_protocol_run.py'
        mpi = 1
    protStrId = protocol.strId()
    params = '%s %s' % (protocol.dbPath, protStrId)
    command = buildRunCommand(None, program, params, 
                              mpi, protocol.numberOfThreads.get(), bg)
    # Check if need to submit to queue
    hostConfig = protocol.getHostConfig()
    submitToQueue = hostConfig.isQueueMandatory() or protocol.useQueue()
    
    if submitToQueue:        
        submitDict = protocol.getSubmitDict()
        submitDict['JOB_COMMAND'] = command
        jobId = _submitProtocol(hostConfig, submitDict)
    else:
        jobId = _runProtocol(command, wait)

    return jobId

    
def _launchRemoteProtocol(protocol, wait):
    raise Exception("Not yet implemented")


def _submitProtocol(hostConfig, submitDict):
    """ Submit a protocol to a queue system. 
    """
    # Create forst the submission script to be launched
    # formatting using the template
    template = hostConfig.getSubmitTemplate() % submitDict
    #FIXME: CREATE THE PATH FIRST
    scripPath = submitDict['JOB_SCRIPT']
    f = open(scripPath, 'w')
    #Ensure the path exists
    makeFilePath(scripPath)
    f.write(template)
    f.close()
    # This should format the command using a template like: 
    # "qsub %(JOB_SCRIPT)s"
    command = hostConfig.getSubmitCommand() % submitDict
    gcmd = greenStr(command)
    print "** Submiting to queue: '%s'" % gcmd
    p = Popen(command, shell=True, stdout=PIPE)
    out = p.communicate()[0]
    # Try to parse the result of qsub, searching for a number (jobId)
    jobId = UNKNOWN_JOBID
    s = re.search('(\d+)', out)
    if s:
        jobId = int(s.group(0))
    else:
        print "** Couldn't parse %s ouput: %s" % (gcmd, redStr(out)) 
        
    return jobId 
    
    
def _runProtocol(command, wait):
    """Directly execute the protocol.
    mpiComm is only used when the protocol steps are execute
    with several MPI nodes.
    """
    gcmd = greenStr(command)
    print "** Running protocol: '%s'" % gcmd
    p = Popen(command, shell=True)
    jobId = p.pid
    if wait:
        p.wait()
        
    return jobId

