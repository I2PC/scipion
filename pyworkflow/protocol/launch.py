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
import pyworkflow as pw
from pyworkflow.utils import buildRunCommand, redStr, greenStr, makeFilePath, join
from pyworkflow.utils import process
from pyworkflow.protocol import STEPS_PARALLEL

UNKNOWN_JOBID = -1
LOCALHOST = 'localhost'


# ******************************************************************
# *         Public functions provided by the module
# ******************************************************************

def launch(protocol, wait=False):
    """ This function should be used to launch a protocol
    This function will decide wich case, A or B will be used.
    """
    if _isLocal(protocol):
        jobId = _launchLocal(protocol, wait)
    else:
        jobId = _launchRemote(protocol, wait)
    
    protocol.setJobId(jobId)
    
    
def stop(protocol):
    """ 
    """
    if _isLocal(protocol):
        return _stopLocal(protocol)
    else:
        return _stopRemote(protocol)



# ******************************************************************
# *         Internal utility functions
# ******************************************************************

def _isLocal(protocol):
    return protocol.getHostName() == LOCALHOST
    
# ******************************************************************
# *                 Function related to LAUNCH
# ******************************************************************
def _getAppsProgram(prog):
    """ Get a command to launch a program under the apps folder.
    And also using a different python if configured in SCIPION_PYTHON var.
    """
    return pw.join('apps', prog)
    
def _launchLocal(protocol, wait):
    bg = not wait
    # Check first if we need to launch with MPI or not
    if (protocol.stepsExecutionMode == STEPS_PARALLEL and
        protocol.numberOfMpi > 1):
        script = _getAppsProgram('pw_protocol_mpirun.py')
        mpi = protocol.numberOfMpi.get() + 1
    else:
        script = _getAppsProgram('pw_protocol_run.py')
        mpi = 1
    protStrId = protocol.strId()
    threads = protocol.numberOfThreads.get()
    params = '%s %s' % (protocol.getDbPath(), protStrId)
    hostConfig = protocol.getHostConfig()
    command = process.buildRunCommand(None, script, params, mpi, threads, bg, hostConfig)
    # Check if need to submit to queue
    useQueue = hostConfig.isQueueMandatory() or protocol.useQueue()
    
    if useQueue:        
        submitDict = protocol.getSubmitDict()
        submitDict['JOB_COMMAND'] = command
        jobId = _submit(hostConfig, submitDict)
    else:
        jobId = _run(command, wait)

    return jobId
    
def _launchRemote(protocol, wait):
    from pyworkflow.utils.remote import sshConnectFromHost, RemotePath
    # Establish connection
    host = protocol.getHostConfig()
    ssh = sshConnectFromHost(host)
    rpath = RemotePath(ssh)
    # Copy protocol files
    _copyFiles(protocol, rpath)
    # Run remote program to launch the protocol
    cmd  = 'cd %s; pw_protocol_launch.py %s %s %s' % (host.getHostPath(), 
                                                   protocol.getDbPath(), 
                                                   protocol.strId(), str(wait))
    print "Running remote %s" % greenStr(cmd)
    stdin, stdout, stderr = ssh.exec_command(cmd)
    for l in stdout.readlines():
        if l.startswith('OUTPUT jobId'):
            jobId = int(l.split()[2])
            break
    return jobId

def _copyFiles(protocol, rpath):
    """ Copy all required files for protocol to run
    in a remote execution host.
    NOTE: this function should always be execute with 
    the current working dir pointing to the project dir.
    And the remotePath is assumed to be in protocol.getHostConfig().getHostPath()
    Params:
        protocol: protocol to copy files
        ssh: an ssh connection to copy the files.
    """
    remotePath = protocol.getHostConfig().getHostPath()
    
    
    for f in protocol.getFiles():
        remoteFile = join(remotePath, f)
        rpath.putFile(f, remoteFile)

def _submit(hostConfig, submitDict):
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
    
def _run(command, wait):
    """ Execute a command in a subprocess and return the pid. """
    gcmd = greenStr(command)
    print "** Running command: '%s'" % gcmd
    p = Popen(command, shell=True)
    jobId = p.pid
    if wait:
        p.wait()
        
    return jobId

# ******************************************************************
# *                 Function related to STOP
# ******************************************************************

def _stopLocal(protocol):
    
    if protocol.useQueue():     
        jobId = protocol.getJobId()
        raise Exception("_stopLocal for queue not implemented yet")   
        host = protocol.getHostConfig()
        cancelCmd = host.getCancelCommand() % {'JOB_ID': jobId}
        _run(cancelCmd, wait=True)
    else:
        process.killWithChilds(protocol.getPid())
    
def _stopRemote(protocol):
    raise Exception("_stopRemote not implemented yet")
    from pyworkflow.utils.remote import sshConnectFromHost, RemotePath
    # Establish connection
    host = protocol.getHostConfig()
    ssh = sshConnectFromHost(host)
    rpath = RemotePath(ssh)
    # Copy protocol files
    _copyFiles(protocol, rpath)
    # Run remote program to launch the protocol
    cmd  = 'cd %s; pw_protocol_launch.py %s %s %s' % (host.getHostPath(), 
                                                   protocol.getDbPath(), 
                                                   protocol.strId(), str(wait))
    print "Running remote %s" % greenStr(cmd)
    stdin, stdout, stderr = ssh.exec_command(cmd)
    for l in stdout.readlines():
        if l.startswith('OUTPUT jobId'):
            jobId = int(l.split()[2])
            break
    return jobId