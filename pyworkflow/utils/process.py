# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              Laura del Cano         (ldelcano@cnb.csic.es)
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
This module handles process execution
"""

import sys
import resource
from subprocess import call

from utils import greenStr, envVarOn


# The job should be launched from the working directory!
def runJob(log, programname, params,           
           numberOfMpi=1, numberOfThreads=1, 
           runInBackground=False, hostConfig=None, env=None, cwd=None):

    command = buildRunCommand(log, programname, params,
                              numberOfMpi, numberOfThreads, runInBackground,
                              hostConfig)
    
    if log is None:
        print "** Running command: %s" % greenStr(command)
    else:
        log.info(greenStr(command), True)

    return runCommand(command, env, cwd)
        

def runCommand(command, env=None, cwd=None):
    """ Execute command with given environment env and directory cwd """

    # First let us create core dumps if in debug mode
    if envVarOn('SCIPION_DEBUG', env):
        resource.setrlimit(resource.RLIMIT_CORE,
                           (resource.RLIM_INFINITY, resource.RLIM_INFINITY))
        # This is like "ulimit -u 99999999", so we can create core dumps

    try:
        retcode = call(command, shell=True, stdout=sys.stdout, stderr=sys.stderr, env=env, cwd=cwd)
        if retcode != 0:
            raise Exception("Process returned with code %d, command: %s" % (retcode,command))
    except OSError, e:
        raise Exception("Execution failed %s, command: %s" % (e, command))

    return retcode
    
    
def buildRunCommand(log, programname, params,
                    numberOfMpi, numberOfThreads, runInBackground,
                    hostConfig=None):
    if numberOfMpi <= 1:
        command = programname + ' ' + params
    else:
        if programname.startswith('xmipp'):
            programname = programname.replace('xmipp', 'xmipp_mpi')
        paramsDict = {'JOB_NODES': numberOfMpi,
                      'COMMAND': "`which %s` %s" % (programname, params)
                      }
        if hostConfig is None:
            raise Exception('buildRunCommand: hostConfig is needed to launch MPI processes.')
            #hostConfig = loadHostConfig()
        command = hostConfig.mpiCommand.get() % paramsDict

    if runInBackground:
        command += " &"

    return command


def loadHostConfig(host='localhost'):
    """ This function will load the execution host configuration.
    In there will be information to know how to launch MPI processes
    and how to submit jobs to the queue system if exists.
    """
    from pyworkflow.hosts import HostMapper
    from pyworkflow.apps.config import getConfigPath
    fn = getConfigPath()
    mapper = HostMapper(fn)
    return mapper.selectByLabel(host)


def killWithChilds(pid):
    """ Kill the process with given pid and all children processes.
    Params:
     pid: the process id to terminate
    """
    import psutil
    proc = psutil.Process(pid)
    for c in proc.get_children(recursive=True):
        print "Terminating child pid: %d" % c.pid
        c.kill()
    print "Terminating process pid: %d" % pid
    proc.kill()
