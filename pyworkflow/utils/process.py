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
import os.path
import resource
from subprocess import check_call

from utils import greenStr, envVarOn


# The job should be launched from the working directory!
def runJob(log, programname, params,           
           numberOfMpi=1, numberOfThreads=1, 
           hostConfig=None, env=None, cwd=None):

    command = buildRunCommand(programname, params, numberOfMpi, hostConfig, env)
    
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

    # TODO: maybe have to set PBS_NODEFILE in case it is used by "command"
    # (useful for example with gnu parallel)
    check_call(command, shell=True, stdout=sys.stdout, stderr=sys.stderr, env=env, cwd=cwd)
    # It would be nice to avoid shell=True and calling buildRunCommand()...

    
def buildRunCommand(programname, params, numberOfMpi, hostConfig=None, env=None):
    """ Return a string with the command line to run """

    # Convert our list of params to a string, with each element escaped
    # with "" in case there are spaces.
    if not isinstance(params, basestring):
        params = ' '.join('"%s"' % p for p in params)

    if numberOfMpi <= 1:
        return '%s %s' % (programname, params)
    else:
        assert hostConfig is not None, 'hostConfig needed to launch MPI processes.'

        if programname.startswith('xmipp'):
            programname = programname.replace('xmipp', 'xmipp_mpi')
            
        mpiFlags = '' if env is None else env.get('SCIPION_MPI_FLAGS', '') 

        return hostConfig.mpiCommand.get() % {
            'JOB_NODES': numberOfMpi,
            'COMMAND': "%s `which %s` %s" % (mpiFlags, programname, params),
        }


def loadHostConfig(host='localhost'):
    """ This function will load the execution host configuration.
    In there will be information to know how to launch MPI processes
    and how to submit jobs to the queue system if exists.
    """
    from pyworkflow.hosts import HostMapper
    from pyworkflow import HOME
    mapper = HostMapper(os.path.join(HOME, '..', 'settings'))
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

def isProcessAlive(pid):

    import psutil
    try:
        proc = psutil.Process(pid.get())
        return proc.is_running()
    except:
        return False

