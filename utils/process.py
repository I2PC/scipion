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

# The job should be launched from the working directory!
def runJob(log, 
           programname,
           params,           
           NumberOfMpi = 1,
           NumberOfThreads = 1,
           RunInBackground=False):

    command = buildRunCommand(log,
               programname,
               params,
               NumberOfMpi,
               NumberOfThreads,
               RunInBackground)
    if log:
        #TODO: printLog("Running command: %s" % greenStr(command),log)
        pass

    from subprocess import call
    retcode = 1000
    try:
        retcode = call(command, shell=True, stdout=sys.stdout, stderr=sys.stderr)
        if log:
            #TODO: printLog("Process returned with code %d" % retcode,log)
            if retcode != 0:
                raise Exception("Process returned with code %d, command: %s" % (retcode,command))
    except OSError, e:
        raise Exception("Execution failed %s, command: %s" % (e, command))

    return retcode

def buildRunCommand(
               log,
               programname,
               params,
               NumberOfMpi,
               NumberOfThreads,
               RunInBackground):

    DoParallel = NumberOfMpi > 1
    paramsDict={}
    
    if not DoParallel:
        command = programname + ' ' + params
    else:
        paramsDict['nodes'] = NumberOfMpi
        paramsDict['command'] = "`which %(programname)s` %(params)s" % locals()
        launch = loadLaunchModule()
        command = launch.MpiProgram + " " + launch.MpiArgsTemplate % paramsDict

    if RunInBackground:
        command+=" &"

    return command

#TODO: Adapt all this from xmipp to general
# taking into account different configuration for each Execution Host

def loadLaunchModule():
    ''' Load the launch module containing queue and mpi related parameters
    the actual configuration should be in [parallel] section of the XMIPP/.xmipp.cfg file
    '''
    pass
#    launchModuleName = os.environ['XMIPP_PARALLEL_LAUNCH']
#    return loadModule(launchModuleName)

#def loadModule(modulePath, report=True):
#    directory , moduleName = os.path.split(modulePath)
#    moduleName = moduleName.replace('.py', '')
#    if directory=='':
#        sys.path.insert(0, '.')
#    else:
#        sys.path.insert(0, directory)
#    try:
#        if moduleName in sys.modules:
#            module = sys.modules[moduleName]
#            reload(module)
#        else:
#            module = __import__(moduleName)
#    except ImportError, e:
#        if report:
#            reportError(str(e))
#        module = None
#    del sys.path[0]
#    return module