'''
/***************************************************************************
 * Authors:     J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
 *
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/
 '''

import os
import sys
import tkMessageBox

#---------------------------------------------------------------------------
# Logging utilities
#---------------------------------------------------------------------------
class XmippLog:
    '''This class is a simple wrapper around the most rich Python logging system
    Also providing a basic file logging for older python versions
    '''    
    def __init__(self, filename):
        self.is_basic = False
        try:
            import logging
            mylog = logging.getLogger(filename)
            hdlr = logging.FileHandler(filename)
            formatter = logging.Formatter('(%(asctime)s) %(levelname)s (%(lineno)4d) %(message)s')
            hdlr.setFormatter(formatter)
            mylog.addHandler(hdlr) 
            mylog.setLevel(logging.INFO)
            self._log = mylog
        except ImportError:
            self._logfile = open(filename, 'a')
            self.is_basic = True
        # append a line with user, machine and date
        import socket, pwd
        myusername = pwd.getpwuid( os.getuid() )[ 0 ] # str(os.environ.get('USERNAME'))
        myhost = str(socket.gethostname())
        mypwd = str(os.environ.get('PWD'))
        event = "\n"
        event += "NEW LOG SESSION\n" 
        event += "===============\n" 
        event += "Running at: "+myusername + '@' 
        event += myhost + ':' 
        event += mypwd 
        self.info(event)
        
    def error(self, message):
        if self.is_basic:
            self.info("ERROR: " + message)
        else:
            self._log.error(message)

    def debug(self, message):
        if self.is_basic:

            self.info("DEBUG: " + message)
        else:
            self._log.debug(message)

    def info(self, message):
        if self.is_basic:
            self.fh_log.write(getLogMessage(message))
            self.fh_log.flush()
        else:
            self._log.info(message)
            
    def cat(self, filename):
        '''Cat a file contents into a log'''
        fh = open(filename, 'r')
        self.info("# Content of " + filename)
        for line in fh:
            self.info("#       " + line.strip())
        fh.close()       

    def __del__(self):
        if self.is_basic:
            self.fh_log.close()

def getLogMessage(message):
    return "%s:   %s" % (getCurrentTimeStr(), message)
#---------------------------------------------------------------------------
# Time helper functions
#---------------------------------------------------------------------------

def getCurrentTimeStr():
    import time
    return time.asctime(time.localtime(time.time()))

#---------------------------------------------------------------------------
# Naming conventions
#---------------------------------------------------------------------------
def getScriptPrefix(script):
    '''This function will extract the root of the protocol filename
    By example:
    xmipp_protocol_ml2d_00001.py -> xmipp_protocol_ml2d
    xmipp_protocol_ml2d.py       -> xmipp_protocol_ml2d
    
    script - script filename, could also contains path
    '''
    import re
    script = os.path.basename(script)
    #all protocols script should start by 'xmipp_protocol', the numbering part is optional
    s = re.match('((?:\w*[a-zA-Z])+)(?:_(\d+))?(?:\.py)?', script)
    if not s:
        raise Exception('script %s doesn\'t conform Xmipp protocol name convention' % script)
    return s.groups()

#---------------------------------------------------------------------------
# Parsing of arguments
#---------------------------------------------------------------------------
def getRangeValuesFromString(rangeStr):
    import re
    elements=re.compile(r'[, ]').split(rangeStr)
    values=[]
    for element in elements:
        if element.isdigit():
            values.append(int(element))
        else:
            limits=element.split('-')
            values+=range(int(limits[0]),int(limits[1])+1)
    return values

def getComponentFromVector(__vector, _iteration):
    ''' Convert a string to a vector of parameters'''
    _vector = __vector.strip()
    listValues = getListFromVector(_vector)
    if _iteration < 0: _iteration = 0
    if _iteration < len(listValues): return listValues[_iteration]
    else:                          return listValues[len(listValues) - 1]


def getListFromVector(_vector,numberIteration=None):
    ''' Convert a string to a list '''
    import string
    intervalos = string.split(_vector)
    if len(intervalos) == 0:
        raise RuntimeError, "Empty vector"
    listValues = []
    for i in range(len(intervalos)):
        intervalo = intervalos[i]
        listaIntervalo = string.split(intervalo, 'x')
        if len(listaIntervalo) == 1:
            listValues += listaIntervalo
        elif len(listaIntervalo) == 2:
            listValues += [ listaIntervalo[1] ] * string.atoi(listaIntervalo[0])
        else:
            raise RuntimeError, "Unknown syntax: " + intervalos
    #fill with last value the iterations
    if( numberIteration):
        for i in range(len(listValues),numberIteration):
            listValues.append(listValues[-1])
        
    return listValues

def getBoolListFromVector(_vector,numberIteration=None):
    ''' Convert a string to a list of booleans'''
    listValues = getListFromVector(_vector,numberIteration)
    listValuesBool = []
    for i in range(len(listValues)):
        if listValues[i]=='0':
            listValuesBool.append(False)
        else:
            listValuesBool.append(True)
    return listValuesBool

#---------------------------------------------------------------------------
# Error handling
#---------------------------------------------------------------------------
def reportError(msg):
    from protlib_xmipp import failStr
    '''Function to write error message to stderr and raise an Exception'''
    print >> sys.stderr, failStr("ERROR: %s" %  msg)
    raise Exception(msg)

def showWarnings(warningList, notConfirm=False):
    '''Function to write error message to log, to screen and exit'''
    if not warningList or len(warningList) == 0:
        return True
    for warning in warningList:
        from protlib_xmipp import warnStr 
        print >> sys.stderr, warnStr("WARNING: %s"% warning)
    if notConfirm:
        return True
    answer = raw_input('Do you want to proceed? [y/N]:')
    if not answer or answer.lower() == 'n':
        return False
    return True
        
def printLog(msg, log=None, out=True, err=False, isError=False):
    '''Just print a msg in the log'''
    from protlib_xmipp import failStr, findColor
    if not log is None:
        if isError: log.error(failStr(msg))
        else:       log.info(msg)
        
    if out or err:
        if isError:
            color_tuple = findColor(msg)
            if not color_tuple is None:
                msg = color_tuple[3]
            msg = failStr("ERROR: "+ msg)
        msg = getLogMessage(msg)
        if out:
            print msg
            sys.stdout.flush()
        if err:
            print >> sys.stderr, msg
            sys.stderr.flush()
    
#---------------------------------------------------------------------------
# Jobs launching and management
#---------------------------------------------------------------------------    
from subprocess import Popen, PIPE

class Process():
    keys = ['pid','ppid','cputime','etime','state','pcpu','pmem','args']
    def __init__(self, values):
        self.info = dict(zip(Process.keys, values))
        self.__dict__.update(self.info) 
        self.type = 0
        self.host = 'localhost'
    
    ''' Terminate process '''
    def terminate(self):
        p = Popen('kill %s > "/dev/null" 2>&1' % self.pid, shell=True, stdout=PIPE)
        os.waitpid(p.pid, 0)
        
class ProcessManager():
    def __init__(self, run):
        self.run = run
               
    def getProcessFromCmd(self, cmd):
        procs = []
        self.hostfile = self.run['script'].replace('.py', '.nodes')
        
        if os.path.exists(self.hostfile):
            hosts = {}
            f = open(self.hostfile)
            for line in f:                
                hosts[line.strip()] = True
            for h in hosts.keys():
                newProcs = self.__getProcessFromCmd("ssh %(h)s '%(cmd)s'" % locals() )
                for p in newProcs:
                    p.host = h
                procs += newProcs
            return procs
        else:
            procs = self.__getProcessFromCmd(cmd)
        
        if len(procs):
            return procs
        return None
    
    def __getProcessFromCmd(self, cmd):
        ''' Return process data from previous built command'''
        ps = Popen(cmd, shell=True, stdout=PIPE)
        out = ps.communicate()[0]
        if out:
            # return list of processes
            return [Process(l.split()) for l in out.splitlines()]
        return []
            
    def getUniqueProcessFromCmd(self, cmd):
        ''' Return process data from previous built command'''
        procList = self.getProcessFromCmd(cmd)
        if not procList:
            return None
        if len(procList) > 1:
            msg = [str(p) for p in procList]
            reportError("More than one process match query, only one expected\n" + "\n".join(msg))
        return procList[0]

    def getProcessFromPid(self):
        ''' Return the process data, using its arguments to match'''
        pid = self.run['pid']
        return self.getUniqueProcessFromCmd('ps -p %(pid)s -o pid,ppid,cputime,etime,state,pcpu,pmem,args| grep %(pid)s' % locals())
    
    def getProcessGroup(self):
        '''Return a list of process using the same working dir'''
        return self.getProcessFromCmd('ps -A -o pid,ppid,cputime,etime,state,pcpu,pmem,args| grep "%s" ' % self.run['script'])

    def stopProcessGroup(self):
        if os.path.exists(self.hostfile):
            launch = loadLaunchModule()
            cmd = launch.StopCommand + " " + launch.StopArgsTemplate
            p = Popen(cmd % self.run, shell=True, stdout=PIPE)
            os.waitpid(p.pid, 0)
        else:
            childs = self.getProcessGroup()
            for c in childs:
                c.terminate()        
    
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
        from protlib_xmipp import greenStr
        printLog("Running command: %s" % greenStr(command),log)

    from subprocess import call
    retcode = 1000
    try:
        retcode = call(command, shell=True, stdout=sys.stdout, stderr=sys.stderr)
        if log:
            printLog("Process returned with code %d" % retcode,log)
            if retcode!=0:
                raise Exception("Process returned with code %d" % retcode)
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
    if programname.startswith("xmipp_") and os.environ.has_key('PROTOCOL_SCRIPT'):
        params += ' --xmipp_protocol_script ' + os.environ['PROTOCOL_SCRIPT']
    
    if not DoParallel:
        command = programname + ' ' + params
    else:
        launch = loadLaunchModule()
        SystemFlavour = launch.SystemFlavour
        paramsDict['prog'] = programname.replace('xmipp', 'xmipp_mpi')
        paramsDict['jobs'] = NumberOfMpi
        paramsDict['params'] = params
        
        if (SystemFlavour == 'SLURM-MPICH'): # like BSCs MareNostrum, LaPalma etc
            mpicommand = 'srun '
        elif (SystemFlavour == 'TORQUE-OPENMPI'): # like our crunchy
            mpicommand = 'mpirun -mca mpi_yield_when_idle 1 -np %(jobs)d'
            if (int(NumberOfThreads) > 1):
                mpicommand += ' --bynode'
        elif (SystemFlavour == 'SGE-OPENMPI'): # like cluster at imp.ac.at (no variable nr_cpus yet...)
            mpicommand = 'mpiexec -n  %(jobs)d' 
        elif (SystemFlavour == 'PBS'): # like in Vermeer and FinisTerrae
            paramsDict['file']  = os.environ.get('PBS_NODEFILE')
            mpicommand = 'mpirun -np  %(jobs)d -hostfile %(file)s'
        elif (SystemFlavour == 'XMIPP_MACHINEFILE'): # environment variable $XMIPP_MACHINEFILE points to machinefile
            paramsDict['file']  = os.environ.get('XMIPP_MACHINEFILE')
            mpicommand = 'mpirun -np  %(jobs)d -machinefile %(file)s'
        elif (SystemFlavour == 'HOME_MACHINEFILE'): # machinefile is called $HOME/machines.dat
            paramsDict['file'] = os.environ.get('HOME') + '/machinefile.dat'
            mpicommand = 'mpirun -np   %(jobs)d -machinefile %(file)s'
        elif (SystemFlavour == ''):
            mpicommand = 'mpirun -mca mpi_yield_when_idle 1 -np %(jobs)d'
        else:
            from protlib_xmipp import failStr
            printLog(failStr('Unrecognized SystemFlavour %s' % SystemFlavour),log,err=True,isError=True)
        command = (mpicommand + ' `which %(prog)s` %(params)s') % paramsDict
    if RunInBackground:
        command+=" &"

    return command

def loadModule(modulePath, report=True):
    directory , moduleName = os.path.split(modulePath)
    moduleName = moduleName.replace('.py', '')
    if directory=='':
        sys.path.insert(0, '.')
    else:
        sys.path.insert(0, directory)
    try:
        if moduleName in sys.modules:
            module = sys.modules[moduleName]
            reload(module)
        else:
            module = __import__(moduleName)
    except ImportError, e:
        if report:
            reportError(str(e))
        module = None
    del sys.path[0]
    return module

def loadLaunchModule():
    ''' Load the launch module containing queue and mpi related parameters'''
    return loadModule('config_launch.py')
        
def submitProtocol(protocolPath, **params):
    '''Launch a protocol, to a queue or executing directly.
    If the queue options are found, it will be launched with 
    configuration (command and file template) found in project settings
    This function should be called from ProjectDir
    '''
    #Load the config module
    launch = loadLaunchModule()
    launchFilename = protocolPath.replace('.py', '.job')
    # This is for make a copy of nodes files
    nodesFile = protocolPath.replace('.py', '.nodes')
    params['pbsNodeBackup']= nodesFile
    params['file'] = launchFilename
    #create launch file
    launchfile = open(launchFilename, 'w')
    launchfile.write(launch.FileTemplate % params)
    launchfile.close()
    command = launch.Program + " " + launch.ArgsTemplate % params
    from protlib_xmipp import greenStr
    print "** Submiting to queue: '%s'" % greenStr(command)
    ps = Popen(command, shell=True, stdout=PIPE)
    out = ps.communicate()[0]
    return int(out.split('.')[0])
    
def getImageJPluginCmd(memory, macro, args, batchMode=False):
    from protlib_filesystem import getXmippPath
    if len(memory) == 0:
        memory = "512m"
        print "No memory size provided. Using default: " + memory
    imagej_home = getXmippPath("external/imagej")
    plugins_dir = os.path.join(imagej_home, "plugins")
    macro = os.path.join(imagej_home, "macros", macro)
    imagej_jar = os.path.join(imagej_home, "ij.jar")
    cmd = """ java -Xmx%s -Dplugins.dir=%s -jar %s -macro %s "%s" """ % (memory, plugins_dir, imagej_jar, macro, args)
    if batchMode:
        cmd += " &"
    print "getImageJPluginCmd:", cmd
    return cmd

def runImageJPlugin(memory, macro, args, batchMode=False):
    os.system(getImageJPluginCmd(memory, macro, args, batchMode))

def runExternalAppWithResponse(cmd):
    #Create a simple socket server to wait for response
    HOST = ''                 # Symbolic name meaning the local host
    PORT = 14321
    if cmd[-1]=="&":
        cmd=cmd.replace("&"," -port %d &" % PORT)
    else:
        cmd += " -port %d" % PORT
    #Launch the plugin that will send the data response over a socket
    print cmd
    os.system(cmd)
    #os.system('java SocketsTest 54321 &')
    #Create the socket server
    msg = ''
    try:
        import  socket
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        s.settimeout(5.0)
        s.bind((HOST, PORT))
        s.listen(1)
        conn = s.accept()[0]
        #Read the awake message
        data = conn.recv(256)
        if not data or data.strip() != '__STARTED__':
            raise Exception('Invalid START message received')
        else:
            conn.close() 
            s.settimeout(None) # reset timeout
            conn = s.accept()[0]
        #Read len of message (max 4 bytes)
        while True:
            data = conn.recv(1024)
            msg += data
            if not data or msg.find('__END__') !=-1: break      
        conn.close()
    except Exception, e:
        tkMessageBox.showerror("Error waiting for response", "No reponse, returning empty string. ERROR: " + str(e))

    return msg.replace('__END__', '')
 
def runImageJPluginWithResponse(memory, macro, args):
    return runExternalAppWithResponse(getImageJPluginCmd(memory, macro, args, True))
    
def getJavaIJappCmd(memory, appName, args, batchMode=False):
    '''Launch an Java application based on ImageJ '''
    from protlib_filesystem import getXmippPath
    if len(memory) == 0:
        memory = "512m"
        print "No memory size provided. Using default: " + memory
    imagej_home = getXmippPath("external", "imagej")
    plugins_dir = os.path.join(imagej_home, "plugins")

    cmd = "java -classpath %(plugins_dir)s/*:%(imagej_home)s/*: %(appName)s %(args)s" % locals()
    if batchMode:
        cmd += " &"
    return cmd
    
def runJavaIJapp(memory, appName, args, batchMode=True):
    os.system(getJavaIJappCmd(memory, appName, args, batchMode))
    
def runJavaJar(memory, jarName, args, batchMode=True):
    from protlib_filesystem import getXmippPath
    jarPath = getXmippPath(jarName)
    runJavaIJapp(memory, '-jar %s' % jarPath, args, batchMode)

def runJavaIJappWithResponse(memory, appName, args):
    return runExternalAppWithResponse(getJavaIJappCmd(memory, appName, args, True))

def runShowJ(inputFiles, memory="512m"):
    runImageJPlugin(memory, "XmippBrowser.txt", "-i "+inputFiles, True)
    
def runChimera(inputFile):
    if which("chimera") and os.path.exists(inputFile):
        os.system("chimera "+inputFile+" &")
    
   
""" Return the machine name """
def getHostname():
    import socket
    return socket.gethostname()


# Copyright (c) 2002-2005 ActiveState Corp.
# See LICENSE.txt for license details.
# Author:
#   Trent Mick (TrentM@ActiveState.com)
# Home:
#   http://trentm.com/projects/which/

"""Find the full path to commands.

which(command, path=None, verbose=0, exts=None)
    Return the full path to the first match of the given command on the
    path.

whichall(command, path=None, verbose=0, exts=None)
    Return a list of full paths to all matches of the given command on
    the path.

whichgen(command, path=None, verbose=0, exts=None)
    Return a generator which will yield full paths to all matches of the
    given command on the path.
    
"""


import stat
#---- exceptions

class WhichError(Exception):
    pass

#---- internal support stuff

def _getRegisteredExecutable(exeName):
    """Windows allow application paths to be registered in the registry."""
    registered = None
    if sys.platform.startswith('win'):
        if os.path.splitext(exeName)[1].lower() != '.exe':
            exeName += '.exe'
        import _winreg
        try:
            key = "SOFTWARE\\Microsoft\\Windows\\CurrentVersion\\App Paths\\" +\
                  exeName
            value = _winreg.QueryValue(_winreg.HKEY_LOCAL_MACHINE, key)
            registered = (value, "from HKLM\\"+key)
        except _winreg.error:
            pass
        if registered and not os.path.exists(registered[0]):
            registered = None
    return registered

def _samefile(fname1, fname2):
    if sys.platform.startswith('win'):
        return ( os.path.normpath(os.path.normcase(fname1)) ==\
            os.path.normpath(os.path.normcase(fname2)) )
    else:
        return os.path.samefile(fname1, fname2)

def _cull(potential, matches, verbose=0):
    """Cull inappropriate matches. Possible reasons:
        - a duplicate of a previous match
        - not a disk file
        - not executable (non-Windows)
    If 'potential' is approved it is returned and added to 'matches'.
    Otherwise, None is returned.
    """
    for match in matches:  # don't yield duplicates
        if _samefile(potential[0], match[0]):
            if verbose:
                sys.stderr.write("duplicate: %s (%s)\n" % potential)
            return None
    else:
        if not stat.S_ISREG(os.stat(potential[0]).st_mode):
            if verbose:
                sys.stderr.write("not a regular file: %s (%s)\n" % potential)
        elif not os.access(potential[0], os.X_OK):
            if verbose:
                sys.stderr.write("no executable access: %s (%s)\n"\
                                 % potential)
        else:
            matches.append(potential)
            return potential

        
#---- module API

def whichgen(command, path=None, verbose=0, exts=None):
    """Return a generator of full paths to the given command.
    
    "command" is a the name of the executable to search for.
    "path" is an optional alternate path list to search. The default it
        to use the PATH environment variable.
    "verbose", if true, will cause a 2-tuple to be returned for each
        match. The second element is a textual description of where the
        match was found.
    "exts" optionally allows one to specify a list of extensions to use
        instead of the standard list for this system. This can
        effectively be used as an optimization to, for example, avoid
        stat's of "foo.vbs" when searching for "foo" and you know it is
        not a VisualBasic script but ".vbs" is on PATHEXT. This option
        is only supported on Windows.

    This method returns a generator which yields either full paths to
    the given command or, if verbose, tuples of the form (<path to
    command>, <where path found>).
    """
    matches = []
    if path is None:
        usingGivenPath = 0
        path = os.environ.get("PATH", "").split(os.pathsep)
        if sys.platform.startswith("win"):
            path.insert(0, os.curdir)  # implied by Windows shell
    else:
        usingGivenPath = 1

    # Windows has the concept of a list of extensions (PATHEXT env var).
    if sys.platform.startswith("win"):
        if exts is None:
            exts = os.environ.get("PATHEXT", "").split(os.pathsep)
            # If '.exe' is not in exts then obviously this is Win9x and
            # or a bogus PATHEXT, then use a reasonable default.
            for ext in exts:
                if ext.lower() == ".exe":
                    break
            else:
                exts = ['.COM', '.EXE', '.BAT']
        elif not isinstance(exts, list):
            raise TypeError("'exts' argument must be a list or None")
    else:
        if exts is not None:
            raise WhichError("'exts' argument is not supported on "\
                             "platform '%s'" % sys.platform)
        exts = []

    # File name cannot have path separators because PATH lookup does not
    # work that way.
    if os.sep in command or os.altsep and os.altsep in command:
        pass
    else:
        for i in range(len(path)):
            dirName = path[i]
            # On windows the dirName *could* be quoted, drop the quotes
            if sys.platform.startswith("win") and len(dirName) >= 2\
               and dirName[0] == '"' and dirName[-1] == '"':
                dirName = dirName[1:-1]
            for ext in ['']+exts:
                absName = os.path.abspath(
                    os.path.normpath(os.path.join(dirName, command+ext)))
                if os.path.isfile(absName):
                    if usingGivenPath:
                        fromWhere = "from given path element %d" % i
                    elif not sys.platform.startswith("win"):
                        fromWhere = "from PATH element %d" % i
                    elif i == 0:
                        fromWhere = "from current directory"
                    else:
                        fromWhere = "from PATH element %d" % (i-1)
                    match = _cull((absName, fromWhere), matches, verbose)
                    if match:
                        if verbose:
                            yield match
                        else:
                            yield match[0]
        match = _getRegisteredExecutable(command)
        if match is not None:
            match = _cull(match, matches, verbose)
            if match:
                if verbose:
                    yield match
                else:
                    yield match[0]


def which(command, path=None, verbose=0, exts=None):
    """Return the full path to the first match of the given command on
    the path.
    
    "command" is a the name of the executable to search for.
    "path" is an optional alternate path list to search. The default is
        to use the PATH environment variable.
    "verbose", if true, will cause a 2-tuple to be returned. The second
        element is a textual description of where the match was found.
    "exts" optionally allows one to specify a list of extensions to use
        instead of the standard list for this system. This can
        effectively be used as an optimization to, for example, avoid
        stat's of "foo.vbs" when searching for "foo" and you know it is
        not a VisualBasic script but ".vbs" is on PATHEXT. This option
        is only supported on Windows.

    If no match is found for the command, a WhichError is raised.
    """
    try:
        match = whichgen(command, path, verbose, exts).next()
    except StopIteration:
        return ''
    return match


def whichall(command, path=None, verbose=0, exts=None):
    """Return a list of full paths to all matches of the given command
    on the path.  

    "command" is a the name of the executable to search for.
    "path" is an optional alternate path list to search. The default it
        to use the PATH environment variable.
    "verbose", if true, will cause a 2-tuple to be returned for each
        match. The second element is a textual description of where the
        match was found.
    "exts" optionally allows one to specify a list of extensions to use
        instead of the standard list for this system. This can
        effectively be used as an optimization to, for example, avoid
        stat's of "foo.vbs" when searching for "foo" and you know it is
        not a VisualBasic script but ".vbs" is on PATHEXT. This option
        is only supported on Windows.
    """
    return list( whichgen(command, path, verbose, exts) )

def pretty_date(time=False):
    """
    Get a datetime object or a int() Epoch timestamp and return a
    pretty string like 'an hour ago', 'Yesterday', '3 months ago',
    'just now', etc
    """
    from datetime import datetime
    now = datetime.now()
    if type(time) is int:
        diff = now - datetime.fromtimestamp(time)
    elif isinstance(time,datetime):
        diff = now - time 
    elif not time:
        diff = now - now
    second_diff = diff.seconds
    day_diff = diff.days

    if day_diff < 0:
        return ''

    if day_diff == 0:
        if second_diff < 10:
            return "just now"
        if second_diff < 60:
            return str(second_diff) + " seconds ago"
        if second_diff < 120:
            return  "a minute ago"
        if second_diff < 3600:
            return str( second_diff / 60 ) + " minutes ago"
        if second_diff < 7200:
            return "an hour ago"
        if second_diff < 86400:
            return str( second_diff / 3600 ) + " hours ago"
    if day_diff == 1:
        return "Yesterday"
    if day_diff < 7:
        return str(day_diff) + " days ago"
    if day_diff < 31:
        return str(day_diff/7) + " weeks ago"
    if day_diff < 365:
        return str(day_diff/30) + " months ago"
    return str(day_diff/365) + " years ago"
