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
import platform

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
def getListFromRangeString(rangeStr):
    ''' Create a list of integer from a string with range definitions
    Examples:
    "1,5-8,10" -> [1,5,6,7,8,10]
    "2,6,9-11" -> [2,6,9,10,11]
    "2 5, 6-8" -> [2,5,6,7,8]
    '''
    elements = rangeStr.split(',')
    values = []
    for e in elements:
        if '-' in e:
            limits = e.split('-')
            values += range(int(limits[0]), int(limits[1])+1)
        else:
            # If values are separated by comma also splitted 
            values += map(int, e.split())
    return values

def getRangeStringFromList(list):
    left = None
    right = None
    ranges = []

    def addRange():
        if left == right: # Single element
            ranges.append("%d" % right)
        else:
            ranges.append("%(left)d-%(right)d" % locals())
    
    for item in list:
        if right is None:
            left = right = item
        else:
            if item == right + 1:
                right += 1
            else:
                addRange()
                left = right = item
    addRange()
    return ','.join(ranges)
    
                

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
            try:
                listValues += [ listaIntervalo[1] ] * string.atoi(listaIntervalo[0])
            except:
                listValues += listaIntervalo
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
        
    def __str__(self):
        return str(self.info)
        
class ProcessManager():
    def __init__(self, run):
        self.run = run
        from protlib_sql import SqliteDb
        self.isBatch = run['jobid'] > SqliteDb.NO_JOBID;
               
    def getProcessFromCmd(self, cmd):
        procs = []
        self.hostfile = self.run['script'].replace('.py', '.nodes')
        
        if os.path.exists(self.hostfile) and self.isBatch:
            hosts = {}
            f = open(self.hostfile)
            for line in f:                
                hosts[line.strip()] = True
            f.close()
            
            def setHost(p, h):
                p.host = h
                
            for h in hosts.keys():
                newProcs = self.__getProcessFromCmd("ssh %(h)s '%(cmd)s'" % locals() )
                for p in newProcs:
                    setHost(p, h) # Python annoying looping problem
                procs += newProcs
            return procs
        else:
            procs = self.__getProcessFromCmd(cmd)
        
        return procs
    
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
        n = len(procList)
        if n == 0:
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
        script = 'xmipp_protocol_script %s' % os.path.abspath(self.run['script'])
        return self.getProcessFromCmd('ps -A -o pid,ppid,cputime,etime,state,pcpu,pmem,args| grep "%(script)s" | grep -v grep ' % locals())

    def stopProcessGroup(self, project=None):
        if project != None:
            from protlib_sql import SqliteDb
            project.projectDb.updateRunState(SqliteDb.RUN_ABORTED, self.run['run_id'])
            
        if self.isBatch:
            launch = loadLaunchModule()
            cmd = launch.StopCommand + " " + launch.StopArgsTemplate
            p = Popen(cmd % self.run, shell=True, stdout=PIPE)
            os.waitpid(p.pid, 0)
        else:
            childs = self.getProcessGroup()
            for c in childs:
                c.terminate()
        
                
    def isAlive(self):
        if self.isBatch:
            launch = loadLaunchModule()
            cmd = launch.QueryCommand + " " + launch.QueryArgsTemplate            
            from subprocess import call
            retcode = 1000
            try:
                fnull = open(os.devnull, 'w')
                retcode = call(cmd % self.run, shell=True, stdout=fnull, stderr=fnull)                
                fnull.close()
            except OSError, e:
                raise Exception("Failed %s, command: %s" % (e, cmd))
            return (retcode == 0)
        else:
            childs = self.getProcessGroup()
            return len(childs) > 0
    
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
    if programname.startswith("xmipp_") and os.environ.has_key('PROTOCOL_SCRIPT'):
        params += ' --xmipp_protocol_script ' + os.environ['PROTOCOL_SCRIPT']
    
    if not DoParallel:
        command = programname + ' ' + params
    else:
        paramsDict['nodes'] = NumberOfMpi
        prog = programname.replace('xmipp', 'xmipp_mpi')
        paramsDict['command'] = "`which %(prog)s` %(params)s" % locals()
        launch = loadLaunchModule()
        command = launch.MpiProgram + " " + launch.MpiArgsTemplate % paramsDict

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
    ''' Load the launch module containing queue and mpi related parameters
    the actual configuration should be in [parallel] section of the XMIPP/.xmipp.cfg file
    '''
    launchModuleName = os.environ['XMIPP_PARALLEL_LAUNCH']
    return loadModule(launchModuleName)
        
def submitProtocol(script, **params):
    '''Launch a protocol, to a queue or executing directly.
    If the queue options are found, it will be launched with 
    configuration (command and file template) found in project settings
    This function should be called from ProjectDir
    '''
    #Load the config module
    launch = loadLaunchModule()
    launchFilename = script.replace('.py', '.job')
    # This is for make a copy of nodes files
    nodesFile = script.replace('.py', '.nodes')
    params['nodesfileBackup'] = nodesFile
    params['file'] = launchFilename
    #create launch file
    launchfile = open(launchFilename, 'w')
    launchfile.write(launch.FileTemplate % params)
    launchfile.close()
    command = launch.Program + " " + launch.ArgsTemplate % params
    from protlib_xmipp import greenStr, redStr
    from protlib_sql import SqliteDb
    print "** Submiting to queue: '%s'" % greenStr(command)
    ps = Popen(command, shell=True, stdout=PIPE)
    out = ps.communicate()[0]
    import re
    s = re.search('(\d+)', out)
    if s:
        return int(s.group(0))
    else:
        print "** Couldn't parse %s ouput: %s" % (greenStr(launch.Program), redStr(out)) 
        return SqliteDb.UNKNOWN_JOBID

def submitProgram(script, **params):
    ''' Same function as submitProtocol but just for single
    programs, not need to be inside a Project '''
        #Load the config module
    launch = loadLaunchModule()
    launchFilename = script.replace('.py', '.job')
    nodesFile = script.replace('.py', '.nodes')
    params['nodesfileBackup'] = nodesFile
    params['file'] = launchFilename
    #create launch file
    launchfile = open(launchFilename, 'w')
    launchfile.write(launch.FileTemplate % params)
    launchfile.close()
    command = launch.Program + " " + launch.ArgsTemplate % params
    from protlib_xmipp import greenStr
    print "** Submiting to queue: '%s'" % greenStr(command)
    os.system(command)
    
def getImageJPluginCmd(memory, macro, args, batchMode=False):
    from protlib_filesystem import getXmippPath
    if len(memory) == 0:
        memory = "1g"
        print "No memory size provided. Using default: " + memory
    imagej_home = getXmippPath("external/imagej")
    plugins_dir = os.path.join(imagej_home, "plugins")
    macro = os.path.join(imagej_home, "macros", macro)
    imagej_jar = os.path.join(imagej_home, "ij.jar")
    cmd = """ java -Xmx%s -Dplugins.dir=%s -jar %s -macro %s "%s" """ % (memory, plugins_dir, imagej_jar, macro, args)
    if batchMode:
        cmd += " &"
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
        from protlib_gui_ext import showError
        showError("Error waiting for response", "No reponse, returning empty string. ERROR: " + str(e))

    return msg.replace('__END__', '')
 
def runImageJPluginWithResponse(memory, macro, args):
    return runExternalAppWithResponse(getImageJPluginCmd(memory, macro, args, True))

def getArchitecture():
    arch = platform.architecture()[0]
    for a in ['32', '64']:
        if a in arch:
            return a
    return 'NO_ARCH' 
    
def getJavaIJappCmd(memory, appName, args, batchMode=False):
    '''Launch an Java application based on ImageJ '''
    from protlib_filesystem import getXmippPath
    if len(memory) == 0:
        memory = "2g"
        print "No memory size provided. Using default: " + memory
    imagej_home = getXmippPath("external", "imagej")
    lib = getXmippPath("lib")
    javaLib = getXmippPath('java', 'lib')
    plugins_dir = os.path.join(imagej_home, "plugins")
    arch = getArchitecture()
    cmd = "java -Xmx%(memory)s -d%(arch)s -Djava.library.path=%(lib)s -Dplugins.dir=%(plugins_dir)s -cp %(imagej_home)s/*:%(javaLib)s/* %(appName)s %(args)s" % locals()
    if batchMode:
        cmd += " &"
    return cmd
    
def runJavaIJapp(memory, appName, args, batchMode=True):
    cmd = getJavaIJappCmd(memory, appName, args, batchMode)
    print cmd
    os.system(cmd)
    
def runJavaJar(memory, jarName, args, batchMode=True):
    from protlib_filesystem import getXmippPath
    jarPath = getXmippPath(jarName)
    runJavaIJapp(memory, '-jar %s' % jarPath, args, batchMode)

def runJavaIJappWithResponse(memory, appName, args):
    return runExternalAppWithResponse(getJavaIJappCmd(memory, appName, args, True))

def runShowJ(inputFiles, memory="1g", extraParams=""):
    runJavaIJapp(memory, "'xmipp.viewer.Viewer'", "-i %s %s" % (inputFiles, extraParams), True)
    
def runMaskToolbar(inputFile, memory="1g", extraParams=""):
    runShowJ(inputFile, memory, extraParams + " --mask_toolbar")
    
def runChimera(inputFile,extraParams=""):
    if which("chimera") and os.path.exists(inputFile):
        from protlib_filesystem import hasSpiderExt
        if hasSpiderExt(inputFile):
            inputFile = 'spider:%s' % inputFile
        os.system("chimera %s %s &" % (inputFile,extraParams))
    else:
        print "Error Chimera not available or inputFile %s does not exits."%inputFile

def runVMD(inputFile,extraParams=""):
    if which("vmd") and os.path.exists(inputFile):
        os.system("vmd %s %s" % (inputFile,extraParams))
    else:
        print "Error VMD not available or inputFile %s does not exits."%inputFile

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

def grepFirst(filename, text):
    """Returns the first line in filename in which text appears"""
    for line in open(filename):
        if text in line:
            return line
    return None

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

def pretty_size(size):
    """Human friendly file size"""
    from math import log
    unit_list = zip(['bytes', 'kB', 'MB', 'GB', 'TB', 'PB'], [0, 0, 1, 2, 2, 2])
    if size > 1:
        exponent = min(int(log(size, 1024)), len(unit_list) - 1)
        quotient = float(size) / 1024**exponent
        unit, num_decimals = unit_list[exponent]
        format_string = '{:.%sf} {}' % (num_decimals)
        return format_string.format(quotient, unit)
    if size == 0:
        return '0 bytes'
    if size == 1:
        return '1 byte'

def createUniqueFileName(fn):
    '''
    This function creates a file name that is similar to the original 
    by adding a unique numeric suffix. check   NamedTemporaryFile
    from tempfile for alternatives
    '''
    if not os.path.exists(fn):
        return fn

    path, name = os.path.split(fn)
    name, ext = os.path.splitext(name)

    make_fn = lambda i: os.path.join(path, '%s__tmp_%d__%s' % (name, i, ext))

    for i in xrange(2, sys.maxint):
        uni_fn = make_fn(i)
        if not os.path.exists(uni_fn):
            return uni_fn

    return None

def createMetaDataFromPattern(pattern, isStack=False, label="image"):
    ''' Create a metadata from files matching pattern'''
    import glob
    files = glob.glob(pattern)
    files.sort()

    from xmipp import MetaData, FileName, getImageSize, MDL_ENABLED, str2Label
    label = str2Label(label) #Check for label value
    
    mD = MetaData()
    inFile = FileName()
    
    nSize = 1
    for file in files:
        fileAux=file
        if isStack:
            if file.endswith(".mrc"):
                fileAux=file+":mrcs"
            x, x, x, nSize = getImageSize(fileAux)
        if nSize != 1:
            counter = 1
            for jj in range(nSize):
                inFile.compose(counter, fileAux)
                objId = mD.addObject()
                mD.setValue(label, inFile, objId)
                mD.setValue(MDL_ENABLED, 1, objId)
                counter += 1
        else:
            objId = mD.addObject()
            mD.setValue(label, fileAux, objId)
            mD.setValue(MDL_ENABLED, 1, objId)
    return mD

def getMemoryAvailable():
    return int(os.popen("free -m").readlines()[1].split()[1])
    