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

#---------------------------------------------------------------------------
# Logging utilities
#---------------------------------------------------------------------------
class XmippLog:
    '''This class is a simple wrapper around the most rich Python logging system
    Also providing a basic file logging for older python versions
    '''    
    def __init__(self, logname, filename):
        self.is_basic = False
        try:
            import logging
            mylog = logging.getLogger(logname)
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
        import os, socket
        myusername = str(os.environ.get('USERNAME'))
        myhost = str(socket.gethostname())
        mypwd = str(os.environ.get('PWD'))
        event = "\n"
        event += "NEW LOG SESSION\n" 
        event += "===============\n" 
        event += myusername + '@' 
        event += myhost + ':' 
        event += mypwd 
        self.info(event)
        
    def debug(self, message):
        if self.is_basic:
            self.info("DEBUG: " + message)
        else:
            self._log.debug(message)

    def info(self, message):
        if self.is_basic:
            import time
            self.fh_log.write("%s %s\n" % (message, time.asctime(time.localtime(time.time()))))
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
    import os
    script = os.path.basename(script)
    #all protocols script should start by 'xmipp_protocol', the numbering part is optional
    s = re.match('((?:\w*[a-zA-Z])+)(?:_\d+)?.py', script)
    if not s:
        raise Exception('script %s doesn\'t conform Xmipp protocol name convention' % script)
    return s.group(1)

#---------------------------------------------------------------------------
# Other utilities
#---------------------------------------------------------------------------
def makeScriptBackup(log, script, WorkingDir):
    '''Make a backup of the script
    This function assumes the execution from ProjectDir
    '''
    log.info("Making backup of script " + script);
    import shutil, os
    try:        
        script_prefix = getScriptPrefix(script)
        script_out = os.path.join(WorkingDir, script_prefix + '_backup.py')
        shutil.copy(script, script_out)
    except shutil.Error, e:
        printLog(log, e.message)     
    
#---------------------------------------------------------------------------
# Parsing of arguments
#---------------------------------------------------------------------------
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
def printLogError(log, msg):
    '''Function to write error message to log, to screen and exit'''
    log.error(msg)
    print "ERROR: ", msg
    exit(1)
    
def printLog(log, msg):
    '''Just print a msg and log'''
    log.info(msg)
    print msg


