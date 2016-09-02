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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import sys
import os
import re
from datetime import datetime
import traceback
import numpy as np


def prettyDate(time=False):
    """
    Get a datetime object or a int() Epoch timestamp and return a
    pretty string like 'an hour ago', 'Yesterday', '3 months ago',
    'just now', etc
    """
    now = datetime.now()
    if type(time) is int:
        diff = now - datetime.fromtimestamp(time)
    elif type(time) is float:
        diff = now - datetime.fromtimestamp(int(time))
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


def dateStr(dt=None, time=True, secs=False, dateFormat=None):
    """ Get a normal string representation of datetime. 
    If dt is None, use NOW.
    """
    if dt is None:
        dt = datetime.now()
    elif isinstance(dt, float) or isinstance(dt, int):
        dt = datetime.fromtimestamp(dt)

    if dateFormat is None:
        dateFormat = '%d-%m-%Y'
        if time:
            dateFormat += ' %H:%M'
            if secs:
                dateFormat += ':%S'

    return dt.strftime(dateFormat)

prettyTime = dateStr


def prettyTimestamp(dt=None, format='%Y-%m-%d_%H%M%S'):
    if dt is None:
        dt = datetime.now()

    return dt.strftime(format)


def prettySize(size):
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
    

def prettyDelta(timedelta):
    """ Remove the milliseconds of the timedelta. """
    return str(timedelta).split('.')[0]


def prettyLog(msg):
    print cyan(prettyTime(datetime.now(), secs=True)), msg


class Timer(object):
    """ Simple Timer base in datetime.now and timedelta. """
    def tic(self):
        self._dt = datetime.now()
        
    def getToc(self):
        return prettyDelta(datetime.now()-self._dt)
        
    def toc(self, message='Elapsed:'):
        print message, self.getToc()
        

def timeit(func):
    """ Decorator function to have a simple measurement
    of the execution time of a given function.
    Just use:
    @timeit
    def func(...)
        ...
    to use it.
    """
    def timedFunc(*args, **kwargs):
        t = Timer()
        t.tic()
        result = func(*args, **kwargs)
        t.toc("Function '%s' took" % func)
        
        return result
        
    return timedFunc


def trace(nlevels, separator=' --> ', stream=sys.stdout):
    # Example:
    #   @trace(3)
    #   def doRefresh(...
    # gives as output whenever doRefresh is called lines like:
    #   text.py:486 _addFileTab --> text.py:330 __init__ --> doRefresh

    def realTrace(f):
        """ Decorator function to print stack call in a human-readable way.
        """
        def tracedFunc(*args, **kwargs):
            stack = traceback.extract_stack()[-nlevels-1:-1]
            fmt = lambda x: '%s:%d %s' % (os.path.basename(x[0]), x[1], x[2])
            stream.write(separator.join(map(fmt, stack)+[f.__name__]) + '\n')
            return f(*args, **kwargs)

        return tracedFunc
    return realTrace

    
def prettyDict(d):
    print "{"
    for k, v in d.iteritems():
        print "    %s: %s" % (k, v)
    print "}"


def prettyXml(elem, level=0):
    """ Add indentation for XML elements for more human readable text. """
    i = "\n" + level*"  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for _elem in elem:
            prettyXml(_elem, level+1)
        if not _elem.tail or not _elem.tail.strip():
            _elem.tail = i
    
    
def getUniqueItems(originalList):
    """ Method to remove repeated items from one list 
    originalList -- Original list with repeated items, or not.
    returns -- New list with the content of original list without repeated items
    """  
    auxDict = {}
    resultList = [auxDict.setdefault(x,x) for x in originalList if x not in auxDict]
    return resultList


def executeRemoteX (command, hostName, userName, password):
    """ Execute a remote command with X11 forwarding.
    Params:
        command: Command to execute.
        hostName: Remote host name.
        userName: User name.
        password: Password.
    Returns: 
        Tuple with standard output and error output.
    """
    scriptPath = os.path.abspath(os.path.join(os.path.dirname( __file__ ), "sshAskPass.sh"))
    pswCommand = "echo '" + password + "' | " + scriptPath + " ssh -X " + userName + "@" + hostName + " " + command
    import subprocess
    p = subprocess.Popen(pswCommand, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    return stdout, stderr

def executeRemote (command, hostName, userName, password):
    """ Execute a remote command.
    Params:
        command: Command to execute.
        hostName: Remote host name.
        userName: User name.
        password: Password.
    Returns: 
        Tuple with standard input, standard output and error output.
    """
    import paramiko
    ssh = paramiko.SSHClient()
    ssh.load_system_host_keys()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect(hostName, 22, userName, password)
    stdin, stdout, stderr = ssh.exec_command(command)
    ssh.close()
    return stdin, stdout, stderr
    
    
def executeLongRemote (command, hostName, userName, password):
    """ Execute a remote command.
    Params:
        command: Command to execute.
        hostName: Remote host name.
        userName: User name.
        password: Password.
    Returns: 
        Tuple with standard input, standard output and error output.
    """
    import paramiko
    import select
    ssh = paramiko.SSHClient()
    ssh.load_system_host_keys()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect(hostName, 22, userName, password)
    transport = ssh.get_transport()
    channel = transport.open_session()
    channel.exec_command(command)
    while True:
        if channel.exit_status_ready():
            break
        rl, wl, xl = select.select([channel], [], [], 0.0)
        if len(rl) > 0:
            print channel.recv(1024)


def getLocalUserName():
    """ Recover local machine user name.
    returns: Local machine user name.
    """
    import getpass
    return getpass.getuser()


def getLocalHostName():
    """ Recover local machine name.
    returns: Local machine name.
    """
    import socket
    return socket.gethostname()


def isInFile(text, filePath):
    """ Checks if given text is in the given file.
    params:
        text: Text to check.
        filePath : File path to check.
    returns: True if the given text is in the given file, 
             False if it is not in the file.
    """
    return any(text in line for line in open(filePath))


def getLineInFile(text, fileName):
    """ Find the line where the given text is located in the given file.
    params:
       text: Text to check.
       filePath : File path to check.
    returns: File number where the text was located.
    """
    with open(fileName) as f:
        for i, line in enumerate(f):
            if text in line:
                return i + 1
    return None


#------------- Colored message strings -------------------------------
def getColorStr(text, color, bold=False):
    """ Add ANSI color codes to the string if there is a terminal sys.stdout.
    Params:
     text: text to be colored
     color: red or green
     bold: bold the text
    """
    if envVarOn('SCIPION_SAFE_COLORS') and not sys.stdout.isatty():
        return text
    
    colors = {'gray': 30, 'red': 31, 'green': 32, 'yellow': 33, 'blue': 34, 'magenta': 35, 'cyan': 36}
    attr = [str(colors[color])]
    
    if bold:
        attr.append('1')
    return '\x1b[%sm%s\x1b[0m' % (';'.join(attr), text)

def greenStr(text):
    return getColorStr(text, color='green')

def redStr(text):
    return getColorStr(text, color='red')

def magentaStr(text):
    return getColorStr(text, color='magenta')

def ansi(n, bold=False):
    """Return function that escapes text with ANSI color n."""
    return lambda txt: '\x1b[%d%sm%s\x1b[0m' % (n, ';1' if bold else '', txt)

black, red, green, yellow, blue, magenta, cyan, white = map(ansi, range(30, 38))
blackB, redB, greenB, yellowB, blueB, magentaB, cyanB, whiteB = [
    ansi(i, bold=True) for i in range(30, 38)]


#-------------- Hyper text highlighting ----------------------------
"""
We use a subset of TWiki hyper text conventions.
In particular:
    *some_text* will display some_text in bold
    _some_text_ will display some_text in italic
    Links:
        http://www.link-page.com  -> hyperlink using the url as label
        [[http://www.link-page.com][Link page]] -> hyperlink using "Link page" as label
"""
# Types of recognized styles
HYPER_BOLD = 'bold'
HYPER_ITALIC = 'italic'
HYPER_LINK1 = 'link1'
HYPER_SCIPION_OPEN = 'sci-open'
HYPER_LINK2 = 'link2'
HYPER_ALL = 'all'

# Associated regular expressions
PATTERN_BOLD = "(^|[\s])[*](?P<bold>[^\s*][^*]*[^\s*]|[^\s*])[*]"
#PATTERN_BOLD = r"[\s]+[*]([^\s][^*]+[^\s])[*][\s]+"
PATTERN_ITALIC = "(^|[\s])[_](?P<italic>[^\s_][^_]*[^\s_]|[^\s_])[_]"
#PATTERN_ITALIC = r"[\s]+[_]([^\s][^_]+[^\s])[_][\s]+"
PATTERN_LINK1 = '(?P<link1>http[s]?://(?:[a-zA-Z]|[0-9]|[$-_@.&+]|[!*\(\),]|(?:%[0-9a-fA-F][0-9a-fA-F]))+)'
PATTERN_LINK2 = "[\[]{2}(?P<link2>[^\s][^\]]+[^\s])[\]][\[](?P<link2_label>[^\s][^\]]+[^\s])[\]]{2}"
# __PATTERN_LINK2 should be first since it could contain __PATTERN_LINK1
PATTERN_ALL = '|'.join([PATTERN_BOLD, PATTERN_ITALIC, PATTERN_LINK2, PATTERN_LINK1])

# Compiled regex
# Not need now, each pattern compiled separately
#HYPER_REGEX = {
#               HYPER_BOLD: re.compile(PATTERN_BOLD),
#               HYPER_ITALIC: re.compile(PATTERN_ITALIC),
#               HYPER_LINK1: re.compile(PATTERN_LINK1),
#               HYPER_LINK2: re.compile(PATTERN_LINK1),
#               }
HYPER_ALL_RE = re.compile(PATTERN_ALL)

def parseHyperText(text, matchCallback):
    """ Parse the text recognizing Hyper definitions below.
    Params:
        matchCallback: a callback function to processing each matching,
                       it should accept the type of match (HYPER_BOLD, ITALIC or LINK)
    Return:
        The input text with the replacements made by matchCallback
    """
    def _match(match):
        """ Call the proper matchCallback with some extra info. """
        m = match.group().strip()
        if m.startswith('*'):
            tag = HYPER_BOLD
        elif m.startswith('_'):
            tag = HYPER_ITALIC
        elif m.startswith('http'):
            tag = HYPER_LINK1
        elif m.startswith('[['):
            tag = HYPER_LINK2
        else:
            raise Exception("Bad prefix for HyperText match")
        return matchCallback(match, tag)
        
    return HYPER_ALL_RE.sub(_match, text)
#    for hyperMode, hyperRegex in HYPER_REGEX.iteritems():
#        text = hyperRegex.sub(lambda match: matchCallback(match, hyperMode), text)
#
#    return text

def parseBibTex(bibtexStr):
    """ Parse a bibtex file and return a dictionary. """
    import bibtexparser

    if hasattr(bibtexparser, 'loads'):
        return bibtexparser.loads(bibtexStr).entries_dict

    # For older bibtexparser version 0.5
    from bibtexparser.bparser import BibTexParser
    from StringIO import StringIO

    f = StringIO()
    f.write(bibtexStr)
    f.seek(0, 0)
    parser = BibTexParser(f)

    return parser.get_entry_dict()


def isPower2(num):
    """ Return True if 'num' is a power of 2. """
    return num != 0 and ((num & (num - 1)) == 0)

#---------------------------------------------------------------------------
# Parsing of arguments
#---------------------------------------------------------------------------
def getListFromRangeString(rangeStr):
    """ Create a list of integers from a string with range definitions.
    Examples:
    "1,5-8,10" -> [1,5,6,7,8,10]
    "2,6,9-11" -> [2,6,9,10,11]
    "2 5, 6-8" -> [2,5,6,7,8]
    """
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


def getListFromValues(valuesStr, length=None):
    """ Convert a string representing list items into a list.
    The items should be separated by spaces and a multiplier 'x' can be used.
    If length is not None, then the last element will be repeated
    until the desired length is reached.
    Examples:
    '1 1 2x2 4 4' -> ['1', '1', '2', '2', '4', '4']
    '2x3, 3x4, 1' -> ['3', '3', '4', '4', '4', '1']
    """
    result = []
    
    for chunk in valuesStr.split():
        values = chunk.split('x')
        n = len(values)
        if n == 1: # 'x' is not present in the chunk, single value
            result += values
        elif n == 2: # multiple the values by the number after 'x'
            result += [values[1]] * int(values[0])
        else:
            raise Exception("More than one 'x' is not allowed in list string value.")
            
    # If length is passed, we fill the list with 
    # the last element until length is reached
    if length is not None and length > len(result):
        item = result[-1]
        result += [item] * (length - len(result))
        
    return result
        
    
def getFloatListFromValues(valuesStr, length=None):
    ''' Convert a string to a list of floats'''
    return [float(v) for v in getListFromValues(valuesStr, length)]


def getBoolListFromValues(valuesStr, length=None):
    ''' Convert a string to a list of booleans'''
    from pyworkflow.object import Boolean
    return [Boolean(value=v).get() for v in getListFromValues(valuesStr, length)]


def getStringListFromValues(valuesStr, length=None):
    ''' Convert a string to a list of booleans'''
    from pyworkflow.object import String
    return [String(value=v).get() for v in getListFromValues(valuesStr, length)]


class Environ(dict):
    """ Some utilities to handle environment settings. """
    REPLACE = 0
    BEGIN = 1
    END = 2
    
    def set(self, varName, varValue, position=REPLACE):
        """ Modify the value for some variable.
        Params:
            varName: for example LD_LIBRARY_PATH
            varValue: the value to add or replace.
            position: controls how the value will be changed.
                If REPLACE, it will overwrite the value of
                the var.
                BEGIN or END will preserve the current value
                and add (at begin or end) the new value.
        """
        if varName in self and position != self.REPLACE:
            if position == self.BEGIN:
                self[varName] = varValue + os.pathsep + self[varName]
            elif position == self.END:
                self[varName] = self[varName] + os.pathsep + varValue
        else:
            self[varName] = varValue
                
            
    def update(self, valuesDict, position=REPLACE):
        """ Use set for each key, value pair in valuesDict. """
        for k, v in valuesDict.iteritems():
            self.set(k, v, position)
            
            
def environAdd(varName, newValue, valueFirst=False):
    """ Add a new value to some environ variable.
    If valueFirst is true, the new value will be at the beginning.
    """
    varList = [os.environ[varName]]
    i = 1
    if valueFirst:
        i = 0
    varList.insert(i, newValue)
    os.environ[varName] = os.pathsep.join(varList)


def envVarOn(varName, env=None):
    """ Is variable set to True in the environment? """
    v = env.get(varName) if env else os.environ.get(varName)
    return v is not None and v.lower() in ['true', 'yes', 'on', '1']


def getMemoryAvailable():
    """ Return the total memory of the system in MB """
    from psutil import virtual_memory
    return virtual_memory().total // 1024**2


def startDebugger(mode='SCIPION_DEBUG', password='a'):
    if mode != 'SCIPION_DEBUG' or envVarOn('SCIPION_DEBUG'):
        try:
            from rpdb2 import start_embedded_debugger
            print "Starting debugger..."
            start_embedded_debugger(password)
        except Exception:
            print "Error importing rpdb2 debugging module, consider installing winpdb."


def getFreePort(basePort=0,host=''):
        import socket
        port=0
        try:
            s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            s.bind((host, basePort))
            ipaddr, port = s.getsockname()
            s.close()
        except Exception, e:
            print e
            return 0
        return port
    
    
def readProperties(propsFile):
    myprops = {}
    with open(propsFile, 'r') as f:
        for line in f:
            line = line.rstrip() #removes trailing whitespace and '\n' chars
    
            if "=" not in line: continue #skips blanks and comments w/o =
            if line.startswith("#"): continue #skips comments which contain =
    
            k, v = line.split("=", 1)
            myprops[k] = v
    return myprops


# ---------------------Color utils --------------------------
def hex_to_rgb(value):
    value = value.lstrip('#')
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))


def rgb_to_hex(rgb):
    return '#%02x%02x%02x' % rgb


def lighter(color, percent):
    '''assumes color is rgb between (0, 0, 0) and (255, 255, 255)'''
    color = np.array(color)
    white = np.array([255, 255, 255])
    vector = white - color
    return tuple(np.around(color + vector * percent))

