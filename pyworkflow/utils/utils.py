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
This module contains utilities functions and classes.
"""

import os, sys, re


def prettyDate(time=False):
    """
    Get a datetime object or a int() Epoch timestamp and return a
    pretty string like 'an hour ago', 'Yesterday', '3 months ago',
    'just now', etc
    """
    from datetime import datetime
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
    f = open(filePath, 'r')
    for line in f:
        if text in line:
            f.close()
            return True
    f.close()
    return False

def getLineInFile(text, fileName):
    """ Find the line where the given text is located in the given file.
    params:
       text: Text to check.
       filePath : File path to check.
    returns: File number where the text was located.
    """
    f = open(fileName, 'r')
    cont = 0
    for line in f:
        cont += 1
        if text in line:
            f.close()
            return cont
    f.close()
    return None

#------------- Colored message strings -------------------------------
def getColorStr(text, color, bold=False):
    """ Add ANSI color codes to the string if there is a terminal sys.stdout.
    Params:
     text: text to be colored
     color: red or green
     bold: bold the text
    """
    if not sys.stdout.isatty():
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


#-------------- Hyper text highlighting ----------------------------
"""
We use a subset of TWiki hyper text convetions.
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
HYPER_LINK2 = 'link2'
HYPER_ALL = 'all'

# Associated regular expressions
PATTERN_BOLD = "(^|[\s])[*](?P<bold>[^\s*][^*]*[^\s*]|[^\s*])[*]"
#PATTERN_BOLD = r"[\s]+[*]([^\s][^*]+[^\s])[*][\s]+"
PATTERN_ITALIC = "(^|[\s])[_](?P<italic>[^\s_][^_]*[^\s_]|[^\s_])[_]"
#PATTERN_ITALIC = r"[\s]+[_]([^\s][^_]+[^\s])[_][\s]+"
PATTERN_LINK1 = '(?P<link1>http[s]?://(?:[a-zA-Z]|[0-9]|[$-_@.&+]|[!*\(\),]|(?:%[0-9a-fA-F][0-9a-fA-F]))+)'
PATTERN_LINK2 = "[\[]{2}(?P<link2>[^\s][^\]]+[^\s])[\]][\[](?P<link2_label>[^\s][^\]]+[^\s])[\]]{2}"
# __PATTERN_LINK2 should be first since could contains __PATTERN_LINK1
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
