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
This sub-package will contains Spider protocols
"""
import os
from os.path import join, dirname, abspath
import subprocess
import re

from pyworkflow.object import String
from pyworkflow.utils import runJob, Environ
from pyworkflow.utils.path import replaceBaseExt, removeBaseExt
from pyworkflow.em.data import EMObject


END_HEADER = 'END BATCH HEADER'
SPIDER = 'spider_linux_mp_intel64'

PATH = abspath(dirname(__file__))
TEMPLATE_DIR = 'templates'
SCRIPTS_DIR = 'scripts'

# Regular expressions for parsing vars in scripts header

# Match strings of the type
# key = value ; some comment
REGEX_KEYVALUE = re.compile("(?P<var>\[?[a-zA-Z0-9_-]+\]?)(?P<s1>\s*)=(?P<s2>\s*)(?P<value>\S+)(?P<rest>\s+.*)")
# Match strings of the type [key]value
# just before a 'fr l' line
REGEX_KEYFRL = re.compile("(?P<var>\[?[a-zA-Z0-9_-]+\]?)(?P<value>\S+)(?P<rest>\s+.*)")


def getEnviron():
    """ Load the environment variables needed for Spider.
    If SPIDER_DIR is defined, the bin, man and proc folders will be 
    defined from it. If not, each of them should be defined separately. 
    """
    global SPIDER
    env = Environ(os.environ)
    SPIDER_DIR = env.get('SPIDER_DIR', None) # Scipion definition
    
    if SPIDER_DIR is None:
        errors = ''
        for var in ['SPBIN_DIR', 'SPMAN_DIR', 'SPPROC_DIR']:
            if not var in env:
                errors += "\n   Missing SPIDER variable: '%s'" % var
        if len(errors):
            print "ERRORS: " + errors
    else: 
        env.update({'SPBIN_DIR': join(SPIDER_DIR, 'bin', ''),
                    'SPMAN_DIR': join(SPIDER_DIR, 'man', ''),
                    'SPPROC_DIR': join(SPIDER_DIR, 'proc', '')
                    })
    
    # Get the executable or 'spider' by default
    SPIDER = join(env['SPBIN_DIR'], env.get('SPIDER', 'spider_linux_mp_intel64'))
    # expand ~ and vars
    SPIDER = abspath(os.path.expanduser(os.path.expandvars(SPIDER)))
    # Check that executable exists
    if not os.path.exists(SPIDER):
        msg = "SPIDER executable not found at:\n   '%s'" % SPIDER
        msg += "\nPlease create a link inside the bin folder: \n   '%s'" % env['SPBIN_DIR']
        msg += "\n named 'spider' or define the SPIDER environment variable"
        print msg
        return None
        
    env.set('PATH', env['SPBIN_DIR'], env.END)
    
    return env
    

environment = getEnviron()


def _getFile(*paths):
    return join(PATH, *paths)


def getScript(*paths):
    return _getFile(SCRIPTS_DIR, *paths)


def __substituteVar(match, paramsDict, lineTemplate):
    if match and match.groupdict()['var'] in paramsDict:
        d = match.groupdict()
        d['value'] = paramsDict[d['var']]
        return lineTemplate % d
    return None
    
    
def runScript(inputScript, ext, paramsDict, log=None, cwd=None):
    """ This function will create a valid Spider script
    by copying the template and replacing the values in dictionary.
    After the new file is read, the Spider interpreter is invoked.
    Usually the execution should be done where the results will
    be left.
    """
    outputScript = replaceBaseExt(inputScript, ext)
    
    if cwd is not None:
        outputScript = join(cwd, outputScript)

    fIn = open(getScript(inputScript), 'r')
    fOut = open(outputScript, 'w')
    inHeader = True # After the end of header, not more value replacement
    inFrL = False
    
    for i, line in enumerate(fIn):
        if END_HEADER in line:
            inHeader = False
        if inHeader:
            try:
                newLine = __substituteVar(REGEX_KEYVALUE.match(line), paramsDict, 
                                          "%(var)s%(s1)s=%(s2)s%(value)s%(rest)s\n")
                if newLine is None and inFrL:
                    newLine = __substituteVar(REGEX_KEYFRL.match(line), paramsDict, 
                                              "%(var)s%(value)s%(rest)s\n")
                if newLine:
                    line = newLine
            except Exception, ex:
                print ex, "on line (%d): %s" % (i+1, line)
                raise ex
            inFrL = line.lower().startswith("fr ")
        fOut.write(line)
    fIn.close()
    fOut.close()    

    scriptName = removeBaseExt(outputScript)
    args = " %s @%s" % (ext, scriptName)
    
    runJob(log, SPIDER, args, env=dict(environment), cwd=cwd)
    

def runCustomMaskScript(filterRadius1, sdFactor,
                        filterRadius2, maskThreshold,
                        workingDir, ext='stk',
                        inputImage='input_image',
                        outputMask='stkmask'):
    """ Utility function to run the custommask.msa script.
    This function will be called from the custom mask protocol
    and from the wizards to create the mask.
    """
    params = {'[filter-radius1]': filterRadius1,
              '[sd-factor]': sdFactor,
              '[filter-radius2]': filterRadius2,
              '[mask-threshold2]': maskThreshold,
              '[input_image]': inputImage,
              '[output_mask]': outputMask,
              } 
    # Run the script with the given parameters
    runScript('mda/custommask.msa', ext, params, cwd=workingDir)
    
    
class SpiderShell(object):
    """ This class will open a child process running Spider interpreter
    and will keep conection to send commands. 
    """
    def __init__(self, ext='spi', **kwargs):
        self._debug = kwargs.get('debug', True)
        self._log = kwargs.get('log', None)
        cwd = kwargs.get('cwd', None)
        
        FNULL = open(os.devnull, 'w')
        self._proc = subprocess.Popen(SPIDER, shell=True, 
                                      stdin=subprocess.PIPE,
                                      stdout=FNULL, stderr=FNULL,
                                      env=environment,
                                      cwd=cwd)
        if self._debug and self._log:
            self._log = open(self._log, 'w+')
            
        self.runCmd(ext)
        
    def runFunction(self, funcName, *args):
        cmd = funcName
        for a in args:
            cmd += '\n' + str(a)
        self.runCmd(cmd)
    
    def runCmd(self, cmd):
        if self._debug:
            #print "SPIDER: ", cmd
            print >> self._log, cmd
        print >> self._proc.stdin, cmd
        self._proc.stdin.flush()
        
    def close(self, end=True):
        if end:
            self.runCmd("end")
        self._proc.wait()
        # self._proc.kill() TODO: Check if necesary
        

class SpiderDocFile(object):
    """ Handler class to read/write spider docfile. """
    def __init__(self, filename, mode='r'):
        self._file = open(filename, mode)
        self._count = 0
        
    def writeValues(self, *values):
        """ Write values in spider docfile. """
        self._count += 1
            # write data lines
        line = "%5d %2d" % (self._count, len(values))
        for v in values:
            line += " %11g" % float(v)
            
        print >> self._file, line
        
    def iterValues(self):
        for line in self._file:
            line = line.strip()
            if not line.startswith(';'):
                values = [float(s) for s in line.split()[2:]]
                yield values

    def close(self):
        self._file.close()
        
        
class PcaFile(EMObject):
    """ This is a container of files produced by CA PCA Spider protocol.
    It is possible to use the cas_IMC or cas_SEQ files.
    """
    def __init__(self, **args):
        EMObject.__init__(self, **args)
        
        self.filename = String()
        
    def getFileName(self):
        return self.filename.get()
        
     
def getDocsLink(op, label):
    """ Return a label for documentation url of a given command. """
    from constants import SPIDER_DOCS
    return '[[%(SPIDER_DOCS)s/%(op)s.html][%(label)s]]' % locals()
    