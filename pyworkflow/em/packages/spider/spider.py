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
from os.path import join, dirname, abspath, exists, basename
from pyworkflow.utils.path import copyFile, removeExt, replaceExt
from pyworkflow.utils import runJob
import subprocess

END_HEADER = 'END BATCH HEADER'


def loadEnvironment():
    """ Load the environment variables needed for use EMAN2 tools. """
    SPIDER_DIR = os.environ['SPIDER_DIR']
    
    os.environ['SPBIN_DIR'] = join(SPIDER_DIR, 'bin', '')
    os.environ['SPMAN_DIR'] = join(SPIDER_DIR, 'man', '')
    os.environ['SPPROC_DIR'] = join(SPIDER_DIR, 'proc', '')
    
    os.environ['PATH'] = os.environ['PATH'] + os.pathsep + os.environ['SPBIN_DIR']
    
    

TEMPLATE_DIR = abspath(join(dirname(__file__), 'templates'))
        
def getTemplate(templateName):
    """ Return the path to the template file given its name. """
    templateFile = join(TEMPLATE_DIR, templateName)
    if not exists(templateFile):
        raise Exception("getTemplate: template '%s' was not found in templates directory" % templateName)
    
    return templateFile

def copyTemplate(templateName, destDir):
    """ Copy a template file to a diretory """
    template = getTemplate(templateName)
    templateDest = join(destDir, basename(template))
    copyFile(template, templateDest)
    
def runSpiderTemplate(templateName, ext, paramsDict):
    """ This function will create a valid Spider script
    by copying the template and replacing the values in dictionary.
    After the new file is read, the Spider interpreter is invoked.
    """
    loadEnvironment()
    copyTemplate(templateName, '.')
    scriptName = replaceExt(templateName, ext)
    print "scriptName:", scriptName
    
    fIn = open(templateName, 'r')
    fOut = open(scriptName, 'w')
    replace = True # After the end of header, not more value replacement
    
    for line in fIn:
        if END_HEADER in line:
            replace = False
        if replace:
            line = line % paramsDict
        fOut.write(line)
    fIn.close()
    fOut.close()    

    scriptName = removeExt(scriptName)  
    runJob(None, "spider", "%(ext)s @%(scriptName)s" % locals())


class SpiderShell(object):
    """ This class will open a child process running Spider interpreter
    and will keep conection to send commands. 
    """
    def __init__(self, ext='spi', **args):
        self._debug = args.get('debug', True)
        self._log = args.get('log', None)
        
        loadEnvironment()
        FNULL = open(os.devnull, 'w')
        self._proc = subprocess.Popen("spider", shell=True, 
                                      stdin=subprocess.PIPE,
                                      stdout=FNULL, stderr=FNULL)
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
            print "SPIDER: ", cmd
            print >> self._log, cmd
        print >> self._proc.stdin, cmd
        self._proc.stdin.flush()
        
    def runScript(self, templateName, paramsDict):
        """ Run all lines in the template script after replace the params
        with their values.
        """
        templateFile = getTemplate(templateName)        
        f = open(templateFile, 'r')
        replace = True # After the end of header, not more value replacement
        
        for line in f:
            line = line.strip()
            
            if END_HEADER in line:
                replace = False
            
            if not line.startswith(';'): # Skip comment lines
                try:
                    if replace:
                        line = line % paramsDict
                except Exception, ex:
                    print ex, "on line: ", line
            self.runCmd(line)
        
        f.close()
        if self._debug and self._log:
            self._log.close()
        
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

    def close(self):
        self._file.close()
        
        
