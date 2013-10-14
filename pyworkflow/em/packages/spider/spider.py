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
from os.path import join
import subprocess


def loadEnvironment():
    """ Load the environment variables needed for use EMAN2 tools. """
    # TODO: Read SPIDER_HOME from the host config.
    SPIDER_HOME = os.environ['SPIDER_HOME']
    SPIDER_DIR = join(SPIDER_HOME, "spider")
    
    os.environ['SPBIN_DIR'] = join(SPIDER_DIR, 'bin', '')
    os.environ['SPMAN_DIR'] = join(SPIDER_DIR, 'man', '')
    os.environ['SPPROC_DIR'] = join(SPIDER_DIR, 'proc', '')
    
    os.environ['PATH'] = os.environ['PATH'] + os.pathsep + os.environ['SPBIN_DIR']
    
    
class SpiderShell(object):
    """ This class will open a child process running Spider interpreter
    and will keep conection to send commands. 
    """
    def __init__(self, ext='spi', **args):
        self._debug = args.get('debug', True)
        loadEnvironment()
        FNULL = open(os.devnull, 'w')
        self._proc = subprocess.Popen("spider", shell=True, 
                                      stdin=subprocess.PIPE,
                                      stdout=FNULL, stderr=FNULL)
        self.runCmd(ext)
        
    def runFunction(self, funcName, *args):
        cmd = funcName
        for a in args:
            cmd += '\n' + str(a)
        self.runCmd(cmd)
    
    def runCmd(self, cmd):
        if self._debug:
            print "SPIDER: ", cmd
        print >> self._proc.stdin, cmd
        self._proc.stdin.flush()
        
    def close(self):
        self.runCmd("end")
        # self._proc.kill() TODO: Check if necesary
