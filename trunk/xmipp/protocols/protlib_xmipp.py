#!/usr/bin/env python
'''
#/***************************************************************************
# * Authors:     J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *
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
# *  e-mail address 'xmipp@cnb.csic.es'
# ***************************************************************************
 '''

from xmipp import Program
from protlib_utils import runImageJPlugin


class XmippScript():
    ''' This class will serve as wrapper around the XmippProgram class
    to have same facilities from Python scripts'''
    def __init__(self, runWithoutArgs=False):
        self._prog = Program(runWithoutArgs)
        
    def defineParams(self):
        ''' This function should be overwrited by subclasses for 
        define its own parameters'''
        pass
    
    def readParams(self):
        ''' This function should be overwrited by subclasses for 
        and take desired params from command line'''
        pass
    
    def checkParam(self, param):
        return self._prog.checkParam(param)
    
    def getParam(self, param, index=0):
        return self._prog.getParam(param, index)
    
    def getIntParam(self, param, index=0):
        return int(self._prog.getParam(param, index))
    
    def getDoubleParam(self, param, index=0):
        return float(self._prog.getParam(param, index))
    
    def getListParam(self, param):
        return self._prog.getListParam(param)
    
    def addUsageLine(self, line, verbatim=False):
        self._prog.addUsageLine(line, verbatim)

    def addExampleLine(self, line, verbatim=True):
        self._prog.addExampleLine(line, verbatim)
        
    def addParamsLine(self, line):
        self._prog.addParamsLine(line)
    
    def run(self):
        ''' This function should be overwrited by subclasses and
        it the main body of the script'''   
        pass
     
    def tryRun(self):
        ''' This function should be overwrited by subclasses and
        it the main body of the script'''
        try:
            import sys
            self.defineParams()
            doRun = self._prog.read(sys.argv)
            if doRun:
                self.readParams()
                self.run()
        except Exception, e:
            import traceback
            traceback.print_exc(file=sys.stdout)
            
class ScriptPluginIJ(XmippScript):
    def __init__(self, macro):
        XmippScript.__init__(self)
        self.macro = macro
        
    def defineOtherParams(self):
        pass
    
    def readOtherParams(self):
        pass
    
    def defineParams(self):
        self.addParamsLine('  --input <...>                         : Input files to show');
        self.addParamsLine('         alias -i;');
        self.addParamsLine('  [--memory <mem="512m">]              : Memory ammount for JVM');
        self.addParamsLine('         alias -m;');
        self.defineOtherParams()
            
    def readParams(self):
        self.memory = self.getParam('--memory')
        if self.memory == "512m":
            print "No memory size provided. Using default: " + self.memory

        self.args = "-i %s" % ' '.join(self.getListParam('-i'))
        self.readOtherParams()
        
    def run(self):
        runImageJPlugin(self.memory, self.macro, self.args)
            
def getXmippPrograms():
    import os
    from glob import glob
    from protlib_filesystem import getXmippPath
    programs = [os.path.basename(p) for p in glob(os.path.join(getXmippPath(), 'bin', 'xmipp_*'))]
    programs.sort()
    return programs


        