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

import os
from xmipp import Program
from protlib_utils import runImageJPlugin, failStr
from protlib_filesystem import getXmippPath


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
        except Exception:
            import traceback
            traceback.print_exc(file=sys.stderr)
            
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
     
#------------- FUNCTION TO WORK WITH PROGRAMS META-INFORMATION -----------------
class LabelData():
    def __init__(self):
        pass
    
def getXmippLabels():
    labelHeader = getXmippPath(os.path.join('libraries', 'data', 'metadata_label.h'))
    f = open(labelHeader)
    labels = []
    for line in f:
        line = line.strip()
        if line.startswith('MDL::addLabel(MDL_'):
            l = line.find('(')
            r = line.find(')')
            parts = line[l+1:r].split(',')
            labels.append({'name':parts[2].replace('"', ''), 'type':parts[1], 'enum':parts[0]})
    return labels

'''Return the list of Xmipp's programs, taken from from bin/ folder'''     
def getXmippPrograms():
    from glob import glob
    programs = [os.path.basename(p) for p in glob(os.path.join(getXmippPath(), 'bin', 'xmipp_*'))]
    programs.sort()
    return programs

#FIXME: this is only while development
def skipProgram(programName):
    if programName in ['xmipp_sqlite3','xmipp_mpi_steps_runner',
                       'xmipp_angular_commonline',
                    'xmipp_transform_threshold']:
        return True
    for p in ['xmipp_test', 'xmipp_template']:
        if programName.find(p) != -1:
            return True
    return False

def getProgramsDbName():
    return os.path.join(getXmippPath(), 'programs.sqlite')

#Some helper functions
def createProgramsDb(dbName=None):
    from protlib_sql import ProgramDb
    if not dbName:
        dbName = getProgramsDbName()
    db = ProgramDb(dbName)
    print 'Created db with name: %(dbName)s' % locals()
    db.create()
    #Create categories dictionary to classify programs
    #looking in program name and category prefixes
    categories = db.selectCategories()
    categoryDict = {}
    for c in categories:
        prefixes = c['prefixes'].split()
        for p in prefixes:
            categoryDict[p] = c
            
    programs = getXmippPrograms()
    for p in programs:
        p = os.path.basename(p)
        try:
            if not skipProgram(p):
                cmd = [p, "--xmipp_write_definition"]
                if p.find('_mpi') != -1:                    
                    cmd = ['mpirun', '-np', '1'] + cmd
                print ' '.join(cmd)
                from subprocess import Popen, PIPE
                ps = Popen(cmd, stdout=PIPE, stderr=PIPE)
                stderrdata = ps.communicate()[1]
                if stderrdata != '':
                    raise Exception(stderrdata)
                for prefix, category in categoryDict.iteritems():
                    if prefix in p:
                        db.updateProgramCategory(p, category)
                        break
        except Exception, e:
            print failStr("PROGRAM: " + p)
            print failStr("ERROR: " + str(e)) 
    return db

class ProgramKeywordsRank():
    def __init__(self, keywords=None):
        self.keywords = keywords
        self.weights = {'name': 5, 'keywords': 3, 'usage': 1}
        
    def getRank(self, program):
        if not self.keywords:
            return 1 
        rank = 0
        for k in self.keywords:
            for wkey, wvalue in self.weights.iteritems():
                if program[wkey].find(k) != -1:
                    rank += wvalue
        return rank


        