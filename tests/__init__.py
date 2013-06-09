#!/usr/bin/env python
# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
#                Laura del Cano         (ldelcano@cnb.csic.es)
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
import sys
import os
from os.path import join, exists, isdir, relpath
from unittest import TestResult
from pyworkflow.utils.utils import getColorStr
import time

try:
   from unittest.runner import _WritelnDecorator # Python 2.7+
except ImportError:
   from unittest import _WritelnDecorator # Python <2.6
from pyworkflow.utils.path import cleanPath, makePath
from pyworkflow.manager import Manager


if "SCIPION_HOME" not in os.environ:
    raise Exception("SCIPION_HOME is not defined as environment variable")

TESTS_HOME = join(os.environ['SCIPION_HOME'], 'tests')

def getInputPath(*filenames):
    """Return the path to the SCIPION_HOME/tests/input dir
    joined with filename"""
    return join(TESTS_HOME, "input", *filenames)

def getGoldPath(*filenames):
    """Return the path to the SCIPION_HOME/tests/gold dir
    joined with filename"""
    return join(TESTS_HOME, "gold", *filenames)

def getOutputPath(*filenames):
    """Return the path to the SCIPION_HOME/tests/output dir
    joined with filename"""
    return join(TESTS_HOME, "output", *filenames)

def getRelPath(filename):
    """Return the path relative to SCIPION_HOME/tests"""
    return relpath(filename, TESTS_HOME)

def setupOutput(test, outputDir):
    """ Define the output path for the calling test and 
    define a function to retrieve output path from this root. 
    """
    test.outputPath = getOutputPath(outputDir)
    cleanPath(test.outputPath)
    
def setupProject(testClass):
    """ Create and setup a project for this test. """
    projName = testClass.__name__
    proj = Manager().createProject(projName) # Now it will be loaded if exists
    # Check that exists hosts for execution
    hosts = proj.getHosts()
    if len(hosts) <= 0:
        raise Exception("Project: %s can't load host configuration." % projName)
    
    testClass.projName = projName
    testClass.proj = proj
    
    
def greenStr(msg):
    return getColorStr(msg, 'green')

def failStr(msg):
    return getColorStr(msg, 'red')


class GTestResult(TestResult):
    """ Subclass TestResult to ouput tests results with colors (green for success and red for failure)
    and write a report on an .xml file. 
    """
    xml = None
    testFailed = 0
    numberTests = 0
    
    def __init__(self):
        TestResult.__init__(self)
        self.startTimeAll = time.time()
    
    def openXmlReport(self, classname, filename):
        #self.xml = open(filename, 'w')
        #self.xml.write('<testsuite name="%s">\n' % classname)
        pass
        
    def doReport(self):
        secs = time.time() - self.startTimeAll
        print >> sys.stderr, greenStr("\n[==========]") + " run %d tests (%0.3f secs)" % (self.numberTests, secs)
        if self.testFailed:
            print >> sys.stderr, failStr("[  FAILED  ]") + " %d tests" % self.testFailed
        print >> sys.stdout, greenStr("[  PASSED  ]") + " %d tests" % (self.numberTests - self.testFailed)
        #self.xml.write('</testsuite>\n')
        #self.xml.close()
             
    def tic(self):
        self.startTime = time.time()
        
    def toc(self):
        return time.time() - self.startTime
        
    def startTest(self, test):
        self.tic()
        self.numberTests += 1         
    
    def getTestName(self, test):
        parts = str(test).split()
        name = parts[0]
        parts = parts[1].split('.')
        classname = parts[-1].replace(")", "")
        return "%s.%s" % (classname, name)
    
    def addSuccess(self, test):
        secs = self.toc()
        sys.stderr.write("%s %s (%0.3f secs)\n" % (greenStr('[ RUN   OK ]'), self.getTestName(test), secs))
    
    def reportError(self, test, err):
        sys.stderr.write("%s %s\n" % (failStr('[   FAILED ]'), self.getTestName(test)))
        self.testFailed += 1
                
    def addError(self, test, err):
        self.reportError(test, err)
        
    def addFailure(self, test, err):
        self.reportError(test, err)
