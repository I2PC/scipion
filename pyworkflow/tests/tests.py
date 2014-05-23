


import unittest
import os, sys, time
from os.path import join, exists, isdir, relpath
from unittest import TestResult
from pyworkflow.utils.path import cleanPath, makePath
from pyworkflow.manager import Manager
from pyworkflow.utils.utils import getColorStr
from pyworkflow.object import *
from pyworkflow.protocol import *

TESTS_INPUT = join(os.environ['SCIPION_HOME'], 'data', 'tests')
TESTS_OUTPUT = join(os.environ['SCIPION_USER_DATA'], 'Tests')


class DataSet:

    _datasetDict = {} # store all created datasets

    def __init__(self, name, folder, files):
        """ 
        Params:
            
        #filesDict is dict with key, value pairs for each file
        """
        self._datasetDict[name] = self
        self.folder = folder
        self.path = join(TESTS_INPUT, folder)
        self.filesDict = files
        
        
    def getFile(self, key):
        return join(self.path, self.filesDict[key])
    
    def getPath(self):
        return self.path
    
    @classmethod
    def getDataSet(cls, name):
        """
        This method is called everytime the dataset want to be retreived
        """
        assert name in cls._datasetDict, "Dataset: %s dataset doesn't exist." % name
        folder = cls._datasetDict[name].folder
        if not 'SCIPION_DATA_NOSYNC' in os.environ:
            command = "%s %s/scipion testdata %s download" % (os.environ['SCIPION_PYTHON'],
                                                              os.environ['SCIPION_HOME'],
                                                              folder)
            print ">>>> " + command
            os.system(command)
        return cls._datasetDict[name]


class BaseTest(unittest.TestCase):
    
    @classmethod
    def getOutputPath(cls, *filenames):
        """Return the path to the SCIPION_HOME/tests/output dir
        joined with filename"""
        return join(cls.outputPath, *filenames)   
    
    @classmethod
    def getRelPath(cls, basedir, filename):
        """Return the path relative to SCIPION_HOME/tests"""
        return relpath(filename, basedir)
    
    @classmethod
    def launchProtocol(cls, prot):
        """ Launch a given protocol using cls.proj and the
        flag wait=True.
        """
        cls.proj.launchProtocol(prot, wait=True)
        
        if prot.isFailed():
            print "\n>>> ERROR running protocol %s" % prot.getRunName()
            print "    FAILED with error: %s\n" % prot.getErrorMessage()
            raise Exception("ERROR launching protocol.")
    
    @classmethod    
    def saveProtocol(cls, prot):
        """ Save protocol using cls.proj """
        cls.proj.saveProtocol(prot)   
        
    @classmethod
    def newProtocol(cls, protocolClass, **kwargs):
        """ Create new protocols instances throught the project
        and return a newly created protocol of the given class
        """
        return cls.proj.newProtocol(protocolClass, **kwargs)

 
def setupTestOutput(cls):
    """ Create the output folder for a give Test class. """
    cls.outputPath = join(TESTS_OUTPUT, cls.__name__)
    cleanPath(cls.outputPath)
    makePath(cls.outputPath)
       

def setupTestProject(cls):
    """ Create and setup a Project for a give Test class. """
    projName = cls.__name__
    proj = Manager().createProject(projName) # Now it will be loaded if exists
    # Check that exists hosts for execution
    hosts = proj.getSettings().getHosts()
    if len(hosts) <= 0:
        raise Exception("Project: %s can't load host configuration." % projName)
    
    cls.projName = projName
    cls.proj = proj
        
        
#class for tests
class Complex(Object):
    
    cGold = complex(1.0, 1.0)
    
    def __init__(self, imag=0., real=0., **args):
        Object.__init__(self, **args)
        self.imag = Float(imag)
        self.real = Float(real)
        # Create reference complex values
        
    def __str__(self):
        return '(%s, %s)' % (self.imag, self.real)
    
    def __eq__(self, other):
        return (self.imag == other.imag and 
                self.real == other.real)
            
    def hasValue(self):
        return True
    
    @classmethod
    def createComplex(cls):
        """Create a Complex object and set
        values with cls.cGold standard"""
        c = Complex() # Create Complex object and set values
        c.imag.set(cls.cGold.imag)
        c.real.set(cls.cGold.real)
        return c
    
        
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
        

def greenStr(msg):
    return getColorStr(msg, 'green')


def failStr(msg):
    return getColorStr(msg, 'red')
