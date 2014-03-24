


import unittest
import os, sys, time
from os.path import join, exists, isdir, relpath
from unittest import TestResult
from pyworkflow.utils.path import cleanPath, makePath
from pyworkflow.manager import Manager


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
    
    @classmethod
    def getDataSet(cls, name):
        return cls._datasetDict[name]

class BaseTest(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        cls.outputPath = join(TESTS_OUTPUT, cls.__name__)
        cleanPath(cls.outputPath)
        makePath(cls.outputPath)
        
        
    def setUp(self):
        pass
    
   
    
    
    @classmethod
    def getOutputPath(cls, *filenames):
        """Return the path to the SCIPION_HOME/tests/output dir
        joined with filename"""
        return join(cls.outputPath, *filenames)   
    
    @classmethod
    def getRelPath(cls, filename):
        """Return the path relative to SCIPION_HOME/tests"""
        return relpath(filename, TESTS_OUTPUT)

    
    def setupOutput(self, outputDir):
        """ Define the output path for the calling test and 
        define a function to retrieve output path from this root. 
        """
        self.outputPath = self.getOutputPath(outputDir)
        cleanPath(self.outputPath)
        
    
    def setupProject(self):
        """ Create and setup a project for this test. """
        projName = self.__name__
        proj = Manager().createProject(projName) # Now it will be loaded if exists
        # Check that exists hosts for execution
        hosts = proj.getSettings().getHosts()
        if len(hosts) <= 0:
            raise Exception("Project: %s can't load host configuration." % projName)
        
        self.projName = projName
        self.proj = proj
        
    def getTmpPath(self, *filenames):
        """Return the filename in /tmp/ folder.
        If the file exists, it will be deleted"""
        path = self.getTestPath('tmp', *filenames)
        if os.path.exists(path) and not os.path.isdir(path):
            os.remove(path)
        return path
        
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