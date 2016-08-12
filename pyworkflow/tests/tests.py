
import sys
import os
import time
from traceback import format_exception
import unittest
from os.path import join, relpath
from itertools import izip

from pyworkflow.utils.path import cleanPath, makePath
from pyworkflow.manager import Manager
from pyworkflow.utils.utils import envVarOn, redStr, greenStr
from pyworkflow.object import Object, Float
from pyworkflow.protocol import MODE_RESTART

TESTS_INPUT = join(os.environ['SCIPION_HOME'], 'data', 'tests')
TESTS_OUTPUT = join(os.environ['SCIPION_USER_DATA'], 'Tests')


SMALL = 'small'
PULL_REQUEST = 'pull'
DAILY = 'daily'
WEEKLY = 'weekly'


# Procedure to check if a test class has an attribute called _labels and if so
# then it checks if the class test matches any of the labels in input label parameter.
def hasLabel(TestClass, labels):
    
    # Get _labels attributes in class if any.
    classLabels = getattr(TestClass, '_labels', None)

    # Check if no label in test class.    
    return classLabels is not None and any(l in classLabels for l in labels)


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
        if key in self.filesDict:
            return join(self.path, self.filesDict[key])
        return join(self.path, key)
    
    def getPath(self):
        return self.path
    
    @classmethod
    def getDataSet(cls, name):
        """
        This method is called every time the dataset want to be retreived
        """
        assert name in cls._datasetDict, "Dataset: %s dataset doesn't exist." % name
        folder = cls._datasetDict[name].folder
        if not envVarOn('SCIPION_TEST_NOSYNC'):
            scipion = "%s %s/scipion" % (os.environ['SCIPION_PYTHON'],
                                         os.environ['SCIPION_HOME'])
            command = scipion + " testdata --download " + folder
            print ">>>> " + command
            os.system(command)
        return cls._datasetDict[name]


class BaseTest(unittest.TestCase):
    
    _labels = [WEEKLY]
     
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
        if getattr(prot, '_run', True):
            cls.proj.launchProtocol(prot, wait=True)
        
        if prot.isFailed():
            print "\n>>> ERROR running protocol %s" % prot.getRunName()
            print "    FAILED with error: %s\n" % prot.getErrorMessage()
            raise Exception("ERROR launching protocol.")
        
        if not prot.isFinished():
            print "\n>>> ERROR running protocol %s" % prot.getRunName()
            raise Exception("ERROR: Protocol not finished")
    
    @classmethod    
    def saveProtocol(cls, prot):
        """ Save protocol using cls.proj """
        cls.proj.saveProtocol(prot)   
        
    @classmethod
    def newProtocol(cls, protocolClass, **kwargs):
        """ Create new protocols instances through the project
        and return a newly created protocol of the given class
        """
        # Try to continue from previous execution
        if envVarOn('SCIPION_TEST_CONTINUE'):
            candidates = cls.proj.mapper.selectByClass(protocolClass.__name__)
            if candidates:
                c = candidates[0]
                if c.isFinished():
                    setattr(c, '_run', False)
                else:
                    c.runMode.set(MODE_RESTART)
                return c
        return cls.proj.newProtocol(protocolClass, **kwargs)

    @classmethod
    def compareSets(cls, test, set1, set2):
        """ Iterate the elements of boths sets and check
        that all elements have equal attributes. """
        for item1, item2 in izip(set1, set2):
            areEqual = item1.equalAttributes(item2)
            if not areEqual:
                print "item 1 and item2 are different: "
                item1.printAll()
                item2.printAll()
            test.assertTrue(areEqual)
 
def setupTestOutput(cls):
    """ Create the output folder for a give Test class. """
    cls.outputPath = join(TESTS_OUTPUT, cls.__name__)
    cleanPath(cls.outputPath)
    makePath(cls.outputPath)
       

def setupTestProject(cls):
    """ Create and setup a Project for a give Test class. """
    projName = cls.__name__
    if os.environ.get('SCIPION_TEST_CONTINUE', None) == '1':
        proj = Manager().loadProject(projName)
    else:
        proj = Manager().createProject(projName) # Now it will be loaded if exists
    
    cls.outputPath = proj.path
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
       
    
        
class GTestResult(unittest.TestResult):
    """ Subclass TestResult to output tests results with colors (green for success and red for failure)
    and write a report on an .xml file. 
    """
    xml = None
    testFailed = 0
    numberTests = 0
    
    def __init__(self):
        unittest.TestResult.__init__(self)
        self.startTimeAll = time.time()
    
    def openXmlReport(self, classname, filename):
        #self.xml = open(filename, 'w')
        #self.xml.write('<testsuite name="%s">\n' % classname)
        pass
        
    def doReport(self):
        secs = time.time() - self.startTimeAll
        sys.stderr.write("\n%s run %d tests (%0.3f secs)\n" %
                         (greenStr("[==========]"), self.numberTests, secs))
        if self.testFailed:
            print >> sys.stderr, redStr("[  FAILED  ]") + " %d tests" % self.testFailed
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
        sys.stderr.write("%s %s\n" % (redStr('[   FAILED ]'),
                                      self.getTestName(test)))
        sys.stderr.write("\n%s" % redStr("".join(format_exception(*err))))
        self.testFailed += 1
                
    def addError(self, test, err):
        self.reportError(test, err)
        
    def addFailure(self, test, err):
        self.reportError(test, err)
