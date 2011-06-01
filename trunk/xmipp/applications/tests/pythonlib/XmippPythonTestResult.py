from unittest import TestResult, _TextTestResult
try:
   from unittest.runner import _WritelnDecorator # Python 2.7+
except ImportError:
   from unittest import _WritelnDecorator # Python <2.6

import sys
from bcolors import *
class XmippPythonTestResult(TestResult):
    xml = None
    testFailed = 0
    numberTests = 0
    def openXmlReport(self, classname, filename):
        self.xml = open(filename, 'w')
        self.xml.write('<testsuite name="%s">\n' % classname)
        
    def closeXmlReport(self):
         self.xml.write('</testsuite>\n')
         self.xml.close()
         sys.stderr.write("%s[==========]%s run %d tests\n" % (bcolors.OKGREEN, bcolors.ENDC, self.numberTests))
         if (self.testFailed):
             sys.stderr.write("%s[  FAILED  ] %s %d tests\n" % (bcolors.FAIL, bcolors.ENDC, self.testFailed))
         else:
             sys.stderr.write("%s[  PASSED  ] %s %d tests\n" % (bcolors.OKGREEN, bcolors.ENDC, self.numberTests - self.testFailed))
             
    
    def getTestNames(self, test):
        parts = str(test).split()
        name = parts[0];
        parts = parts[1].split('.')
        classname = parts[1].replace(")", "")
        self.numberTests += 1 
        return (name, classname)
    
    def addSuccess(self, test):
        name, classname = self.getTestNames(test)
        self.xml.write('   <testcase name="%s" classname="%s"/>\n' % (name, classname))
        sys.stderr.write("%s[ RUN      ]%s %s.%s\n" % (bcolors.OKGREEN, bcolors.ENDC, classname, name))
        sys.stderr.write("%s[      OK  ]%s %s.%s\n" % (bcolors.OKGREEN, bcolors.ENDC, classname, name))
    
    def reportError(self, test, err):
        name, classname = self.getTestNames(test)
        self.xml.write('   <testcase name="%s" classname="%s">\n' % (name, classname))
        self.xml.write('      <failure message=" "/>\n')
        self.xml.write('   </testcase>\n')
        sys.stderr.write("%s[   FAILED ]%s %s.%s\n" % (bcolors.FAIL, bcolors.ENDC, classname, name))
        self.testFailed += 1
                
    def addError(self, test, err):
        self.reportError(test, err)
        
    def addFailure(self, test, err):
        self.reportError(test, err)
                
