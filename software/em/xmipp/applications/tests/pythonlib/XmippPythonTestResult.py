'''
/***************************************************************************
 *
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/
'''
from unittest import TestResult, _TextTestResult
from protlib_xmipp import greenLowStr, failStr

try:
   from unittest.runner import _WritelnDecorator # Python 2.7+
except ImportError:
   from unittest import _WritelnDecorator # Python <2.6

import sys

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
         print >> sys.stderr, greenLowStr("[==========]") + " run %d tests\n" % self.numberTests
         if (self.testFailed):
             print >> sys.stderr, failStr("[  FAILED  ]") + " %d tests\n" % self.testFailed
         else:
             print >> sys.stderr, greenLowStr("[  PASSED  ]") + " %d tests\n" % (self.numberTests - self.testFailed)
             
    
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
        sys.stderr.write("%s %s.%s\n" % (greenLowStr('[ RUN      ]'), classname, name))
        sys.stderr.write("%s %s.%s\n" % (greenLowStr('[       OK ]'), classname, name))
    
    def reportError(self, test, err):
        name, classname = self.getTestNames(test)
        self.xml.write('   <testcase name="%s" classname="%s">\n' % (name, classname))
        self.xml.write('      <failure message=" "/>\n')
        self.xml.write('   </testcase>\n')
        sys.stderr.write("%s %s.%s\n" % (failStr('[   FAILED ]'), classname, name))
        self.testFailed += 1
                
    def addError(self, test, err):
        self.reportError(test, err)
        
    def addFailure(self, test, err):
        self.reportError(test, err)
                
