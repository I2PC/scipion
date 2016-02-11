#!/usr/bin/env python
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
Run or show the selected tests. Tests can be selected by giving
the "case", or by giving the paths and file pattern to use for
searching them.
"""

import os
import sys
from os.path import join, basename
import argparse
import unittest

import pyworkflow as pw
import pyworkflow.utils as pwutils
from pyworkflow.tests import GTestResult


PATH_PATTERN = {'model': ('model em/data', 'test*.py'),
                'xmipp': ('em/workflows', 'test*xmipp*.py'),
                'mixed': ('em/workflows', 'test*mixed*.py'),
                'protocols': ('em/protocols', 'test_protocols*.py'),
                'programs': ('em/programs', 'test_program*py')}

MODULE = 0
CLASS = 1
TEST = 2


class Tester():
    def main(self):
        parser = argparse.ArgumentParser(description=__doc__)
        g = parser.add_mutually_exclusive_group()
        g.add_argument('--run', action='store_true', help='run the selected tests')
        g.add_argument('--show', action='store_true', help='show available tests')
        
        add = parser.add_argument  # shortcut
        
        add('--case', choices=['model', 'xmipp', 'mixed', 'protocols'],
            help='get pre-defined paths and patterns for the specified case')
        add('--paths', default='.',
            help='space-separated list of paths where the tests are searched')
        add('--pattern', default='test*.py',
            help='pattern for the files that will be used in the tests')
        add('--grep', default=None, nargs='+',
            help='only show/run tests containing the provided words')
        add('--skip', default=None, nargs='+',
            help='skip tests that contains these words')
        add('--log', default=None, nargs='?',
            help="Generate logs files with the output of each test.")
        add('--mode', default='classes', choices=['modules', 'classes', 'all'],
            help='how much detail to give in show mode')
        add('tests', metavar='TEST', nargs='*',
            help='test case from string identifier (module, class or callable)')
        args = parser.parse_args()
    
        if not args.run and not args.show and not args.tests:
            sys.exit(parser.format_help())
    
        # If 'case' is given, 'path' and 'pattern' are (re-)written from it
        if args.case:
            print 'Using paths and pattern for given case "%s"' % args.case
            args.paths, args.pattern = PATH_PATTERN[args.case]
    
        if args.tests:
            tests = unittest.TestSuite()
            for t in args.tests:
                prefix = ('' if t.startswith('tests.') else 'tests.')
                try:
                    tests.addTests(unittest.defaultTestLoader.loadTestsFromName(
                        '%s%s%s' % ('pyworkflow.', prefix, t)))
                except Exception as e:
                    print 'Cannot find test %s -- skipping' % t
                    print 'error: ', e
        else:
            tests = self.discoverTests(args.paths.split(), args.pattern)
    
        self.grep = args.grep
        self.skip = args.skip
        self.mode = args.mode
        self.log = args.log
        
        if args.show:
            self.printTests(tests)
            
        elif args.run:
            self.runTests(tests)
        
        elif args.tests:
            self.runSingleTest(tests)
    
    def discoverTests(self, paths, pattern):
        """ Return tests discovered in paths that follow the given pattern """
        tests = unittest.TestSuite()
        for path in [join('pyworkflow', 'tests', x) for x in paths]:
            print "Discovering tests in '%s'" % path
            tests.addTests(unittest.defaultTestLoader.discover(
                path, pattern=pattern, top_level_dir=pw.HOME))
        return tests
    
    def _match(self, itemName):
        itemLower = itemName.lower()
        grep = (not self.grep or
                all(g.lower() in itemLower for g in self.grep))
        skip = (self.skip and
                any(g.lower() in itemLower for g in self.skip))
        
        return (grep and not skip)
        
    
    def _visitTests(self, tests, newItemCallback):
        """ Show the list of tests available """
        mode = self.mode
    
        assert mode in ['modules', 'classes', 'all'], 'Unknown mode %s' % mode
    
        # First flatten the list of tests.
        testsFlat = []
        toCheck = [t for t in tests]
        while toCheck:
            test = toCheck.pop()
            if isinstance(test, unittest.TestSuite):
                toCheck += [t for t in test]
            else:
                testsFlat.append(test)
    
        # Follow the flattened list of tests and show the module, class
        # and name, in a nice way.
        lastClass = None
        lastModule = None
        
        for t in testsFlat:
            moduleName, className, testName = t.id().rsplit('.', 2)
    
            # If there is a failure loading the test, show it
            if moduleName.startswith('unittest.loader.ModuleImportFailure'):
                print pwutils.red(moduleName), "  test:", t.id()
                continue
    
            if moduleName != lastModule:
                lastModule = moduleName
                newItemCallback(MODULE, moduleName)
                #print "scipion test %s" % moduleName
                
            if mode in ['classes', 'all'] and className != lastClass:
                lastClass = className
                newItemCallback(CLASS, "%s.%s" % (moduleName, className))
                #print "  scipion test %s.%s" % (moduleName, className)
                
            if mode == 'all':
                newItemCallback(TEST, "%s.%s.%s" % (moduleName, className, testName))
                #print "    scipion test %s.%s.%s" % (moduleName, className, testName)
    
    def _printNewItem(self, itemType, itemName):
        if self._match(itemName):
            spaces = (itemType * 2) * ' '
            print "%s scipion test %s" % (spaces, itemName)
            
    def printTests(self, tests):
        self._visitTests(tests, self._printNewItem)
        
    def _logTest(self, cmd, runTime, result, logFile):
        with open(self.testLog, "r+") as f:
            lines = f.readlines()
            f.seek(0)
            for l in lines:
                if '<!-- LAST_ROW -->' in l:
                    rowStr = "<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>"
                    if result: # != 0 means failed in os.system
                        resultStr = '<font color="red">[FAILED]</font>'
                    else:
                        resultStr = '<font color="green">[SUCCEED]</font>'
                    logStr = '<a href="file://%s">%s</a>' % (logFile, basename(logFile))

                    print >> f, rowStr % (self.testCount, cmd, runTime, resultStr, logStr)
                if self.headerPrefix  in l:
                    f.write(self.headerPrefix + self.testTimer.getToc() + '</h3>\n')
                else:
                    f.write(l)                
            f.close()
        
    def _runNewItem(self, itemType, itemName):
        if self._match(itemName):
            spaces = (itemType * 2) * ' '
            scipion = join(os.environ['SCIPION_HOME'], 'scipion')
            cmd = "%s %s test %s" % (spaces, scipion, itemName)
            run = ((itemType == MODULE and self.mode == 'module') or
                   (itemType == CLASS and self.mode == 'classes') or
                   (itemType == TEST and self.mode == 'all'))
            if run:
                if self.log:
                    logFile = join(self.testsDir, '%s.txt' % itemName)
                    cmdFull = cmd + " > %s 2>&1" % logFile
                else:
                    logFile = ''
                    cmdFull = cmd
                
                print pwutils.green(cmdFull)
                t = pwutils.Timer()
                t.tic()
                self.testCount += 1
                result = os.system(cmdFull)
                if self.log:
                    self._logTest(cmd.replace(scipion, 'scipion'), 
                                  t.getToc(), result, logFile)

    def runTests(self, tests):
        self.testCount = 0

        if self.log:
            self.testsDir = join(os.environ['SCIPION_USER_DATA'], 'Tests', self.log)
            pwutils.cleanPath(self.testsDir)
            pwutils.makePath(self.testsDir)
            self.testLog = join(self.testsDir, 'tests.html')
            self.testTimer = pwutils.Timer()
            self.testTimer.tic()
            self.headerPrefix = '<h3>Test results (%s) Duration: ' % pwutils.prettyTime()
            f = open(self.testLog, 'w')
            f.write("""<!DOCTYPE html>
    <html>
    <body>
    """)
            f.write(self.headerPrefix + '</h3>')
            f.write("""    
     <table style="width:100%" border="1">
      <tr>
        <th>#</th>
        <th>Command</th>
        <th>Time</th>
        <th>Result</th>
        <th>Log file</th>
      </tr>
    <!-- LAST_ROW -->
    </table> 
    
    </body>
    </html>""")

            f.close()
        self._visitTests(tests, self._runNewItem)

        if self.log:
            print "\n\nOpen results in your browser: \nfile:///%s" % self.testLog
        
    def runSingleTest(self, tests):
        result = GTestResult()
        tests.run(result)
        result.doReport()
        if result.testFailed > 0:
            sys.exit(1)
        else:
            sys.exit(0)



if __name__ == '__main__':
    print "Running tests...."
    Tester().main()
