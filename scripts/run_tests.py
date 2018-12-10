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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

"""
Run or show the selected tests. Tests can be selected by giving
the "case", or by giving the paths and file pattern to use for
searching them.
"""

from __future__ import print_function

import os
import sys
from os.path import basename
import argparse
import unittest
from collections import OrderedDict

import pyworkflow as pw
import pyworkflow.utils as pwutils
import pyworkflow.tests as pwtests
import pyworkflow.em as pwem


from pyworkflow.tests import *


MODULE = 0
CLASS = 1
TEST = 2


class Tester():
    def main(self):

        parser = argparse.ArgumentParser(description=__doc__)
        g = parser.add_mutually_exclusive_group()
        g.add_argument('--run', action='store_true', help='run the selected tests')
        g.add_argument('--show', action='store_true', help='show available tests',
                       default=True)

        add = parser.add_argument  # shortcut

        add('--pattern', default='test*.py',
            help='pattern for the files that will be used in the tests')
        add('--grep', default=None, nargs='+',
            help='only show/run tests containing the provided words')
        add('--skip', default=None, nargs='+',
            help='skip tests that contains these words')
        add('--log', default=None, nargs='?',
            help="Generate logs files with the output of each test.")
        add('--mode', default='classes', choices=['modules', 'classes', 'onlyclasses', 'all'],
            help='how much detail to give in show mode')
        add('tests', metavar='TEST', nargs='*',
            help='test case from string identifier (module, class or callable)')
        args = parser.parse_args()

        if not args.run and not args.show and not args.tests:
            sys.exit(parser.format_help())

        testsDict = OrderedDict()
        testLoader = unittest.defaultTestLoader

        if args.tests:
            # In this case the tests are passed as argument.
            # The full name of the test will be used to load it.
            testsDict['tests'] = unittest.TestSuite()
            for t in args.tests:
                try:
                    testsDict['tests'].addTests(testLoader.loadTestsFromName(t))
                except Exception as e:
                    print('Cannot load test %s -- skipping' % t)
                    print('error: ', e)
        else:
            # In this other case, we will load the test available
            # from pyworkflow and the other plugins
            self.paths = [('pyworkflow', '.')]
            for name, plugin in pwem.Domain.getPlugins().iteritems():
                self.paths.append((name, os.path.dirname(plugin.__path__[0])))
            for k, p in self.paths:
                testsDict[k] = testLoader.discover(os.path.join(p, k),
                                                   pattern=args.pattern,
                                                   top_level_dir=p)

        self.grep = [g.lower() for g in args.grep] if args.grep else []
        self.skip = args.skip
        self.mode = args.mode
        self.log = args.log

        if args.tests:
            self.runSingleTest(testsDict['tests'])

        elif args.run:
            for moduleName, tests in testsDict.iteritems():
                print(">>>> %s" % moduleName)
                self.runTests(moduleName, tests)

        elif args.grep:
            pattern = args.grep[0]
            for moduleName, tests in testsDict.iteritems():
                self.printTests(pattern, tests)

        else:
            for moduleName, tests in testsDict.iteritems():
                if self._match(moduleName):
                    print(">>>> %s" % moduleName)
                    self.printTests(moduleName, tests)

    def _match(self, itemName):
        itemLower = itemName.lower()
        grep = (not self.grep or
                all(g.lower() in itemLower for g in self.grep))
        skip = (self.skip and
                any(g.lower() in itemLower for g in self.skip))

        return (grep and not skip)

    def __iterTests(self, test):
        """ Recursively iterate over a testsuite. """
        print("__iterTests: %s, is-suite: %s" % (str(test.__class__),
                                                 isinstance(test, unittest.TestSuite)))

        if isinstance(test, unittest.TestSuite):
            for t in test:
                self.__iterTests(t)
        else:
            yield test

    def _visitTests(self, moduleName, tests, newItemCallback):
        """ Show the list of tests available """
        mode = self.mode

        assert mode in ['modules', 'classes', 'onlyclasses', 'all'], 'Unknown mode %s' % mode

        # First flatten the list of tests.
        # testsFlat = list(iter(self.__iterTests(tests)))

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
        if testsFlat:
            for t in testsFlat:

                testModuleName, className, testName = t.id().rsplit('.', 2)

                # If there is a failure loading the test, show it
                errorStr = 'unittest.loader.ModuleImportFailure.'
                if testModuleName.startswith(errorStr):
                    newName = t.id().replace(errorStr, '')
                    if self._match(newName):
                        print(pwutils.red('Error loading the test. Please, run the test for more information:'), newName)
                    continue

                if testModuleName != lastModule:
                    lastModule = testModuleName
                    if mode != 'onlyclasses':
                        newItemCallback(MODULE, "%s" % testModuleName)

                if mode in ['classes', 'onlyclasses', 'all'] and className != lastClass:
                    lastClass = className
                    newItemCallback(CLASS, "%s.%s"
                                    % (testModuleName, className))

                if mode == 'all':
                    newItemCallback(TEST, "%s.%s.%s"
                                    % (testModuleName, className, testName))
        else:
            if not self.grep:
                print(pwutils.green(' The plugin does not have any test'))

    def _printNewItem(self, itemType, itemName):
        if self._match(itemName):
            spaces = (itemType * 2) * ' '
            print("%s scipion test %s" % (spaces, itemName))

    def printTests(self, moduleName, tests):
        self._visitTests(moduleName, tests, self._printNewItem)

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
            scipion = pw.getScipionScript()
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

                print(pwutils.green(cmdFull))
                t = pwutils.Timer()
                t.tic()
                self.testCount += 1
                result = os.system(cmdFull)
                if self.log:
                    self._logTest(cmd.replace(scipion, 'scipion'),
                                  t.getToc(), result, logFile)

    def runTests(self, moduleName, tests):
        self.testCount = 0

        if self.log:
            self.testsDir = join(pw.Config.SCIPION_USER_DATA, 'Tests', self.log)
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
        self._visitTests(moduleName, tests, self._runNewItem)

        if self.log:
            print("\n\nOpen results in your browser: \nfile:///%s"
                  % self.testLog)

    def runSingleTest(self, tests):
        result = pwtests.GTestResult()
        tests.run(result)
        result.doReport()
        sys.exit(1 if result.testFailed > 0 else 0)


if __name__ == '__main__':
    print("Running tests....")
    Tester().main()
