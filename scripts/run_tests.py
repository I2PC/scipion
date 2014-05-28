#!/usr/bin/env python

import sys
import os
from os.path import join, dirname, exists
import unittest
import pyworkflow as pw
from pyworkflow.tests import *
# TODO: import only what is necessary from pyworkflow.tests, not *


def discoverTests(pathList, pattern):
    tests = unittest.TestSuite()
    for path in pathList:
        testPath = join('pyworkflow', 'tests', path)
        print "discovering '%s'" % testPath
        tests.addTests(unittest.defaultTestLoader.discover(testPath, pattern=pattern, top_level_dir=pw.HOME))
    return tests


lastClass = None
lastModule = None


def printTests(tests, mode='modules'):
    global lastClass
    global lastModule
    
    for t in tests:
        if isinstance(t, unittest.TestSuite):
#             if t.countTestCases():
#                 print indent, "Suite, count=", t.countTestCases(), "class", t.__class__
            printTests(t, mode)
        elif isinstance(t, unittest.TestCase):
            parts = t.id().split('.')
            testName = parts[-1]
            className = parts[-2]
            moduleName = '.'.join(parts[:-2])  
            
            # Check of Failure loading tests
            if moduleName.startswith('unittest.loader.ModuleImportFailure'):
                print failStr(moduleName)
                print "  test: ", t.id()
            else:
                if moduleName != lastModule:
                    lastModule = moduleName
                    print moduleName
                if mode in ['classes', 'all']:
                    if className != lastClass:
                        lastClass = className
                        print "  ", className
                if mode == 'all':
                    print "    ", testName

       
def runTests(tests):
    result = GTestResult()
    tests.run(result)
    result.doReport()
    if result.testFailed > 0:
        sys.exit(1)
    else:
        sys.exit(0)
       
  
PATH = 'path'
CASE = 'case'
PRINT = 'print'
PATTERN = 'pattern'

CASE_DICT = {'model': {PATH:'model em/data', PATTERN: 'test*.py'},
             'xmipp': {PATH:'em/workflows', PATTERN: 'test*xmipp*.py'},
             'mixed': {PATH:'em/workflows', PATTERN: 'test*mixed*.py'},
             'protocols': {PATH:'em/protocols', PATTERN: 'test_protocols*.py'}}
  
def parseArgs(args):
    """ Parse arguments in the form of:
        path="em/data classes"
        case="short"
        print
        And store in a dictionary
    """
    d = {}
    for a in args:
        if '=' in a:
            k, v = a.split('=')
        else:
            k, v = a, True
        d[k.lower()] = v

    return d


if __name__ == '__main__':
    if len(sys.argv) > 1:
        argsDict = parseArgs(sys.argv[1:])
    else:
        argsDict = {}
      
    # CASE and PATH are excusive
    if CASE in argsDict:
        case = argsDict[CASE]
        argsDict.update(CASE_DICT[case])
    
    if PRINT in argsDict:
        defaultPath = '.'
    else:
        defaultPath = 'model'
    path = argsDict.get(PATH, defaultPath)
    pattern = argsDict.get(PATTERN, 'test*.py')    
    pathList = path.split()
    
    tests = discoverTests(pathList, pattern)
    
    if PRINT in argsDict:
        printTests(tests, argsDict.get('mode', 'modules'))
    else:
        runTests(tests)
