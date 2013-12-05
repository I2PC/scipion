#!/usr/bin/env python

import os, sys
from os.path import join, dirname, exists
import unittest
import pyworkflow as pw
from pyworkflow.tests import GTestResult


def discoverTests(pathList, pattern):
    tests = unittest.TestSuite()
    for path in pathList:
        testPath = join('pyworkflow','tests', path)
        tests.addTests(unittest.defaultTestLoader.discover(testPath, pattern=pattern, top_level_dir=pw.HOME))
    return tests

def printTests(tests):
    for t in tests:
        print '-'*100, '\n', t
        
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

CASE_DICT = {'classes': {PATH:'classes em/data', PATTERN: 'test*.py'},
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
        argsDict = {PATH: 'classes'}
      
    # CASE and PATH are excusive
    if CASE in argsDict:
        case = argsDict[CASE]
        argsDict.update(CASE_DICT[case])
    
    path = argsDict[PATH]
    pattern = argsDict.get(PATTERN, 'test*.py')    
    pathList = path.split()
    
    tests = discoverTests(pathList, pattern)
    
    if PRINT in argsDict:
        printTests(tests)
    else:
        runTests(tests)
