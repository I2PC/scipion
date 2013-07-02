#!/usr/bin/env python

import os, sys
from os.path import join, dirname, exists
import unittest
import pyworkflow as pw
from pyworkflow.tests import GTestResult


def discoverTests(path='.'):
    return unittest.defaultTestLoader.discover(join('pyworkflow','tests', path), top_level_dir=pw.HOME)

def printTests(tests):
    for t in tests:
        print '-'*100, '\n', t
        
def runTests(tests):
    result = GTestResult()
    tests.run(result)
    result.doReport()
       
       
if __name__ == '__main__':
    
    arg = None
    tests = []
    doPrint = False
    
    if len(sys.argv) > 1:
        arg = sys.argv[1].lower()
      
    if arg is None:
        pass  
    elif arg == 'print':
        tests = discoverTests()
        doPrint = True
    else:
        tests = discoverTests(arg)
        if 'print' in sys.argv:
            doPrint = True

    if doPrint:
        printTests(tests)
    else:
        runTests(tests)
