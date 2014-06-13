#!/usr/bin/env python

import sys
from os.path import join
import unittest

import pyworkflow as pw
from pyworkflow.tests import failStr, GTestResult


CASE = {
    'model': {'path': 'model em/data', 'pattern': 'test*.py'},
    'xmipp': {'path': 'em/workflows', 'pattern': 'test*xmipp*.py'},
    'mixed': {'path': 'em/workflows', 'pattern': 'test*mixed*.py'},
    'protocols': {'path': 'em/protocols', 'pattern': 'test_protocols*.py'}}


def main():
    args = parseArgs(sys.argv[1:])

    # If 'case' is given, 'path' and 'pattern' are (re-)written from it
    if 'case' in args:
        args.update(CASE[args['case']])

    defaultPath = 'model' if not 'print' in args else '.'

    tests = discoverTests(args.get('path', defaultPath).split(),
                          args.get('pattern', 'test*.py'))

    if 'print' in args:
        printTests(tests, args.get('mode', 'modules'))
    else:
        runTests(tests)


def discoverTests(pathList, pattern):
    """ Return tests discovered in pathList that follow the given pattern """

    tests = unittest.TestSuite()
    for path in pathList:
        testPath = join('pyworkflow', 'tests', path)
        print "Discovering tests in '%s'" % testPath
        tests.addTests(unittest.defaultTestLoader.discover(
            testPath, pattern=pattern, top_level_dir=pw.HOME))
    return tests


def printTests(tests, mode='modules'):
    """ Show the list of tests available """

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
            print failStr(moduleName), "  test: ", t.id()
            continue

        if moduleName != lastModule:
            lastModule = moduleName
            print moduleName
        if mode in ['classes', 'all'] and className != lastClass:
            lastClass = className
            print "  %s" % className
        if mode == 'all':
            print "    %s" % testName


def runTests(tests):
    result = GTestResult()
    tests.run(result)
    result.doReport()
    if result.testFailed > 0:
        sys.exit(1)
    else:
        sys.exit(0)


def parseArgs(args):
    """ Parse arguments in the form of:
          path="em/data classes"
          case="short"
          print
        And store them in a dictionary
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
    main()
