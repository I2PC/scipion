#!/usr/bin/env python

"""
Launch tests.
"""

import sys
from os.path import join
import argparse
import unittest

import pyworkflow as pw
from pyworkflow.tests import failStr, GTestResult


PATH_PATTERN = {'model': ('model em/data', 'test*.py'),
                'xmipp': ('em/workflows', 'test*xmipp*.py'),
                'mixed': ('em/workflows', 'test*mixed*.py'),
                'protocols': ('em/protocols', 'test_protocols*.py')}


def main():
    # First a hack (add --) so we use argparse in spite of our peculiar syntax.
    sys.argv[1:] = [x if x[0] == '-' else '--'+x for x in sys.argv[1:]]
    # Also, if no parameters given at all, do as before, run only paths=model
    if not sys.argv[1:]:
        print 'No explicit paths given, running only tests in paths=model'
        sys.argv[1:] = ['--paths=model']

    parser = argparse.ArgumentParser(description=__doc__)
    add = parser.add_argument  # shortcut
    add('--case', choices=['model', 'xmipp', 'mixed', 'protocols'],
        help='get pre-defined paths and patterns for the specified case')
    add('--paths', default='.',
        help='space-separated list of paths where the tests are searched')
    add('--pattern', default='test*.py',
        help='pattern for the files that will be used in the tests')
    add('--show', action='store_true', help='show the tests available')
    add('--mode', default='modules', choices=['modules', 'classes', 'all'],
        help='how much detail to give in show mode')
    args = parser.parse_args()

    # If 'case' is given, 'path' and 'pattern' are (re-)written from it
    if args.case:
        print 'Using paths and pattern for given case "%s"' % args.case
        args.paths, args.pattern = PATH_PATTERN[args.case]

    tests = discoverTests(args.paths.split(), args.pattern)

    if args.show:
        printTests(tests, args.mode)
    else:
        runTests(tests)


def discoverTests(paths, pattern):
    """ Return tests discovered in paths that follow the given pattern """

    tests = unittest.TestSuite()
    for path in [join('pyworkflow', 'tests', x) for x in paths]:
        print "Discovering tests in '%s'" % path
        tests.addTests(unittest.defaultTestLoader.discover(
            path, pattern=pattern, top_level_dir=pw.HOME))
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
            print failStr(moduleName), "  test:", t.id()
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



if __name__ == '__main__':
    main()
