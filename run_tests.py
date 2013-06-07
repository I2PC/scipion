#!/usr/bin/env python

import os, sys
from os.path import join, dirname, exists
import unittest
import pyworkflow as pw


        
if __name__ == '__main__':
    tests = unittest.defaultTestLoader.discover('tests', top_level_dir=pw.HOME)
    
    if len(sys.argv) > 1:
        arg = sys.argv[1]
        
    arg
    for t in tests:
        print '-'*100, '\n', t
    
    