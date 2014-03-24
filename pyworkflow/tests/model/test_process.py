#!/usr/bin/env python

import os, time
import unittest
from subprocess import Popen
from pyworkflow.utils.process import killWithChilds
from pyworkflow.tests import *


#FIXME:Jose Miguel
class TestProccess(BaseTest):
    """ Some tests for utils.process module. """

    @classmethod
    def setUpClass(cls):
        setupTestOutput(cls)
        
    def test_Process(self):
        p = Popen('pw_sleep.py 500', shell=True)
        print "pid: ", p.pid
        time.sleep(5)
        killWithChilds(p.pid)

        

 
        
if __name__ == '__main__':
    unittest.main()
    