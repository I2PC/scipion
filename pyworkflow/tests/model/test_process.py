#!/usr/bin/env python

import os, time
import unittest
from subprocess import Popen
import psutil
from pyworkflow.utils.process import killWithChilds


class TestProccess(unittest.TestCase):
    """ Some tests for utils.process module. """

    def setUp(self):
        pass
        
    def test_Object(self):
        p = Popen('pw_sleep.py 500', shell=True)
        print "pid: ", p.pid
        time.sleep(5)
        killWithChilds(p.pid)

        

 
        
if __name__ == '__main__':
    unittest.main()
    