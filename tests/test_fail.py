'''
Created on Jun 12, 2013

@author: antonio
'''
import unittest


class TestFail(unittest.TestCase):


    def testFail(self):
        self.assertTrue(False)


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()