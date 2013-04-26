"""
Created on Apr 9, 2013

@author: antonio
"""
import unittest
import logging
import logging.config


class TestScipionLog(unittest.TestCase):
        
    def testSimpleScipionLog(self):
        logging.config.fileConfig('/home/antonio/Desarrollo/Projects/EclipseProjects/Scipion/pyworkflow/settings/log.config')
        log = logging.getLogger('pyworkflow.tests.utils.testScipionLog')        
        log.info('INNNNFOOOOOO!!!!!!')
        log.debug('DEEBBBUUUUGGGGG!!!!!!')
        log.warning("WARRNINNNGGG!!!")
        self.assertTrue(True)
    


if __name__ == "__main__":
    #import sys;sys.argv = ["", "Test.testName"]
    unittest.main()