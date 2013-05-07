"""
Created on Apr 9, 2013

@author: antonio
"""
import unittest
from utils.log import *

class TestLog(unittest.TestCase):
      
    def testSimpleFileLog(self):
        log = getGeneralLogger('pyworkflow.test.log.test_scipon_log')
        log.info('General info')
        log.debug('General debug')
        log.warning("General warning")
        
        log = getFileLogger('/home/antonio/scipionLog/fileLog.log')
        log.info('File info!!!!!!')
        log.debug('File debug!!!!!!')
        log.warning("File warning!!!")
        
        log = getGeneralLogger('pyworkflow.test.log.test_scipon_log')
        log.error('General error')
        
        log = getFileLogger('/home/antonio/scipionLog/fileLog.log')
        log.error('File error!!!!!!')
        
        log = getClassAndFileLogger('pyworkflow.test.log.test_scipon_log', '/home/antonio/scipionLog/classAndFileLog.log')
        log.info('Class and File info!!!!!!')
        log.debug('Class and File debug!!!!!!')
        log.warning("Class and File warning!!!")
        
        self.assertTrue(True)  
        
    """
    config = {  'version': 1,              
                'disable_existing_loggers': False,
                'formatters': {
                    'standard': {
                        'format': '%(asctime)s [%(levelname)s] %(name)s: %(message)s'
                    },
                },
                'handlers': {
                    'default': {
                        'level':'INFO',    
                        'class':'logging.handlers.RotatingFileHandler',
                        'filename' : "/home/antonio/Documents/testScipion.log",
                        'maxBytes' : 100000,
                    },  
                },
                'loggers': {
                    '': {                  
                        'handlers': ['default'],        
                        'level': 'INFO',  
                        'propagate': True  
                    },
                }
            }
    logging.config.dictConfig(config)    
    """

if __name__ == "__main__":
    #import sys;sys.argv = ["", "Test.testName"]
    unittest.main()