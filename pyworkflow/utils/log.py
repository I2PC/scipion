'''
Created on Apr 30, 2013

@author: antonio
'''
import logging
import logging.config
import os
from pyworkflow.utils.path import *
import sys

SCIPION_PATH = 'Scipion'
LOG_PATH = 'logs'

""" Get general log file path """
logPath = join (getHomePath(), SCIPION_PATH, LOG_PATH, 'scipionLog.log')

""" Create the folders path if it does not exist """
makeFilePath(logPath)

""" Config the log """
config = {  'version': 1,              
            'disable_existing_loggers': False,
            'formatters': {
                'standard': {
                    'format': '%(asctime)s %(levelname)s %(name)s (%(lineno)d):  %(message)s'
                },
                'fileFormat': {
                    'format': '%(asctime)s %(levelname)s:  %(message)s'
                },
            },
            'handlers': {
                'fileHandler': {
                    'level': 'INFO',    
                    'class': 'logging.handlers.RotatingFileHandler',
                    'formatter': 'standard',
                    'filename': logPath,
                    'maxBytes': 100000,
                },
                'consoleHandler': {
                    'level': 'INFO',    
                    'class': 'logging.StreamHandler',
                    'formatter': 'standard',
                },
            },
            'loggers': {
                '': {                  
                    'handlers': ['consoleHandler','fileHandler'],        
                    'level': 'INFO',  
                    'propagate': False,
                    'qualname': 'pyworkflow',
                },
            }
        }

logging.config.dictConfig(config)


def getGeneralLogger(classPath):
    """ Method that returns the general log """
    return logging.getLogger(classPath)
    

class ScipionLogger():
    def __init__(self, filePath):
        makeFilePath(filePath)
        self._filePath = filePath

        if filePath not in config['loggers']:
            config['handlers'][filePath] = {'level': 'INFO',    
                                            'class': 'logging.handlers.RotatingFileHandler',
                                            'formatter': 'fileFormat',
                                            'filename': filePath,
                                            'maxBytes': 100000,}
            config['loggers'][filePath] = {'handlers': ['consoleHandler', filePath],        
                                           'level': 'INFO',  
                                           'propagate': False,}
            logging.config.dictConfig(config)
            
        self._log = logging.getLogger(filePath) 
        
    def getLog(self):
        return self._log    
        
    def info(self, message, redirectStandard = False, *args, **kwargs):
        if redirectStandard:
            print message
        self._log.info(message, *args, **kwargs)
    
    def warning(self, message, redirectStandard = False, *args, **kwargs):
        if redirectStandard:
            print message
        self._log.warning(message, *args, **kwargs)
        
    def error(self, message, redirectStandard = False, *args, **kwargs):
        if redirectStandard:
            print >> sys.stderr, message
        self._log.error(message, *args, **kwargs)    
        
    def close(self):
        if self._filePath in config['loggers']:
            del config['handlers'][self._filePath]
            del config['loggers'][self._filePath]   

#def closeFileLogger(filePath):
#    """ This method should be called to un-register a previous acquired
#    file logger with the method getFileLogger, the same filePath should
#    be used.
#    """
#    if filePath in config['loggers']:
#        del config['handlers'][filePath]
#        del config['loggers'][filePath]