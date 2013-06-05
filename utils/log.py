'''
Created on Apr 30, 2013

@author: antonio
'''
import logging
import logging.config
import os
from pyworkflow.utils.path import *

SCIPION_PATH = 'Scipion'
LOG_PATH = 'logs'

""" Get general log file path """
logPath = join (getHomePath(), SCIPION_PATH, LOG_PATH, 'scipionLog.log')
""" Create the folders path if it does not exist """
createFolderForFile(logPath)
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
    

def getFileLogger(filePath):
    """ Method that creates and returns a log for the given file """
    # Create the folders path if it does not exist 
    createFolderForFile(filePath)
    
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
    return logging.getLogger(filePath)