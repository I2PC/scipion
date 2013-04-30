'''
Created on Apr 30, 2013

@author: antonio
'''
import logging
import logging.config
import os
from pyworkflow import HOME

logging.config.fileConfig(os.path.join(HOME, 'settings/log.config'))

def getGeneralLogger(classPath):
    return logging.getLogger(classPath)
    
def getFileLogger(fileName):
    mylog = logging.getLogger(fileName)
    for handler in mylog.handlers:
        if handler.get_name() == fileName:
            return mylog
    hdlr = logging.FileHandler(fileName)
    hdlr.set_name(fileName)
    formatter = logging.Formatter('%(asctime)s %(levelname)s %(name)s (%(lineno)d):  %(message)s')
    hdlr.setFormatter(formatter)
    mylog.addHandler(hdlr)
    mylog.setLevel(logging.INFO)
    return mylog