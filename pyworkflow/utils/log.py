#!/usr/bin/env python
# **************************************************************************
# *
# * Authors:     Antonio Poza (Apr 30, 2013)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import sys
import os
import logging
import logging.config
from pyworkflow.utils.path import makeFilePath



# Get general log file path
LOG_FILE = os.path.join(os.environ['SCIPION_LOGS'], 'scipion.log')


# Log configuration
config = {  'version': 1,              
            'disable_existing_loggers': False,
            'formatters': {
                'standard': {
                    'format': '%(asctime)s %(levelname)s:  %(message)s'
                    # TODO: use formattime to show the time less verbose
                },
                'fileFormat': {
                    'format': '%(asctime)s %(levelname)s:  %(message)s'
                },
            },
            'handlers': {
                'fileHandler': {
                    'level': 'NOTSET',
                    'class': 'logging.handlers.RotatingFileHandler',
                    'formatter': 'standard',
                    'filename': LOG_FILE,
                    'maxBytes': 100000,
                },
                'consoleHandler': {
                    'level': 'NOTSET',    
                    'class': 'logging.StreamHandler',
                    'formatter': 'standard',
                },
            },
            'loggers': {
                '': {                  
                    'handlers': ['consoleHandler', 'fileHandler'],
                    'level': 'INFO',  
                    'propagate': False,
                    'qualname': 'pyworkflow',
                },
            }
        }

logging.config.dictConfig(config)


class ScipionLogger():
    def __init__(self, filePath=''):
        """ If filePath is empty string, the general logger is used. """
        self._filePath = filePath
        makeFilePath(self._filePath)

        if self._filePath not in config['loggers']:
            config['handlers'][self._filePath] = {
                'level': 'NOTSET',
                'class': 'logging.handlers.RotatingFileHandler',
                'formatter': 'fileFormat',
                'filename': self._filePath,
                'maxBytes': 100000}

            config['loggers'][self._filePath] = {
                'handlers': [self._filePath],
                'level': 'NOTSET',
                'propagate': False}
            # Note: if we want to see in the console what we also have in
            # run.log, add 'consoleHandler' to the list of 'handlers'.

            logging.config.dictConfig(config)
            
        self._log = logging.getLogger(self._filePath) 
        
    def getLog(self):
        return self._log  
    
    def getLogString(self):
        return open(self._filePath, 'r').readlines()  
        
    def info(self, message, redirectStandard=False, *args, **kwargs):
        if redirectStandard:
            print message
            sys.stdout.flush()
        self._log.info(message, *args, **kwargs)

    def warning(self, message, redirectStandard=False, *args, **kwargs):
        if redirectStandard:
            print message
            sys.stdout.flush()
        self._log.warning(message, *args, **kwargs)
        
    def error(self, message, redirectStandard=False, *args, **kwargs):
        if redirectStandard:
            print >> sys.stderr, message
            sys.stderr.flush()
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
