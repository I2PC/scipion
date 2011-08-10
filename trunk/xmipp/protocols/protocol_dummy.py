#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# Protocol for Xmipp-based 2D alignment and classification,
# using maximum-likelihood principles
#
#   Author:  Sjors Scheres, January 2008
#  Updated:  J. M. de la Rosa Trevin July 2011
#

from protlib_base import XmippProtocol, protocolMain
from config_protocols import protDict
import os
from protlib_filesystem import copyFile
from xmipp import MetaData

class ProtDummy(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.dummy.key, scriptname, project)
        self.Import = 'from protocol_dummy import *'

    def validate(self):
        return []
        #return ["Protocol not implemented yet..."]
    
    def summary(self):
        return ["This is a test summary",
                "Need a real summary here..."
                ]
        
    def defineSteps(self):
        filename = self.InputMd
        backup = os.path.join(self.WorkingDir, os.path.basename(filename) + ".backup")
        self.Db.insertStep('backupMetaData', [backup], filename=filename, backup=backup)
        self.Db.insertStep('splitMetaData', filename=filename, parts=self.NumberOfParts)

def backupMetaData(log, filename, backup):    
    copyFile(log, filename, backup)
    
def splitMetaData(log, filename, parts):
    md = MetaData(filename)
    print md

