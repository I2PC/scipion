#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
#  DUMMY PROTOCOL FOR TESTING
#

from protlib_base import XmippProtocol, protocolMain
from config_protocols import protDict
import os
from protlib_filesystem import copyFile
from xmipp import MetaData
from protlib_sql import SqliteDb

class ProtDummy(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.dummy.name, scriptname, project)
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
        backup = self.workingDirPath( os.path.basename(filename) + ".backup")
        self.Db.insertStep('backupMetaData', [backup], filename=filename, backup=backup)
        id = self.Db.insertStep('splitMetaData', filename=filename,workingDir=self.WorkingDir,  
                           parts=self.NumberOfParts)
        self.Db.insertStep('runStepGapsMpi',passDb=True, script=self.scriptName, 
                                 NumberOfMpi=2)
        for i in range(self.NumberOfParts):
            self.Db.insertStep('sortMetaDataPart', execution_mode=SqliteDb.EXEC_GAP,
                               parent_step_id=id,
                               filename="%s/part%03d_%s" % (self.WorkingDir, i, filename))
        backup += ".2"
        self.Db.insertStep('backupMetaData', [backup], filename=filename, backup=backup)
        
def backupMetaData(log, filename, backup):    
    copyFile(log, filename, backup)
    
def splitMetaData(log, workingDir, filename, parts):
    md1 = MetaData(filename)
    md2 = MetaData()
    md2.randomize(md1)
    chunk = md2.size() / parts
    for i in range(parts):
        md1.selectPart(md2, i*chunk, chunk)
        md1.write("%(workingDir)s/part%(i)03d_%(filename)s" % locals())

def sortMetaDataPart(log, filename):
    md = MetaData(filename)
    md.sort()
    md.write(filename)