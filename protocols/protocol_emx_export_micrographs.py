#!/usr/bin/env xmipp_python

#------------------------------------------------------------------------------------------------
#   General script for exporting micrographs to EMX format
#
#   Authors: J.M. de la Rosa Trevin, Sept 2013
#            Roberto Marabini
# 

from os.path import relpath, join, abspath, dirname, exists
from glob import glob
import math

import xmipp
from emx import *
from emx.emxmapper import *

from protlib_base import *
from protlib_filesystem import replaceBasenameExt, renameFile, copyFile
from protlib_utils import runJob
from protlib_xmipp import redStr, RowMetaData
from protlib_emx import *


class ProtEmxExportMicrographs(XmippProtocol):
    
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.emx_export_micrographs.name, scriptname, project)
        self.Import = "from protocol_emx_export_micrographs import *"

            
    def createFilenameTemplates(self):
        return {
                 'emxDir': join('%(WorkingDir)s','emxData')
            } 
        
    def defineSteps(self):
        emxDir = self.getFilename('emxDir')
        self.insertStep("createDir", verifyfiles=[emxDir], 
                        path=emxDir)
        
        emxMicrographs = join(emxDir, 'micrographs.emx')
        self.insertStep("exportMicrographs", verifyfiles=[emxMicrographs],
                        emxDir=emxDir, inputMd=self.MicrographsMd,
                        )
        
    def validate(self):
        errors = []
        
        return errors

    def summary(self):
        summary = ['Images to export: [%s]' % self.MicrographsMd,
                   ]
        return summary
    
    def visualize(self):
        pass



def exportMicrographs(log, emxDir, inputMd):
    emxData = EmxData()
    xmippMicrographsToEmx(inputMd, emxData, emxDir)

