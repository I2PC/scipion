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
        emxFile = join(emxDir, 'micrographs.emx')
        self.insertStep("createDir", 
                        verifyfiles=[emxDir], 
                        path=emxDir)
        
        emxMicrographs = join(emxDir, 'micrographs.emx')
        self.insertStep("exportMicrographs", 
                        verifyfiles=[emxMicrographs],
                        emxDir=emxDir, 
                        inputMd=self.MicrographsMd,
                        )
        
    def validate(self):
        errors = []
        
        return errors

    def summary(self):
        emxDir = self.getFilename('emxDir')
        emxFile = join(emxDir, 'micrographs.emx')
        message=[]
        message.append ('Images to export: [%s]' % self.MicrographsMd)
        if os.path.exists(emxFile):
            message.append ('EMX file: [%s]' % emxFile)
            message.append ('Directory with Binary files: [%s]' % emxDir)
        return message
    
    def visualize(self):
        emxDir = self.getFilename('emxDir')
        emxFile = join(emxDir, 'micrographs.emx')
        if os.path.exists(emxFile):
            from protlib_gui_ext import showTextfileViewer
            showTextfileViewer(emxFile,[emxFile])



def exportMicrographs(log, emxDir, inputMd):
    #TODO we get sampling form CTF, that is not correct
    emxData = EmxData()
    #get sampling rate
    xmippMicrographsToEmx(inputMd, emxData, emxDir)

