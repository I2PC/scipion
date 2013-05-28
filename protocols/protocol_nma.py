#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Protocol for Normal Mode analysis of atomic and EM structures
# Author: Carlos Oscar Sanchez Sorzano, May 2013
#         Slavica Jonic
#

import glob,os,re,sys,shutil,time
from protlib_base import *
from config_protocols import protDict
from protlib_utils import runJob

class ProtNMA(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.nma.name, scriptname, project)
        self.Import = 'from protocol_nma import *'    

    def createFilenameTemplates(self):
        return {
            'normal_modes': "normal_modes.xmd"
            }

    def defineSteps(self):
        self.insertStep('createDir',path=self.ExtraDir)
    
    def summary(self):
        message=[]
        message.append('NMA of [%s]'%self.InputStructure)
        return message
    
    def validate(self):
        errors = []
        return errors
    
    def visualize(self):
        pass
