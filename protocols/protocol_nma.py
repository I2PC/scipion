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
from protlib_filesystem import createLink, getExt

class ProtNMA(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.nma.name, scriptname, project)
        self.Import = 'from protocol_nma import *'    

    def createFilenameTemplates(self):
        return {
            'normal_modes': "%(WorkingDir)s/normal_modes.xmd",
            'pseudoatoms':  '%(WorkingDir)s/pseudoatoms'
            }

    def defineSteps(self):
        self.insertStep('createDir',path=self.ExtraDir)
        if self.StructureType=="EM":
            # Mask or link input structure
            fnMask=""
            if self.MaskMode=="Threshold":
                fnMask=self.extraPath('mask.vol')
                self.insertRunJobStep("xmipp_transform_threshold", params="-i %s -o %s --select below %f --substitute binarize"%\
                                      (self.InputStructure,fnMask,self.Threshold),verifyFiles=[fnMask])
            else:
                fnMask=self.MaskFile
            
            # Convert to pseudoatoms
            fnOut=self.getFileName("pseudoatoms")
            params="-i %s -o %s --sigma %f --targetError %f --sampling_rate %f"%\
                (self.InputStructure,fnOut,self.PseudoAtomRadius*self.Sampling,self.PseudoAtomTarget,self.Sampling)
            if fnMask!="":
                params+=" --mask binary_file %s"%fnMask
            self.insertRunJobStep("xmipp_volume_to_pseudoatoms", params=params,verifyFiles=[fnOut])
    
    def summary(self):
        message=[]
        message.append('NMA of [%s]'%self.InputStructure)
        return message
    
    def validate(self):
        errors = []
        if self.MaskMode=="Binary mask" and not os.path.exists(self.MaskFile):
            errors.append(self.MaskFile+" does not exist")
        return errors
    
    def visualize(self):
        pass
