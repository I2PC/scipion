#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
#
# General script for Xmipp-based pre-processing of single-particles: 
#  - phase flipping
#  - extraction of particles
#  - normalization
#  - sort_junk

# Author: Carlos Oscar, August 2011
#
from protlib_base import *
import glob
import os
from protlib_utils import runJob, runShowJ
from protlib_filesystem import createLink, linkAcquisitionInfo
from xmipp import MetaData, MDL_ITEMID

class ProtSubsetParticles(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.subset_particles.name, scriptname, project)
        self.Import = 'from protocol_subset_particles import *'
        self.outputFile=self.workingDirPath('images.xmd')

    def defineSteps(self):
        self.Db.insertStep("linkAcquisitionInfo",InputFile=self.InputFile,dirDest=self.WorkingDir)
        self.Db.insertStep('createSubset',inputFile=self.InputFile,subsetFile=self.subsetFile,outputFile=self.outputFile)
                
    def summary(self):
        message=[]
        message.append("Input  file: [%s]" % self.InputFile)
        message.append("Subset file: [%s]" % self.Subsetfile)
        return message
    
    def validate(self):
        errors=[]
        mdIn=MetaData(self.InputFile)
        if not mdIn.containsLabel(MDL_ITEMID):
            errors.append(self.InputFile+" does not contain an item_id column")
        mdSubset=MetaData(self.SubsetFile)
        if not mdSubset.containsLabel(MDL_ITEMID):
            errors.append(self.SubsetFile+" does not contain an item_id column")
        return errors

    def visualize(self):
        if not os.path.exists(self.outputFile):
            from protlib_gui_ext import showWarning
            showWarning("Error", "There is no result yet",self.master)
        else:   
            runShowJ(self.outputFile)                                     

def createSubset(log, inputFile, subsetFile, outputFile):
    runJob(log,"xmipp_image_sort_by_statistics","-i "+inputFile+"--set join "+subsetFile+" itemId -o "+outputFile)
