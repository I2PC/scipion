#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
#
# General script for Xmipp-based merging of different particle sets 
#
# Author: Roberto Marabini        (roberto@cnb.csic.es)     July 2013
#         J. M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
#
#
from protlib_base import *
import xmipp
import glob
import os
from os.path import relpath, dirname, exists
from protlib_utils import runJob
from protlib_filesystem import findAcquisitionInfo


class ProtMergeParticles(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.merge_particles.name, scriptname, project)
        self.Import += 'from protocol_merge_particles import *'
        self.acquisionInfo1 = self.findAcquisitionInfo(self.InputImages1)
        self.acquisionInfo2 = self.findAcquisitionInfo(self.InputImages2)

    def defineSteps(self):
        self.insertCopyFile(self.acquisionInfo1, self.getFilename('acquisition'))
        fnOut = self.getFilename('images')
        self.insertStep('createOutputMd', verifyfiles=[fnOut], 
                        input1=self.InputImages1, input2=self.InputImages2, fnOut=fnOut)
          
    def getOutputImages(self):
        return self.getFilename('images')
    
    def validate(self):
        errors = []
        if self.acquisionInfo1 is None or self.acquisionInfo2 is None:
            errors.append("One of the <acquisition_info.xmd> file not found for input images.")
        else:
            md1 = xmipp.MetaData(self.acquisionInfo1)
            md2 = xmipp.MetaData(self.acquisionInfo2)
            samplingRate1 = md1.getValue(xmipp.MDL_SAMPLINGRATE, md1.firstObject())
            samplingRate2 = md2.getValue(xmipp.MDL_SAMPLINGRATE, md2.firstObject())
            if samplingRate1 != samplingRate2:
                errors.append('Image sets have different sampling rates:\n   <%f> and <%f>, respectively.' % (samplingRate1, samplingRate2))
            xdim1, ydim1, zdim1, _, _ = xmipp.MetaDataInfo(self.InputImages1)
            xdim2, ydim2, zdim2, _, _ = xmipp.MetaDataInfo(self.InputImages2)
            dim1 = (xdim1, ydim1, zdim1)
            dim2 = (xdim2, ydim2, zdim2)
            if  dim1 != dim2:
                errors.append("Image sets have different dimensions:\n   <%s> and <%s>, respectively" % (dim1, dim2)) 
        return errors

    def summary(self):
        messages = []
        messages.append("Import images set 1 from: [%s]" % self.InputImages1)
        messages.append("Import images set 2 from: [%s]" % self.InputImages2)
        
        if exists(self.getOutputImages()):
            messages.append("Output images: [%s]" % self.getOutputImages())
            
        return messages
    
    def visualize(self):
        fn = self.getFilename('images')        
        if exists(fn):
            from protlib_utils import runShowJ
            runShowJ(fn)
	    
def createAcquisitionMd(log, samplingRate, fnOut):
    md = MetaData()
    md.setValue(MDL_SAMPLINGRATE, float(samplingRate), md.addObject())
    md.write(fnOut)

def createOutputMd(log, input1, input2, fnOut):
    """ Create the output metadata from set1 and set2. """
    md = xmipp.MetaData(input1)
    md2 = xmipp.MetaData(input2)
    md.unionAll(md2)
    md.write("images@" + fnOut)
