#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
#
# General script for Xmipp-based pre-processing of volumes 

# Author: Carlos Oscar, August 2013
#
from protlib_base import *
import os
from xmipp import MetaData
from os.path import exists, split, splitext
from protlib_utils import runJob, runShowJ
from protlib_filesystem import linkAcquisitionInfo
import glob
from protlib_gui_ext import showError

class ProtMetadataSplit(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.metadata_split.name, scriptname, project)
        self.Import = 'from protocol_metadata_split import *'

    def defineSteps(self):
        self.insertStep('linkAcquisitionInfo',InputFile=self.InMetadata, dirDest=self.WorkingDir)
        self.insertStep('splitMetadata',InMetadata=self.InMetadata,WorkingDir=self.WorkingDir,Nparts=self.Nparts,
                        SortBy=self.SortBy,RemoveDisabled=self.RemoveDisabled,RandomSplit=self.RandomSplit)
        
    def validate(self):
        errors = []
        return errors
        
    def summary(self):
        messages = []      
        messages.append("Input metadata: [%s]" % self.InMetadata)
        messages.append("Number of parts: %d"%self.Nparts)
        if self.SortBy!="Do not sort":
            messages.append("Sorted by %s"%self.SortBy)
        if self.RemoveDisabled:
            messages.append("Disabled images have been removed")
        if self.RandomSplit:
            messages.append("Random split")
        else:
            messages.append("Deterministic split")
        return messages

    def visualize(self):
        from protlib_utils import runShowJ
        import glob
        files=glob.glob(self.workingDirPath('images*.xmd'))
        if files:
            runShowJ(" ".join(files))

def splitMetadata(log,InMetadata,WorkingDir,Nparts,SortBy,RemoveDisabled,RandomSplit):
    fnRoot=os.path.join(WorkingDir,'images')
    args="-i %s --oroot %s -n %d"%(InMetadata,fnRoot,Nparts)
    if SortBy=="Do not sort":
        args+=" --dont_sort"
    elif SortBy=="image name":
        args+=" --l image"
    elif SortBy=="micrograph name":
        args+=" --l micrograph"
    if not RandomSplit:
        args+=" --dont_randomize"
    runJob(log,"xmipp_metadata_split",args)
