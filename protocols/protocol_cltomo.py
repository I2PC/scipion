#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------

# Author: Carlos Oscar, September 2013
#
from protlib_base import *
from xmipp import MetaData
import glob
import os
from protlib_filesystem import createLink
from protlib_utils import getListFromVector, runJob

class ProtCLTomo(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.cltomo.name, scriptname, project)
        self.Import = 'from protocol_cltomo import *'

    def defineSteps(self):
        self.insertStep('createDir',path=self.ExtraDir)
        params= ' -i '            + self.VolumeList + \
                ' --oroot '       + self.extraPath("results") + \
                ' --iter '        + str(self.NumberOfIterations) + \
                ' --nref '        + str(self.NumberOfReferences) + \
                ' --sym '         + self.Symmetry + \
                ' --maxFreq '     + str(self.MaximumResolution) + \
                ' --sparsity '    + str(self.Sparsity/100.0) + \
                ' --DWTsparsity ' + str(self.DWTSparsity/100.0) + \
                ' --maxShiftX '   + str(self.MaxShiftX) + \
                ' --maxShiftY '   + str(self.MaxShiftY) + \
                ' --maxShiftZ '   + str(self.MaxShiftZ) + \
                ' --maxRot '      + str(self.MaxRot) + \
                ' --maxTilt '     + str(self.MaxTilt) + \
                ' --maxPsi '      + str(self.MaxPsi)
        if self.DoGenerateInitial:
            params+=' --nref0 '+str(self.NumberOfReferences0)
            if self.RandomizeOrientation:
                params+=' --randomizeStartingOrientation'
        else:
            params+=' --ref0 '+self.RefMd
        if self.Mask!="":
            params+=' --mask '+self.Mask
        if self.GenerateAligned:
            params+=" --generateAlignedVolumes"
        if self.DontAlign:
            params+=" --dontAlign"
        self.insertStep("runCLTomo", WorkingDir=self.WorkingDir, params=params, nproc=self.NumberOfMpi)
                
    def validate(self):
        errors = []
        if not os.path.exists(self.VolumeList):
            errors.append(self.VolumeList+" does not exist")
        if self.MaximumResolution>0.5:
            errors.append("Maximum resolution must be smaller than 0.5")
        if self.Sparsity>100 or self.Sparsity<0:
            errors.append("Sparsity must be a value between 0 and 100");
        if self.DWTSparsity>100 or self.DWTSparsity<0:
            errors.append("Wavelet sparsity must be a value between 0 and 100");
        return errors

    def summary(self):
        message=[]
        message.append("CLTomo on [%s] with %d references" % (self.VolumeList, self.NumberOfReferences))
        if not self.DoGenerateInitial:
            message.append("Initial classes: [%s]"%(self.RefMd))
        return message

    def visualize(self):
        from protlib_utils import runShowJ
        fnClasses=self.workingDirPath('classes.xmd')
        if os.path.exists(fnClasses):
            runShowJ(fnClasses)

def runCLTomo(log, WorkingDir, params, nproc):
    runJob(log,'xmipp_mpi_classify_CLTomo','%d %s'%(nproc,params))
    levelFiles=glob.glob(WorkingDir+"/extra/results_classes_level*.xmd")
    if levelFiles:
        levelFiles.sort()
        lastLevelFile=levelFiles[-1]
        createLink(log,lastLevelFile,os.path.join(WorkingDir,'results.xmd'))
