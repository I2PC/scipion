#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------

# Author: Carlos Oscar, December 2011
#
from protlib_base import *
from xmipp import MetaData
import glob
import os
from protlib_filesystem import deleteFile, deleteDir, createLink, copyFile
from protlib_utils import getListFromVector

class ProtMLTomo(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.mltomo.name, scriptname, project)
        self.Import = 'from protocol_mltomo import *'

    def createFilenameTemplates(self):
        return {
                'missingXmd': self.workingDirPath('missingRegions.xmd'),
                }

    def defineSteps(self):
        if self.MissingMd!="":
            self.Db.insertStep("createLink",source=self.MissingMd,dest=self.getFilename("missingXmd"))

        AngSamplingList=getListFromVector(self.AngSampling)
        AngLimitList=getListFromVector(self.AngLimit)
        ShiftLimitList=getListFromVector(self.ShiftLimit)
        NumberOfIterationsList=getListFromVector(self.NumberOfIterations)

        N=len(AngSamplingList)
        for n in range(0,N):
            workingDir=os.path.join(self.WorkingDir,"refinement_%02d"%n)
            self.Db.insertStep("createDir",verifyfiles=[workingDir],path=workingDir)
            fnIn=os.path.join(workingDir,'inputVolumeList.xmd')
            fnRoot=os.path.join(workingDir,'mltomo')
            if n==0:
                self.Db.insertStep("createLink",verifyfiles=[fnIn],source=self.VolumeList,dest=fnIn)
            else:
                workingDir_1=os.path.join(self.WorkingDir,"refinement_%02d"%(n-1))
                self.Db.insertStep("createLink",verifyfiles=[fnIn],source=os.path.join(workingDir_1,'mltomo_img.xmd'),dest=fnIn)
                
            params= ' -i '       + fnIn + \
                    ' --oroot '  + fnRoot + \
                    ' --iter '   + NumberOfIterationsList[n] + \
                    ' --sym '    + self.Symmetry +\
                    ' --maxres ' + str(self.MaximumResolution)
            # getComponentFromVector, getListFromVector
            if self.MissingMd!="":
                params+= ' --missing ' + self.getFilename("missingXmd")
            if self.Dimension > 0:
                params+= ' --dim ' + str(self.Dimension)
    
            # Alignment parameters
            params+= ' --ang ' + AngSamplingList[n]
            if int(AngLimitList[n])<360:
                params+= ' --ang_search ' + AngLimitList[n]
            if self.DoPerturb:
                params+= ' --perturb '
            if int(ShiftLimitList[n])>0:
                params+=" --limit_trans " + ShiftLimitList[n]
    
            # Classification parameters
            if self.DoGenerateReferences and n==0:
                params+=' --nref ' + str(self.NumberOfReferences)
            else:
                if n==0:
                    params+=' --ref ' + self.RefMd
                else:
                    params+=' --ref ' + os.path.join(workingDir_1,"mltomo_ref.xmd")
            if self.InitialRegularization>0:
                params+= ' --reg0 ' + str(self.InitialRegularization)
                params+= ' --reg_steps ' + str(self.NumberRegularizationSteps)
    
            # Extra parameters
            if (len(self.ExtraParamsMLtomo) > 0):
                params+=' ' + str (self.ExtraParamsMLtomo)
    
            # Thread parallelization
            if self.NumberOfThreads > 1:
                params+=' --thr ' + str(self.NumberOfThreads)
            verifyFiles=[fnRoot+"_img.xmd"]
            self.insertRunJobStep("xmipp_ml_tomo", params, verifyFiles)
        lastDir=os.path.join(self.WorkingDir,"refinement_%02d"%(N-1))
        fnOut=self.workingDirPath("mltomo.fsc")
        self.Db.insertStep("createLink",verifyfiles=[fnOut],source=os.path.join(lastDir,"mltomo.fsc"),dest=fnOut)
        fnOut=self.workingDirPath("mltomo_img.xmd")
        self.Db.insertStep("createLink",verifyfiles=[fnOut],source=os.path.join(lastDir,"mltomo_img.xmd"),dest=fnOut)
        fnOut=self.workingDirPath("mltomo_ref.xmd")
        self.Db.insertStep("createLink",verifyfiles=[fnOut],source=os.path.join(lastDir,"mltomo_ref.xmd"),dest=fnOut)
                
    def validate(self):
        errors = []
        if not os.path.exists(self.VolumeList):
            errors.append(self.VolumeList+" does not exist")
        if self.MissingMd!="" and not os.path.exists(self.MissingMd):
            errors.append(self.MissingMd+" does not exist")
        if self.MaximumResolution>0.5:
            errors.append("Maximum resolution must be smaller than 0.5")

        AngSamplingList=getListFromVector(self.AngSampling)
        AngLimitList=getListFromVector(self.AngLimit)
        ShiftLimitList=getListFromVector(self.ShiftLimit)
        NumberOfIterationsList=getListFromVector(self.NumberOfIterations)
        if len(AngSamplingList)!=len(AngLimitList) or \
           len(AngSamplingList)!=len(ShiftLimitList) or \
           len(AngSamplingList)!=len(NumberOfIterationsList):
           errors.append("You must supply the same number of parameters for Angular Sampling, Angular Limit, Shift Limit and Number of Iterations")

        return errors

    def summary(self):
        message=[]
        if self.DoGenerateReferences:
            N=self.NumberOfReferences
        else:
            from protlib_xmipp import getMdSize
            N = getMdSize(self.RefMd)
        message.append("MLTomo on [%s] with %d references" % (self.VolumeList, N))

        return message

    def visualize(self):
        pass
