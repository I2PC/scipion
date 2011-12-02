#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------

# Author: Carlos Oscar, December 2011
#
from protlib_base import *
from xmipp import MetaData
import glob
import os
from protlib_filesystem import deleteFile, deleteDir, createLink, copyFile

class ProtMLTomo(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.mltomo.name, scriptname, project)
        self.Import = 'from protocol_mltomo import *'

    def createFilenameTemplates(self):
        return {
                'inputXmd': self.workingDirPath('inputVolumeList.xmd'),
                'missingXmd': self.workingDirPath('missingRegions.xmd'),
                'oroot': self.workingDirPath('mltomo')
                }

    def defineSteps(self):
        self.Db.insertStep("createLink",source=self.VolumeList,dest=self.getFilename("inputXmd"))
        if self.MissingMd!="":
            self.Db.insertStep("createLink",source=self.MissingMd,dest=self.getFilename("missingXmd"))
        
        params= ' -i '       + self.getFilename("inputXmd") + \
                ' --oroot '  + self.getFilename("oroot") + \
                ' --iter '   + str(self.NumberOfIterations) + \
                ' --sym '    + str(self.Symmetry) +\
                ' --maxres ' + str(self.MaximumResolution)
        if self.MissingMd!="":
            params+= ' --missing ' + str(self.MissingDocFile)
        if self.Dimension > 0:
            params+= ' --dim ' + str(self.Dimension)

        # Alignment parameters
        params+= ' --ang ' + str(self.AngSampling)
        if self.AngSampling<360:
            params+= ' -ang_search ' + str(self.AngLimit)
        if self.DoPerturb:
            params+= ' --perturb '
        if self.ShiftLimit>0:
            params+=" --limit_trans " + str(self.ShiftLimit)

        # Classification parameters
        if self.DoGenerateReferences:
            params+= ' --nref ' + str(self.NumberOfReferences)
        else:
            params+= ' --ref ' + str(self.RefMd)
        if self.InitialRegularization>0:
            params+= ' --reg0 ' + str(self.InitialRegularization)
            params+= ' --reg_steps ' + str(self.NumberRegularizationSteps)

        # Extra parameters
        if (len(self.ExtraParamsMLtomo) > 0):
            params+=' ' + str (self.ExtraParamsMLtomo)

        # Thread parallelization
        if self.NumberOfThreads > 1:
            params+=' --thr ' + str(self.NumberOfThreads)
        verifyFiles=[]
        self.insertRunJobStep("xmipp_ml_tomo", params, verifyFiles)
                
    def validate(self):
        errors = []
        if not os.path.exists(self.VolumeList):
            errors.append(self.VolumeList+" does not exist")
        if self.MissingMd!="" and not os.path.exists(self.MissingMd):
            errors.append(self.MissingMd+" does not exist")
        if self.MaximumResolution>0.5:
            errors.append("Maximum resolution must be smaller than 0.5")
        return errors

    def summary(self):
        message=[]
        if self.DoGenerateReferences:
            N=self.NumberOfReferences
        else:
            mD=MetaData(self.RefMd)
            N=mD.size()
        message.append("MLTomo on "+self.VolumeList+" with "+str(N)+" references")
        return message

    def visualize(self):
        pass
