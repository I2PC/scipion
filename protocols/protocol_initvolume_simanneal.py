#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Author: Carlos Oscar Oct 2013 
#

from protlib_base import *
from protlib_filesystem import createLink

class ProtInitVolSimAnneal(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.initvolume_simanneal.name, scriptname, project)
        self.Import = 'from protocol_initvolume_simanneal import *'
        
    def defineSteps(self):
        self.insertStep('createDir',path=self.ExtraDir)
        args="-i %s --oroot %s --sym %s --randomIter %d --greedyIter %d --rejection %f --keepIntermediateVolumes --T0 %f --angularSampling %f"\
            %(self.Classes,self.extraPath("proposedVolume"),self.SymmetryGroup,self.NIterRandom,self.NIterGreedy,self.Rejection,self.T0,
              self.AngularSampling)
        if self.DontApplyPositiveConstraint:
            args+=" --dontApplyPositive"
        if self.InitialVolume!="":
            args+=" --initial "+self.InitialVolume
        self.insertRunJobStep("xmipp_volume_initial_simulated_annealing",args)
        self.insertStep("createLink",source=self.extraPath("proposedVolume.vol"),dest=self.workingDirPath("proposedVolume.vol"))
        
    def summary(self):
        message=[]
        message.append("Input images: [%s]"%self.Classes)
        message.append("Simulated annealing iterations: %d"%self.NIterRandom)
        if self.NIterGreedy>0:
            message.append("Greedy iterations: %d"%self.NIterGreedy)
        if self.InitialVolume!="":
            message.append("Initial volume: %s"%self.InitialVolume)
        message.append("Symmetry: %s"%self.SymmetryGroup)
        return message
    
    def visualize(self):
        fnVolume = self.workingDirPath('proposedVolume.vol')
        if os.path.exist(fnVolume):
            runShowJ(fnVolume)
