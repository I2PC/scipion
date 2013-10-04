#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Author: Carlos Oscar October 2013 
#

from protlib_initvolume import *

class ProtInitVolH(ProtInitVolumeBase):
    def __init__(self, scriptname, project):
        ProtInitVolumeBase.__init__(self, protDict.initvolume_heuristic.name, scriptname, project)
        self.Import += 'from protocol_initvolume_heuristic import *'
        self.fnRoot=self.workingDirPath('proposedVolume')
        
    def defineSteps(self):
        ProtInitVolumeBase.defineSteps(self)
        args="-i %s --oroot %s --sym %s --randomIter %d --greedyIter %d --rejection %f --thr %d"%\
            (self.Classes,self.fnRoot,self.SymmetryGroup,self.NIterRandom, self.NIterGreedy, self.Rejection, self.NumberOfThreads)
        if self.InitialVolume!="":
            args+=" --initial "+self.InitialVolume
        if not self.Positive:
            args+=" --dontApplyPositive"
        if self.KeepIntermediate:
            args+=" --keepIntermediateVolumes"
        
        self.insertRunJobStep("xmipp_volume_initial_H", args, [self.fnRoot+".vol"])
        if self.Xdim!=self.Xdim2:
            self.insertRunJobStep("xmipp_image_resize","-i %s --dim %d"%(self.fnRoot+".vol",self.Xdim))

    def visualize(self):
        os.system("xmipp_chimera_client -i %s.vol --mode projector 256 --angulardist %s.xmd"%(self.fnRoot,self.fnRoot))
