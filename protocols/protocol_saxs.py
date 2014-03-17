#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Protocol for Xmipp-based 2D alignment and classification,
# using hierarchical clustering principles
# Author: Carlos Oscar Sanchez Sorzano, August 2011
#

import glob,os,re,sys,shutil,time
from protlib_base import *
from config_protocols import protDict
from protlib_utils import runJob, which
from protlib_filesystem import createLink, changeDir
from protlib_gui_figure import XmippPlotter

class ProtSaxs(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.saxs.name, scriptname, project)
        self.Import = 'from protocol_saxs import *'

    def createFilenameTemplates(self):
        return {
            'pseudoatoms':  '%(WorkingDir)s/pseudoatoms',
            'extra_pseudoatoms':  '%(ExtraDir)s/pseudoatoms',
            }

    def defineSteps(self):
        self.insertStep('createDir',path=self.ExtraDir)

        # Link the input
        fnLocalInputStructure=self.workingDirPath(os.path.basename(self.InputStructure))
        self.insertStep("createLink",source=self.InputStructure,dest=fnLocalInputStructure)
        
        # Mask or link input structure
        fnMask=""
        if self.MaskMode=="Threshold":
            fnMask=self.extraPath('mask.vol')
            self.insertRunJobStep("xmipp_transform_threshold", params="-i %s -o %s --select below %f --substitute binarize"%\
                                  (self.InputStructure,fnMask,self.Threshold),verifyFiles=[fnMask])
        elif self.MaskMode=="Binary mask":
            fnMask=self.MaskFile
        
        # Convert to pseudoatoms
        fnOut=self.getFilename("pseudoatoms")
        params="-i %s -o %s --sigma %f --targetError %f --sampling_rate %f -v 2 --intensityColumn Bfactor"%\
            (self.InputStructure,fnOut,self.PseudoAtomRadius*self.Sampling,self.PseudoAtomTarget,self.Sampling)
        if fnMask!="":
            params+=" --mask binary_file %s"%fnMask
        self.insertRunJobStep("xmipp_volume_to_pseudoatoms", params=params,verifyFiles=[fnOut+".pdb"])
        self.insertStep('moveFile',source=fnOut+"_approximation.vol",dest=self.getFilename("extra_pseudoatoms")+"_approximation.vol")
        self.insertStep('moveFile',source=fnOut+"_distance.hist",dest=self.getFilename("extra_pseudoatoms")+"_distance.hist")
        self.insertRunJobStep("rm", params=fnOut+"_*")
        
        # Compute SAXS curve
        self.insertStep("computeSAXS",WorkingDir=self.WorkingDir,NumberOfHarmonics=self.NumberOfHarmonics,
                        MaxFreq=self.MaxFreq, NumberOfSamples=self.NumberOfSamples, OtherCrysol=self.OtherCrysol,
                        ExperimentalSAXS=self.ExperimentalSAXS)
    
    def summary(self):
        message=[]
        message.append("Input volume: [%s] "%self.InputStructure)
        return message
    
    def validate(self):
        errors = []
        return errors
    
    def visualize(self):
        if self.ExperimentalSAXS=="":
            fnInt=self.workingDirPath("pseudoatoms00.int")
        else:
            fnInt=self.workingDirPath("pseudoatoms00.fit")
        
        import numpy
        x=numpy.loadtxt(fnInt,skiprows=1)
        xplotter = XmippPlotter(*[1,1],windowTitle="SAXS Curves")
        a = xplotter.createSubPlot('SAXS curves', 'Armstrongs^-1', 'log(SAXS)', yformat=False)
        a.plot(x[:,0], numpy.log(x[:,1]))
        a.plot(x[:,0], numpy.log(x[:,2]))
        if self.ExperimentalSAXS=="":
            xplotter.showLegend(['SAXS in solution','SAXS in vacuo'])
        else:
            xplotter.showLegend(['Experimental SAXS','SAXS from volume'])
        xplotter.draw()
        xplotter.show()

def computeSAXS(log, WorkingDir,NumberOfHarmonics,MaxFreq, NumberOfSamples, OtherCrysol, ExperimentalSAXS):
    if ExperimentalSAXS!="":
        ExperimentalSAXS=os.path.abspath(ExperimentalSAXS)
    currentDir=os.getcwd()
    changeDir(log,WorkingDir)
    if ExperimentalSAXS!="":
        ExperimentalSAXS=os.path.relpath(ExperimentalSAXS)
    runJob(log,"crysol","pseudoatoms.pdb %s /lm %d /sm %f /ns %d %s"%(ExperimentalSAXS,NumberOfHarmonics,MaxFreq,NumberOfSamples,OtherCrysol))
    runJob(log,"mv","*log *txt extra")
    if ExperimentalSAXS=="":
        runJob(log,"mv","*alm extra")
    changeDir(log,currentDir)
