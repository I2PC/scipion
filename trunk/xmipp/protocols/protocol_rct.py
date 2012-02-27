#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Protocol for random conical tilt reconstruction
#
# Example use:
# ./xmipp_protocol_rct.py
#
# Author: Sjors Scheres, March 2007
#
from protlib_base import *
from protlib_utils import getListFromRangeString, runJob

class ProtRCT(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.rct.name, scriptname, project)
        self.Import = 'from protocol_rct import *'

    def defineSteps(self):
        extractionProt = self.getProtocolFromRunName(self.ExtractionRun)
        pickingDir = getWorkingDirFromRunName(extractionProt.PickingRun)

        classNumbers = getListFromRangeString(self.SelectedClasses)
        extractRootName = os.path.join(getWorkingDirFromRunName(self.ExtractionRun),"Default")
        self.insertStep('gatherPairs',verifyfiles=[os.path.join(self.WorkingDir,"classes.xmd")],
                               WorkingDir=self.WorkingDir,ClassNumbers=classNumbers,
                               ClassFile=self.ClassifMd,ExtractRootName=extractRootName,
                               PickingDir=pickingDir)
        for classNo in classNumbers:
            classNameIn="class_%06d@%s"%(classNo,self.workingDirPath("classes.xmd"))
            classNameOut=self.workingDirPath("rct_%06d.xmd"%classNo)
            self.insertParallelStep('reconstructClass',verifyfiles=[classNameOut],
                               WorkingDir=self.WorkingDir,
                               ClassNameIn=classNameIn, ClassNameOut=classNameOut,
                               CenterMaxShift=self.CenterMaxShift,
                               AlignTiltPairsAdditionalParams=self.AlignTiltPairsAdditionalParams)
                
    def execute_reconstruction(self):
        import os
        import launch_job
        for ref in self.untiltclasslist:
            til_selfile=self.untiltclasslist[ref][2]
            outname=til_selfile.replace('.sel','')
            if (self.ReconstructMethod=='fourier'):
                program='xmipp_reconstruct_fourier'
                command=' -i ' + til_selfile + \
                        ' -o ' + outname+'.vol'
            if not self.ReconstructAdditionalParams=="":
                command += ' ' + self.ReconstructAdditionalParams

            launchJob(program,
                                  command,
                                  self.log,
                                  False,1,1,'')

    def execute_filter(self):
        import os
        import launch_job
        for ref in self.untiltclasslist:
            til_selfile=self.untiltclasslist[ref][2]
            volname=til_selfile.replace('.sel','.vol')
            if os.path.exists(volname):
                filname=volname.replace('.vol','_filtered.vol')
                command=' -o ' + filname + \
                 ' -i ' + volname  + \
                 ' -sampling ' + str(self.PixelSize) + \
                 ' -low_pass ' + str(self.LowPassFilter)
                launchJob("xmipp_fourier_filter",
                                      command,
                                      self.log,
                                      False,1,1,'')

def gatherPairs(log,WorkingDir,ClassNumbers,ClassFile,ExtractRootName,PickingDir):
    from xmipp import MetaData, MDL_IMAGE, MDL_IMAGE_TILTED, \
                      MDL_MICROGRAPH, MDL_MICROGRAPH_TILTED, LEFT_JOIN, MD_APPEND
    MDpairs=MetaData(ExtractRootName+".xmd")
    
    MDtiltAngles=MetaData(os.path.join(PickingDir,"tilted_pairs.xmd"))

    MDuntilted=MetaData(ExtractRootName+"_untilted.xmd")
    MDuntiltedAux=MetaData()
    MDuntiltedAux.setColumnValues(MDL_IMAGE,MDuntilted.getColumnValues(MDL_IMAGE))
    MDuntiltedAux.setColumnValues(MDL_MICROGRAPH,MDuntilted.getColumnValues(MDL_MICROGRAPH))

    MDtilted=MetaData(ExtractRootName+"_tilted.xmd")
    MDtiltedAux=MetaData()
    MDtiltedAux.setColumnValues(MDL_IMAGE_TILTED,MDtilted.getColumnValues(MDL_IMAGE))
    MDtiltedAux.setColumnValues(MDL_MICROGRAPH_TILTED,MDtilted.getColumnValues(MDL_MICROGRAPH))

    MDjoin1=MetaData()
    MDjoin2=MetaData()
    MDjoin3=MetaData()
    MDjoin4=MetaData()
    fnOut=os.path.join(WorkingDir,"classes.xmd")
    for classNo in ClassNumbers:
        MDclass=MetaData("class_%06d@%s"%(classNo,ClassFile))
        MDjoin1.join(MDclass,MDpairs,MDL_IMAGE,MDL_IMAGE,LEFT_JOIN)
        MDjoin2.join(MDjoin1,MDuntiltedAux,MDL_IMAGE,MDL_IMAGE,LEFT_JOIN)
        MDjoin3.join(MDjoin2,MDtiltedAux,MDL_IMAGE_TILTED,MDL_IMAGE_TILTED,LEFT_JOIN)
        MDjoin4.join(MDjoin3,MDtiltAngles,MDL_MICROGRAPH,MDL_MICROGRAPH,LEFT_JOIN)
        MDjoin4.write("class_%06d@%s"%(classNo,fnOut),MD_APPEND)

def reconstructClass(log,WorkingDir,ClassNameIn,ClassNameOut,CenterMaxShift,AlignTiltPairsAdditionalParams):
    runJob(log,"xmipp_image_align_tilt_pairs","-i %s -o %s --max_shift %d %s"%\
           (ClassNameIn,ClassNameOut,CenterMaxShift,AlignTiltPairsAdditionalParams))
