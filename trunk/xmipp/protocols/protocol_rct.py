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
from protlib_utils import getListFromRangeString, runJob, runShowJ

class ProtRCT(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.rct.name, scriptname, project)
        self.Import = 'from protocol_rct import *'

    def defineSteps(self):
        extractionProt = self.getProtocolFromRunName(self.ExtractionRun)
        pickingDir = getWorkingDirFromRunName(extractionProt.PickingRun)
        self.insertStep("createLink2",filename="acquisition_info.xmd",dirSrc=pickingDir,destDir=self.WorkingDir)

        classNumbers = getListFromRangeString(self.SelectedClasses)
        extractRootName = os.path.join(getWorkingDirFromRunName(self.ExtractionRun),"Default")
        self.insertStep('gatherPairs',verifyfiles=[os.path.join(self.WorkingDir,"classes.xmd")],
                               WorkingDir=self.WorkingDir,ClassNumbers=classNumbers,
                               ClassFile=self.ClassifMd,ExtractRootName=extractRootName,
                               PickingDir=pickingDir)
        from xmipp import MetaData, MDL_REF, MDL_IMAGE
        MDrepresentatives=MetaData("classes@"+self.ClassifMd)
        for classNo in classNumbers:
            # Locate class representative
            classRepresentative=""
            for id in MDrepresentatives:
                if MDrepresentatives.getValue(MDL_REF,id)==classNo:
                    classRepresentative=MDrepresentatives.getValue(MDL_IMAGE,id)
            if classRepresentative!="":
                classNameIn="class_%06d@%s"%(classNo,self.workingDirPath("classes.xmd"))
                classNameOut=self.workingDirPath("rct_%06d.xmd"%classNo)
                classVolumeOut=self.workingDirPath("rct_%06d.vol"%classNo)
                self.insertParallelStep('reconstructClass',verifyfiles=[classNameOut,classVolumeOut],
                                   WorkingDir=self.WorkingDir,
                                   ClassNameIn=classNameIn, ClassNameOut=classNameOut,
                                   ClassRepresentative=classRepresentative,
                                   CenterMaxShift=self.CenterMaxShift,
                                   ThinObject=self.ThinObject,
                                   SkipTiltedTranslations=self.SkipTiltedTranslations,
                                   ClassVolumeOut=classVolumeOut,
                                   ReconstructAdditionalParams=self.ReconstructAdditionalParams,
                                   DoLowPassFilter=self.DoLowPassFilter,LowPassFilter=self.LowPassFilter)
            else:
                self.Log.error("Cannot find representative for class %d"%classNo)
                
    def summary(self):
        message=[]
        classNumbers = getListFromRangeString(self.SelectedClasses)
        message.append("Random conical tilt reconstruction of "+str(len(classNumbers))+" classes of "+self.ClassifMd)
        import glob
        reconstructedFiles=glob.glob(self.WorkingDir+"/rct_??????.vol")
        if not reconstructedFiles:
            message.append("No reconstruction has been generated")
        else:
            message.append(str(len(reconstructedFiles))+" reconstructions have finished")
        return message
    
    def validate(self):
        errors = []
        classNumbers = getListFromRangeString(self.SelectedClasses)
        from xmipp import MetaData
        for classNo in classNumbers:
            try:
                blockName="class_%06d"%classNo
                mD=MetaData(blockName+"@"+self.ClassifMd)
            except:
                errors.append(blockName+" cannot be found at "+self.ClassifMd)
        if self.CenterMaxShift>100 or self.CenterMaxShift<0:
            errors.append("Maximum shift must be in the range 0-100")
        if self.LowPassFilter>0.5 or self.LowPassFilter<0:
            errors.append("Low pass filter must be in the range 0-0.5")
        return errors    

    def visualize(self):
        classNumbers = getListFromRangeString(self.SelectedClasses)
        for classNo in classNumbers:
            fnVol=self.workingDirPath("rct_%06d.vol"%classNo)
            if os.path.exists(fnVol):
                runShowJ(fnVol)

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

def reconstructClass(log,WorkingDir,ClassNameIn,ClassNameOut,ClassRepresentative,
                     CenterMaxShift,ThinObject,SkipTiltedTranslations,ClassVolumeOut,
                     ReconstructAdditionalParams,DoLowPassFilter,LowPassFilter):
    alignmentParameters="-i %s -o %s --ref %s --max_shift %d"%(ClassNameIn,ClassNameOut,ClassRepresentative,CenterMaxShift)
    if ThinObject:
        alignmentParameters+=" --do_stretch"
    if SkipTiltedTranslations:
        alignmentParameters+=" --do_not_align_tilted"
    runJob(log,"xmipp_image_align_tilt_pairs",alignmentParameters)
    runJob(log,"xmipp_reconstruct_fourier","-i %s -o %s %s"%(ClassNameOut,ClassVolumeOut,ReconstructAdditionalParams))
    if DoLowPassFilter:
        if os.path.exists(ClassVolumeOut):
            filteredVolume=ClassVolumeOut.replace('.vol','_filtered.vol')
            runJob(log,"xmipp_transform_filter","-i %s -o %s --fourier low_pass %f"%(ClassVolumeOut,filteredVolume,LowPassFilter))
