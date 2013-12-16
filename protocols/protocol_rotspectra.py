#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Xmipp protocol for rotational spectra classification
# using self-organizing maps
#
# Author:Roberto Marabini, March 2007
#        Carlos Oscar Sorzano, January 2011

from os.path import join
from protlib_base import *
from config_protocols import protDict
from protlib_utils import runJob, runShowJ
from xmipp import MetaData, MDL_X, MDL_Y

class ProtRotSpectra(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.rotspectra.name, scriptname, project)
        self.Import = 'from protocol_rotspectra import *'    
       
    def defineSteps(self):
        self.insertStep('createDir',path=self.ExtraDir)
        self.Db.insertStep("linkAcquisitionInfo",InputFile=self.InSelFile,dirDest=self.WorkingDir)
        self.Db.insertStep('findCenter', [join(self.ExtraDir, "center2d_center.xmd")],
                           HowCenter=self.HowCenter,
                           Selfile=self.InSelFile, ExtraDir=self.ExtraDir,
                           SpectraInnerRadius=self.SpectraInnerRadius,
                           SpectraOuterRadius=self.SpectraOuterRadius)
        self.Db.insertStep('calculateSpectra', [join(self.ExtraDir, "rotSpectra.xmd")],
                           Selfile=self.InSelFile, ExtraDir=self.ExtraDir,
                           SpectraInnerRadius=self.SpectraInnerRadius,
                           SpectraOuterRadius=self.SpectraOuterRadius,
                           SpectraLowHarmonic=self.SpectraLowHarmonic,
                           SpectraHighHarmonic=self.SpectraHighHarmonic)
        self.Db.insertStep('kerdensom', [self.extraPath("kerdensom_vectors.xmd"),
                                        self.extraPath("kerdensom_classes.xmd"),
                                        self.extraPath("kerdensom_images.xmd")],
                           ExtraDir=self.ExtraDir, SomXdim=self.SomXdim, SomYdim=self.SomYdim,
                           SomReg0=self.SomReg0, SomReg1=self.SomReg1, SomSteps=self.SomSteps,
                           KerdensomExtraCommand=self.KerdensomExtraCommand)
        fnClasses=self.workingDirPath('classes.xmd')
        self.Db.insertStep('createLink',verifyfiles=[fnClasses],source=self.extraPath("kerdensom_classes.xmd"),dest=fnClasses)
        fnImages=self.workingDirPath('images.xmd')
        self.Db.insertStep('createLink',verifyfiles=[fnImages],source=self.extraPath("kerdensom_images.xmd"),dest=fnImages)

    def validate(self):
        errors = []
        if self.SomReg0 < self.SomReg1:
            errors.append("Regularization must decrease over iterations:")
            errors.append("    Initial regularization must be larger than final")
        return errors
    
    def summary(self):
        message = ["Classification of the rotational spectra"]
        message.append("  Input classes: [%s]" % self.InSelFile)
        message.append("  Map size: <%(SomYdim)d> x <%(SomXdim)d>" % self.ParamsDict)
        
        if self.getRunState() == SqliteDb.RUN_STARTED:
            lines = []
            for line in open(self.LogPrefix + ".err").readlines():
                if "Training Deterministic Annealing" in line:
                    lines.append(line)
            message.append("Currently at iteration <%d> out of <%d>" % (len(lines), self.SomSteps))
        return message
    
    def papers(self):
        papers=[]
        papers.append('Pascual-Montano, Ultramic (2000) [http://www.ncbi.nlm.nih.gov/pubmed/10896143]')
        papers.append('Pascual-Montano, JSB (2001) [http://www.ncbi.nlm.nih.gov/pubmed/11472094]')
        papers.append('Pascual-Montano, JSB (2002) [http://www.ncbi.nlm.nih.gov/pubmed/12160707]')
        return papers

    def visualize(self):
        runShowJ("classes@%s" % self.extraPath("kerdensom_classes.xmd"), extraParams="--mode rotspectra --columns %d" % self.SomXdim)
  
def findCenter(log, HowCenter, Selfile, ExtraDir, SpectraInnerRadius, SpectraOuterRadius):
    if HowCenter == 'Minimize first harmonic':
        if SpectraOuterRadius + 20 > 100:
            R3 = SpectraOuterRadius + (100 - SpectraOuterRadius) / 2
            R4 = 100
        else:
            R3 = SpectraOuterRadius + 10
            R4 = SpectraOuterRadius + 20
        runJob(log, 'xmipp_image_find_center',
                '-i ' + Selfile + \
                ' --oroot ' + ExtraDir + "/center2d" + \
                ' --r1 ' + str(SpectraInnerRadius) + \
                ' --r2 ' + str(SpectraOuterRadius) + \
                ' --r3 ' + str(R3) + ' --r4 ' + str(R4))
    elif HowCenter == 'Use the middle of the image':
        from xmipp import MetaDataInfo
        dims = MetaDataInfo(Selfile)
        MD = MetaData()
        id = MD.addObject()
        MD.setValue(MDL_X, float(dims[0] / 2), id)
        MD.setValue(MDL_Y, float(dims[1] / 2), id)
        MD.write(join(ExtraDir, "center2d_center.xmd"))

def calculateSpectra(log, Selfile, ExtraDir, SpectraInnerRadius, SpectraOuterRadius,
                     SpectraLowHarmonic, SpectraHighHarmonic):
    MD = MetaData(join(ExtraDir, "center2d_center.xmd"))
    id = MD.firstObject()
    xOffset=MD.getValue(MDL_X, id)
    yOffset=MD.getValue(MDL_Y, id)
    
    runJob(log, "xmipp_image_rotational_spectra",
         ' -i ' + Selfile + \
         ' -o ' + join(ExtraDir,"rotSpectra.xmd") + \
         ' --x0 ' + str(xOffset) + \
         ' --y0 ' + str(yOffset) + \
         ' --r1 ' + str(SpectraInnerRadius) + \
         ' --r2 ' + str(SpectraOuterRadius) + \
         ' --low ' + str(SpectraLowHarmonic) + \
         ' --high ' + str(SpectraHighHarmonic))

def kerdensom(log,ExtraDir,SomXdim,SomYdim,SomReg0,SomReg1,SomSteps,KerdensomExtraCommand):
    args='-i '+join(ExtraDir,"rotSpectra.xmd")+\
         ' --oroot '+join(ExtraDir,"kerdensom")+\
         ' --xdim ' + str(SomXdim) + \
         ' --ydim ' + str(SomYdim) + \
         ' --deterministic_annealing %f %f %f'%(SomSteps,SomReg0,SomReg1) + \
         ' '+ str(KerdensomExtraCommand)
    runJob(log,"xmipp_classify_kerdensom",args)
    deleteFiles(log, [os.path.join(ExtraDir,"rotSpectra.xmd"),os.path.join(ExtraDir,"rotSpectra.vec")], True)
    
