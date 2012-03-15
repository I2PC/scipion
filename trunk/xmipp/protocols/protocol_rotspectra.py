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
from protlib_utils import runJob
from xmipp import MetaData, MDL_X, MDL_Y

class ProtRotSpectra(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.rotspectra.name, scriptname, project)
        self.Import = 'from protocol_rotspectra import *'    
       
    def defineSteps(self):
        self.Db.insertStep('findCenter', [join(self.WorkingDir, "center2d_center.xmd")],
                           HowCenter=self.HowCenter,
                           Selfile=self.InSelFile, WorkingDir=self.WorkingDir,
                           SpectraInnerRadius=self.SpectraInnerRadius,
                           SpectraOuterRadius=self.SpectraOuterRadius)
        self.Db.insertStep('calculateSpectra', [join(self.WorkingDir, "rotSpectra.xmd")],
                           Selfile=self.InSelFile, WorkingDir=self.WorkingDir,
                           SpectraInnerRadius=self.SpectraInnerRadius,
                           SpectraOuterRadius=self.SpectraOuterRadius,
                           SpectraLowHarmonic=self.SpectraLowHarmonic,
                           SpectraHighHarmonic=self.SpectraHighHarmonic)
        self.Db.insertStep('kerdensom', [self.workingDirPath("results_vectors.xmd"),
                                        self.workingDirPath("results_classes.xmd"),
                                        self.workingDirPath("results_images.xmd")],
                           WorkingDir=self.WorkingDir, SomXdim=self.SomXdim, SomYdim=self.SomYdim,
                           SomReg0=self.SomReg0, SomReg1=self.SomReg1, SomSteps=self.SomSteps,
                           KerdensomExtraCommand=self.KerdensomExtraCommand)

    def validate(self):
        errors = []
        if self.SomReg0 < self.SomReg1:
            errors.append("Regularization must decrease over iterations: Initial regularization must be larger than final")
        return errors
    
    def summary(self):
        message = []
        message.append("Classification of the rotational spectra of " + self.InSelFile + " into a map of size " + str(self.SomYdim) + "x" + str(self.SomXdim))
        if self.getRunState() == SqliteDb.RUN_STARTED:
            lines = []
            for line in open(self.LogPrefix + ".err").readlines():
                if "Training Deterministic Annealing" in line:
                    lines.append(line)
            message.append("Currently at iteration " + str(len(lines)) + " out of " + str(self.SomSteps))
        return message
  
def findCenter(log, HowCenter, Selfile, WorkingDir, SpectraInnerRadius, SpectraOuterRadius):
    if HowCenter == 'Minimize first harmonic':
        if SpectraOuterRadius + 20 > 100:
            R3 = SpectraOuterRadius + (100 - SpectraOuterRadius) / 2
            R4 = 100
        else:
            R3 = SpectraOuterRadius + 10
            R4 = SpectraOuterRadius + 20
        runJob(log, 'xmipp_image_find_center',
                '-i ' + Selfile + \
                ' --oroot ' + WorkingDir + "/center2d" + \
                ' --r1 ' + str(SpectraInnerRadius) + \
                ' --r2 ' + str(SpectraOuterRadius) + \
                ' --r3 ' + str(R3) + ' --r4 ' + str(R4))
    elif HowCenter == 'Use the middle of the image':
        from xmipp import ImgSize
        dims = ImgSize(Selfile)
        MD = MetaData()
        id = MD.addObject()
        MD.setValue(MDL_X, float(dims[0] / 2), id)
        MD.setValue(MDL_Y, float(dims[1] / 2), id)
        MD.write(join(WorkingDir, "center2d_center.xmd"))

def calculateSpectra(log, Selfile, WorkingDir, SpectraInnerRadius, SpectraOuterRadius,
                     SpectraLowHarmonic, SpectraHighHarmonic):
    MD = MetaData(join(WorkingDir, "center2d_center.xmd"))
    id = MD.firstObject()
    xOffset=MD.getValue(MDL_X, id)
    yOffset=MD.getValue(MDL_Y, id)
    
    runJob(log, "xmipp_image_rotational_spectra",
         ' -i ' + Selfile + \
         ' -o ' + join(WorkingDir,"rotSpectra.xmd") + \
         ' --x0 ' + str(xOffset) + \
         ' --y0 ' + str(yOffset) + \
         ' --r1 ' + str(SpectraInnerRadius) + \
         ' --r2 ' + str(SpectraOuterRadius) + \
         ' --low ' + str(SpectraLowHarmonic) + \
         ' --high ' + str(SpectraHighHarmonic))

def kerdensom(log,WorkingDir,SomXdim,SomYdim,SomReg0,SomReg1,SomSteps,KerdensomExtraCommand):
    args='-i '+join(WorkingDir,"rotSpectra.xmd")+\
         ' --oroot '+join(WorkingDir,"results")+\
         ' --xdim ' + str(SomXdim) + \
         ' --ydim ' + str(SomYdim) + \
         ' --deterministic_annealing %f %f %f'%(SomSteps,SomReg0,SomReg1) + \
         ' '+ str(KerdensomExtraCommand)
    runJob(log,"xmipp_classify_kerdensom",args)
    deleteFiles(log, [join(WorkingDir,"rotSpectra.xmd"),join(WorkingDir,"rotSpectra.xmd.raw")], True)
   
