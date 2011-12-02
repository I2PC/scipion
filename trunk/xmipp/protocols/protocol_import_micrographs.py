#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# General script for Xmipp-based pre-processing of micrographs
# Author: Carlos Oscar Sorzano, July 2011

from glob import glob
from protlib_base import *
import xmipp
from protlib_filesystem import replaceFilenameExt

class ProtImportMicrographs(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.import_micrographs.name, scriptname, project)
        self.Import = "from protocol_import_micrographs import *"
        self.PatternMicrographs = join(self.DirMicrographs, self.ExtMicrographs)
        
    def defineSteps(self):
        # Create microscope
        self.insertCreateMicroscope()

        # Decide name after preprocessing
        doPreprocess = self.DoPreprocess and (self.DoCrop or self.DoRemoveBadPixelsStddev or self.DoDownsample)
        micrographs = self.getMicrographs()
        micrographs.sort()
        if doPreprocess:
            func = self.insertPreprocessStep
        elif self.CopyMicrographs:
            func = self.insertCopyMicrograph
        else:
            func = lambda i, o: i #Do nothing
        
        filenameDict = {}
        for m in micrographs:
            output = self.workingDirPath(replaceFilenameExt(os.path.basename(m), '.mrc'))
            filenameDict[m] = output
            func(m, output)
        
        self.insertCreateResults(filenameDict)
#        # Create metadata with 
#        summaryFile = self.getFilename('micrographs')
#        self.insertStep('gatherResults', verifyfiles=[summaryFile],
#                        WorkingDir=self.WorkingDir, summaryFile=summaryFile)
#        # Copy tilt pairs description
#        if self.TiltPairs:
#            self.insertStep("createLink", verifyfiles=[self.getFilename('tiltedPairs')],
#                            source=self.PairDescr, dest=self.getFilename('tiltedPairs'))

    def validate(self):
        errors = []
        # Check that there are any micrograph to process
        if len(self.getMicrographs()) == 0:
            errors.append("There are no micrographs to process in " + self.PatternMicrographs)
        return errors

    def summary(self):
        message = []
        message.append("Import of <%d> micrographs from <%s>" % (len(self.getMicrographs()), self.DirMicrographs))
        if self.TiltPairs:
            message.append("Micrographs are in tilt pairs")
        return message
    
    def getMicrographs(self):
        ''' Return a list with micrographs in WorkingDir'''
        return glob(self.PatternMicrographs)
    
    def visualize(self):
        if self.TiltPairs:
            summaryFile = self.getFilename('tiltedPairs')
        else:
            summaryFile = self.getFilename('micrographs')
        if os.path.exists(summaryFile):
            from protlib_utils import runShowJ
            runShowJ(summaryFile)

    def insertCreateMicroscope(self):    
        if self.SamplingRateMode == "From image":
            AngPix = self.SamplingRate
        else:
            AngPix = (10000. * self.ScannedPixelSize * self.Down) / self.Magnification
        fnOut = self.workingDirPath("microscope.xmd")
        self.insertStep("createMicroscope", verifyfiles=[fnOut], fnOut=fnOut, Voltage=self.Voltage,
                        SphericalAberration=self.SphericalAberration,SamplingRate=AngPix,
                        Magnification=self.Magnification)
            
    def insertCreateResults(self, filenameDict):
        micFn = self.getFilename('micrographs')
        print micFn
        vf = [micFn]
        pairMd = ''
        tilted = ''
        if self.TiltPairs:
            tilted = self.getFilename('tiltedPairs')
            vf.append(tilted)
            pairMd = self.PairDescr
        self.insertStep('createResults', verifyfiles=vf, WorkingDir=self.WorkingDir, PairsMd=pairMd, 
                        FilenameDict=filenameDict, MicrographFn=vf[0], TiltedFn=tilted)
        
    def insertCopyMicrograph(self, inputMic, outputMic):
        self.insertStep('copyFile', source=inputMic, dest=outputMic)
        if inputMic.endswith('.raw'):
            self.insertStep('copyFile', source=replaceFilenameExt(inputMic, '.inf'), 
                            dest=replaceFilenameExt(outputMic, '.inf'))
        
    def insertPreprocessStep(self, inputMic, outputMic):
        previousId = XmippProjectDb.FIRST_STEP
        # Crop
        iname = inputMic
        if self.DoCrop:
            previousId = self.insertParallelRunJobStep("xmipp_transform_window", " -i %s -o %s --crop %d -v 0" %(iname,outputMic, self.Crop),
                                                       verifyfiles=[outputMic])
            iname = outputMic
        # Remove bad pixels
        if self.Stddev != -1:
            params = " -i %s --bad_pixels outliers %f -v 0" % (iname, self.Stddev)
            if iname != outputMic:
                params += " -o " + outputMic
                iname = outputMic
            previousId = self.insertParallelRunJobStep("xmipp_transform_filter", params, verifyfiles=[outputMic], parent_step_id=previousId)
        
        # Downsample
        if self.Down != 1:
            tmpFile = outputMic + "_tmp.mrc"
            previousId = self.insertParallelRunJobStep("xmipp_transform_downsample", "-i %s -o %s --step %f --method fourier" % (iname,tmpFile,self.Down),
                                       verifyfiles = [tmpFile], parent_step_id=previousId)
            self.insertParallelStep("renameFile", verifyfiles=[outputMic], parent_step_id=previousId, 
                            source=tmpFile, dest=outputMic)

def createMicroscope(log,fnOut,Voltage,SphericalAberration,SamplingRate,Magnification):
    md = xmipp.MetaData()
    md.setColumnFormat(False)
    objId = md.addObject()
    md.setValue(xmipp.MDL_CTF_VOLTAGE,float(Voltage),objId)    
    md.setValue(xmipp.MDL_CTF_CS,float(SphericalAberration),objId)    
    md.setValue(xmipp.MDL_CTF_SAMPLING_RATE,float(SamplingRate),objId)
    md.setValue(xmipp.MDL_MAGNIFICATION,float(Magnification),objId)
    md.write(fnOut)    

def createResults(log, WorkingDir, PairsMd, FilenameDict, MicrographFn, TiltedFn):
    ''' Create a metadata micrographs.xmd with all micrographs
    and if tilted pairs another one tilted_pairs.xmd'''
    from xmipp import MetaData, MDL_MICROGRAPH, MDL_MICROGRAPH_TILTED
    md = MetaData() 
    for outFn in FilenameDict.values():
        md.setValue(MDL_MICROGRAPH, outFn, md.addObject())
    md.write(MicrographFn)
    
    if len(PairsMd):
        md = MetaData() 
        mdTilted = MetaData(PairsMd)
        for objId in mdTilted:
            u = mdTilted.getValue(MDL_MICROGRAPH, objId)
            t = mdTilted.getValue(MDL_MICROGRAPH_TILTED, objId)
            id2 = md.addObject()
            md.setValue(MDL_MICROGRAPH, u, id2)
            md.setValue(MDL_MICROGRAPH_TILTED, t, id2)
        mdTilted.write(TiltedFn)
