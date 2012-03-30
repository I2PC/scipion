#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# General script for Xmipp-based pre-processing of micrographs
# Author: Carlos Oscar Sorzano, July 2011

from glob import glob
from protlib_base import *
import xmipp
from protlib_filesystem import replaceFilenameExt, renameFile
from protlib_utils import runJob

class ProtImportMicrographs(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.import_micrographs.name, scriptname, project)
        self.Import = "from protocol_import_micrographs import *"
        self.PatternMicrographs = join(self.DirMicrographs, self.ExtMicrographs)
        
    def defineSteps(self):
        # Create microscope
        self.insertCreateMicroscope()
        # Decide name after preprocessing
        doPreprocess = self.DoPreprocess and (self.DoCrop or self.DoRemoveBadPixels or self.DoDownsample)
        micrographs = self.getMicrographs()
        if doPreprocess:
            func = self.insertPreprocessStep
            funcOutput = lambda i: self.workingDirPath(replaceFilenameExt(os.path.basename(i), '.mrc'))
        elif self.CopyMicrographs:
            func = self.insertCopyMicrograph
            funcOutput = lambda i: self.workingDirPath(os.path.basename(i))
        else:
            func = lambda i, o: i #Do nothing
            funcOutput = lambda i: i 
        filenameDict = {}
        for m in micrographs:
            output = funcOutput(m)
            filenameDict[m] = output
            func(m, output)
        # Insert step for result metadatas creation     
        self.insertCreateResults(filenameDict)

    def createFilenameTemplates(self):
        return {
                'microscope': '%(WorkingDir)s/microscope.xmd'
                }
        
    def validate(self):
        errors = []
        # Check that there are any micrograph to process
        if len(self.getMicrographs()) == 0:
            errors.append("There are no micrographs to process in " + self.PatternMicrographs)
        if self.SamplingRateMode == "From image":
            try:
                AngPix = float(self.SamplingRate)
            except:
                errors.append("Sampling rate is not correctly set")
        else:
            try:
                scannedPixelSize=float(self.ScannedPixelSize)
            except:
                errors.append("Sampling rate is not correctly set")
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
            AngPix = float(self.SamplingRate)
        else:
            scannedPixelSize=float(self.ScannedPixelSize)
            if self.DoDownsample:
                AngPix = (10000. * scannedPixelSize * self.DownsampleFactor) / self.Magnification
            else:
                AngPix = (10000. * scannedPixelSize) / self.Magnification
        fnOut = self.getFilename('microscope')
        self.insertStep("createMicroscope", verifyfiles=[fnOut], fnOut=fnOut, Voltage=self.Voltage,
                        SphericalAberration=self.SphericalAberration,SamplingRate=AngPix,
                        Magnification=self.Magnification)
            
    def insertCreateResults(self, filenameDict):
        micFn = self.getFilename('micrographs')
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
        if self.DoRemoveBadPixels:
            params = " -i %s --bad_pixels outliers %f -v 0" % (iname, self.Stddev)
            if iname != outputMic:
                params += " -o " + outputMic
                iname = outputMic
            previousId = self.insertParallelRunJobStep("xmipp_transform_filter", params, verifyfiles=[outputMic], parent_step_id=previousId)
        
        # Downsample
        if self.DoDownsample:
            self.insertParallelStep("doDownsample", verifyfiles=[outputMic], parent_step_id=previousId, 
                            iname=iname, outputMic=outputMic, downsampleFactor=self.DownsampleFactor)

def doDownsample(log,iname,outputMic,downsampleFactor):
    from protlib_filesystem import renameFile
    tmpFile = outputMic + "_tmp.mrc"
    runJob(log,"xmipp_transform_downsample", "-i %s -o %s --step %f --method fourier" % (iname,tmpFile,downsampleFactor))
    renameFile(log,source=tmpFile, dest=outputMic)

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
    micrographs = FilenameDict.values()
    micrographs.sort()
    for m in micrographs:
        md.setValue(MDL_MICROGRAPH, m, md.addObject())
    md.write(MicrographFn)
    
    if len(PairsMd):
        md = MetaData() 
        mdTilted = MetaData(PairsMd)
        for objId in mdTilted:
            u = mdTilted.getValue(MDL_MICROGRAPH, objId)
            t = mdTilted.getValue(MDL_MICROGRAPH_TILTED, objId)
            id2 = md.addObject()
            md.setValue(MDL_MICROGRAPH, FilenameDict[u], id2)
            md.setValue(MDL_MICROGRAPH_TILTED, FilenameDict[t], id2)
        md.write(TiltedFn)
