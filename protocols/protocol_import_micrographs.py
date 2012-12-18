#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# General script for Xmipp-based pre-processing of micrographs
# Author: Carlos Oscar Sorzano, July 2011

from glob import glob
from protlib_base import *
from xmipp import MetaData, MDL_MICROGRAPH, MDL_MICROGRAPH_TILTED, MDL_SAMPLINGRATE, MDL_CTF_VOLTAGE, \
    MDL_CTF_CS, MDL_CTF_SAMPLING_RATE, MDL_MAGNIFICATION, checkImageFileSize, checkImageCorners
from protlib_filesystem import replaceBasenameExt, renameFile
from protlib_utils import runJob
from protlib_xmipp import redStr
import math

class ProtImportMicrographs(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.import_micrographs.name, scriptname, project)
        self.Import = "from protocol_import_micrographs import *"
        self.PatternMicrographs = join(self.DirMicrographs, self.ExtMicrographs)
        if self.TiltPairs:
            self.MicrographsMd = self.getFilename('tilted_pairs')
        else:
            self.MicrographsMd = self.getFilename('micrographs')
        
    def defineSteps(self):
        # Create microscope
        self.insertCreateMicroscope()
        # Decide name after preprocessing
        doPreprocess = self.DoPreprocess and (self.DoCrop or self.DoRemoveBadPixels or self.DoLog)
        micrographs = self.getMicrographs()
        if doPreprocess:
            self.insertStep("createDir",verifyfiles=[self.ExtraDir],path=self.ExtraDir)
            func = self.insertPreprocessStep
            funcOutput = lambda i: os.path.join(self.ExtraDir,replaceBasenameExt(i, '.mrc'))
        elif self.CopyMicrographs:
            self.insertStep("createDir",verifyfiles=[self.ExtraDir],path=self.ExtraDir)
            func = self.insertCopyMicrograph
            funcOutput = lambda i: os.path.join(self.ExtraDir,os.path.basename(i))
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

        if self.DoCheckBorders:
            self.insertStep("checkBorders",MicrographFn=self.getFilename('micrographs'),WarningFn=self.workingDirPath('warnings.xmd'))
        
    def validate(self):
        errors = []
        # Check that there are any micrograph to process
        micrographList = self.getMicrographs()
        if len(micrographList) == 0:
            errors.append("There are no micrographs to process in " + self.PatternMicrographs)
        else:
            for micrograph in micrographList:
                try:
                    if not checkImageFileSize(micrograph):
                        errors.append(micrograph+" seems to be corrupted")
                except Exception:
                    errors.append(micrograph+" seems to be corrupted")

        return errors

    def summary(self):
        message = []
        message.append("Import of <%d> micrographs from [%s]" % (len(self.getMicrographs()), self.DirMicrographs))
        if self.TiltPairs:
            message.append("Micrographs are in tilt pairs")
        fnAcquisition=self.getFilename('acquisition')
        if os.path.exists(fnAcquisition):
            MD=MetaData(fnAcquisition)
            message.append("Sampling rate=<"+str(MD.getValue(MDL_SAMPLINGRATE,MD.firstObject()))+" A/pix>")
        fnWarning=self.workingDirPath('warnings.xmd')
        if os.path.exists(fnWarning):
            message.append(redStr("WARNINGS ON THE BORDERS OF SOME MICROGRAPHS:"))
            message.append("[%s]"%fnWarning)
        return message
    
    def getMicrographs(self):
        ''' Return a list with micrographs in WorkingDir'''
        return glob(self.PatternMicrographs)
    
    def visualize(self):
        if self.TiltPairs:
            summaryFile = self.getFilename('tilted_pairs')
        else:
            summaryFile = self.getFilename('micrographs')
        if os.path.exists(summaryFile):
            from protlib_utils import runShowJ
            runShowJ(summaryFile)

    def insertCreateMicroscope(self):    
        if self.SamplingRateMode == "From image":
            self.AngPix = float(self.SamplingRate)
        else:
            scannedPixelSize=float(self.ScannedPixelSize)
            self.AngPix = (10000. * scannedPixelSize) / self.Magnification
        fnOut = self.getFilename('microscope')
        self.insertStep("createMicroscope", verifyfiles=[fnOut], fnOut=fnOut, Voltage=self.Voltage,
                        SphericalAberration=self.SphericalAberration,SamplingRate=self.AngPix,
                        Magnification=self.Magnification)
            
    def insertCreateResults(self, filenameDict):
        micFn = self.getFilename('micrographs')
        vf = [micFn]
        pairMd = ''
        tilted = ''
        if self.TiltPairs:
            tilted = self.getFilename('tilted_pairs')
            vf.append(tilted)
            pairMd = self.PairDescr
        self.insertStep('createResults', verifyfiles=vf, WorkingDir=self.WorkingDir, PairsMd=pairMd, 
                        FilenameDict=filenameDict, MicrographFn=vf[0], TiltedFn=tilted,
                        PixelSize=self.AngPix)
        
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
        # Take logarithm
        if self.DoLog:
            params = " -i %s --log --fa %f --fb %f --fc %f" % (iname,
                                                               self.log_a,
                                                               self.log_b,
                                                               self.log_c)
            if iname != outputMic:
                params += " -o " + outputMic
                iname = outputMic
            previousId = self.insertParallelRunJobStep("xmipp_transform_filter", params, 
                                                       verifyfiles=[outputMic], 
                                                       parent_step_id=previousId)
        # Remove bad pixels
        if self.DoRemoveBadPixels:
            params = " -i %s --bad_pixels outliers %f -v 0" % (iname, self.Stddev)
            if iname != outputMic:
                params += " -o " + outputMic
                iname = outputMic
            previousId = self.insertParallelRunJobStep("xmipp_transform_filter", params, verifyfiles=[outputMic], parent_step_id=previousId)
        
    def merge(self, PrevRun1, PrevRun2):

        try:
	    #if tilt pair report error since it is not implemented
	    file_name1 = PrevRun1.getFilename('TiltPairs')
	    if exists(file_name1):
        	raise Exception('Error (%s): Merging for "%s" is not yet implemented' 
                        	 % (self.scriptName,'tilt pairs'))

	    #check acquisition is identical in both runs 
	    file_name1 = PrevRun1.getFilename('acquisition')
	    file_name2 = PrevRun2.getFilename('acquisition')
	    mdAcquisition1 = MetaData(file_name1)
	    mdAcquisition2 = MetaData(file_name2)
	    id1=mdAcquisition1.firstObject()
	    id2=mdAcquisition2.firstObject()
	    for label in md.getActiveLabels() : 
        	check1 = mdAcquisition1.getValue(label,id1)
        	check2 = mdAcquisition2.getValue(label,id2)
        	if math.fabs(check1-check2) > EQUAL_ACCURACY:
        	    raise Exception('Error (%s): %s is not the same for both runs. run1=%f while run2=%f' 
                        	 % (self.scriptName,label2Str(label),check1,check2))
	    #check Microscope is identical in both runs 
	    file_name1 = PrevRun1.getFilename('microscope')
	    file_name2 = PrevRun2.getFilename('microscope')
	    mdMicroscope1 = MetaData(file_name1)
	    mdMicroscope2 = MetaData(file_name2)
	    id1=mdMicroscope1.firstObject()
	    id2=mdMicroscope2.firstObject()
	    for label in md.getActiveLabels() : 
        	check1 = mdMicroscope1.getValue(label,id1)
        	check2 = mdMicroscope2.getValue(label,id2)
        	if math.fabs(check1-check2) > EQUAL_ACCURACY:
        	    raise Exception('Error (%s): %s is not the same for both runs. run1=%f while run2=%f' 
                        	 % (self.scriptName,label2Str(label),check1,check2))

	    #open micrographs and union them
	    file_name1 = PrevRun1.getFilename('micrographs')
	    file_name2 = PrevRun2.getFilename('micrographs')
	    mdMicrographs1 = MetaData(file_name1)
	    mdMicrographs2 = MetaData(file_name2)
	    mdMicrographs1.union (mdMicrographs2, MDL_MICROGRAPH)

	    #now save everything in the right place
	    mdMicroscope1.write(self.getFilename('microscope'))
	    mdAcquisition1.write(self.getFilename('acquisition'))
	    mdMicrographs1.write(self.getFilename('micrographs'))
        except Exception, ex:
               return str(ex)
	return ''
def createMicroscope(log,fnOut,Voltage,SphericalAberration,SamplingRate,Magnification):
    md = MetaData()
    md.setColumnFormat(False)
    objId = md.addObject()
    md.setValue(MDL_CTF_VOLTAGE,float(Voltage),objId)    
    md.setValue(MDL_CTF_CS,float(SphericalAberration),objId)    
    md.setValue(MDL_CTF_SAMPLING_RATE,float(SamplingRate),objId)
    md.setValue(MDL_MAGNIFICATION,float(Magnification),objId)
    md.write(fnOut)    

def createResults(log, WorkingDir, PairsMd, FilenameDict, MicrographFn, TiltedFn, PixelSize):
    ''' Create a metadata micrographs.xmd with all micrographs
    and if tilted pairs another one tilted_pairs.xmd'''
    md = MetaData()
    micrographs = FilenameDict.values()
    micrographs.sort()
    for m in micrographs:
        md.setValue(MDL_MICROGRAPH, m, md.addObject())
    md.write("micrographs@"+MicrographFn)
    mdAcquisition = MetaData()
    mdAcquisition.setValue(MDL_SAMPLINGRATE,float(PixelSize),mdAcquisition.addObject())
    mdAcquisition.write(os.path.join(WorkingDir,"acquisition_info.xmd"))
    
    if len(PairsMd):
        md.clear()
        mdTilted = MetaData(PairsMd)
        for objId in mdTilted:
            u = mdTilted.getValue(MDL_MICROGRAPH, objId)
            t = mdTilted.getValue(MDL_MICROGRAPH_TILTED, objId)
            id2 = md.addObject()
            md.setValue(MDL_MICROGRAPH, FilenameDict[u], id2)
            md.setValue(MDL_MICROGRAPH_TILTED, FilenameDict[t], id2)
        md.write("micrographPairs@"+TiltedFn)

def checkBorders(log,MicrographFn,WarningFn):
    md = MetaData()
    md.read("micrographs@"+MicrographFn)
    
    mdOut = MetaData()
    for objId in md:
        micrograph=md.getValue(MDL_MICROGRAPH,objId)
        if not checkImageCorners(micrograph):
            mdOut.setValue(MDL_MICROGRAPH,micrograph,mdOut.addObject())
    if mdOut.size()>0:
        mdOut.write(WarningFn)
