#!/usr/bin/env xmipp_python

#------------------------------------------------------------------------------------------------
#   General script for importing micrographs from EMX format
#
#   Authors: J.M. de la Rosa Trevin, Sept 2013
#            Roberto Marabini
# 

from glob import glob
from protlib_base import *
from xmipp import MetaData, MDL_MICROGRAPH, MDL_MICROGRAPH_TILTED, MDL_SAMPLINGRATE, MDL_CTF_VOLTAGE, \
    MDL_CTF_CS, MDL_CTF_SAMPLING_RATE, MDL_MAGNIFICATION, checkImageFileSize, checkImageCorners
from protlib_filesystem import replaceBasenameExt, renameFile
from protlib_utils import runJob
from protlib_xmipp import redStr, RowMetaData
import math
from emx import *
from emx.emxmapper import *
from protlib_emx import *
from os.path import relpath, join, abspath, dirname, exists


class ProtEmxImportMicrographs(XmippProtocol):
    
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.emx_import_micrographs.name, scriptname, project)
        self.Import = "from protocol_emx_import_micrographs import *"
        #read emx file
        #Class to group EMX objects
        self.emxData = EmxData()
        self.xmlMapper = XmlMapper(self.emxData)
        self.TiltPairs = False
        self.MicrographsMd = self.workingDirPath('micrographs.xmd')
        
        self.propDict = {'SamplingRate': 'pixelSpacing__X',
                         'SphericalAberration': 'cs',
                         'Voltage': 'acceleratingVoltage'}
            
    def createFilenameTemplates(self):
        return {
                 'pos': join('%(ExtraDir)s', '%(micrograph)s.pos'),
                 'config': join('%(ExtraDir)s','config.xmd')
            } 
        
    def defineSteps(self):
        self._loadInfo()
        self.insertStep("validateSchema", verifyfiles=[], 
                        emxFileName=self.EmxFileName
                        )
        self.insertStep("createDir", verifyfiles=[self.ExtraDir], 
                        path=self.ExtraDir)
        _verifyFiles=[]
        #createMicrographs
        micsFn       = self.getFilename('micrographs')
        _verifyFiles.append(micsFn)
        #createAcquisition
        acqFn        = self.getFilename('acquisition')
        _verifyFiles.append(acqFn)
        #createMicroscope
        microscopeFn = self.getFilename('microscope')
        _verifyFiles.append(microscopeFn)
        
        self.insertStep("createMicrographs", verifyfiles=_verifyFiles,
                        emxFileName=self.EmxFileName,
                        binaryFilename=self.binaryFile,
                        micsFileName=micsFn,
                        projectDir=self.projectDir,
                        ctfDir=self.ExtraDir,
                        acqFn=acqFn, 
                        SamplingRate=self.SamplingRate,
                        microscopeFn=microscopeFn, 
                        Voltage=self.Voltage, 
                        SphericalAberration=self.SphericalAberration
                        )
        
        if PARTICLE in self.objDict:
            part = self.objDict[PARTICLE]
            if part.has('centerCoord__X'):
                self.insertStep('createCoordinates', verifyfiles=[],
                                emxFileName=self.EmxFileName,
                                oroot=self.ExtraDir)
            
    def _loadInfo(self):
        #What kind of elements are in the binary file
        emxDir = dirname(self.EmxFileName)
        self.classElement = None
        self.binaryFile = None
        self.object = None
        self.objDict = {}
        
        for classElement in CLASSLIST:
            obj = self.xmlMapper.firstObject(classElement, self.EmxFileName)
            if obj is not None:
                self.objDict[classElement] = obj
                #is the binary file of this type
                binaryFile = join(emxDir, obj.get(FILENAME))
                
                if exists(binaryFile):
                    self.object = obj
                    self.binaryFile = binaryFile
                    self.classElement = classElement
                    
                    for k, v in self.propDict.iteritems():
                        if len(getattr(self, k)) == 0:
                            if self.object.get(v) is not None:
                                setattr(self, k, str(self.object.get(v)))
                    #break
        
    def validate(self):
        errors = []
        # Check that input EMX file exist
        if not exists(self.EmxFileName):
                errors.append("Input EMX file <%s> doesn't exists" % self.EmxFileName)
        else:
            self._loadInfo()   
            
            if self.object is None:
                errors.append('Cannot find any object in EMX file <%s>' % self.EmxFileName)
            else:
                if self.binaryFile is None:
                    errors.append('Cannot find binary data <%s> associated with EMX metadata file' % self.binaryFile)
                for k, v in self.propDict.iteritems():
                        if len(getattr(self, k)) == 0:
                            errors.append('<%s> was left empty and <%s> does not have this property' % (k, self.classElement))
        
        return errors

    def summary(self):
        self._loadInfo()
        summary = ['Input EMX file: [%s]' % self.EmxFileName,
                   'Main class: <%s>' % self.classElement,
                   'Binary file: [%s]' % self.binaryFile]            
        return summary
    
    def visualize(self):
        micsFn = self.getFilename('micrographs')
        if exists(micsFn):
            from protlib_utils import runShowJ
            runShowJ(micsFn)

def loadEmxData(emxFileName):
    """ Given an EMX filename, load data. """
    emxData = EmxData()
    xmlMapper = XmlMapper(emxData)
    xmlMapper.readEMXFile(emxFileName)
    
    return emxData
    
    
def validateSchema(log, emxFileName):
    return
    code, out, err = validateSchema(emxFileName)
    
    if code:
        raise Exception(err)
    
def createMicrographs(log, emxFileName, binaryFilename, micsFileName, projectDir, ctfDir):
    filesPrefix = dirname(emxFileName)
    emxData = loadEmxData(emxFileName)
    Voltage, SphericalAberration, SamplingRate = emxMicsToXmipp(emxData, micsFileName, filesPrefix, ctfDir)
    createAcquisition(log, acqFn, samplingRate )
    createMicroscope(log, microscopeFn, Voltage, SphericalAberration, SamplingRate)
    
def createAcquisition(log, fnacq, SamplingRate):
        # Create the acquisition info file
    mdAcq = RowMetaData()
    mdAcq.setValue(MDL_SAMPLINGRATE, float(SamplingRate))
    mdAcq.write(acqFn)
    
def createMicroscope(log, microscopeFn, Voltage, SphericalAberration, SamplingRate):
    md = RowMetaData()
    md.setValue(MDL_CTF_VOLTAGE, float(Voltage))    
    md.setValue(MDL_CTF_CS, float(SphericalAberration))    
    md.setValue(MDL_CTF_SAMPLING_RATE, float(SamplingRate))
    #md.setValue(MDL_MAGNIFICATION, 60000.0)
    md.write(microscopeFn)
    
def createCoordinates(log, emxFileName, oroot):
    emxData = loadEmxData(emxFileName)
    emxCoordsToXmipp(emxData, oroot)
    
################### Old functions ########################3

def getWdFile(wd, key):
    return getProtocolFilename(key, WorkingDir=wd)

def validateSameMd(fnOut, fn1, fn2):
    EQUAL_ACCURACY = 1.e-8
    '''Check fn1 and fn2 metadata are identical and write 
    a 3rd one with the union of both '''
    #check md1 and m2 are identical  
    print "validateSameMd", fnOut, fn1, fn2
    md1 = MetaData(fn1)
    md2 = MetaData(fn2)
    id1=md1.firstObject()
    id2=md2.firstObject()
    for label in md1.getActiveLabels(): 
        check1 = md1.getValue(label,id1)
        check2 = md2.getValue(label,id2)
        if math.fabs(check1-check2) > EQUAL_ACCURACY:
            raise Exception('Error: %s is not the same for both runs. run1=%f while run2=%f' 
                         % (label2Str(label), check1, check2))
    print "WRINTING: ", fnOut
    md1.write(fnOut)
                
def merge(log, OutWd, InWd1, InWd2):
    ''' Generate import micrograph files from two existing runs.
    Respective working directories are passed as argument '''
    
        #if tilt pair report error since it is not implemented
#        fn1 = PrevRun1.getFilename('TiltPairs')
#        if exists(fn1):
#        	raise Exception('Error (%s): Merging for "%s" is not yet implemented' 
#                        	 % (self.scriptName,'tilt pairs'))
        #check acquisition and microscopy are identical in both runs
    for key in ['acquisition', 'microscope']:
        validateSameMd(getWdFile(OutWd, key), getWdFile(InWd1, key), getWdFile(InWd2, key))

    #open micrographs and union them
    fn1 = getWdFile(InWd1, 'micrographs')
    fn2 = getWdFile(InWd2, 'micrographs')
    md1 = MetaData(fn1)
    md2 = MetaData(fn2)
    #md1.union(md2, MDL_MICROGRAPH)
    md1.unionAll(md2)
    md1.write(getWdFile(OutWd, 'micrographs'))



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
            mdOut.setValue(MDL_MICROGRAPH,micrograph, mdOut.addObject())
    if mdOut.size() > 0:
        mdOut.write(WarningFn)
