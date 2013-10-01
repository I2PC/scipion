#!/usr/bin/env xmipp_python

#------------------------------------------------------------------------------------------------
#   General script for importing particles from EMX format
#
#   Authors: J.M. de la Rosa Trevin, Sept 2013
#            Roberto Marabini
# 

from os.path import relpath, join, abspath, dirname, exists
from glob import glob
import math

import xmipp
from emx import *
from emx.emxmapper import *

from protlib_base import *
from protlib_filesystem import replaceBasenameExt, renameFile
from protlib_utils import runJob
from protlib_xmipp import redStr, RowMetaData
from protlib_emx import *


class ProtEmxImportParticles(XmippProtocol):
    
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.emx_import_particles.name, scriptname, project)
        self.Import = "from protocol_emx_import_particles import *"
        #read emx file
        #Class to group EMX objects
        self.emxData = EmxData()
        self.xmlMapper = XmlMapper(self.emxData)
        self.TiltPairs = False
        self.MicrographsMd = self.workingDirPath('micrographs.xmd')
        
        self.propDict = {'SamplingRate': 'pixelSpacing__X',
                         }
            
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
        
        imgsFn = self.getFilename('images')
        self.insertStep("createParticles", verifyfiles=[imgsFn],
                        emxFileName=self.EmxFileName,
                        binaryFilename=self.binaryFile,
                        micsFileName=imgsFn,
                        projectDir=self.projectDir,
                        ctfDir=self.ExtraDir,
                        )
        
        acqFn = self.getFilename('acquisition')
        self.insertStep('createAcquisition', verifyfiles=[acqFn],
                        fnOut=acqFn, SamplingRate=self.SamplingRate)

            
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
                errors.append('Canot find any object in EMX file <%s>' % self.EmxFileName)
            else:
                if self.binaryFile is None:
                    errors.append('Canot find binary data <%s> associated with EMX metadata file' % self.binaryFile)
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
        pass



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
    
    
def createParticles(log, emxFileName, binaryFilename, micsFileName, projectDir, ctfDir):
    filesPrefix = dirname(emxFileName)
    emxData = loadEmxData(emxFileName)
    emxParticlesToXmipp(emxData, micsFileName, filesPrefix, ctfDir)
    
    
def createAcquisition(log, fnOut, SamplingRate):
        # Create the acquisition info file
    mdAcq = RowMetaData()
    mdAcq.setValue(xmipp.MDL_SAMPLINGRATE, float(SamplingRate))
    mdAcq.write(fnOut)
    
    
def createMicroscope(log, fnOut, Voltage, SphericalAberration, SamplingRate):
    md = RowMetaData()
    md.setValue(xmipp.MDL_CTF_VOLTAGE, float(Voltage))    
    md.setValue(xmipp.MDL_CTF_CS, float(SphericalAberration))    
    md.setValue(xmipp.MDL_CTF_SAMPLING_RATE, float(SamplingRate))
    md.setValue(xmipp.MDL_MAGNIFICATION, 60000.0)
    md.write(fnOut)
    
def createCoordinates(log, emxFileName, oroot):
    emxData = loadEmxData(emxFileName)
    emxCoordsToXmipp(emxData, oroot)
    
