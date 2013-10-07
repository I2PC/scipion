#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# General script for Xmipp-based pre-processing of micrographs
# Author: Joaquin Oton, Oct 2013

from glob import glob
from protlib_base import *
from xmipp import MetaData, MDL_MICROGRAPH, MDL_MICROGRAPH_TILTED, MDL_SAMPLINGRATE, MDL_CTF_VOLTAGE, \
    MDL_CTF_CS, MDL_CTF_SAMPLING_RATE, MDL_MAGNIFICATION, checkImageFileSize, checkImageCorners, MD_APPEND
from protlib_filesystem import replaceBasenameExt, renameFile, join
from protlib_utils import runJob
from protlib_xmipp import redStr, RowMetaData
import math


class ProtXrayImport(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.xray_import.name, scriptname, project)
        self.Import = "from protocol_xray_import import *"
        self.PatternMicrographs = join(self.DirMicrographs, self.ExtMicrographs)
        self.MicrographsMd = self.getFilename('micrographs')
            
    def defineSteps(self):
        
        for t in self.getTomograms():
            self.insertTomogramsStep(t)


    def validate(self):
        errors = []

        return errors

    def summary(self):
        message = []
            
        return message
    
    def getTomograms(self):
        ''' Return a list with micrographs in WorkingDir'''
        return glob(join(self.DirTomograms, '*.hdf*'))
    
    def visualize(self):
        pass
        
    def insertTomogramStep(self, tomogram):
        params = '-i ' + tomogram
        self.insertRunJobStep("xmipp_image_header", params)
        


################### Old functions #####################################3  

def createResults(log, WorkingDir, PairsMd, FilenameDict, MicrographFn, PixelSize):
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
        md.write("micrographPairs@"+MicrographFn,MD_APPEND)

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
