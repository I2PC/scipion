#!/usr/bin/env xmipp_python
'''
#/***************************************************************************
# * Authors:     J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'xmipp@cnb.csic.es'
# ***************************************************************************
 '''

import os
from glob import glob
from subprocess import Popen
from os.path import join, exists, splitext, split
import Tkinter as tk
import tkFont
from xmipp import MetaData
from protlib_base import getWorkingDirFromRunName, getExtendedRunName,\
    XmippProject
from protlib_utils import loadModule, which, runShowJ,\
    runImageJPluginWithResponse, runMaskToolbar
from protlib_gui_ext import centerWindows, changeFontSize, askYesNo, Fonts, registerCommonFonts, \
    showError, showInfo, showBrowseDialog, showWarning, AutoScrollbar, FlashMessage,\
    TaggedText
from protlib_filesystem import getXmippPath, xmippExists, xmippRelpath
from config_protocols import protDict
from config_protocols import FontName, FontSize, MaxHeight, MaxWidth, WrapLenght
from config_protocols import LabelTextColor, SectionTextColor, CitationTextColor
from config_protocols import BgColor, EntryBgColor, SectionBgColor, LabelBgColor, ButtonActiveBgColor, ButtonBgColor                         
from protlib_sql import SqliteDb
from protlib_include import *
from protlib_parser import ProtocolParser
from protlib_xmipp import redStr, greenStr


# This group of function are called Wizards and should help
# to set some parameters  in the GUI, they will receive as parameters
# ProtocolGUI instance and the variable to setup 
# wizard functions should usually set the variable value
# if not, can be used as viewers

def wizardSelectFromList(master, frame, list):
    '''Helper function to select elements from a list '''
    L=len(list)
    if L == 0:
        showWarning("Warning", "No elements to select", parent=master)
        return
    if L == 1:
        return list[0]
    from protlib_gui_ext import ListboxDialog
    d = ListboxDialog(frame, list, selectmode=tk.SINGLE)
    if len(d.result) > 0:
        index = d.result[0]
        return(list[index])
    else:
        return None

def wizardNotFound(self, var):
    showError("Wizard not found", "The wizard <%s> for this parameter has not been found" % var.tags['wizard']
                           , parent=self.master)
    
def wizardDummy(self, var):
    showInfo("Wizard test", "This is only a test on wizards setup", parent=self.master)
    
def wizardShowJ(self, var):
    value = var.getTkValue().strip()
    if len(value):
        runShowJ(var.getTkValue())
    else:
        showWarning("Empty file", "Please select a file to visualize", parent=self.master)
    
def wizardBrowse(self, var):
    if 'file' in var.tags.keys():
        seltype="file"
        filterExt = var.tags['file']
    else:
        seltype="folder"
        filterExt = ''
    files = showBrowseDialog(parent=self.master, seltype=seltype, filter=filterExt)
    if files:
        var.setTkValue(', '.join([xmippRelpath(f) for f in files]))

#Helper function to select Downsampling wizards
def wizardHelperSetDownsampling(self, var, path, filterExt, value, freqs=None, md=None):  
    from protlib_gui_ext import XmippBrowserCTF
    results = showBrowseDialog(path=path, parent=self.master, browser=XmippBrowserCTF,title="Select Downsampling", 
                                    seltype="file", selmode="browse", filter=filterExt, previewDim=256, 
                                    extra={'freqs':freqs, 'downsampling':value, 'previewLabel': 'Micrograph', \
                                           'computingMessage': 'Estimating PSD...', 'md':md}) # a list is returned
    if results:
        self.setVarValue('DownsampleFactor', results[0])
    return results
      
#This wizard is specific for import_micrographs protocol
def wizardBrowseCTF(self, var):    
    importRunName = self.getVarValue('ImportRun')
    downsample = self.getVarValue('DownsampleFactor')
    prot = self.project.getProtocolFromRunName(importRunName)
    path = prot.WorkingDir
    md = MetaData()
    fnMicrographs = join(path, "micrographs.xmd")
    if exists(fnMicrographs):
        md.read(fnMicrographs)
    fnMicrographs = join(path, "tilted_pairs.xmd")
    if exists(fnMicrographs):
        md.read(fnMicrographs)
    wizardHelperSetDownsampling(self, var, '.', None, downsample, md=md)
    
#This wizard is specific for screen_micrographs protocol
#it will help to select downsampling, and frequencies cutoffs
def wizardBrowseCTF2(self, var):
    error = None
    vList = ['LowResolCutoff', 'HighResolCutoff']
    freqs = self.getVarlistValue(vList)
    importRunName = self.getVarValue('ImportRun')
    prot = self.project.getProtocolFromRunName(importRunName)
    path = prot.WorkingDir
    if path and exists(path):
        mdPath = prot.getFilename('micrographs')
        if exists(mdPath):
            from xmipp import MetaData, MDL_MICROGRAPH
            md = MetaData(mdPath)
            if md.size():
                image = md.getValue(MDL_MICROGRAPH, md.firstObject())     
                if image:         
                    filterExt = "*" + splitext(image)[1]
                    value = self.getVarValue('DownsampleFactor')
                    results = wizardHelperSetDownsampling(self, var, path, filterExt, value, freqs, md)
                    if results:
                        self.setVarlistValue(vList, results[1:])
                else:
                    error = "Not micrograph found in metadata <%s>" % mdPath
                    #self.setVarValue('LowResolCutoff', results[1])
                    #self.setVarValue('HighResolCutoff', results[2])
            else:
                error = "Micrograph metadata <%s> is empty" % mdPath
        else:
            error = "Micrograph metadata <%s> doesn't exists" % mdPath
    else:
        error = "Import run <%s> doesn't exists" % str(path)
    if error:
        showWarning("Select Downsampling Wizard", error, self.master)
        return None
    else:
        return results
            
#Select family from extraction run
def wizardChooseFamily(self, var):
    extractionDir = getWorkingDirFromRunName(self.getVarValue('PreviousRun'))
    if not extractionDir:
        showWarning("Warning", "No previous Run has been found", parent=self.master)
        return
    familyList = []
    for file in glob(join(extractionDir, "*_sorted.sel")):
        familyList.append(split(file)[1].replace("_sorted.sel",""))
    if len(familyList)==1:
        var.setTkValue(familyList[0])
    else:
        self.selectFromList(var, familyList)        

def wizardHelperFilter(self, browser, title, **args):
    extra = {'previewLabel': 'Image', 'computingMessage': 'Applying filter...'}
    extra.update(args)
    selfile = self.getVarValue('InSelFile')
    path, filename = split(selfile)
    if not exists(selfile):
        showWarning("Warning", "The input selfile is not a valid file", parent=self.master)
        return
    return showBrowseDialog(path=path, parent=self.master, browser=browser,title=title, 
                            seltype="file", selmode="browse", filter=filename, previewDim=256, extra=extra)        
    
def wizardChooseBandPassFilter(self, var):
    '''Wizard dialog to help choosing Bandpass filter parameters (used in protocol_preprocess_particles) '''
    vList = ['Freq_low','Freq_high','Freq_decay']
    from protlib_gui_ext import XmippBrowserBandpassFilter
    results = wizardHelperFilter(self, XmippBrowserBandpassFilter, "Bandpass Filter", freqs=self.getVarlistValue(vList))
    if results:
        self.setVarlistValue(vList, results)
        
#Choose Gaussian Filter
def wizardChooseGaussianFilter(self, var):
    '''Wizard dialog to help choosing Gaussian filter(in Fourier space) parameters (used in protocol_preprocess_particles) '''
    from protlib_gui_ext import XmippBrowserGaussianFilter
    results = wizardHelperFilter(self, XmippBrowserGaussianFilter, "Gaussian Filter", freqSigma=self.getVarValue('Freq_sigma'))
    if results:
        var.setTkValue(results) #expecting single result            

#Choose Bad pixels wizard
def wizardChooseBadPixelsFilter(self, var):
    from protlib_gui_ext import XmippBrowserBadpixelFilter
    results = wizardHelperFilter(self, XmippBrowserBadpixelFilter, "Gaussian Filter", dustRemovalThreshold=self.getVarValue('DustRemovalThreshold'))
    if results:
        var.setTkValue(results) #expecting single result            

#Design mask wizard
def wizardDesignMask(self, var):
    selfile = self.getVarValue('InSelFile')
    ##workingDir = getWorkingDirFromRunName(self.getVarValue('RunName'))
    from xmipp import MetaData, MDL_IMAGE
    md = MetaData(selfile)
    fnImg = md.getValue(MDL_IMAGE, md.firstObject())
    #runShowJ(fnImg, extraParams="--mask_toolbar")
    runMaskToolbar(fnImg)
    #fnMask=os.path.join(workingDir,"mask.xmp")
#    fnMask = self.project.projectTmpPath("mask.xmp")
#    from protlib_utils import runJavaIJapp
#    msg = runImageJPluginWithResponse("1g", "Masks Tool Bar", "-i %(selfile)s -mask %(fnMask)s" % locals())
#    msg = msg.strip().splitlines()
#    if len(msg)>0:
#        var.setTkValue(fnMask)            

#Select micrograph extension
def wizardMicrographExtension(self,var):
    import fnmatch
    imgExt=['.raw','.tif','.tiff','.mrc','.dm3','.em','.ser','.spi', '.xmp']
    files = []
    currentDir=self.getVarValue('DirMicrographs')
    if currentDir=="":
        currentDir="."

    possibleLocations=[]
    for ext in imgExt:
        for root, dirnames, filenames in os.walk(currentDir):
            if len(fnmatch.filter(filenames, '*'+ext))>0:
                possibleLocations.append(join(root, '*'+ext))
    selected=wizardSelectFromList(self.master, self.frame, possibleLocations)
    if selected is not None:
        dir,ext = split(selected)
        self.setVarValue('DirMicrographs',dir)
        self.setVarValue('ExtMicrographs',ext)
                
#Select Tilt pairs
def wizardTiltPairs(self, var):
    dirMicrographs = self.getVarValue('DirMicrographs')
    extMicrographs = self.getVarValue('ExtMicrographs')
    resultFilename = var.getTkValue()
    uList = []
    tList = []
    from os.path import basename, dirname
    from xmipp import MDL_MICROGRAPH, MDL_MICROGRAPH_TILTED
    
    if exists(resultFilename):
        md = MetaData(resultFilename)
        for id in md:
            tList.append(md.getValue(MDL_MICROGRAPH_TILTED, id))
            path = md.getValue(MDL_MICROGRAPH, id)
            uList.append(basename(path))
            prefix = dirname(path) # This assumes that all micrograph are in the same folder         
    else:
        if len(resultFilename) == 0:
            resultFilename = "tilted_pairs.xmd"
        micrographs = glob(join(dirMicrographs, extMicrographs))
        micrographs.sort()
        for i, m in enumerate(micrographs):
            m = basename(m)
            if i % 2 == 0:
                tList.append(m)
            else:
                uList.append(m)
        prefix = dirMicrographs
    
    from protlib_gui_ext import showTiltPairsDialog
    results = showTiltPairsDialog((uList, tList), self.master)
    if results:
        var.setTkValue(resultFilename)
        uList, tList = results
        md = MetaData()
        for u, t in zip(uList, tList):
            id = md.addObject()
            md.setValue(MDL_MICROGRAPH, join(prefix,u), id)
            md.setValue(MDL_MICROGRAPH_TILTED, join(prefix,t), id)
        md.write(resultFilename)

#Select family from extraction run
def wizardChooseFamilyToExtract(self, var):
    from xmipp import MDL_PICKING_FAMILY, MDL_PICKING_PARTICLE_SIZE, MDL_CTF_MODEL, MDL_SAMPLINGRATE, MDL_SAMPLINGRATE_ORIGINAL
    from protlib_gui_ext import ListboxDialog
    pickingRun = self.getVarValue('PickingRun')
    pickingProt = self.project.getProtocolFromRunName(pickingRun)
    fnFamilies = pickingProt.getFilename("families")  
    if not exists(fnFamilies):
        showWarning("Warning", "No elements to select", parent=self.master)
        return
    md = MetaData(fnFamilies)
    families = [md.getValue(MDL_PICKING_FAMILY, objId) for objId in md]
    if len(families) == 1:
        d = 0
    else:  
        d = ListboxDialog(self.frame, families, selectmode=tk.SINGLE)
        if len(d.result) > 0:
            d = d.result[0]
        else:
            d = None
    if d is not None:
        selectedFamily = families[d]
        var.setValue(selectedFamily)
        for objId in md:
            if md.getValue(MDL_PICKING_FAMILY, objId) == selectedFamily:
                particleSize = md.getValue(MDL_PICKING_PARTICLE_SIZE, objId)
                type = self.getVarValue("DownsampleType")
                mdAcquisition = MetaData(pickingProt.getFilename('acquisition'))
                objId = mdAcquisition.firstObject()
                tsOriginal = tsPicking = mdAcquisition.getValue(MDL_SAMPLINGRATE, objId)
                
                if mdAcquisition.containsLabel(MDL_SAMPLINGRATE_ORIGINAL):
                    tsOriginal = mdAcquisition.getValue(MDL_SAMPLINGRATE_ORIGINAL, objId)
                
                if type == "same as picking":
                    factor = 1
                else:
                    factor = tsPicking / tsOriginal;
                    if type == "other":
                        try:
                            factor /= float(self.getVarValue("DownsampleFactor"))
                        except Exception, e:
                            showWarning("Warning", "Please select valid downsample factor", parent=self.master)
                            return
                particleSize *= factor 
                self.setVarValue("ParticleSize", str(int(particleSize)))
                self.setVarValue("Family", selectedFamily)
        if getattr(pickingProt,'TiltPairs', False):
            self.setVarValue("DoFlip", str(False))
        else:
            md = MetaData(pickingProt.getFilename("micrographs"))
            self.setVarValue("DoFlip", str(md.containsLabel(MDL_CTF_MODEL)))

#Select family from extraction run
def wizardChooseFamilyToExtractSupervised(self, var):
    from xmipp import MDL_PICKING_FAMILY, MDL_PICKING_PARTICLE_SIZE, MDL_CTF_MODEL, MDL_SAMPLINGRATE, MDL_SAMPLINGRATE_ORIGINAL
    from protlib_gui_ext import ListboxDialog
    pickingRun = self.getVarValue('PickingRun')
    pickingProt = self.project.getProtocolFromRunName(pickingRun)
    fnFamilies = pickingProt.getFilename("families")    
    if not exists(fnFamilies):
        showWarning("Warning", "No elements to select", parent=self.master)
        return
    md = MetaData(fnFamilies)
    families = [md.getValue(MDL_PICKING_FAMILY, objId) for objId in md]
    if len(families) == 1:
        d = 0
    else:  
        d = ListboxDialog(self.frame, families, selectmode=tk.SINGLE)
        if len(d.result) > 0:
            d = d.result[0]
        else:
            d = None
    if d is not None:
        selectedFamily = families[d]
        var.setValue(selectedFamily)
        self.setVarValue("Family", selectedFamily)

#This wizard is specific for cl2d protocol
def wizardCL2DNumberOfClasses(self, var):
    fnSel = self.getVarValue('InSelFile')
    if exists(fnSel):
        md = MetaData(fnSel)
        self.setVarValue("NumberOfReferences", int(round(md.size()/200.0)))

#Select micrograph extension
def wizardHelperSetRadii(self, inputVarName, outerVarName, innerVarName=None, ):
    fileList = self.getVarValue(inputVarName).split()
    showInner = innerVarName is not None
    innerRadius = 0
    try:
        if showInner:
            innerRadius = int(self.getVarValue(innerVarName))
        outerRadius = int(self.getVarValue(outerVarName))
    except Exception, e:
        showError("Conversion error", "Error trying to parse integer value from \ninnerRadius: <%(innerVarName)s> or\n outerRadius: <%(outerVarName)s>" % locals() 
                  ,parent=self.master)
        return
    
    from xmipp import Image, FileName, HEADER, MDL_IMAGE
    
    def getFilename():
        if len(fileList) == 0:
            showError("Input error", "File list is empty", parent=self.master)
            return None
        fn = FileName(fileList[0])  
        if not xmippExists(fn):
            showError("Input error", "Filename <%s> doesn't exists" % str(fn),parent=self.master)
            return None
        return fn
    
    fn = getFilename()
    if fn is None:
        return 
    
    if fn.isMetaData():
        md = MetaData(fn)
        fileList = []
        for objId in md:
            fileList.append(md.getValue(MDL_IMAGE, objId))
            if len(fileList) > 10:
                break                      
        fn = getFilename()
        if fn is None:
            return 
            
    if outerRadius < 0:
        if fn.isImage():
            img = Image()
            img.read(fn, HEADER)
            xdim = img.getDimensions()[0]
            outerRadius = xdim / 2 
            self.setVarValue(outerVarName, outerRadius)
    from protlib_gui_ext import XmippBrowserMask
    results = showBrowseDialog(parent=self.master, browser=XmippBrowserMask, title="Select mask radius", allowFilter=False, 
                                    extra={'fileList': fileList, 'outerRadius': outerRadius, 
                                           'innerRadius': innerRadius, 'showInner': showInner})
    if results:
        self.setVarValue(outerVarName, int(results[1]))
        if showInner:
            self.setVarValue(innerVarName, int(results[0]))
        
def wizardSetMaskRadius(self, var):
    wizardHelperSetRadii(self, 'ReferenceFileNames', 'MaskRadius')
    
def wizardSetAlignRadii(self, var):
    wizardHelperSetRadii(self, 'ReferenceFileNames', 'OuterRadius', 'InnerRadius')
    
def wizardSetBackgroundRadius(self, var):
    wizardHelperSetRadii(self, 'InSelFile', var.name)

# This group of functions are called Validator, and should serve
# for validation of user input for each variable
# The return value of each validator should be an error message string 
# or None if not error
def validatorNonEmpty(var):
    if len(var.getTkValue().strip()) == 0:
        return "Input for <%s> is empty" % var.comment
    return None
    
def validatorPathExists(var):
    if not var.satisfiesCondition():
        return None
    err = validatorNonEmpty(var)
    if not err:
        pathList = var.getTkValue().split()
        err = ''
        for p in pathList:
            if not xmippExists(p):
                err += "\n<%s>" % p
        if len(err):
            err = "Following path: %s\ndoesn't exist\nFor input <%s>" % (err, var.comment)
        else: 
            err = None
    return err  

def validatorIsFloat(var):
    err = validatorNonEmpty(var)
    if not err:
        try:
            value = var.getTkValue()
            float(value)
        except ValueError:
            err = "Input value: <%s> for <%s> isn't a valid number" % (value, var.comment)
    return err    

def validatorIsInt(var):
    err = validatorNonEmpty(var)
    if not err:
        value = var.getTkValue()
        try:
            int(value)
        except ValueError:
            err = "Input value: <%s> for <%s> isn't a valid integer" % (value, var.comment)
    return err

def validatorValidRun(var):
    ''' Check if the variable value is a valid posible run '''
    project = XmippProject()
    project.load()
    runList = project.getFinishedRunList(var.getTagValues('run'))
    run = var.getTkValue()
    if not run in runList:
        err = "Run <%s> is not valid" % run
    else:
        err = None
    return err
