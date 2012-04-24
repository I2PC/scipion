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
from os.path import join, relpath, exists
import Tkinter as tk
import tkFont

from protlib_base import getWorkingDirFromRunName, getExtendedRunName,\
    XmippProject
from protlib_utils import loadModule, which, runShowJ
from protlib_gui_ext import centerWindows, changeFontSize, askYesNo, Fonts, registerCommonFonts, \
    showError, showInfo, showBrowseDialog, showWarning, AutoScrollbar, FlashMessage,\
    TaggedText
from protlib_filesystem import getXmippPath, xmippExists
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
        var.setTkValue(', '.join([relpath(f) for f in files]))

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
    from xmipp import MetaData
    importRunName = self.getVarValue('ImportRun')
    downsample = self.getVarValue('DownsampleFactor')
    prot = self.project.getProtocolFromRunName(importRunName)
    path = prot.WorkingDir
    MD=MetaData()
    fnMicrographs=os.path.join(path,"micrographs.xmd")
    if os.path.exists(fnMicrographs):
        MD.read(fnMicrographs)
    fnMicrographs=os.path.join(path,"tilted_pairs.xmd")
    if os.path.exists(fnMicrographs):
        MD.read(fnMicrographs)
    wizardHelperSetDownsampling(self, var, None, None, downsample, md=MD)
    
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
                    filterExt = "*" + os.path.splitext(image)[1]
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
        familyList.append(os.path.split(file)[1].replace("_sorted.sel",""))
    if len(familyList)==1:
        var.setTkValue(familyList[0])
    else:
        self.selectFromList(var, familyList)        

def wizardHelperFilter(self, browser, title, **args):
    extra = {'previewLabel': 'Image', 'computingMessage': 'Applying filter...'}
    extra.update(args)
    selfile = self.getVarValue('InSelFile')
    path, filename = os.path.split(selfile)
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
def wizardChooseCropSizeFilter(self, var):
    from protlib_gui_ext import XmippBrowserCropSizeFilter
    results = wizardHelperFilter(self, XmippBrowserCropSizeFilter, "Crop Size", cropSize=self.getVarValue('CropSize'))
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
    workingDir = getWorkingDirFromRunName(self.getVarValue('RunName'))
    #fnMask=os.path.join(workingDir,"mask.xmp")
    fnMask = self.project.projectTmpPath("mask.xmp")
    from protlib_utils import runJavaIJapp
    msg = runJavaIJapp("1g", "XmippMaskDesignWizard", "-i %(selfile)s -mask %(fnMask)s" % locals())
    msg = msg.strip().splitlines()
    if len(msg)>0:
        var.setTkValue(fnMask)            

#Select micrograph extension
def wizardMicrographExtension(self,var):
    import fnmatch
    imgExt=['.raw','.tif','.mrc','.dm3','.ser','.spi']
    files = []
    currentDir=self.getVarValue('DirMicrographs')
    if currentDir=="":
        currentDir="."

    possibleLocations=[]
    for ext in imgExt:
        for root, dirnames, filenames in os.walk(currentDir):
            if len(fnmatch.filter(filenames, '*'+ext))>0:
                possibleLocations.append(os.path.join(root, '*'+ext))
    selected=wizardSelectFromList(self.master, self.frame, possibleLocations)
    if selected is not None:
        dir,ext=os.path.split(selected)
        self.setVarValue('DirMicrographs',dir)
        self.setVarValue('ExtMicrographs',ext)
                
#Select Tilt pairs
def wizardTiltPairs(self, var):
    dirMicrographs = self.getVarValue('DirMicrographs')
    extMicrographs = self.getVarValue('ExtMicrographs')
    resultFilename = var.getTkValue()
    uList = []
    tList = []
    from xmipp import MetaData, MDL_MICROGRAPH, MDL_MICROGRAPH_TILTED
    
    if exists(resultFilename):
        md = MetaData(resultFilename)
        for id in md:
            tList.append(md.getValue(MDL_MICROGRAPH_TILTED, id))
            uList.append(md.getValue(MDL_MICROGRAPH, id))                
    else:
        if len(resultFilename) == 0:
            resultFilename = "tilted_pairs.xmd"
        micrographs = glob(join(dirMicrographs, extMicrographs))
        micrographs.sort()
        for i, m in enumerate(micrographs):
            m = os.path.basename(m)
            if i % 2 == 0:
                tList.append(m)
            else:
                uList.append(m)
    
    from protlib_gui_ext import showTiltPairsDialog
    results = showTiltPairsDialog((uList, tList), self.master)
    if results:
        var.setTkValue(resultFilename)
        uList, tList = results
        md = MetaData()
        for u, t in zip(uList, tList):
            id = md.addObject()
            md.setValue(MDL_MICROGRAPH, join(dirMicrographs,u), id)
            md.setValue(MDL_MICROGRAPH_TILTED, join(dirMicrographs,t), id)
        md.write(resultFilename)

#Select family from extraction run
def wizardChooseFamilyToExtract(self, var):
    from xmipp import MetaData, MDL_PICKING_FAMILY, MDL_PICKING_PARTICLE_SIZE, MDL_CTFMODEL
    from protlib_gui_ext import ListboxDialog
    pickingRun = self.getVarValue('PickingRun')
    pickingProt = self.project.getProtocolFromRunName(pickingRun)
    fnFamilies = pickingProt.workingDirPath("families.xmd")    
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
                self.setVarValue("ParticleSize", str(particleSize))
                self.setVarValue("Family",selectedFamily)
        if hasattr(pickingProt,'TiltPairs'):
            if pickingProt.TiltPairs:
                self.setVarValue("DoFlip", str(False))
        else:
            md=MetaData(pickingProt.getFilename("micrographs"))
            self.setVarValue("DoFlip", str(md.containsLabel(MDL_CTFMODEL)))

#This wizard is specific for cl2d protocol
def wizardCL2DNumberOfClasses(self, var):
    from xmipp import MetaData
    fnSel = self.getVarValue('InSelFile')
    if os.path.exists(fnSel):
        MD=MetaData(fnSel)
        self.setVarValue("NumberOfReferences", int(round(MD.size()/200.0)))

#Select micrograph extension
def wizardProjMatchRadius(self,var):
    from xmipp import SingleImgSize, FileName
    volumeList = self.getVarValue('ReferenceFileNames')
    fnVol = FileName(volumeList.split(' ')[0])
    if fnVol.exists() and fnVol.isImage():
        (Xdim, Ydim, Zdim, Ndim) = SingleImgSize(fnVol)
        self.setVarValue("MaskRadius", str(Xdim/2))

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
