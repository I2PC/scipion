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
    runImageJPluginWithResponse, runMaskToolbar, \
    getComponentFromVector
from protlib_gui_ext import centerWindows, changeFontSize, askYesNo, Fonts, registerCommonFonts, \
    showError, showInfo, showBrowseDialog, showWarning, AutoScrollbar, FlashMessage,\
    TaggedText, openFile
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

def postSelectRunProjmatch(gui):
    #get projection matching dir
    runName = gui.getVarValue('ImportRun')
    prot = gui.project.getProtocolFromRunName(runName)
    #get last iteration in projection matching
    iterationNo       = int(prot.NumberOfIterations)
    gui.setVarValue('iterationNo', prot.NumberOfIterations)
    #gui.setVarValue('doMask',prot.DoMask)
    gui.setVarValue('InnerRadius',prot.InnerRadius)
    gui.setVarValue('OuterRadius',prot.OuterRadius)
    angSamplingRateDeg = getComponentFromVector(prot.AngSamplingRateDeg,iterationNo - 1)
    gui.setVarValue('AngSamplingRateDeg',angSamplingRateDeg)
    MaxChangeInAngles = getComponentFromVector(prot.MaxChangeInAngles,iterationNo - 1)
    gui.setVarValue('MaxChangeInAngles',MaxChangeInAngles)
    gui.setVarValue('SymmetryGroup',prot.SymmetryGroup)
    gui.setVarValue('CTFDatName',prot.CTFDatName)
    #gui.setVarValue('doCTFCorrection',prot.DoCtfCorrection)
    

def wizardNotFound(gui, var):
    showError("Wizard not found", "The wizard <%s> for this parameter has not been found" % var.tags['wizard']
                           , parent=gui.master)
    
def wizardDummy(gui, var):
    showInfo("Wizard test", "This is only a test on wizards setup", parent=gui.master)
    
def wizardShowJ(gui, var):
    value = var.getTkValue().strip()
    if len(value):
        openFile(var.getTkValue())
    else:
        showWarning("Empty file", "Please select a file to visualize", parent=gui.master)
    
def wizardBrowse(gui, var):
    if 'file' in var.tags.keys():
        seltype="file"
        filterExt = var.tags['file']
    else:
        seltype="folder"
        filterExt = ''
    files = showBrowseDialog(parent=gui.master, seltype=seltype, filter=filterExt)
    if files:
        var.setTkValue(', '.join([xmippRelpath(f) for f in files]))

#Helper function to select Downsampling wizards
def wizardHelperSetDownsampling(gui, var, path, filterExt, value, freqs=None, md=None):  
    from protlib_gui_ext import XmippBrowserCTF
    results = showBrowseDialog(path=path, parent=gui.master, browser=XmippBrowserCTF,title="Select Downsampling", 
                                    seltype="file", selmode="browse", filter=filterExt, previewDim=256, 
                                    extra={'freqs':freqs, 'downsampling':value, 'previewLabel': 'Micrograph', \
                                           'computingMessage': 'Estimating PSD...', 'md':md}) # a list is returned
    if results:
        gui.setVarValue('DownsampleFactor', results[0])
    return results
      
#This wizard is specific for import_micrographs protocol
def wizardBrowseCTF(gui, var):    
    importRunName = gui.getVarValue('ImportRun')
    downsample = gui.getVarValue('DownsampleFactor')
    prot = gui.project.getProtocolFromRunName(importRunName)
    path = prot.WorkingDir
    md = MetaData()
    fnMicrographs = join(path, "micrographs.xmd")
    if exists(fnMicrographs):
        md.read(fnMicrographs)
    fnMicrographs = join(path, "tilted_pairs.xmd")
    if exists(fnMicrographs):
        md.read(fnMicrographs)
    wizardHelperSetDownsampling(gui, var, '.', None, downsample, md=md)
    
#This wizard is specific for screen_micrographs protocol
#it will help to select downsampling, and frequencies cutoffs
def wizardBrowseCTF2(gui, var):
    error = None
    vList = ['LowResolCutoff', 'HighResolCutoff']
    freqs = gui.getVarlistValue(vList)
    importRunName = gui.getVarValue('ImportRun')
    prot = gui.project.getProtocolFromRunName(importRunName)
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
                    value = gui.getVarValue('DownsampleFactor')
                    results = wizardHelperSetDownsampling(gui, var, path, filterExt, value, freqs, md)
                    if results:
                        gui.setVarlistValue(vList, results[1:])
                else:
                    error = "Not micrograph found in metadata <%s>" % mdPath
                    #gui.setVarValue('LowResolCutoff', results[1])
                    #gui.setVarValue('HighResolCutoff', results[2])
            else:
                error = "Micrograph metadata <%s> is empty" % mdPath
        else:
            error = "Micrograph metadata <%s> doesn't exists" % mdPath
    else:
        error = "Import run <%s> doesn't exists" % str(path)
    if error:
        showWarning("Select Downsampling Wizard", error, gui.master)
        return None
    else:
        return results
            
def wizardHelperFilter(gui, browser, title, **args):
    extra = {'previewLabel': 'Image', 'computingMessage': 'Applying filter...'}
    extra.update(args)
    selfile = gui.getVarValue('InSelFile')
    path, filename = split(selfile)
    if not exists(selfile):
        showWarning("Warning", "The input selfile is not a valid file", parent=gui.master)
        return
    return showBrowseDialog(path=path, parent=gui.master, browser=browser,title=title, 
                            seltype="file", selmode="browse", filter=filename, previewDim=256, extra=extra)        
    
def wizardChooseBandPassFilter(gui, var):
    '''Wizard dialog to help choosing Bandpass filter parameters (used in protocol_preprocess_particles) '''
    vList = ['Freq_low','Freq_high','Freq_decay']
    from protlib_gui_ext import XmippBrowserBandpassFilter
    results = wizardHelperFilter(gui, XmippBrowserBandpassFilter, "Bandpass Filter", freqs=gui.getVarlistValue(vList))
    if results:
        gui.setVarlistValue(vList, results)
        
#Choose Gaussian Filter
def wizardChooseGaussianFilter(gui, var):
    '''Wizard dialog to help choosing Gaussian filter (in Fourier space) parameters (used in protocol_preprocess_particles) '''
    from protlib_gui_ext import XmippBrowserGaussianFilter
    results = wizardHelperFilter(gui, XmippBrowserGaussianFilter, "Gaussian Filter", freqSigma=gui.getVarValue('Freq_sigma'))
    if results:
        var.setTkValue(results) #expecting single result            

#Choose Gaussian Filter
def wizardChooseRealGaussianFilter(gui, var):
    '''Wizard dialog to help choosing Real Gaussian filter (in real space) parameters (used in protocol_preprocess_particles) '''
    from protlib_gui_ext import XmippBrowserRealGaussianFilter
    results = wizardHelperFilter(gui, XmippBrowserRealGaussianFilter, "Real Gaussian Filter", realSigma=gui.getVarValue('Real_sigma'))
    if results:
        var.setTkValue(results) #expecting single result            

#Choose Bad pixels wizard
def wizardChooseBadPixelsFilter(gui, var):
    from protlib_gui_ext import XmippBrowserBadpixelFilter
    results = wizardHelperFilter(gui, XmippBrowserBadpixelFilter, "Bad pixels Filter", dustRemovalThreshold=gui.getVarValue('DustRemovalThreshold'))
    if results:
        var.setTkValue(results) #expecting single result            

#Design mask wizard
def wizardDesignMask(gui, var):
    selfile = gui.getVarValue('InSelFile')
    ##workingDir = getWorkingDirFromRunName(gui.getVarValue('RunName'))
    from xmipp import MetaData, MDL_IMAGE
    md = MetaData(selfile)
    fnImg = md.getValue(MDL_IMAGE, md.firstObject())
    #runShowJ(fnImg, extraParams="--mask_toolbar")
    runMaskToolbar(fnImg)
    #fnMask=os.path.join(workingDir,"mask.xmp")
#    fnMask = gui.project.projectTmpPath("mask.xmp")
#    from protlib_utils import runJavaIJapp
#    msg = runImageJPluginWithResponse("1g", "Masks Tool Bar", "-i %(selfile)s -mask %(fnMask)s" % locals())
#    msg = msg.strip().splitlines()
#    if len(msg)>0:
#        var.setTkValue(fnMask)            

#Select micrograph extension
def wizardMicrographExtension(gui,var):
    import fnmatch
    imgExt=['.raw','.tif','.tiff','.mrc','.dm3','.em','.ser','.spi', '.xmp']
    files = []
    currentDir=gui.getVarValue('DirMicrographs')
    if currentDir=="":
        currentDir="."

    possibleLocations=[]
    for ext in imgExt:
        for root, dirnames, filenames in os.walk(currentDir):
            if len(fnmatch.filter(filenames, '*'+ext))>0:
                possibleLocations.append(join(root, '*'+ext))
    selected=wizardSelectFromList(gui.master, gui.frame, possibleLocations)
    if selected is not None:
        dir,ext = split(selected)
        gui.setVarValue('DirMicrographs',dir)
        gui.setVarValue('ExtMicrographs',ext)
                
#Select Tilt pairs
def wizardTiltPairs(gui, var):
    dirMicrographs = gui.getVarValue('DirMicrographs')
    extMicrographs = gui.getVarValue('ExtMicrographs')
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
    results = showTiltPairsDialog((uList, tList), gui.master)
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
def wizardChooseSizeToExtract(gui, var):
    from xmipp import MDL_PICKING_PARTICLE_SIZE, MDL_CTF_MODEL, MDL_SAMPLINGRATE, MDL_SAMPLINGRATE_ORIGINAL
    from protlib_gui_ext import ListboxDialog
    pickingRun = gui.getVarValue('PickingRun')
    pickingProt = gui.project.getProtocolFromRunName(pickingRun)
    fnConfig = pickingProt.getFilename("config")  
    if not exists(fnConfig):
        showWarning("Warning", "No elements to select", parent=gui.master)
        return
    md = MetaData(fnConfig)
    particleSize = md.getValue(MDL_PICKING_PARTICLE_SIZE, md.firstObject())
    type = gui.getVarValue("DownsampleType")
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
                factor /= float(gui.getVarValue("DownsampleFactor"))
            except Exception, e:
                showWarning("Warning", "Please select valid downsample factor", parent=gui.master)
                return
    particleSize *= factor 
    gui.setVarValue("ParticleSize", str(int(particleSize)))
    if getattr(pickingProt,'TiltPairs', False):
        gui.setVarValue("DoFlip", str(False))
    else:
        md = MetaData(pickingProt.getFilename("micrographs"))
        gui.setVarValue("DoFlip", str(md.containsLabel(MDL_CTF_MODEL)))

#This wizard is specific for cl2d protocol
def wizardCL2DNumberOfClasses(gui, var):
    fnSel = gui.getVarValue('InSelFile')
    if exists(fnSel):
        md = MetaData(fnSel)
        gui.setVarValue("NumberOfReferences", int(round(md.size()/200.0)))

#Select micrograph extension
def wizardHelperSetRadii(gui, inputVarName, outerVarName, innerVarName=None, ):
    fileList = gui.getVarValue(inputVarName).split()
    showInner = innerVarName is not None
    innerRadius = 0
    try:
        if showInner:
            innerRadius = int(gui.getVarValue(innerVarName))
        outerRadius = int(gui.getVarValue(outerVarName))
    except Exception, e:
        showError("Conversion error", "Error trying to parse integer value from \ninnerRadius: <%(innerVarName)s> or\n outerRadius: <%(outerVarName)s>" % locals() 
                  ,parent=gui.master)
        return
    from xmipp import Image, FileName, HEADER, MDL_IMAGE
    
    def getFilename():
        if len(fileList) == 0:
            showError("Input error", "File list is empty", parent=gui.master)
            return None
        fn = FileName(fileList[0])  
        if not xmippExists(fn):
            showError("Input error", "Filename <%s> doesn't exists" % str(fn),parent=gui.master)
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
            #ROB: must convert to str first or PyString_AsString will complain
            img.read(str(fn), HEADER)
            xdim = img.getDimensions()[0]
            outerRadius = xdim / 2 
            gui.setVarValue(outerVarName, outerRadius)
    from protlib_gui_ext import XmippBrowserMask
    results = showBrowseDialog(parent=gui.master, browser=XmippBrowserMask, title="Select mask radius", allowFilter=False, 
                                    extra={'fileList': fileList, 'outerRadius': outerRadius, 
                                           'innerRadius': innerRadius, 'showInner': showInner})
    if results:
        gui.setVarValue(outerVarName, int(results[1]))
        if showInner:
            gui.setVarValue(innerVarName, int(results[0]))
        
def wizardSetMaskRadius(gui, var):
    wizardHelperSetRadii(gui, 'ReferenceFileNames', 'MaskRadius')
    
def wizardSetMaskRadiusAlign(gui, var):
    wizardHelperSetRadii(gui, 'ReferenceVolume', 'MaskRadius')
    
def wizardSetMaskRadiusPreprocess(gui, var):
    wizardHelperSetRadii(gui, 'InModel', 'MaskRadius')
    
def wizardSetMaskRadiusPreprocess2(gui, var):
    wizardHelperSetRadii(gui, 'InModel', 'MaskRadiusNormalize')
    
def wizardSetAlignRadii(gui, var):
    wizardHelperSetRadii(gui, 'ReferenceFileNames', 'OuterRadius', 'InnerRadius')

def wizardSetBackgroundRadius(gui, var):
    wizardHelperSetRadii(gui, 'InSelFile', var.name)

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

def validatorIsFloatOrEmpty(var):
    """ Same as validator float, but allows empty values. """
    err = None
    try:
        value = var.getTkValue()
        if len(value):
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
