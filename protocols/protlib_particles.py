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
# This library contains some common utilities 
# for all particles related protocols: Extract, Import
from protlib_base import *
from xmipp import MetaData, ImgSize, MDL_ZSCORE
from protlib_utils import runJob, runShowJ
#MDL_CTF_SAMPLING_RATE, MDL_CTF_VOLTAGE, MDL_CTF_DEFOCUSU, MDL_CTF_DEFOCUSV, \
#MDL_CTF_DEFOCUS_ANGLE, MDL_CTF_CS, MDL_CTF_CA, MDL_CTF_Q0, MDL_CTF_K, label2Str, MetaData,\
#MDL_XCOOR, MDL_YCOOR, MDL_PICKING_FAMILY, MDL_PICKING_MICROGRAPH_FAMILY_STATE, MD_APPEND

class ProtParticlesBase(XmippProtocol):
    '''This class will serve as base for Particles related protocols'''
    def __init__(self, protocolName, scriptname, project):
        XmippProtocol.__init__(self, protocolName, scriptname, project)        
        self.Import = 'from protlib_particles import *; '
        
    def insertFilterMaskSteps(self, stackFn):
        if getattr(self, 'DoScale', False):
            self.insertStep('doScale',stack=stackFn,new_size=self.NewSize,Nproc=self.NumberOfMpi)
        if self.DoFourier:
            self.insertStep('doFourier',stack=stackFn,freq_low=self.Freq_low,freq_high=self.Freq_high,freq_decay=self.Freq_decay,Nproc=self.NumberOfMpi)
        if self.DoGaussian:
            self.insertStep('doGaussian',stack=stackFn,freq_sigma=self.Freq_sigma, Nproc=self.NumberOfMpi)
        if getattr(self, 'DoCrop', False):
            self.insertStep('doCrop',stack=stackFn, cropSize=self.CropSize, tmpStack=self.tmpPath('tmpCrop.stk'))
        if self.DoRemoveDust:
            self.insertStep('doRemoveDust',stack=stackFn,threshold=self.DustRemovalThreshold,Nproc=self.NumberOfMpi)
        if self.DoNorm:
            self.insertStep('doNorm',stack=stackFn,normType=self.NormType,bgRadius=self.BackGroundRadius,Nproc=self.NumberOfMpi)
        if self.DoMask:
            if self.Substitute == "value":
                self.Substitute = str(self.SubstituteValue)
            params = "-i %s --substitute %s --mask %s " % (stackFn, self.Substitute, self.MaskType)
            if self.MaskType == 'raised_cosine':
                params += "-%d -%d" % (self.MaskRadius, self.MaskRadius + self.MaskRadiusOuter)
            elif self.MaskType == 'circular':
                params += '-%d' % self.MaskRadius
            else: # from file:
                params += self.MaskFile
            self.insertRunJobStep("xmipp_transform_mask", params)
        
    def visualize(self):
        fn = self.getFilename('images')        
        if not exists(fn):
            showError("Error", "There is no result 'images.xmd' yet")
        else:
            from protlib_utils import runShowJ
            if getattr(self, 'TiltPairs', False):
                runShowJ(fn,extraParams="--mode metadata --render first")
            else:
                runShowJ(fn)

            from protlib_gui_figure import XmippPlotter
            from xmipp import MDL_ZSCORE
            md = MetaData(fn)
            if md.containsLabel(MDL_ZSCORE):
                #MD.sort(MDL_ZSCORE)
                xplotter = XmippPlotter(windowTitle="Zscore particles sorting")
                xplotter.createSubPlot("Particle sorting", "Particle number", "Zscore")
                xplotter.plotMd(md, False, mdLabelY=MDL_ZSCORE)
                xplotter.show()   

def runFourierFilter(log,stack,freq_low,freq_high,freq_decay,Nproc):
    program = "xmipp_transform_filter"
    args = "-i %(stack)s --fourier "
    if freq_low == 0:
        args += "low_pass %(freq_high)f %(freq_decay)f"
    elif freq_high == 0.5:
        args += "high_pass %(freq_low)f %(freq_decay)f"
    else:
        args += "band_pass %(freq_low)f %(freq_high)f %(freq_decay)f"
    runJob(log, program, args % locals())

def runGaussianFilter(log, stack, freq_sigma, Nproc):
    runJob(log,"xmipp_transform_filter","-i %(stack)s --fourier gaussian %(freq_sigma)f" % locals(),Nproc)

def runCrop(log, stack, cropSize, tmpStack):
    runJob(log,"xmipp_transform_window","-i %(stack)s --size %(cropSize)d -o %(tmpStack)s" % locals())
    moveFile(log, tmpStack, stack)

def runResize(log,stack,new_size,Nproc):
    runJob(log,"xmipp_image_resize","-i %(stack)s --fourier %(new_size)d" % locals(),Nproc)

def invert(log,ImagesMd,Nproc):
    runJob(log,'xmipp_image_operate','-i %(ImagesMd)s --mult -1' % locals(),Nproc)

def removeDust(log,ImagesMd,threshold,Nproc):
    runJob(log,'xmipp_transform_filter','-i %(ImagesMd)s --bad_pixels outliers %(threshold)f' % locals(),Nproc)

def runNormalize(log,stack,normType,bgRadius,Nproc):
    program = "xmipp_transform_normalize"
    args = "-i %(stack)s "
    
    if bgRadius <= 0:
        particleSize = ImgSize(stack)[0]
        bgRadius = int(particleSize/2)
    
    if normType=="OldXmipp":
        args += "--method OldXmipp"
    elif normType=="NewXmipp":
        args += "--method NewXmipp --background circle %(bgRadius)d"
    else:
        args += "--method Ramp --background circle %(bgRadius)d"
    runJob(log, program, args % locals(), Nproc)

def doMask(log,stack,maskFile,substitute,Nproc):
    runJob(log,"xmipp_transform_mask","-i %(stack)s --mask binary_file %(maskFile)s %(substitute)s" % locals(),Nproc)

def sortImages(log, ImagesFn):    
    md = MetaData(ImagesFn)
    if not md.isEmpty():
        runJob(log, "xmipp_image_sort_by_statistics","-i %(ImagesFn)s --multivariate --addToInput" % locals())
        md.read(ImagesFn) # Should have ZScore label after runJob
        md.sort(MDL_ZSCORE)
        md.write(ImagesFn)        

