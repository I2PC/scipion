# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     J.L. Vilas (jlvilas@cnb.csic.es)
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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
"""
Local Resolution
"""

from pyworkflow import VERSION_1_1
import pyworkflow.protocol.params as params
# from pyworkflow.utils import getExt
from pyworkflow.em.protocol.protocol_3d import ProtAnalysis3D
from pyworkflow.em import ImageHandler
from pyworkflow.object import Float
# from collections import OrderedDict
from pyworkflow.em.data import Volume
import numpy as np
# import pyworkflow.em.metadata as md
       
"""
Bsoft program: blocres
It calculates the local resolution map from to half maps.
The method is based on a local measurement inside a mobile window.
"""

FN_HALF1 = 'half1'
FN_HALF2 = 'half2'
FN_MASKVOL = 'maskvol'
FN_RESOLMAP = 'resolutionMap'

class BsoftProtBlocres(ProtAnalysis3D):
    """ Wrapper around blocres program.
    """
    _label = 'blocres'
    _version = VERSION_1_1

    def __init__(self, **args):
        ProtAnalysis3D.__init__(self, **args)
        self.min_res_init = Float() 
        self.max_res_init = Float()
        self.halfVolumes = True
        
    #--------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputVolume', params.PointerParam,
                      pointerClass='Volume',
                      label="Half 1",
                      help="""Select first half volume to compute the
                      local resolution.
                      """)
        form.addParam('inputVolume2', params.PointerParam,
                      pointerClass='Volume',
                      label="Half 2",
                      help="""Select first half volume to compute the
                      local resolution.
                      """)
        form.addParam('mask', params.PointerParam, allowsNull=True,
                      pointerClass='Volume',
                      label='Mask',
                      help="""Mask file to use for limiting the analysis to a 
                      defined region and level value to use 
                      (optional, default: all but zero).""")
        
        form.addParam('method', params.BooleanParam, default=True,
                      label='Use Box',
                      help="""The local (box) and shell (shell) resolution 
                      calculations are mutually exclusive.""")
        
        form.addParam('box', params.IntParam, default=20,
                      condition = '(method)',
                      label='Box',
                      help="""Kernel size for determining 
                      local resolution (pixels/voxels)""")
        form.addParam('shell', params.IntParam, default=20, 
                      condition ='(not method)',
                      label='Shell',
                      help="""Shell width for determining 
                      radial resolution (pixels/voxels).""")
        
        line = form.addLine('Resolution Criterion')
        line.addParam('resolutionCriterion', params.EnumParam, 
                      choices=['FSC', 'DPR', 'SSNR', 'RAB'],
                      default=0,
                      help="""Resolution criterion: 
                            FSC = Fourier Shell Correlation (default),
        	                DPR = Differential Phase Residual, 
        	                SSNR = Spectral Signal-to-Noise Ratio,
        	                RAB = R(A+B) figure of merit.""")

        line.addParam('cutoff', params.FloatParam,
                      default=0.5,
                      label='Cutoff',
                      help="""Resolution cutoff for the criterion chosen
            	                        (default: FSC: 0.5, DPR: 45, 
            	                        SSNR: 1, RAB: 0.5).""")

        form.addSection(label='Parameters')
        form.addParam('step', params.IntParam, default=1,
                      label='Step',
                      help="""Interval between voxel samples or shells for 
                      resolution analysis (pixels, default: 1)""")
        form.addParam('maxresolution', params.FloatParam, default=2,
                      label='Maximum Resolution',
                      help="""Maximum frequency available in the data (angstrom)""")
        form.addParam('fill', params.IntParam, allowsNull=True,
                      label='Fill',
                      help="""Value to fill the background 
                      (non-masked regions; default 0).
                          """)
        form.addParam('pad', params.EnumParam, 
                      choices=['None', 'Box', 'Shell'],
                      default=1,
                      label='Padding Factor',
                      help="""Resolution box padding factor
                       (0 = none, default: 1 (box) and 0 (shell)).
                          """)
        form.addParam('symmetry', params.StringParam, allowsNull=True,
                      default='',
                      label='Symmetry',
                      help="""Point group symmetry.""")
        form.addParam('smooth', params.BooleanParam, default=True,
                      label='Smooth',
                      help="""Smooth the shell edge.""")

    def _createFilenameTemplates(self):
        """ Centralize how files are called """
        myDict = {
                 FN_HALF1: self._getTmpPath("half1.map"),
                 FN_HALF2: self._getTmpPath("half2.map"),
                 FN_MASKVOL: self._getTmpPath("mask.map"),
                 FN_RESOLMAP: self._getExtraPath("resolutionMap.map")
                 }
        self._updateFilenamesDict(myDict)
        
    #--------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._createFilenameTemplates()
        self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('resolutionStep')
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions -------------------------------   
    def convertInputStep(self):
        # blocres works with .map
        vol1Fn = self.inputVolume.get().getFileName()
        vol2Fn = self.inputVolume2.get().getFileName()
        if self.mask.get().getFileName() != '':
            maskFn = self.mask.get().getFileName()

        self.fnvol1 = self._getFileName("half1")
        self.fnvol2 = self._getFileName("half2")
        ImageHandler().convert(vol1Fn, self.fnvol1)
        ImageHandler().convert(vol2Fn, self.fnvol2)
        if self.mask.get().getFileName() != '':
            self.fnmask = self._getFileName("maskvol")
            ImageHandler().convert(maskFn, self.fnmask)     
        else:      
            self.fnmask = self.maks.get().getFileName()
 
 
    def resolutionStep(self):
        """ blocres parameters. """
        sampling = self.inputVolume.get().getSamplingRate()
        #Actions
        params =  ' -v 1'  # No Verbose
        params += ' -criterion %s' % self.getEnumText("resolutionCriterion")
        if (self.method):
            params += ' -box %i' %self.box.get()
        else:
            params += ' -shell %i' % self.shell.get()

        #Parameters
        params += ' -sampling %f,%f,%f' % (sampling, sampling, sampling)
        params += ' -step %f' % self.step.get()
        params += ' -maxresolution %f' % self.maxresolution.get()
        params += ' -cutoff %f' % self.cutoff.get()
        if self.fill.get() != '':
            params += ' -fill %i' % self.fill.get()

        #Parameters for local resolution
        params += ' -pad %f' % self.pad.get()
        if self.symmetry.get() !='':
            params += ' -symmetry %s' % self.symmetry.get()
        if self.smooth.get():
            params += ' -smooth'
        if self.mask.get().getFileName() !='':
            params += ' -Mask %s' % self.fnmask

        # Input
        # input halves and output map
        params += ' %s %s %s' % (self.fnvol1, self.fnvol2, 
                                 self._getFileName(FN_RESOLMAP))

        self.runJob('blocres', params)


    def createOutputStep(self):
        volume = Volume()
        volume.setFileName(self._getFileName(FN_RESOLMAP))
        volume.setSamplingRate(self.inputVolume.get().getSamplingRate())
        self._defineOutputs(resolution_Volume=volume)
        self._defineSourceRelation(self.inputVolume, volume)
        
        imageFile = self._getFileName(FN_RESOLMAP)
        min_, max_ = self.getMinMax(imageFile)
        self.min_res_init.set(round(min_*100)/100)
        self.max_res_init.set(round(max_*100)/100)
        self._store(self.min_res_init)
        self._store(self.max_res_init)
    
    def getMinMax(self, imageFile):
        img = ImageHandler().read(imageFile)
        imgData = img.getData()
        min_res = round(np.amin(imgData) * 100) / 100
        max_res = round(np.amax(imgData) * 100) / 100
        return min_res, max_res

#--------------------------- INFO functions -----------------------------------
    def _validate(self):
        errors = []
        return errors
    
    def _citations(self):
        cites = ['Cardone2013']
        return cites
    
    def _summary(self):
        summary = []
        summary.append("Highest resolution %.2f Å,   "
                       "Lowest resolution %.2f Å. \n" % (self.min_res_init,
                                                         self.max_res_init))
    
    def _methods(self):
        methods = []
        return methods
