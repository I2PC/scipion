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
from pyworkflow.utils import getExt
from pyworkflow.em.protocol.protocol_3d import ProtAnalysis3D
from pyworkflow.em import ImageHandler
from collections import OrderedDict
from pyworkflow.em.packages.xmipp3 import readSetOfVolumes
import numpy as np
import pyworkflow.em.metadata as md


       
"""
Bsoft program: blocres

Usage: blocres [options] input1.img input2.img output.img
The local (-box) and shell (-shell) resolution calculations are mutually exclusive.

Actions:
-box 20	Kernel size for determining local resolution (pixels/voxels).
-shell 14	Shell width for determining radial resolution (pixels/voxels).
-criterion DPR	Resolution criterion.
	FSC = Fourier Shell Correlation (default)
	DPR = Differential Phase Residual
	SSNR = Spectral Signal-to-Noise Ratio
	RAB = R(A+B) figure of merit


Parameters:
-verbose 7	Verbosity of output.
-sampling 1.5,1.5,1.5	Sampling (A/pixel; a single value can be given).
-origin 0.8,-10,15.7	Set the origin (default from input image).
-step 5	Interval between voxel samples or shells for resolution analysis (pixels, default: 1).
-maxresolution 5	Maximum frequency available in the data (angstrom).
-cutoff 0.3	Resolution cutoff for the criterion chosen
	(default: FSC: 0.5, DPR: 45, SSNR: 1, RAB: 0.5).
	(multiple values can be given for FSC).
-fill 200	Value to fill the background (non-masked regions; default 0).


Parameters for local resolution:
-pad 2	Resolution box padding factor (0 = none, default: 1 (box) and 0 (shell)).
-edge 10	Pixels at the edge of the whole map to exclude from analysis (default: half the size of the resolution box).
-taper edge	Apply a windowing aperture to box volumes before resolution analysis.
	none = do not taper
	edge = circular mask with smooth edges
	hann = hanning window (default)
-nofill	Do not fill values in between (if step > 1).
-symmetry D5	Point group symmetry.


Parameters for shell resolution:
-radii 20,130	Radial limits for shell resolution (default 0,).
-smooth	Smooth the shell edge.


Input:
-Mask mask.tif,2	Mask file to use for limiting the analysis to a defined region
	and level value to use (optional, default: all but zero).


Bsoft version 1.9.0-20150223 (64 bit)
"""

OUTPUT_RESOLUTION_FILE = 'resolutionMap.map'

# RESOLUTION CRITERIA
CRITERION_FSC  = 0
CRITERION_DPR  = 1
CRITERION_SSNR = 2
CRITERION_RAB  = 3

CRITERION_CHOICES = OrderedDict()

CRITERION_CHOICES[CRITERION_FSC]  = 'FSC'
CRITERION_CHOICES[CRITERION_DPR]  = 'DPR'
CRITERION_CHOICES[CRITERION_SSNR] = 'SSNR'
CRITERION_CHOICES[CRITERION_RAB]  = 'RAB'

PAD_NONE  = 0
PAD_BOX  = 1
PAD_SHELL = 2


PAD_CHOICES = OrderedDict()

PAD_CHOICES[PAD_NONE]  = 'None'
PAD_CHOICES[PAD_BOX]  = 'Box'
PAD_CHOICES[PAD_SHELL] = 'Shell'

METHOD_CHOICES = OrderedDict()

METHOD_BOX = 0
METHOD_SHELL = 1

METHOD_CHOICES[METHOD_BOX]  = 'Box'
METHOD_CHOICES[METHOD_SHELL]  = 'Shell'




class BsoftProtBlocres(ProtAnalysis3D):
    """ Wrapper around blocres program.
    """
    _label = 'blocres'
    _version = VERSION_1_1

    def __init__(self, **args):
        ProtAnalysis3D.__init__(self, **args)
        self.halfVolumes = True

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputVolume', params.PointerParam,
                      pointerClass='Volume',
                      label="Half 1",
                      help="""Select first half volume to compute the local resolution.
                      """)
        form.addParam('inputVolume2', params.PointerParam,
                      pointerClass='Volume',
                      label="Half 2",
                      help="""Select first half volume to compute the local resolution.
                      """)
        form.addParam('mask', params.PointerParam, allowsNull=True,
                      pointerClass='Volume',
                      label='Mask',
                      help="""Mask file to use for limiting the analysis to a defined region
                      and level value to use (optional, default: all but zero).""")
        
        form.addParam('method', params.BooleanParam, default=True,
                      label='Use Box',
                      help="""The local (box) and shell (shell) resolution 
                      calculations are mutually exclusive.""")
        
        form.addParam('box', params.IntParam, default=20, condition = '(method)',
                      label='Box',
                      help="""Kernel size for determining local resolution (pixels/voxels)""")
        form.addParam('shell', params.IntParam, default=20, condition ='(not method)',
                      label='Shell',
                      help="""Shell width for determining radial resolution (pixels/voxels).""")
        
        line = form.addLine('Resolution Criterion')
        line.addParam('resolutionCriterion', params.EnumParam, choices=CRITERION_CHOICES.values(),
                      default=CRITERION_FSC,
                      help="""Resolution criterion: FSC = Fourier Shell Correlation (default),
        	                DPR = Differential Phase Residual, SSNR = Spectral Signal-to-Noise Ratio,
        	                RAB = R(A+B) figure of merit.""")

        line.addParam('cutoff', params.FloatParam,
                      default=0.5,
                      label='Cutoff',
                      help="""Resolution cutoff for the criterion chosen
            	                        (default: FSC: 0.5, DPR: 45, SSNR: 1, RAB: 0.5).""")

        form.addSection(label='Parameters')
        form.addParam('step', params.IntParam, default=1,
                      label='Step',
                      help="""Interval between voxel samples or shells for resolution analysis
                           (pixels, default: 1)""")
        form.addParam('maxresolution', params.FloatParam, default=2,
                      label='Maximum Resolution',
                      help="""Maximum frequency available in the data (angstrom)""")
        form.addParam('fill', params.IntParam, allowsNull=True,
                      label='Fill',
                      help="""Value to fill the background (non-masked regions; default 0).
                          """)
        form.addParam('pad', params.EnumParam, choices=PAD_CHOICES.values(),
                      default=PAD_BOX,
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


    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self.fnvol1 = self.inputVolume.get().getFileName()
        self.fnvol2 = self.inputVolume2.get().getFileName()
        if self.mask.get().getFileName() != '':
            self.fnmask = self.mask.get().getFileName()

        convertId = self._insertFunctionStep('convertInputStep')
        ResId = self._insertFunctionStep('resolutionStep', prerequisites=[convertId])
        self._insertFunctionStep('createOutputStep', prerequisites=[ResId])

    #--------------------------- STEPS functions --------------------------------------------   
    def convertInputStep(self):
        # blocres works with .map
        vol1Fn = self.inputVolume.get().getFileName()
        vol2Fn = self.inputVolume2.get().getFileName()
        if self.mask.get().getFileName() != '':
            maskFn = self.mask.get().getFileName()

        extVol1 = getExt(vol1Fn)
        extVol2 = getExt(vol2Fn)

        if self.mask.get().getFileName() != '':
            extMask = getExt(maskFn)

        if (extVol1 != '.map') or (extVol2 != '.map') or (extMask != '.map'):
            self.fnvol1 = self._getTmpPath('half1.map')
            self.fnvol2 = self._getTmpPath('half2.map')
            ImageHandler().convert(vol1Fn, self.fnvol1)
            ImageHandler().convert(vol2Fn, self.fnvol2)
            if self.mask.get().getFileName() != '':
                self.fnmask = self._getTmpPath('mask.map')
                ImageHandler().convert(maskFn, self.fnmask)


    def resolutionStep(self):
        """ blocres parameters. """
        sampling = self.inputVolume.get().getSamplingRate()
        #Actions
        params =  ' -v 1'  # No Verbose
        params += ' -criterion %s' % CRITERION_CHOICES[self.resolutionCriterion.get()]
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
        params += ' %s %s %s' % (self.fnvol1, self.fnvol2, self._getExtraPath(OUTPUT_RESOLUTION_FILE))

        self.runJob('blocres', params)


    def createOutputStep(self):
        fnResolutionMap = self._getExtraPath(OUTPUT_RESOLUTION_FILE)
        self.volumesSet = self._createSetOfVolumes('resolutionVol')
        self.volumesSet.setSamplingRate(self.inputVolume.get().getSamplingRate())
        readSetOfVolumes(fnResolutionMap, self.volumesSet)
        self._defineOutputs(outputVolume=self.volumesSet)
        self._defineSourceRelation(self.inputVolume, self.volumesSet)

#--------------------------- INFO functions -------------------------------------------- 
    def _validate(self):
        errors = []
        return errors
    
    def _citations(self):
        cites = ['Cardone2013']
        return cites
    
    def _summary(self):
        summary = []
        return summary
    
    def _methods(self):
        methods = []
        return methods
