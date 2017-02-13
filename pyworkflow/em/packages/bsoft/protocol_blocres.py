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
from .convert import ImageHandler
from collections import OrderedDict
from pyworkflow.em.packages.xmipp3 import readSetOfVolumes

       
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


class BsoftProtBlocres(ProtAnalysis3D):
    """ Wrapper around blocres program.
    """
    _label = 'blocres'
    _version = VERSION_1_1

    #def __init__(self, **args):
        #ProtAnalysis3D.__init__(self, **args)

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('half1', params.PointerParam,
                      pointerClass='Volume',
                      label="Half 1",
                      help="""Select first half volume to compute the local resolution.
                      """)
        form.addParam('half2', params.PointerParam,
                      pointerClass='Volume',
                      label="Half 2",
                      help="""Select first half volume to compute the local resolution.
                      """)
        form.addParam('mask', params.PointerParam,
                      pointerClass='Volume',
                      label='Mask',
                      help="""Mask file to use for limiting the analysis to a defined region
                      and level value to use (optional, default: all but zero).""")
        form.addParam('box', params.IntParam, default=20,
                      label='Box',
                      help="""Kernel size for determining local resolution (pixels/voxels)""")
        form.addParam('shell', params.IntParam, default=20,
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
        form.addParam('fill', params.FloatParam, default=0,
                      label='Fill',
                      help="""Value to fill the background (non-masked regions; default 0).
                          """)
        form.addParam('pad', params.EnumParam, choices=PAD_CHOICES.values(),
                      default=PAD_BOX,
                      label='Padding Factor',
                      help="""Resolution box padding factor (0 = none, default: 1 (box) and 0 (shell)).
                          """)
        form.addParam('symmetry', params.StringParam, default='C1',
                      label='Symmetry',
                      help="""Point group symmetry.""")
        form.addParam('smooth', params.BooleanParam, default='true',
                      label='Smooth',
                      help="""Smooth the shell edge.""")


    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        # Insert processing steps

        self.vol1Fn = self.half1.get().getFileName()
        self.vol2Fn = self.half2.get().getFileName()
        self.maskFn = self.mask.get().getFileName()

        self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('resolutionStep')
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions --------------------------------------------   
    def convertInputStep(self):
        # we need to put all images into a single stack to ease the call of bsoft programs
        extVol1 = getExt(self.vol1Fn)
        extVol2 = getExt(self.vol2Fn)
        extMask = getExt(self.maskFn)

        if (extVol1 != '.map') or (extVol2 != '.map') or (extMask != '.map'):
            fnvol1 = 'half1.map'
            fnvol2 = 'half2.map'
            fnmask = 'mask.map'

            ImageHandler().convert(self.vol1Fn, self._getTmpPath(fnvol1))
            ImageHandler().convert(self.vol2Fn, self._getTmpPath(fnvol2))
            ImageHandler().convert(self.maskFn, self._getTmpPath(fnmask))

    def resolutionStep(self):
        """ Util function used by wizard. """
        print self.vol1Fn
        print self.vol2Fn
        print self._getExtraPath(OUTPUT_RESOLUTION_FILE)
        sampling = self.half1.get().getSamplingRate()
        #Actions
        params =  ' -v 1'  # No Verbose
        params += ' -criterion %s' % CRITERION_CHOICES[self.resolutionCriterion.get()]
        params += ' -box %f' %self.box.get()
        params += ' -shell %f' % self.shell.get()

        #Parameters
        params += ' -sampling %f,%f,%f' % (sampling, sampling, sampling)
        params += ' -step %f' % self.step.get()
        params += ' -maxresolution %f' % self.maxresolution.get()
        params += ' -cutoff %f' % self.cutoff.get()
        params += ' -fill %f' % self.fill.get()

        #Parameters for local resolution
        params += ' -pad %f' % self.pad.get()
        params += ' -symmetry %s' % self.symmetry.get()
        if self.smooth.get():
            params += ' -smooth'

        # Input
        # input halves and output map
        params += ' %s %s %s' % (self.vol1Fn, self.vol2Fn, self._getExtraPath(OUTPUT_RESOLUTION_FILE))

        self.runJob('blocres', params)
        
    def createOutputStep(self):
        fnResolutionMap = self._getExtraPath(OUTPUT_RESOLUTION_FILE)
        self.volumesSet = self._createSetOfVolumes('resolutionVol')
        self.volumesSet.setSamplingRate(self.half1.get().getSamplingRate())
        readSetOfVolumes(fnResolutionMap, self.volumesSet)
        self._defineOutputs(outputVolume=self.volumesSet)
        self._defineSourceRelation(self.half1, self.volumesSet)

#--------------------------- INFO functions -------------------------------------------- 
    def _validate(self):
        errors = []
        return errors
    
    def _citations(self):
        cites = []
        return cites
    
    def _summary(self):
        summary = []
        return summary
    
    def _methods(self):
        methods = []
        return methods
