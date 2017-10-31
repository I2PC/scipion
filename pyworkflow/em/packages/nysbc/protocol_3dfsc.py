# **************************************************************************
# *
# * Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk)
# *
# * MRC Laboratory of Molecular Biology (MRC-LMB)
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
Protocol wrapper around the 3D FSC tool for directional
resolution estimation
"""

import os

import pyworkflow.protocol.params as params
from pyworkflow.em.protocol import ProtAnalysis3D
from pyworkflow.em.convert import ImageHandler
from pyworkflow.em.data import Volume
from pyworkflow.utils import join, basename
from convert import getEnviron, findSphericity


class Prot3DFSC(ProtAnalysis3D):
    """
    3D FSC is software tool for quantifying directional
    resolution using 3D Fourier shell correlation volumes.
     
    Find more information at https://github.com/nysbc/Anisotropy
    """
    _label = '3D FSC'

    INPUT_HELP = """ Required input volumes for 3D FSC:
        1. First half map of 3D reconstruction. Can be masked or unmasked.
        2. Second half map of 3D reconstruction. Can be masked or unmasked.
        3. Full map of 3D reconstruction. Can be masked or unmasked, sharpened or unsharpened.
    """
    
    def __init__(self, **kwargs):
        ProtAnalysis3D.__init__(self, **kwargs)

    def _initialize(self):
        """ This function is mean to be called after the
        working dir for the protocol have been set. (maybe after recovery from mapper)
        """
        self._createFilenameTemplates()

    def _createFilenameTemplates(self):
        """ Centralize how files are called for iterations and references. """
        myDict = {
                  'input_volFn': self._getExtraPath('volume_full.mrc'),
                  'input_half1Fn': self._getExtraPath('volume_half1.mrc'),
                  'input_half2Fn': self._getExtraPath('volume_half2.mrc'),
                  'input_maskFn': self._getExtraPath('mask.mrc'),
                  'out_histogram': self._getExtraPath('Results_3D-FSC/histogram.png'),
                  'out_plot3DFSC': self._getExtraPath('Results_3D-FSC/Plots3D-FSC.jpg'),
                  'out_plotFT': self._getExtraPath('Results_3D-FSC/FTPlot3D-FSC.jpg'),
                  'out_vol3DFSC': self._getExtraPath('Results_3D-FSC/3D-FSC.mrc'),
                  'out_vol3DFSC-th': self._getExtraPath('Results_3D-FSC/3D-FSC_Thresholded.mrc'),
                  'out_vol3DFSC-thbin': self._getExtraPath('Results_3D-FSC/3D-FSC_ThresholdedBinarized.mrc'),
                  'out_cmdChimera': self._getExtraPath('Results_3D-FSC/Chimera/3DFSCPlot_Chimera.cmd'),
                  'out_globalFSC': self._getExtraPath('Results_3D-FSC/ResEM3D-FSCOutglobalFSC.csv')
                  }

        self._updateFilenamesDict(myDict)

    #--------------------------- DEFINE param functions ------------------------

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputVolume', params.PointerParam,
                      pointerClass='Volume',
                      label="Input volume", important=True,
                      help=self.INPUT_HELP)
        form.addParam('volumeHalf1', params.PointerParam,
                      label="Volume half 1", important=True,
                      pointerClass='Volume',
                      help=self.INPUT_HELP)
        form.addParam('volumeHalf2', params.PointerParam,
                      pointerClass='Volume',
                      label="Volume half 2", important=True,
                      help=self.INPUT_HELP)

        form.addParam('applyMask', params.BooleanParam, default=False,
                      label="Mask input volume?",
                      help='If given, it would be used to mask the half maps '
                           'during 3DFSC generation and analysis.')
        form.addParam('maskVolume', params.PointerParam, label="Mask volume",
                      pointerClass='VolumeMask', condition="applyMask",
                      help='Select a volume to apply as a mask.')

        group = form.addGroup('Extra parameters')
        group.addParam('dTheta', params.FloatParam, default=20,
                       label='Angle of cone (deg)',
                       help='Angle of cone to be used for 3D FSC sampling in '
                            'degrees. Default is 20 degrees.')
        group.addParam('fscCutoff', params.FloatParam, default=0.143,
                       label='FSC cutoff',
                       help='FSC cutoff criterion. 0.143 is default.')
        group.addParam('thrSph', params.FloatParam, default=0.5,
                       label='Sphericity threshold',
                       help='Threshold value for 3DFSC volume for calculating '
                            'sphericity. 0.5 is default.')
        group.addParam('hpFilter', params.FloatParam, default=200,
                       label='High-pass filter (A)',
                       help='High-pass filter for thresholding in Angstrom. '
                            'Prevents small dips in directional FSCs at low spatial '
                            'frequency due to noise from messing up the '
                            'thresholding step. Decrease if you see a huge wedge '
                            'missing from your thresholded 3DFSC volume. 200 '
                            'Angstroms is default.')
        group.addParam('numThr', params.IntParam, default=1,
                       label='Number of threshold for sphericity',
                       help='Calculate sphericities at different threshold cutoffs '
                            'to determine sphericity deviation across spatial '
                            'frequencies. This can be useful to evaluate possible '
                            'effects of overfitting or improperly assigned '
                            'orientations.')

    #--------------------------- INSERT steps functions ------------------------
    
    def _insertAllSteps(self):
        # Insert processing steps
        self._initialize()
        self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('run3DFSCStep')
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions -------------------------------
    
    def convertInputStep(self):
        """ Convert input volumes to .mrc as expected by 3DFSC."""
        ih = ImageHandler()

        ih.convert(self.volumeHalf1.get().getLocation(),
                   self._getFileName('input_half1Fn'))
        ih.convert(self.volumeHalf2.get().getLocation(),
                   self._getFileName('input_half2Fn'))
        ih.convert(self.inputVolume.get().getLocation(),
                   self._getFileName('input_volFn'))
        if self.maskVolume.get() is not None:
            ih.convert(self.maskVolume.get().getLocation(),
                       self._getFileName('input_maskFn'))

    def run3DFSCStep(self):
        """ Call ResMap.py with the appropriate parameters. """
        args = self._getArgs()
        param = ' '.join(['%s=%s' % (k, str(v)) for k, v in args.iteritems()])
        program = self._getProgram()
        self.info("**Running:** %s %s" % (program, param))

        f = open(self._getExtraPath('script.sh'), "w")
        line = """#!/bin/bash
source ~/eman22.bashrc
unset PYTHONPATH
source activate fsc
python %s %s
source deactivate
""" % (program, param)
        f.write(line)
        f.close()

        self.runJob('/bin/bash ./script.sh', '', cwd=self._getExtraPath(), env=getEnviron())
        #self.runJob('unset PYTHONPATH && source activate fsc && python %s %s && source deactivate' %
        #            (program, param), '', cwd=self._getExtraPath(),
        #            env=getEnviron())

    def createOutputStep(self):
        inputVol = self.inputVolume.get()
        vol = Volume()
        vol.setFileName(self._getFileName('out_vol3DFSC'))
        vol.setSamplingRate(inputVol.getSamplingRate())

        self._defineOutputs(outputVolume=vol)
        self._defineSourceRelation(self.inputVolume, vol)

    #--------------------------- INFO functions --------------------------------
    
    def _summary(self):
        summary = []
        if self.getOutputsSize() > 0:
            logFn = self.getLogPaths()[0]
            sph = findSphericity(logFn)
            summary.append('Sphericity: %0.3f ' % sph)
        else:
            summary.append("Output is not ready yet.")

        return summary
    
    def _validate(self):
        errors = []

        half1 = self.volumeHalf1.get()
        half2 = self.volumeHalf2.get()
        if half1.getSamplingRate() != half2.getSamplingRate():
            errors.append('The selected half volumes have not the same pixel size.')
        if half1.getXDim() != half2.getXDim():
            errors.append('The selected half volumes have not the same dimensions.')
                
        return errors
    
    #--------------------------- UTILS functions -------------------------------
 
    def _getArgs(self):
        """ Prepare the args dictionary."""

        args = {'--halfmap1': basename(self._getFileName('input_half1Fn')),
                '--halfmap2': basename(self._getFileName('input_half2Fn')),
                '--fullmap': basename(self._getFileName('input_volFn')),
                '--apix': self.inputVolume.get().getSamplingRate(),
                '--ThreeDFSC': '3D-FSC',
                '--dthetaInDegrees': self.dTheta.get(),
                '--FSCCutoff': self.fscCutoff.get(),
                '--ThresholdForSphericity': self.thrSph.get(),
                '--HighPassFilter': self.hpFilter.get(),
                '--numThresholdsForSphericityCalcs': self.numThr.get()
                }
        if self.applyMask and self.maskVolume:
            args.update({'--mask': basename(self._getFileName('input_maskFn'))})

        return args

    def _getProgram(self):
        """ Return the program binary that will be used. """
        if 'NYSBC_3DFSC_HOME' not in os.environ:
            return None
        cmd = join(os.environ['NYSBC_3DFSC_HOME'], 'ThreeDFSC', 'ThreeDFSC_Start.py')
        return str(cmd)
