# **************************************************************************
# *
# *  Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk)
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

import os
from pyworkflow.protocol.params import (PointerParam, FloatParam,
                                        LabelParam, EnumParam, StringParam,
                                        BooleanParam, IntParam, LEVEL_ADVANCED)
import pyworkflow.utils as pwutils
from pyworkflow.em.protocol import ProtAnalysis3D
from eman2 import getEmanProgram, validateVersion
from constants import *
from convert import writeSetOfParticles


class EmanProtTiltValidate(ProtAnalysis3D):
    """
    This protocol wraps the *e2tiltvalidate.py* Eman2 program,
    that performs tilt validation using
    the method described in Rosenthal and Henderson, JMB (2003).
    """

    _label = 'tilt validate'

    def __init__(self, **kwargs):
        ProtAnalysis3D.__init__(self, **kwargs)

    def _createFilenameTemplates(self):
        """ Centralize the names of the files. """
        myDict = {
            'untiltPartSet': 'sets/untilted_ptcls.lst',
            'tiltPartSet': 'sets/tilted_ptcls.lst',
            'outputAngles': self._getExtraPath('TiltValidate_01/perparticletilts.json'),
            'outputContourPlot': self._getExtraPath('TiltValidate_01/contour.hdf')
        }
        self._updateFilenamesDict(myDict)

    # --------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputVolume', PointerParam, pointerClass='Volume',
                      label="Input volume",
                      help='Select the input volume that will be validated.')
        form.addParam('inputTiltPair', PointerParam,
                      label="Input tilt pair particles",
                      pointerClass='ParticlesTiltPair',
                      help='Select the input set of tilt pair particles.')
        form.addParam('symmetry', StringParam, default='c1',
                      label='Symmetry group',
                      help='Set the symmetry; if no value is given then '
                           'the model is assumed to have no symmetry. \n'
                           'Choices are: *i, c, d, tet, icos, or oct* \n'
                           'See http://blake.bcm.edu/emanwiki/EMAN2/Symmetry\n'
                           'for a detailed description of symmetry in Eman.')
        form.addParam('maxtilt', FloatParam, default=180.0,
                      label='Max tilt angle',
                      help='Maximum tilt angle permitted when finding tilt '
                           'distances.')
        form.addParam('quaternion', BooleanParam, default=False,
                      label='Use quaternions', expertLevel=LEVEL_ADVANCED,
                      help='Use quaternions for tilt distance computation')
        form.addParam('delta', FloatParam, default=5.0,
                      label='Projection step (deg.)',
                      help='Angular step size for alignment')
        form.addParam('shrink', IntParam, default=1,
                      expertLevel=LEVEL_ADVANCED,
                      label='Shrink particles',
                      help='Optionally shrink the input particles by an integer '
                           'amount prior to computing similarity scores. '
                           'For speed purposes.')
        form.addParam('doContourPlot', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label='Do contour plot?',
                      help='Also make a contour plot similar to fig. 6 '
                           'in Henderson paper')
        form.addParam('tiltRange', IntParam, default=15,
                      expertLevel=LEVEL_ADVANCED,
                      condition='doContourPlot',
                      label='Tilt range',
                      help='The angular tilt range to search')
        form.addParam('verbose', IntParam, default=0,
                      expertLevel=LEVEL_ADVANCED,
                      label='Verbose level',
                      help='Verbose level from 0 to 9. ')

        form.addSection(label='Similarity matrix')
        form.addParam('paramsMsg', LabelParam, default=True,
                      label='These parameters are for advanced users only!\n',
                      help='For help please address to EMAN2 %s or run:\n'
                           '*e2help.py cmp -v 2* or\n'
                           '*e2help.py aligners -v 2*' % WIKI_URL)
        line = form.addLine('simcmp: ',
                            help='The name of a cmp to be used in comparing '
                                 'the aligned images (default=ccc)')
        line.addParam('simcmpType', EnumParam,
                      choices=['ccc', 'dot', 'frc', 'lod', 'optsub',
                               'optvariance', 'phase', 'quadmindot',
                               'sqeuclidean', 'vertical', 'None'],
                      label='type', default=CMP_CCC,
                      display=EnumParam.DISPLAY_COMBO)
        line.addParam('simcmpParams', StringParam,
                      default='', label='params')

        group = form.addGroup('First stage aligner')
        line = group.addLine('simalign: ')
        line.addParam('simalignType', EnumParam,
                      choices=['frm2d', 'rotate_flip',
                               'rotate_flip_iterative', 'rotate_precenter',
                               'rotate_trans_flip_scale',
                               'rotate_trans_flip_scale_iter',
                               'rotate_trans_scale_iter',
                               'rotate_translate', 'rotate_translate_flip',
                               'rotate_translate_flip_iterative',
                               'rotate_translate_flip_resample',
                               'rotate_translate_iterative',
                               'rotate_translate_resample',
                               'rotate_translate_scale',
                               'rotate_translate_tree',
                               'rotational', 'rotational_iterative',
                               'rtf_exhaustive',
                               'rtf_slow_exhaustive', 'scale', 'symalign',
                               'symalignquat', 'translational', 'None'],
                      label='type', default=ALN_ROTATE_TRANSLATE_TREE,
                      display=EnumParam.DISPLAY_COMBO)
        line.addParam('simalignParams', StringParam,
                      default='', label='params')
        line = group.addLine('simaligncmp: ')
        line.addParam('simaligncmpType', EnumParam,
                      choices=['ccc', 'dot', 'frc', 'lod', 'optsub',
                               'optvariance', 'phase', 'quadmindot',
                               'sqeuclidean', 'vertical', 'None'],
                      label='type', default=CMP_CCC,
                      display=EnumParam.DISPLAY_COMBO)
        line.addParam('simaligncmpParams', StringParam,
                      default='', label='params')

        group = form.addGroup('Second stage aligner')
        line = group.addLine('simralign: ')
        line.addParam('simralignType', EnumParam,
                      choices=['None', 'refine',
                               'refine_3d', 'refine_3d_grid'],
                      label='type', default=RALN_NONE,
                      display=EnumParam.DISPLAY_COMBO)
        line.addParam('simralignParams', StringParam,
                      default='', label='params')
        line = group.addLine('simraligncmp: ')
        line.addParam('simraligncmpType', EnumParam,
                      choices=['ccc', 'dot', 'frc', 'lod', 'optsub',
                               'optvariance', 'phase', 'quadmindot',
                               'sqeuclidean', 'vertical', 'None'],
                      label='type', default=CMP_DOT,
                      display=EnumParam.DISPLAY_COMBO)
        line.addParam('simraligncmpParams', StringParam,
                      default='', label='params')

        form.addParallelSection(threads=1, mpi=0)

    # --------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._createFilenameTemplates()
        self._insertFunctionStep('convertImagesStep')
        args = self._prepareParams()
        self._insertFunctionStep('runValidateStep', args)
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions ------------------------------
    def convertImagesStep(self):
        part = self.inputTiltPair.get()
        partUnt = part.getUntilted()
        partTilt = part.getTilted()
        storePath = self._getExtraPath("particles")
        pwutils.makePath(storePath)
        print "Converting input particle set.."

        for partSet, suffix in zip([partUnt, partTilt],
                                   ['_untilted_ptcls', '_tilted_ptcls']):
            partAlign = partSet.getAlignment()
            writeSetOfParticles(partSet, storePath,
                                alignType=partAlign, suffix=suffix)

            setName = suffix.split('_')[1]
            program = getEmanProgram('e2buildsets.py')
            args = " particles/*%s.hdf --setname=%s --minhisnr=-1" % (
                suffix, setName)
            self.runJob(program, args, cwd=self._getExtraPath(),
                        numberOfMpi=1, numberOfThreads=1)

    def runValidateStep(self, args):
        program = getEmanProgram('e2tiltvalidate.py')
        self.runJob(program, args, cwd=self._getExtraPath(), numberOfThreads=1)

    def createOutputStep(self):
        pass

    # --------------------------- INFO functions -------------------------------
    def _validate(self):
        errors = []
        validateVersion(self, errors)
        self._validateDim(self.inputTiltPair.get().getUntilted(),
                          self.inputVolume.get(), errors,
                          'Input tilt pair particles', 'Input volume')

        return errors

    def _summary(self):
        summary = []
        summary.append("Max. tilt angle: *%0.2f*" % self.maxtilt.get())
        summary.append("Projection step: *%d deg.*" % self.delta.get())
        summary.append("Symmetry: *%s*" % self.symmetry.get())

        return summary

    # --------------------------- UTILS functions ------------------------------

    def _prepareParams(self):
        args = " --untiltdata=%(untilt)s --tiltdata=%(tilt)s --volume=%(volume)s"
        args += " --maxtiltangle=%(maxtilt)f --sym=%(sym)s --delta=%(delta)f"
        args += " --verbose=%(verb)d"

        if self.shrink.get() != 1:
            args += " --shrink=%d" % self.shrink.get()
        if self.quaternion:
            args += " --quaternion"
        if self.doContourPlot:
            args += " --docontourplot --tiltrange %d" % self.tiltRange.get()
        if self.numberOfThreads.get() > 1:
            args += " --parallel=thread:%d" % self.numberOfThreads.get()

        params = {'untilt': self._getFileName("untiltPartSet"),
                  'tilt': self._getFileName("tiltPartSet"),
                  'volume': os.path.relpath(self.inputVolume.get().getFileName(),
                                            self._getExtraPath()).replace(":mrc", ""),
                  'maxtilt': self.maxtilt.get(),
                  'sym': self.symmetry.get(),
                  'delta': self.delta.get(),
                  'verb': self.verbose.get()
                  }

        args = args % params

        for param in ['simcmp', 'simalign', 'simaligncmp',
                      'simralign', 'simraligncmp']:
            args += self._getSimmxOpts(param)

        return args

    def _getParticlesStack(self):
        if not (self.inputParticles.get().isPhaseFlipped() and
                self.inputParticles.get().hasCTF()):
            return self._getFileName("partFlipSet")
        else:
            return self._getFileName("partSet")

    def _getSimmxOpts(self, option):
        optionType = "optionType = self.getEnumText('" + option + "Type')"
        optionParams = 'optionParams = self.' + option + 'Params.get()'
        exec (optionType)
        exec (optionParams)

        if optionType == 'None':
            return ''
        if optionParams != '':
            argStr = ' --%s=%s:%s' % (option, optionType, optionParams)
        else:
            argStr = ' --%s=%s' % (option, optionType)

        return argStr
