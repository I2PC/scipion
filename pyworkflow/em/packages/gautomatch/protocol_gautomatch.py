# **************************************************************************
# *
# * Authors:     Grigory Sharov (sharov@igbmc.fr)
# *
# * L'Institut de genetique et de biologie moleculaire et cellulaire (IGBMC)
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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************

import os
from os.path import join, exists, basename

import pyworkflow.utils as pwutils
import pyworkflow.protocol.params as params
import pyworkflow.em as em
from pyworkflow.utils.properties import Message
from pyworkflow.protocol.constants import LEVEL_ADVANCED


class ProtGautomatch(em.ProtParticlePicking):
    """
    Gautomatch is a GPU accelerated program for fully automatic particle picking.

    Note about CTF: It is perfectly fine to use raw micrograph(before CTF correction) for particle picking
    since 30~50A is sufficient to auto-pick the particles. Usually for the micrograph with defocus around 2-5um,
    the first CTF zero is around 20~30A, so it is not useful to do CTF correction in general.
    However, you can determine CTF and flip the phases before picking. Full CTF correction on micrographs or
    applying full CTF on templates is not suggested, because these operation is normally targeting for high resolution.
    Since particle picking is basically a low-resolution operation, higher resolution will only introduce more false picking
    and template-bias, known as the so-called Einstein noise.
    """
    _label = 'auto-picking'

    # --------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):

        em.ProtParticlePicking._defineParams(self, form)
        form.addParam('inputReferences', params.PointerParam,
                      pointerClass='SetOfAverages',
                      label='Input References', important=True,
                      help="Template images (2D class averages or reprojections "
                           "from a reference volume) to be used in picking.")
        form.addParam('invertTemplatesContrast', params.BooleanParam, default=False,
                      label='References have inverted contrast',
                      help='Set to Yes to indicate that the reference have inverted \n'
                           'contrast with respect to the particles in the micrographs.')
        form.addParam('angStep', params.IntParam, default=5,
                      label='Angular step size',
                      help='Angular step size for picking, in degrees')
        form.addParam('threshold', params.FloatParam, default=0.1,
                      label='Threshold',
                      help='Particles with CCC above the threshold will be picked')
        form.addParam('particleSize', params.IntParam, default=250,
                      label='Particle radius (A)',
                      help="Particle radius in Angstrom")

        form.addSection(label='Advanced')

        form.addParam('advanced', params.BooleanParam, default=True,
                      label='Guess advanced parameters?',
                      help="By default, the program will optimize advanced parameters by itself, "
                           "however if you want to modify them, select No")
        form.addParam('boxSize', params.IntParam, default=128,
                      label='Box size (pix)', condition='not advanced',
                      help="Box size, in pixels; a suggested value will be "
                           "automatically calculated using pixel size and particle size")
        form.addParam('maxDist', params.IntParam, default=300,
                      label='Max search distance (A)', condition='not advanced',
                      help='Maximum search distance in Angstrom\n Use value of 0.9~1.1X diameter; '
                           'can be 0.3~0.5X for filament-like particle\n'
                           'For each particle Pi we try to search a local area using this distance option: '
                           'once program finds a better candidate Pj, the particle Pi is rejected. '
                           'If all peaks <= max_dist have been checked and no better candidates were found, '
                           'then Pi is considered as a successfully picked particle.')
        form.addParam('speed', params.IntParam, default=2,
                      label='Speed', condition='not advanced',
                      help='Speed level {0,1,2,3,4}. The bigger the faster, but less accurate.\n'
                           'Suggested values: 2 for >1 MDa complex, 1 for <500 kD complex, 1 or 2 for 500~1000 kD.\n'
                           '0 is not suggested, because the accuracy is simply fitting noise, unless for special noise-free micrographs. '
                           'Use 3 for huge viruses, but 2 is still preferred. Probably do not use 4 at all, it is not accurate in general.')

        group = form.addGroup('Local sigma parameters')
        group.addParam('localSigmaCutoff', params.FloatParam, default=1.2,
                       label='Local sigma cut-off', condition='not advanced',
                       help='Local sigma cut-off (relative value), 1.2~1.5 should be a good range\n'
                            'Normally a value >1.2 will be ice, protein aggregation or contamination')
        group.addParam('localSigmaDiam', params.IntParam, default=100,
                       label='Local sigma diameter (A)', condition='not advanced',
                       help='Diameter for estimation of local sigma, in Angstrom')

        group = form.addGroup('Local average parameters')
        line = group.addLine('Local average range', condition='not advanced',
                             help="Local average cut-off (relative value), "
                                  "any pixel values outside the range will be considered as ice/aggregation/carbon etc.")
        line.addParam('localAvgMin', params.FloatParam, default=-1.0,
                      label='Min')
        line.addParam('localAvgMax', params.FloatParam, default=1.0,
                      label='Max')
        group.addParam('localAvgDiam', params.IntParam, default=100,
                       label='Local average diameter (A)', condition='not advanced',
                       help='Diameter for estimation of local average, in Angstrom. '
                            '1.5~2.0X particle diameter suggested\n'
                            'However, if you have sharp/small ice or any dark/bright dots, '
                            'using a smaller value will be much better to get rid of these areas')

        line = form.addLine('Micrograph band-pass filter range (A)', condition='not advanced',
                            help="Apply band-pass filter on the micrographs:\n"
                                 "low-pass filter to increase the contrast of raw micrographs, suggested range 20~50 A\n"
                                 "high-pass filter to get rid of the global background of raw micrographs, suggested range 200~2000 A")
        line.addParam('lowPass', params.IntParam, default=30,
                      label='Min')
        line.addParam('highPass', params.IntParam, default=1000,
                      label='Max')

        form.addParam('GPUId', params.IntParam, default=0,
                      label='GPU ID', condition='not advanced',
                      help='GPU ID, normally it is 0')


        # --------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('convertInputStep',
                                 self.getInputMicrographs().strId(),
                                 self.inputReferences.get().strId())
        self._insertPickingSteps()
        self._insertFunctionStep('createOutputStep')

    def _insertPickingSteps(self):
        """ Insert the steps to launch the picker with different modes. 
        Prepare the command arguments that will be passed. 
        """
        args = ' --T %s' % self._getTmpPath('references.mrcs')
        args += ' --apixM %0.2f' % self.inputMicrographs.get().getSamplingRate()
        args += ' --apixT %0.2f' % self.inputReferences.get().getSamplingRate()
        args += ' --ang_step %d' % self.angStep
        args += ' --diameter %d' % (2 * self.particleSize.get())
        args += ' --cc_cutoff %0.2f' % self.threshold

        if not self.advanced:
            args += ' --speed %d' % self.speed
            args += ' --boxsize %d' % self.boxSize
            args += ' --max_dist %d' % self.maxDist
            args += ' --gid %d' % self.GPUId
            args += ' --lsigma_cutoff %0.2f' % self.localSigmaCutoff
            args += ' --lsigma_D %d' % self.localSigmaDiam
            args += ' --lave_max %0.2f' % self.localAvgMax
            args += ' --lave_min %0.2f' % self.localAvgMin
            args += ' --lave_D %d' % self.localAvgDiam
            args += ' --lp %d' % self.lowPass
            args += ' --hp %d' % self.highPass

        if not self.invertTemplatesContrast:
            args += ' --dont_invertT'

        args += ' %s' % join(self._getExtraPath('micrographs'), '*.mrc')

        self._insertFunctionStep('runGautomatchStep', args)

    # --------------------------- STEPS functions ---------------------------------------------------
    def convertInputStep(self, micsId, refsId):
        """ This step will take of the convertions from the inputs.
        Micrographs: they will be linked if are in '.mrc' format, converted otherwise.
        References: will always be converted to '.mrcs' format
        """
        ih = em.ImageHandler()
        micDir = self._getExtraPath('micrographs')  # put output and mics in extra dir
        pwutils.makePath(micDir)

        for mic in self.getInputMicrographs():
            # Create micrograph folder
            micName = mic.getFileName()
            # If micrographs are in .mrc format just link it
            # otherwise convert them
            outMic = join(micDir, pwutils.replaceBaseExt(micName, 'mrc'))

            if micName.endswith('.mrc'):
                pwutils.createLink(micName, outMic)
            else:
                ih.convert(mic, outMic)

        # We will always convert the templates to mrcs stack
        inputRefs = self.inputReferences.get()
        inputRefs.writeStack(self._getTmpPath('references.mrcs'))

    def runGautomatchStep(self, args):
        self.runJob(self._getProgram(), args)

    def createOutputStep(self):
        micSet = self.getInputMicrographs()
        ih = em.ImageHandler()
        coordSet = self._createSetOfCoordinates(micSet)
        if self.boxSize and self.boxSize > 0:
            coordSet.setBoxSize(self.boxSize.get())
        else:
            coordSet.setBoxSize(self.inputReferences.get().getDim()[0])

        for mic in micSet:
            fnMic = pwutils.removeExt(mic.getFileName())
            fnCoords = basename(fnMic) + '_automatch.box'
            fn2parse = self._getExtraPath('micrographs', fnCoords)
            print fn2parse
            with open(fn2parse, "r") as source:
                for line in source:
                    tokens = line.split()
                    coord = em.Coordinate()
                    coord.setPosition(int(tokens[0]) + int(tokens[2]) / 2, int(tokens[1]) + int(tokens[3]) / 2)
                    coord.setMicrograph(mic)
                    coordSet.append(coord)
                    # FIXME this should be outside the loop
                    if int(tokens[2]) != coordSet.setBoxSize:
                        coordSet.setBoxSize(int(tokens[2]))

        self._defineOutputs(outputCoordinates=coordSet)
        self._defineSourceRelation(micSet, coordSet)

    # --------------------------- INFO functions --------------------------------------------
    def _validate(self):
        errors = []
        # Check that the program exists
        if not exists(self._getProgram()):
            errors.append("Binary '%s' does not exist. \n"
                          "Check configuration file: ~/.config/scipion/scipion.conf\n"
                          "and set GAUTOMATCH variables properly." % self._getProgram())
            print "os.environ['GAUTOMATCH_HOME']", os.environ['GAUTOMATCH_HOME']
            print "os.environ['GAUTOMATCH']", os.environ['GAUTOMATCH']

        if not self.localAvgMin < self.localAvgMax:
            errors.append('Wrong values of local average cut-off!')
        return errors

    def _summary(self):
        summary = []
        summary.append("Number of input micrographs: %d" % self.getInputMicrographs().getSize())
        if (self.getOutputsSize() > 0):
            summary.append("Number of particles picked: %d" % self.getCoords().getSize())
            summary.append("Particle size: %d px" % self.getCoords().getBoxSize())
            summary.append("Threshold min: %0.2f" % self.threshold)
        else:
            summary.append(Message.TEXT_NO_OUTPUT_CO)
        return summary

    def _methods(self):
        methodsMsgs = []
        if self.getInputMicrographs() is None:
            return ['Input micrographs not available yet.']
        methodsMsgs.append("Input micrographs %s." % (self.getObjectTag(self.getInputMicrographs())))

        if self.getOutputsSize() > 0:
            output = self.getCoords()
            methodsMsgs.append("%s: User picked %d particles with a particle size of %d px and minimal threshold %0.2f."
                               % (self.getObjectTag(output), output.getSize(), output.getBoxSize(),
                                  self.threshold.get()))
        else:
            methodsMsgs.append(Message.TEXT_NO_OUTPUT_CO)

        return methodsMsgs

    def _citations(self):
        return ['Zhang2016']

    # --------------------------- UTILS functions --------------------------------------------------
    def _getProgram(self):
        """ Return the program binary that will be used. """
        binary = os.environ['GAUTOMATCH']
        program = join(os.environ['GAUTOMATCH_HOME'], 'bin', basename(binary))
        return program
