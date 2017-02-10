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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import os

import pyworkflow.utils as pwutils
import pyworkflow.protocol.params as params
import pyworkflow.em as em
from pyworkflow import VERSION_1_1
from pyworkflow.utils.properties import Message

from convert import (readSetOfCoordinates, writeSetOfCoordinates,
                     runGautomatch, getProgram)



class ProtGautomatch(em.ProtParticlePicking):
    """
    Gautomatch is a GPU accelerated program for accurate, fast, flexible and fully
    automatic particle picking from cryo-EM micrographs with or without templates.
    """
    _label = 'auto-picking'
    _version = VERSION_1_1

    # --------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):

        em.ProtParticlePicking._defineParams(self, form)
        form.addParam('inputReferences', params.PointerParam,
                      pointerClass='SetOfAverages',
                      label='Input References', important=True, allowsNull=True,
                      help="Template images (2D class averages or reprojections "
                           "from a reference volume) to be used in picking.\n"
                           "If not provided, references will be auto-generated. This is fine "
                           "for *spherical particles* like virus or ribosome.")
        form.addParam('invertTemplatesContrast', params.BooleanParam, default=False,
                      label='References have inverted contrast',
                      help='Set to Yes to indicate that the reference have inverted '
                           'contrast with respect to the particles in the micrographs.\n'
                           'Keep in mind that auto-generated templates will be WHITE.')
        form.addParam('angStep', params.IntParam, default=5,
                      label='Angular step size',
                      help='Angular step size for picking, in degrees')
        form.addParam('threshold', params.FloatParam, default=0.1,
                      label='Threshold',
                      help='Particles with CCC above the threshold will be picked')
        form.addParam('particleSize', params.IntParam, default=250,
                      label='Particle radius (A)',
                      help="Particle radius in Angstrom")
        form.addParam('GPUId', params.IntParam, default=0,
                      label='GPU ID',
                      help='GPU ID, normally it is 0')

        form.addSection(label='Advanced')
        form.addParam('advanced', params.BooleanParam, default=True,
                      label='Guess advanced parameters?',
                      help="By default, the program will optimize advanced "
                           "parameters by itself, however if you want to "
                           "modify them, select No")
        form.addParam('boxSize', params.IntParam, default=128,
                      label='Box size (pix)', condition='not advanced',
                      help="Box size, in pixels; a suggested value will be "
                           "automatically calculated using pixel size and "
                           "particle size")
        form.addParam('minDist', params.IntParam, default=300,
                      label='Min inter-particle distance (A)',
                      condition='not advanced',
                      help='Minimum distance between particles in Angstrom\n '
                           'Use value of 0.9~1.1X diameter; '
                           'can be 0.3~0.5X for filament-like particle')
        form.addParam('speed', params.IntParam, default=2,
                      label='Speed', condition='not advanced',
                      help='Speed level {0,1,2,3,4}. The bigger the faster, but '
                           'less accurate.\n'
                           'Suggested values: 2 for >1 MDa complex, 1 for <500 '
                           'kD complex, 1 or 2 for 500~1000 kD.\n'
                           '0 is not suggested, because the accuracy is simply '
                           'fitting noise, unless for special noise-free micrographs. '
                           'Use 3 for huge viruses, but 2 is still preferred. '
                           'Probably do not use 4 at all, it is not accurate '
                           'in general.')

        form.addSection(label='Sigma & avg')
        form.addParam('advLabel', params.LabelParam,
                      important=True,
                      label='To adjust these parameters, select "No" for the '
                            '"Guess advanced parameters?" on the Advanced tab.')
        group = form.addGroup('Local sigma parameters',
                              condition='not advanced')
        group.addParam('localSigmaCutoff', params.FloatParam, default=1.2,
                       label='Local sigma cut-off',
                       help='Local sigma cut-off (relative value), 1.2~1.5 '
                            'should be a good range\n'
                            'Normally a value >1.2 will be ice, protein '
                            'aggregation or contamination.\n'
                            'This option is designed to get rid of sharp carbon/ice '
                            'edges or sharp metal particles.')
        group.addParam('localSigmaDiam', params.IntParam, default=100,
                       label='Local sigma diameter (A)',
                       help='Diameter for estimation of local sigma, in Angstrom.\n'
                            'Usually this diameter could be 0.5-2x of your particle '
                            'diameter according to several factors. When using bigger values, '
                            'normally you should decrease *Local sigma cut-off*. '
                            'For smaller and sharper high density contamination/ice/metal particles '
                            'you could use a smaller diameter and larger *Local sigma cut-off*')

        group = form.addGroup('Local average parameters',
                              condition='not advanced')
        line = group.addLine('Local average range',
                             help="Local average cut-off (relative value), "
                                  "any pixel values outside the range will be "
                                  "considered as ice/aggregation/carbon etc.\n"
                                  "Min parameter is used to reject the central parts of ice, "
                                  "carbon etc. which normally have lower density than the particles.\n"
                                  "Max parameter is usually not useful for cryo-EM micrograph with black "
                                  "particles, but might be helpfull to get rid of 'hot' area. "
                                  "For negative stain micrograph, if it rejects most of the true particles, "
                                  "just set Max to very big value like 10.0.")
        line.addParam('localAvgMin', params.FloatParam, default=-1.0,
                      label='Min')
        line.addParam('localAvgMax', params.FloatParam, default=1.0,
                      label='Max')
        group.addParam('localAvgDiam', params.IntParam, default=100,
                       label='Local average diameter (A)',
                       help='Diameter for estimation of local average, in Angstrom. '
                            '1.5~2.0X particle diameter suggested\n'
                            'However, if you have sharp/small ice or any '
                            'dark/bright dots, using a smaller value will be '
                            'much better to get rid of these areas')

        form.addSection(label='Filter')
        line = form.addLine('Micrograph band-pass filter range (A)',
                            help="Apply band-pass filter on the micrographs:\n"
                                 "low-pass filter to increase the contrast of "
                                 "raw micrographs, suggested range 20~50 A\n"
                                 "high-pass filter to get rid of the global "
                                 "background of raw micrographs, suggested "
                                 "range 200~2000 A. This filter is applied after ice/carbon/"
                                 "contamination detection, but before true particle detection")
        line.addParam('lowPass', params.IntParam, default=30,
                      label='Min')
        line.addParam('highPass', params.IntParam, default=1000,
                      label='Max')

        form.addParam('preFilt', params.BooleanParam, default=False,
                      label='Pre-filter micrographs?',
                      help="This band-pass pre-filter is normally not suggested, "
                      "because it can affect ice/carbon detection. "
                      "Use it only if you have a severe ice gradient.")
        line = form.addLine('Pre-filter range (A)', condition='preFilt')
        line.addParam('prelowPass', params.IntParam, default=8,
                      label='Min')
        line.addParam('prehighPass', params.IntParam, default=1000,
                      label='Max')

        form.addSection(label='Exclusive picking')
        form.addParam('exclusive', params.BooleanParam, default=False,
                      label='Exclusive picking?',
                      help='Exclude user-provided areas. This can be useful in the '
                           'following cases:\n\n(a) Another cycle of auto-picking '
                           'after 2D classification: in this case, usually you are '
                           'pretty sure that some of the particles are completely rubbish, '
                           'it will be much better to exclude them during picking.\n'
                           '(b) Picking for partial structure: sometimes, you might have '
                           'two/multiple domain complex, one is severely dominant and affect '
                           'picking of the other (the rest). If you want to focus on another '
                           'domain, it might be quite helpful to exclude such good particles '
                           'from 2D classification.\n(c) Strong orientation preference: if '
                           'your templates were severely biased and mainly picked the '
                           'preferred views, then it might be nice to exclude the preferred '
                           'views and focused on rare views.')
        form.addParam('inputBadCoords', params.PointerParam, allowsNull=True,
                      pointerClass='SetOfCoordinates', condition='exclusive',
                      label='Coordinates to be excluded',
                      help='Coordinates can be imported beforehand or generated from '
                           'particles using scipion - extract coordinates protocol.')
        form.addParam('inputDefects', params.PointerParam, allowsNull=True,
                      pointerClass='SetOfCoordinates', condition='exclusive',
                      label='Detector defects coordinates',
                      help='Occasionally you might have detector defects, e.g. a '
                           'black/white stripe. This will help to get rid of these bad areas. '
                           'The boxes in this case should be overlapping and organized in a stripe/line '
                           'like this: http://www.mrc-lmb.cam.ac.uk/kzhang/Gautomatch/Gautomatch_v0.53/examples/exclusive_picking/global_excluded.png')

        form.addSection(label='Debug')
        form.addParam('writeCC', params.BooleanParam, default=False,
                      label='Write CC files?',
                      help='Specify to write out cross-correlation files in MRC stack')
        form.addParam('writeFilt', params.BooleanParam, default=False,
                      condition='preFilt',
                      label='Write pre-filtered micrographs?',
                      help='Specify to write out pre-filted micrographs')
        form.addParam('writeBg', params.BooleanParam, default=False,
                      label='Write estimated background?',
                      help='Specify to write out estimated background of the micrographs')
        form.addParam('writeBgSub', params.BooleanParam, default=False,
                      label='Write background-subtracted micrographs?',
                      help='Specify to write out background-subtracted micrographs')
        form.addParam('writeSigma', params.BooleanParam, default=False,
                      label='Write local sigma?',
                      help='Specify to write out local sigma micrographs')
        form.addParam('writeMsk', params.BooleanParam, default=False,
                      label='Write detected mask?',
                      help='Specify to write out the auto-detected mask (ice, '
                           'contamination, aggregation, carbon edges etc.)')

        # --------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):
        convId = self._insertFunctionStep('convertInputStep')
        deps = []
        refStack = None

        if self.inputReferences.get():
            refStack = self._getExtraPath('references.mrcs')

        # Insert one picking step per Micrograph
        for mic in self.inputMicrographs.get():
            micName = mic.getFileName()
            pickId = self._insertFunctionStep('runGautomatchStep',
                                              micName, refStack,
                                              self.getArgs(),
                                              prerequisites=[convId])
            deps.append(pickId)
        self._insertFunctionStep('createOutputStep', prerequisites=deps)

    # --------------------------- STEPS functions ---------------------------------------------------

    def convertInputStep(self):
        """ This step will take of the conversions from the inputs.
        Micrographs: they will be linked if are in '.mrc' format, converted otherwise.
        References: will always be converted to '.mrcs' format
        """
        micDir = self.getMicrographsDir()  # put output and mics in extra dir
        pwutils.makePath(micDir)
        # We will always convert the templates to mrcs stack
        self.convertReferences(self._getExtraPath('references.mrcs'))
        # Convert input coords for exclusive picking
        self.convertCoordinates(micDir)

    def runGautomatchStep(self, micName, refStack, args):
        # We convert the input micrograph on demand if not in .mrc
        runGautomatch(micName, refStack, self.getMicrographsDir(), args, env=self._getEnviron())

    def createOutputStep(self):
        micSet = self.getInputMicrographs()
        ih = em.ImageHandler()
        coordSet = self._createSetOfCoordinates(micSet)
        if self.boxSize and self.boxSize > 0:
            coordSet.setBoxSize(self.boxSize.get())
        else:
            coordSet.setBoxSize(self.inputReferences.get().getXDim or 100)

        readSetOfCoordinates(self.getMicrographsDir(), micSet, coordSet)
        coordSetAux = self._createSetOfCoordinates(micSet, suffix='_rejected')
        coordSetAux.setBoxSize(coordSet.getBoxSize())
        readSetOfCoordinates(self.getMicrographsDir(), micSet, coordSetAux, suffix='_rejected.star')
        coordSetAux.write()

        # debug output
        if self.writeCC:
            self.createDebugOutput(suffix='_ccmax')
        if self.writeFilt:
            self.createDebugOutput(suffix='_pref')
        if self.writeBg:
            self.createDebugOutput(suffix='_bg')
        if self.writeBgSub:
            self.createDebugOutput(suffix='_bgfree')
        if self.writeSigma:
            self.createDebugOutput(suffix='_lsigma')
        if self.writeMsk:
            self.createDebugOutput(suffix='_mask')

        self._defineOutputs(outputCoordinates=coordSet)
        self._defineSourceRelation(micSet, coordSet)

    def createDebugOutput(self, suffix):
        micSet = self.getInputMicrographs()
        pixSize = micSet.getSamplingRate()
        outputDebugMics = self._createSetOfMicrographs(suffix=suffix)
        # debug output images are downsampled by a factor of 4
        outputDebugMics.setSamplingRate(float(pixSize * 4))
        for mic in micSet:
            micFn = self.getOutputName(mic.getFileName(), suffix)
            mic.setFileName(micFn)
            outputDebugMics.append(mic)
        outputDebugMics.write()

    # --------------------------- INFO functions --------------------------------------------
    def _validate(self):
        errors = []
        # Check that the program exists
        program = getProgram()
        if program is None:
            errors.append("Missing variables GAUTOMATCH and/or GAUTOMATCH_HOME")
        elif not os.path.exists(program):
            errors.append("Binary '%s' does not exists.\n" % program)

        # If there is any error at this point it is related to config variables
        if errors:
            errors.append("Check configuration file: ~/.config/scipion/scipion.conf")
            errors.append("and set GAUTOMATCH and GAUTOMATCH_HOME variables properly.")
            if program is not None:
                errors.append("Current values:")
                errors.append("GAUTOMATCH_HOME = %s" % os.environ['GAUTOMATCH_HOME'])
                errors.append("GAUTOMATCH = %s" % os.environ['GAUTOMATCH'])

        if not self.localAvgMin < self.localAvgMax:
            errors.append('Wrong values of local average cut-off!')
        if self.exclusive:
            if not self.inputBadCoords.get() and not self.inputDefects.get():
                errors.append("You have to provide at least one set of coordinates ")
                errors.append("for exclusive picking!")

        return errors

    def _summary(self):
        summary = []
        summary.append("Number of input micrographs: %d"
                       % self.getInputMicrographs().getSize())
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
        methodsMsgs.append("Input micrographs %s."
                           % (self.getObjectTag(self.getInputMicrographs())))

        if self.getOutputsSize() > 0:
            output = self.getCoords()
            methodsMsgs.append("%s: User picked %d particles with a particle "
                               "size of %d px and minimal threshold %0.2f."
                               % (self.getObjectTag(output), output.getSize(),
                                  output.getBoxSize(), self.threshold))
        else:
            methodsMsgs.append(Message.TEXT_NO_OUTPUT_CO)

        return methodsMsgs

    def _citations(self):
        return ['Zhang2016']

    # --------------------------- UTILS functions --------------------------------------------------
    def getMicrographsDir(self):
        return self._getExtraPath('micrographs')

    def getArgs(self, threshold=True, mindist=True):
        """ Return the Gautomatch parameters for picking one micrograph.
         The command line will depends on the protocol selected parameters.
         Params:
            micFn: micrograph filename
            refStack: filename with the references stack (.mrcs)
        """
        args = ' --apixM %0.2f' % self.inputMicrographs.get().getSamplingRate()
        args += ' --ang_step %d' % self.angStep
        args += ' --diameter %d' % (2 * self.particleSize.get())
        args += ' --lp %d' % self.lowPass
        args += ' --hp %d' % self.highPass
        args += ' --gid %d' % self.GPUId

        if self.inputReferences.get():
            args += ' --apixT %0.2f' % self.inputReferences.get().getSamplingRate()

        if not self.invertTemplatesContrast:
            args += ' --dont_invertT'

        if threshold:
            args += ' --cc_cutoff %0.2f' % self.threshold

        if not self.advanced:
            args += ' --speed %d' % self.speed
            args += ' --boxsize %d' % self.boxSize
            if mindist:
                args += ' --min_dist %d' % self.minDist
            args += ' --lsigma_cutoff %0.2f' % self.localSigmaCutoff
            args += ' --lsigma_D %d' % self.localSigmaDiam
            args += ' --lave_max %0.2f' % self.localAvgMax
            args += ' --lave_min %0.2f' % self.localAvgMin
            args += ' --lave_D %d' % self.localAvgDiam

        if self.preFilt:
            args += ' --do_pre_filter --pre_lp %d' % self.prelowPass
            args += ' --pre_hp %d' % self.prehighPass

        if self.exclusive:
            if self.inputBadCoords.get():
                args += ' --exclusive_picking --excluded_suffix _rubbish.star'
            if self.inputDefects.get():
                args += ' --global_excluded_box %s' % self._getExtraPath('micrographs/defects.star')

        if self.writeCC:
            args += ' --write_ccmax_mic'
        if self.writeFilt:
            args += ' --write_pref_mic'
        if self.writeBg:
            args += ' --write_bg_mic'
        if self.writeBgSub:
            args += ' --write_bgfree_mic'
        if self.writeSigma:
            args += ' --write_lsigma_mic'
        if self.writeMsk:
            args += ' --write_mic_mask'

        return args

    def convertReferences(self, refStack):
        """ Write input references as an .mrcs stack. """
        inputRefs = self.inputReferences.get()
        if inputRefs:
            self.inputReferences.get().writeStack(refStack)

    def convertCoordinates(self, workingDir):
        if self.exclusive:
            if self.inputBadCoords.get():
                writeSetOfCoordinates(workingDir, self.inputBadCoords.get(), isGlobal=False)
            if self.inputDefects.get():
                writeSetOfCoordinates(workingDir, self.inputDefects.get(), isGlobal=True)

    def getOutputName(self, fn, key):
        """ Give a key, append the mrc extension
        and prefix the protocol working dir.
        """
        template = pwutils.removeBaseExt(fn) + key + '.mrc'

        return pwutils.join(self.getMicrographsDir(), template)
