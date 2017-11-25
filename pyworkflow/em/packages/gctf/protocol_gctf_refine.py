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

import os
from os.path import basename
import pyworkflow.utils as pwutils
import pyworkflow.em as em
import pyworkflow.em.metadata as md
import pyworkflow.protocol.params as params
from pyworkflow.protocol.constants import STEPS_PARALLEL
from pyworkflow.em.protocol import EMProtocol
from pyworkflow.em.data import Coordinate, Particle
from pyworkflow.object import Boolean
from pyworkflow import VERSION_1_2
from convert import getVersion, writeSetOfCoordinates, rowToCtfModel



# Phase shift target type
CCC = 0
MAXRES = 1

# Weighting type for local CTF
WEIGHT_EQUAL = 0
WEIGHT_DIST = 1
WEIGHT_BOTH = 2


class ProtGctfRefine(em.ProtParticles):
    """
    Refines local CTF of a set of particles
    using GPU-accelerated Gctf program.

    To find more information about Gctf go to:
    http://www.mrc-lmb.cam.ac.uk/kzhang
    """
    _label = 'CTF refinement on GPU'
    _lastUpdateVersion = VERSION_1_2

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)
        self.allowMpi = False
        self.allowThreads = False
        self._params = {}
        #self.stepsExecutionMode = STEPS_PARALLEL
        #self.isFirstTime = Boolean(False)

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('recalculate', params.BooleanParam, default=False,
                      condition='recalculate',
                      label="Do recalculate ctf?")
        form.addParam('continueRun', params.PointerParam, allowsNull=True,
                      condition='recalculate', label="Input previous run",
                      pointerClass=self.getClassName())
        form.addHidden('sqliteFile', params.FileParam, condition='recalculate',
                       allowsNull=True)
        form.addParam('inputParticles', params.PointerParam,
                      important=True,
                      label='Input particles',
                      pointerClass='SetOfParticles',
                      help='Provide a set of particles for local CTF refinement.')
        form.addParam('inputMicrographs', params.PointerParam, important=True,
                      condition='not recalculate',
                      label='Input micrographs',
                      pointerClass='SetOfMicrographs',
                      help='Select the SetOfMicrographs related to input particles.')
        form.addParam('ctfDownFactor', params.FloatParam, default=1.,
                      label='CTF Downsampling factor',
                      condition='not recalculate',
                      help='Set to 1 for no downsampling. Non-integer '
                           'downsample factors are possible. This downsampling '
                           'is only used for estimating the CTF and it does not '
                           'affect any further calculation. Ideally the estimation '
                           'of the CTF is optimal when the Thon rings are not too '
                           'concentrated at the origin (too small to be seen) and '
                           'not occupying the whole power spectrum (since this '
                           'downsampling might entail aliasing).')
        line = form.addLine('Resolution', condition='not recalculate',
                            help='Give a value in digital frequency (i.e. between '
                                 '0.0 and 0.5). These cut-offs prevent the typical '
                                 'peak at the center of the PSD and high-resolution '
                                 'terms where only noise exists, to interfere with '
                                 'CTF estimation. The default lowest value is 0.05 '
                                 'but for micrographs with a very fine sampling this '
                                 'may be lowered towards 0. The default highest '
                                 'value is 0.35, but it should be increased for '
                                 'micrographs with signals extending beyond this '
                                 'value. However, if your micrographs extend further '
                                 'than 0.35, you should consider sampling them at a '
                                 'finer rate.')
        line.addParam('lowRes', params.FloatParam, default=0.05,
                      label='Lowest')
        line.addParam('highRes', params.FloatParam, default=0.35,
                      label='Highest')

        line = form.addLine('Defocus search range (microns)',
                            expertLevel=params.LEVEL_ADVANCED,
                            condition='not recalculate',
                            help='Select _minimum_ and _maximum_ values for '
                                 'defocus search range (in microns). '
                                 'Underfocus is represented by a positive '
                                 'number.')
        line.addParam('minDefocus', params.FloatParam, default=0.25,
                      label='Min')
        line.addParam('maxDefocus', params.FloatParam, default=4.,
                      label='Max')
        form.addParam('astigmatism', params.FloatParam, default=100.0,
                      label='Expected (tolerated) astigmatism',
                      help='Estimated astigmatism in Angstroms',
                      expertLevel=params.LEVEL_ADVANCED)
        form.addParam('windowSize', params.IntParam, default=512,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Window size', condition='not recalculate',
                      help='The PSD is estimated from small patches of this '
                           'size. Bigger patches allow identifying more '
                           'details. However, since there are fewer windows, '
                           'estimations are noisier.')
        form.addParam('plotResRing', params.BooleanParam, default=True,
                      label='Plot a resolution ring on a PSD file',
                      help='Whether to plot an estimated resolution ring '
                           'on the power spectrum',
                      expertLevel=params.LEVEL_ADVANCED)
        form.addParam('GPUCore', params.IntParam, default=0,
                      expertLevel=params.LEVEL_ADVANCED,
                      label="Choose GPU core",
                      help='GPU may have several cores. Set it to zero if '
                           'you do not know what we are talking about. '
                           'First core index is 0, second 1 and so on.')

        form.addSection(label='Advanced')
        form.addParam('doEPA', params.BooleanParam, default=False,
                      label="Do EPA",
                      help='Do Equiphase average used for output CTF file. '
                           'Only for nice output, will NOT be used for CTF '
                           'determination.')
        form.addParam('EPAsmp', params.IntParam, default=4,
                      expertLevel=params.LEVEL_ADVANCED,
                      condition='not _oldVersion',
                      label="Over-sampling factor for EPA")
        form.addParam('doBasicRotave', params.BooleanParam, default=False,
                      expertLevel=params.LEVEL_ADVANCED,
                      condition='_oldVersion',
                      label="Do rotational average",
                      help='Do rotational average used for output CTF file. '
                           'Only for nice output, will NOT be used for CTF '
                           'determination.')
        form.addParam('bfactor', params.IntParam, default=150,
                      expertLevel=params.LEVEL_ADVANCED,
                      label="B-factor",
                      help='B-factors used to decrease high resolution '
                           'amplitude, A^2; suggested range 50~300 except '
                           'using REBS method')
        form.addParam('overlap', params.FloatParam, default=0.5,
                      expertLevel=params.LEVEL_ADVANCED,
                      label="Overlap factor",
                      help='Overlapping factor for grid boxes sampling, '
                      'for windowsize=512, 0.5 means 256 pixels overlapping.')
        form.addParam('convsize', params.IntParam, default=85,
                      expertLevel=params.LEVEL_ADVANCED,
                      label="Boxsize for smoothing",
                      help='Boxsize to be used for smoothing, '
                           'suggested 1/5 ~ 1/20 of window size in pixel, '
                           'e.g. 99 for 512 window')

        group = form.addGroup('High-res refinement')
        group.addParam('doHighRes', params.BooleanParam, default=False,
                       label="Do high-resolution refinement",
                       help='Whether to do High-resolution refinement or not, '
                            'very useful for selecting high quality micrographs. '
                            'Especially useful when your data has strong '
                            'low-resolution bias')
        group.addParam('HighResL', params.FloatParam, default=15.0,
                       condition='doHighRes',
                       label="Lowest resolution",
                       help='Lowest resolution  to be used for High-resolution '
                            'refinement, in Angstroms')
        group.addParam('HighResH', params.FloatParam, default=4.0,
                       condition='doHighRes',
                       label="Highest resolution",
                       help='Highest resolution  to be used for High-resolution '
                            'refinement, in Angstroms')
        group.addParam('HighResBf', params.IntParam, default=50,
                       condition='doHighRes',
                       label="B-factor",
                       help='B-factor to be used for High-resolution '
                            'refinement, in Angstroms')

        form.addParam('doValidate', params.BooleanParam, default=False,
                      expertLevel=params.LEVEL_ADVANCED,
                      label="Do validation",
                      help='Whether to validate the CTF determination.')

        form.addSection(label='Phase shift')
        form.addParam('doPhShEst', params.BooleanParam, default=False,
                      label="Estimate phase shift?",
                      help='For micrographs collected with phase-plate. '
                           'It is suggested to import such micrographs with '
                           'amplitude contrast = 0. Also, using smaller '
                           '_lowest resolution_ (e.g. 15A) and smaller '
                           '_boxsize for smoothing_ (e.g. 50 for 1024 '
                           'window size) might be better.')
        line = form.addLine('Phase shift range range (deg)',
                            condition='doPhShEst',
                            help='Select _lowest_ and _highest_ phase shift '
                                 '(in degrees).')
        line.addParam('phaseShiftL', params.FloatParam, default=0.0,
                      condition='doPhShEst',
                      label="Min")
        line.addParam('phaseShiftH', params.FloatParam, default=180.0,
                      condition='doPhShEst',
                      label="Max")
        form.addParam('phaseShiftS', params.FloatParam, default=10.0,
                       condition='doPhShEst',
                       label="Step",
                       help='Phase shift search step. Do not worry about '
                            'the accuracy; this is just the search step, '
                            'Gctf will refine the phase shift anyway.')
        form.addParam('phaseShiftT', params.EnumParam, default=CCC,
                      condition='doPhShEst',
                      label='Target',
                      choices=['CCC', 'Resolution limit'],
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Phase shift target in the search: CCC or '
                           'resolution limit')

        form.addSection(label='Particle refine')
        line = form.addLine('Local resolution (A)',
                            help='Select _lowest_ and _highest_ resolution '
                                 'to be used for local CTF (in Angstrom).')
        line.addParam('locResL', params.IntParam, default=15, label='Low')
        line.addParam('locResH', params.IntParam, default=5, label='High')

        form.addParam('locRad', params.IntParam, default=1024,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Radius for local refinement (px)',
                      help='Radius for local refinement, no weighting '
                           'if the distance is larger than that')
        form.addParam('locAveType', params.EnumParam, default=WEIGHT_BOTH,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Local average type',
                      choices=['Equal weights', 'Distance', 'Distance and freq'],
                      display=params.EnumParam.DISPLAY_COMBO,
                      help='_Equal weights_: equal weights for all local '
                           'areas, neither distance nor frequency is '
                           'weighted\n_Distance_: single weight for each '
                           'local area, only distance is weighted\n'
                           '_Distance and freq_: Guassian weighting for '
                           'both distance and frequency')
        form.addParam('locBoxSize', params.IntParam, default=512,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Boxsize (px)',
                      help='Boxsize for local refinement (in pixels)')
        form.addParam('locOverlap', params.FloatParam, default=0.5,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Overlap',
                      help='Overlapping factor for grid boxes sampling')
        form.addParam('locAstm', params.BooleanParam, default=False,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Refine astigmatism?',
                      help='By default (False) only refine Z-height '
                           'changes in local area (suggested). If True, '
                           'refine local astigmatism (not suggested unless '
                           'SNR is very good).')
        form.addParallelSection(threads=0, mpi=0)

    #--------------------------- STEPS functions -------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('refineCtfStep')
        self._insertFunctionStep('createCtfModelStep')

    def convertInputStep(self):
        inputParticles = self.inputParticles.get()
        firstCoord = inputParticles.getFirstItem().getCoordinate()
        hasMicName = firstCoord.getMicName() is not None
        inputMics = self._getMicrographs()

        coords = self._createSetOfCoordinates(inputMics)
        newCoord = Coordinate()
        scale = inputParticles.getSamplingRate() / inputMics.getSamplingRate()
        print "Scaling coordinates by a factor *%0.2f*" % scale

        # Create the micrographs dict using either micName or micId
        micDict = {}
        micKey2 = None
        for mic in inputMics:
            if hasMicName:
                micKey = mic.getMicName()
                micKey2 = pwutils.removeBaseExt(micKey)
            else:
                micKey = mic.getObjId()
            if micKey in micDict:  # FIXME: maybe also check micKey2?
                print ">>> ERROR: micrograph key %s is duplicated!" % micKey
                print "           Used in micrographs:"
                print "           - %s" % micDict[micKey].getLocation()
                print "           - %s" % mic.getLocation()
                raise Exception("Micrograph key %s is duplicated!" % micKey)
            micDict[micKey] = mic.clone()
            if hasMicName:
                micDict[micKey2] = mic.clone()

        # match the mic from coord with micDict
        for particle in inputParticles:
            coord = particle.getCoordinate() or None
            if coord is None:
                print "Skipping particle, coordinates not found"
                continue

            if hasMicName:
                micKey = coord.getMicName()
                micKey2 = pwutils.removeBaseExt(micKey)
            else:
                micKey = coord.getMicId()
            # find the mapping (either by micName (without ext) or micId)
            mic = micDict.get(micKey, None) or micDict.get(micKey2, None)

            if mic is None:
                print "Skipping particle, key %s not found" % micKey
            else:
                newCoord.copyObjId(particle)
                x, y = coord.getPosition()
                newCoord.setPosition(x * scale, y * scale)
                newCoord.setMicrograph(mic)
                coords.append(newCoord)

        # We convert micrographs if they are not in .mrc
        downFactor = self.ctfDownFactor.get()
        for mic in inputMics:
            # Create micrograph dir
            micName = mic.getFileName()
            micDir = self._getTmpPath(pwutils.removeBaseExt(micName))
            pwutils.makePath(micDir)
            outMic = pwutils.join(micDir, pwutils.replaceBaseExt(micName, 'mrc'))

            if downFactor != 1:
                # Replace extension by 'mrc' cause there are some formats
                # that cannot be written (such as dm3)
                import pyworkflow.em.packages.xmipp3 as xmipp3
                self.runJob("xmipp_transform_downsample",
                            "-i %s -o %s --step %f --method fourier"
                            % (micName, outMic, downFactor),
                            env=xmipp3.getEnviron())
                sps = inputMics.getScannedPixelSize() * downFactor
                self._params['scannedPixelSize'] = sps
            else:
                if micName.endswith('.mrc'):
                    pwutils.createLink(micName, outMic)
                else:
                    em.ImageHandler().convert(micName, outMic)

        # Write out coordinate files
        writeSetOfCoordinates(self._getTmpPath(), coords, inputMics)
        coords.write()
        coords.close()

    def refineCtfStep(self):
        micsSet = self._getMicrographs()
        self._defineValues()
        self._prepareCommand()

        for mic in micsSet:
            micName = mic.getFileName()
            micBase = pwutils.removeBaseExt(micName)
            micDirTmp = self._getTmpPath(pwutils.removeBaseExt(micName))
            outMic = pwutils.join(micDirTmp, pwutils.replaceBaseExt(micName, 'mrc'))
            micFnCtf = pwutils.join(micDirTmp, micBase + '.ctf')
            micFnOut = self._getCtfOutPath(micDirTmp)
            micFnCtfFit = pwutils.join(micDirTmp, micBase + '_EPA.log')
            micFnLocalCtf = pwutils.join(micDirTmp, micBase + '_local.star')

            # Update _params dictionary
            self._params['micFn'] = outMic
            self._params['gctfOut'] = micFnOut

            try:
                self.runJob(self._getProgram(), self._args % self._params,  env=self._getEnviron())
            except:
                print("ERROR: Gctf has failed for micrograph %s" % outMic)

            # move results from tmp to extra folder
            micDir = self._getExtraPath(pwutils.removeBaseExt(micName))
            pwutils.makePath(micDir)
            psdFile = self._getPsdPath(micDir)
            ctfOutFile = self._getCtfOutPath(micDir)
            ctffitFile = self._getCtfFitOutPath(micDir)
            ctflocalFile = self._getCtfLocalPath(micDir, micBase)

            pwutils.moveFile(micFnCtf, psdFile)
            pwutils.moveFile(micFnOut, ctfOutFile)
            pwutils.moveFile(micFnCtfFit, ctffitFile)
            pwutils.moveFile(micFnLocalCtf, ctflocalFile)

            # Let's clean the temporary micrographs
            pwutils.cleanPath(outMic)

    def createCtfModelStep(self):
        newPart = Particle()
        inputSet = self.inputParticles.get()
        partSet = self._createSetOfParticles()
        partSet.copyInfo(inputSet)
        #alignType = inputSet.getAlignment()

        for particle in inputSet:
            coord = particle.getCoordinate()
            pos = coord.getPosition()
            micBase = pwutils.removeBaseExt(coord.getMicName())
            ctfFn = pwutils.join(self._getExtraPath(micBase), micBase + '_local.star')

            if pwutils.exists(ctfFn):
                mdFn = md.MetaData(ctfFn)
                for row in md.iterRows(mdFn):
                    coordX = row.getValue(md.RLN_IMAGE_COORD_X)
                    coordY = row.getValue(md.RLN_IMAGE_COORD_Y)
                    if pos == (coordX, coordY):
                        newPart.copyObjId(particle)
                        newPart.copyLocation(particle)
                        newPart.setCTF(rowToCtfModel(row))
                        # TODO: copy alignment, acquisition, coords etc from particle to newPart
                        partSet.append(newPart)

        self._defineOutputs(outputParticles=partSet)
        self._defineTransformRelation(inputSet, partSet)

    #--------------------------- INFO functions --------------------------------
    def _validate(self):
        errors = []
        # Check that the program exists
        if not pwutils.exists(self._getProgram()):
            errors.append("Binary '%s' does not exits.\n"
                          "Check configuration file: \n"
                          "~/.config/scipion/scipion.conf\n"
                          "and set GCTF variables properly."
                          % self._getProgram())
        return errors

    def _citations(self):
        return ['Zhang2016']

    def _methods(self):
        if self.inputParticles.get() is None:
            return ['Input particles not available yet.']
        methods = "We calculated the CTF of "
        methods += self.getObjectTag('inputParticles')
        methods += " using Gctf [Zhang2016]. "
        methods += self.methodsVar.get('')

        if self.hasAttribute('outputParticles'):
            methods += 'Output particles: %s' % self.getObjectTag('outputParticles')

        return [methods]

    #--------------------------- UTILS functions -------------------------------
    def _defineValues(self):
        """ This function get some parameters of the micrographs"""
        self.inputMics = self._getMicrographs()
        acq = self.inputMics.getAcquisition()

        self._params = {'voltage': acq.getVoltage(),
                        'sphericalAberration': acq.getSphericalAberration(),
                        'magnification': acq.getMagnification(),
                        'ampContrast': acq.getAmplitudeContrast(),
                        'samplingRate': self.inputMics.getSamplingRate(),
                        'scannedPixelSize': self.inputMics.getScannedPixelSize(),
                        'windowSize': self.windowSize.get(),
                        'lowRes': self.lowRes.get(),
                        'highRes': self.highRes.get(),
                        # Convert from microns to Angstroms
                        'minDefocus': self.minDefocus.get() * 1e+4,
                        'maxDefocus': self.maxDefocus.get() * 1e+4
                        }

    def _prepareCommand(self):
        sampling = self._getMicrographs().getSamplingRate() * self.ctfDownFactor.get()
        # Convert digital frequencies to spatial frequencies
        self._params['sampling'] = sampling
        self._params['lowRes'] = sampling / self._params['lowRes']
        if self._params['lowRes'] > 50:
            self._params['lowRes'] = 50
        self._params['highRes'] = sampling / self._params['highRes']
        self._params['step_focus'] = 500.0
        self._argsGctf()

    def _argsGctf(self):
        self._args = " --apix %f " % self._params['sampling']
        self._args += "--kV %f " % self._params['voltage']
        self._args += "--cs %f " % self._params['sphericalAberration']
        self._args += "--ac %f " % self._params['ampContrast']
        self._args += "--dstep %f " % self._params['scannedPixelSize']
        self._args += "--defL %f " % self._params['minDefocus']
        self._args += "--defH %f " % self._params['maxDefocus']
        self._args += "--defS %f " % self._params['step_focus']
        self._args += "--astm %f " % self.astigmatism.get()
        self._args += "--resL %f " % self._params['lowRes']
        self._args += "--resH %f " % self._params['highRes']
        self._args += "--do_EPA %d " % (1 if self.doEPA else 0)
        self._args += "--boxsize %d " % self._params['windowSize']
        self._args += "--plot_res_ring %d " % (1 if self.plotResRing else 0)
        self._args += "--gid %d " % self.GPUCore.get()
        self._args += "--bfac %d " % self.bfactor.get()
        self._args += "--B_resH %f " % (2 * self._params['sampling'])
        self._args += "--overlap %f " % self.overlap.get()
        self._args += "--convsize %d " % self.convsize.get()
        self._args += "--do_Hres_ref %d " % (1 if self.doHighRes else 0)

        # local refine options
        self._args += "--do_local_refine 1 --boxsuffix _coords.star "
        self._args += "--local_radius %d " % self.locRad.get()
        self._args += "--local_avetype %d " % self.locAveType.get()
        self._args += "--local_boxsize %d " % self.locBoxSize.get()
        self._args += "--local_overlap % 0.2f " % self.locOverlap.get()
        self._args += "--local_resL %d " % self.locResL.get()
        self._args += "--local_resH %d " % self.locResH.get()
        self._args += "--refine_local_astm %d " % (1 if self.locAstm else 0)

        if getVersion() == '0.50':
            self._args += "--do_basic_rotave %d " % (1 if self.doBasicRotave else 0)
        else:
            self._args += "--EPA_oversmp %d " % self.EPAsmp.get()

            if self.doPhShEst:
                self._args += "--phase_shift_L %f " % self.phaseShiftL.get()
                self._args += "--phase_shift_H %f " % self.phaseShiftH.get()
                self._args += "--phase_shift_S %f " % self.phaseShiftS.get()
                self._args += "--phase_shift_T %d " % (1 + self.phaseShiftT.get())

        if self.doHighRes:
            self._args += "--Href_resL %d " % self.HighResL.get()
            self._args += "--Href_resH %d " % self.HighResH.get()
            self._args += "--Href_bfac %d " % self.HighResBf.get()

        self._args += "--do_validation %d " % (1 if self.doValidate else 0)
        self._args += "%(micFn)s "
        self._args += "> %(gctfOut)s"

    def _getPsdPath(self, micDir):
        return os.path.join(micDir, 'ctfEstimation.mrc')

    def _getCtfOutPath(self, micDir):
        return os.path.join(micDir, 'ctfEstimation.txt')

    def _getCtfFitOutPath(self, micDir):
        return os.path.join(micDir, 'ctfEstimation_EPA.txt')

    def _getCtfLocalPath(self, micDir, micBase):
        return os.path.join(micDir, micBase + '_local.star')

    def _getProgram(self):
        """ Return the program binary that will be used. """
        binary = os.environ['GCTF']
        program = pwutils.join(os.environ['GCTF_HOME'], 'bin', basename(binary))

        return program

    def _oldVersion(self):
        return True if getVersion() == '0.50' else False

    def _getMicrographs(self):
            return self.inputMicrographs.get()
