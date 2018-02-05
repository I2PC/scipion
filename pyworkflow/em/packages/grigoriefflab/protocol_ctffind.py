# **************************************************************************
# *
# * Authors:     Josue Gomez BLanco (jgomez@cnb.csic.es)
# *              J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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

import os
import sys
import pyworkflow.utils as pwutils
import pyworkflow.em as em
import pyworkflow.protocol.params as params
from grigoriefflab import (CTFFIND_PATH, CTFFINDMP_PATH,
                           CTFFIND4_PATH, getVersion,
                           CTFFIND4_HOME, CTFFIND_HOME)
from convert import (readCtfModel, parseCtffindOutput,
                     parseCtffind4Output)


class ProtCTFFind(em.ProtCTFMicrographs):
    """
    Estimates CTF on a set of micrographs
    using either ctffind3 or ctffind4 program.
    
    To find more information about ctffind4 go to:
    http://grigoriefflab.janelia.org/ctffind4
    """
    _label = 'ctffind'

    @classmethod
    def validateInstallation(cls):
        """ Check if the installation of this protocol is correct.
        Can't rely on package function since this is a "multi package" package
        Returning an empty list means that the installation is correct
        and there are not errors. If some errors are found, a list with
        the error messages will be returned.
        """
        missingPaths = []

        if not os.path.exists(CTFFIND4_PATH) \
                and not os.path.exists(CTFFIND_PATH):
            missingPaths.append("%s, %s : ctffind installation not found"
                                " - %s or %s" % (CTFFIND_HOME, CTFFIND4_HOME,
                                                 CTFFIND_PATH, CTFFIND4_PATH))
        return missingPaths

    def _defineProcessParams(self, form):
        form.addParam('useCtffind4', params.BooleanParam, default=True,
                      label="Use ctffind4 to estimate the CTF?",
                      help='If is true, the protocol will use ctffind4 instead of ctffind3')
        form.addParam('astigmatism', params.FloatParam, default=100.0,
                      label='Expected (tolerated) astigmatism (A)',
                      expertLevel=params.LEVEL_ADVANCED,
                      help='Astigmatism values much larger than this will be penalised '
                           '(Angstroms; set negative to remove this restraint)',
                      condition='useCtffind4')
        form.addParam('findPhaseShift', params.BooleanParam, default=False,
                      label="Find additional phase shift?", condition='useCtffind4',
                      help='If the data was collected with phase plate, this will find '
                           'additional phase shift due to phase plate',
                      expertLevel=params.LEVEL_ADVANCED)

        group = form.addGroup('Phase shift parameters')
        group.addParam('minPhaseShift', params.FloatParam, default=0.0,
                       label="Minimum phase shift (rad)", condition='findPhaseShift',
                       help='Lower bound of the search for additional phase shift. '
                            'Phase shift is of scattered electrons relative to '
                            'unscattered electrons. In radians.',
                       expertLevel=params.LEVEL_ADVANCED)
        group.addParam('maxPhaseShift', params.FloatParam, default=3.15,
                       label="Maximum phase shift (rad)", condition='findPhaseShift',
                       help='Upper bound of the search for additional phase shift. '
                            'Phase shift is of scattered electrons relative to '
                            'unscattered electrons. In radians. '
                            'Please use value between 0.10 and 3.15',
                       expertLevel=params.LEVEL_ADVANCED)
        group.addParam('stepPhaseShift', params.FloatParam, default=0.2,
                       label="Phase shift search step (rad)", condition='findPhaseShift',
                       help='Step size for phase shift search (radians)',
                       expertLevel=params.LEVEL_ADVANCED)

        form.addParam('resamplePix', params.BooleanParam, default=True,
                      label="Resample micrograph if pixel size too small?",
                      condition='useCtffind4 and _isNewCtffind4',
                      help='When the pixel is too small, Thon rings appear very thin '
                           'and near the origin of the spectrum, which can lead to '
                           'suboptimal fitting. This options resamples micrographs to '
                           'a more reasonable pixel size if needed',
                      expertLevel=params.LEVEL_ADVANCED)

        form.addParam('slowSearch', params.BooleanParam, default=True,
                      expertLevel=params.LEVEL_ADVANCED,
                      label="Slower, more exhaustive search?",
                      condition='useCtffind4 and _isNewCtffind4',
                      help="From version 4.1.5 to 4.1.8 the slow (more precise) "
                           "search was activated by default because of reports the "
                           "faster 1D search was significantly less accurate "
                           "(thanks Rado Danev & Tim Grant). "
                           "Set this parameters to *No* to get faster fits.")
    
    #--------------------------- STEPS functions ---------------------------------------------------
    def _estimateCTF(self, micFn, micDir, micName):
        """ Run ctffind, 3 or 4, with required parameters """

        doneFile = os.path.join(micDir, 'done.txt')

        if self.isContinued() and os.path.exists(doneFile):
            return

        try:
            # Create micrograph dir
            pwutils.makePath(micDir)
            downFactor = self.ctfDownFactor.get()
            scannedPixelSize = self.inputMicrographs.get().getScannedPixelSize()
            micFnMrc = self._getTmpPath(pwutils.replaceBaseExt(micFn, 'mrc'))

            if downFactor != 1:
                # Replace extension by 'mrc' because there are some formats
                # that cannot be written (such as dm3)
                import pyworkflow.em.packages.xmipp3 as xmipp3
                args = "-i %s -o %s --step %f --method fourier" % (micFn, micFnMrc, downFactor)
                self.runJob("xmipp_transform_downsample",
                            args, env=xmipp3.getEnviron())
                self._params['scannedPixelSize'] =  scannedPixelSize * downFactor
            else:
                ih = em.ImageHandler()
                if ih.existsLocation(micFn):
                    micFnMrc = self._getTmpPath(pwutils.replaceBaseExt(micFn, "mrc"))
                    ih.convert(micFn, micFnMrc, em.DT_FLOAT)
                else:
                    print >> sys.stderr, "Missing input micrograph %s" % micFn

            # Update _params dictionary
            self._params['micFn'] = micFnMrc
            self._params['micDir'] = micDir
            self._params['ctffindOut'] = self._getCtfOutPath(micDir)
            self._params['ctffindPSD'] = self._getPsdPath(micDir)

        except Exception, ex:
            print >> sys.stderr, "Some error happened: %s" % ex
            import traceback
            traceback.print_exc()

        try:
            self.runJob(self._program, self._args % self._params)
        except Exception, ex:
            print >> sys.stderr, "ctffind has failed with micrograph %s" % micFnMrc

        # Let's notify that this micrograph have been processed
        # just creating an empty file at the end (after success or failure)
        open(doneFile, 'w')
        # Let's clean the temporary mrc micrographs
        pwutils.cleanPath(micFnMrc)

    def _restimateCTF(self, ctfId):
        """ Run ctffind3 with required parameters """

        ctfModel = self.recalculateSet[ctfId]
        mic = ctfModel.getMicrograph()
        micFn = mic.getFileName()
        micDir = self._getMicrographDir(mic)

        out = self._getCtfOutPath(micDir)
        psdFile = self._getPsdPath(micDir)

        pwutils.cleanPath(out)
        micFnMrc = self._getTmpPath(pwutils.replaceBaseExt(micFn, "mrc"))
        em.ImageHandler().convert(micFn, micFnMrc, em.DT_FLOAT)

        # Update _params dictionary
        self._prepareRecalCommand(ctfModel)
        self._params['micFn'] = micFnMrc
        self._params['micDir'] = micDir
        self._params['ctffindOut'] = out
        self._params['ctffindPSD'] = psdFile

        pwutils.cleanPath(psdFile)
        try:
            self.runJob(self._program, self._args % self._params)
        except Exception, ex:
            print >> sys.stderr, "ctffind has failed with micrograph %s" % micFnMrc
        pwutils.cleanPattern(micFnMrc)

    def _createCtfModel(self, mic, updateSampling=True):
        #  When downsample option is used, we need to update the
        # sampling rate of the micrograph associeted with the CTF
        # since it could be downsampled
        if updateSampling:
            newSampling = mic.getSamplingRate() * self.ctfDownFactor.get()
            mic.setSamplingRate(newSampling)

        micDir = self._getMicrographDir(mic)
        out = self._getCtfOutPath(micDir)
        psdFile = self._getPsdPath(micDir)

        ctfModel = em.CTFModel()
        readCtfModel(ctfModel, out, ctf4=self.useCtffind4.get())
        ctfModel.setPsdFile(psdFile)
        ctfModel.setMicrograph(mic)

        return ctfModel

    def _createOutputStep(self):
        pass

    #--------------------------- INFO functions ----------------------------------------------------
    def _validate(self):
        errors = []
        thr = self.numberOfThreads.get()
        ctffind = CTFFIND4_PATH if self.useCtffind4 else CTFFIND_PATH
        if thr > 1 and not self.useCtffind4:
            ctffind = CTFFINDMP_PATH
        if not os.path.exists(ctffind):
            errors.append('Missing %s' % ctffind)

        valueStep = round(self.stepPhaseShift.get(), 2)
        valueMin = round(self.minPhaseShift.get(), 2)
        valueMax = round(self.maxPhaseShift.get(), 2)

        if not (self.minPhaseShift < self.maxPhaseShift and
                valueStep <= (valueMax-valueMin) and
                0.10 <= valueMax <= 3.15):
            errors.append('Wrong values for phase shift search.')

        return errors

    def _citations(self):
        return ['Rohou2015'] if self.useCtffind4 else ['Mindell2003']

    def _methods(self):
        if self.inputMicrographs.get() is None:
            return ['Input micrographs not available yet.']
        methods = "We calculated the CTF of %s using CTFFind. " % self.getObjectTag('inputMicrographs')
        methods += self.methodsVar.get('')
        methods += 'Output CTFs: %s' % self.getObjectTag('outputCTF')

        return [methods]

    #--------------------------- UTILS functions ---------------------------------------------------
    def _isNewCtffind4(self):
        if self.useCtffind4 and getVersion('CTFFIND4') != '4.0.15':
            return True
        else:
            return False

    def _prepareCommand(self):
        sampling = self.inputMics.getSamplingRate() * self.ctfDownFactor.get()
        # Convert digital frequencies to spatial frequencies
        self._params['sampling'] = sampling
        self._params['lowRes'] = sampling / self._params['lowRes']
        if self._params['lowRes'] > 50:
            self._params['lowRes'] = 50
        self._params['highRes'] = sampling / self._params['highRes']
        self._params['step_focus'] = 500.0
        if not self.useCtffind4:
            self._argsCtffind3()
        else:
            self._params['astigmatism'] = self.astigmatism.get()
            if self.findPhaseShift:
                self._params['phaseShift'] = "yes"
                self._params['minPhaseShift'] = self.minPhaseShift.get()
                self._params['maxPhaseShift'] = self.maxPhaseShift.get()
                self._params['stepPhaseShift'] = self.stepPhaseShift.get()
            else:
                self._params['phaseShift'] = "no"

            # ctffind >= v4.1.5
            self._params['resamplePix'] = "yes" if self.resamplePix else "no"

            self._params['slowSearch'] = "yes" if self.slowSearch else "no"

            self._argsCtffind4()

    def _prepareRecalCommand(self, ctfModel):
        line = ctfModel.getObjComment().split()
        self._defineRecalValues(ctfModel)
        # get the size and the image of psd

        imgPsd = ctfModel.getPsdFile()
        imgh = em.ImageHandler()
        size, _, _, _ = imgh.getDimensions(imgPsd)

        mic = ctfModel.getMicrograph()

        # Convert digital frequencies to spatial frequencies
        sampling = mic.getSamplingRate()
        self._params['step_focus'] = 1000.0
        self._params['sampling'] = sampling
        self._params['lowRes'] = sampling / float(line[3])
        self._params['highRes'] = sampling / float(line[4])
        self._params['minDefocus'] = min([float(line[0]), float(line[1])])
        self._params['maxDefocus'] = max([float(line[0]), float(line[1])])
        self._params['windowSize'] = size
        if not self.useCtffind4:
            self._argsCtffind3()
        else:
            self._params['astigmatism'] = self.astigmatism.get()
            if self.findPhaseShift:
                self._params['phaseShift'] = "yes"
                self._params['minPhaseShift'] = self.minPhaseShift.get()
                self._params['maxPhaseShift'] = self.maxPhaseShift.get()
                self._params['stepPhaseShift'] = self.stepPhaseShift.get()
            else:
                self._params['phaseShift'] = "no"
            # ctffind >= v4.1.5
            self._params['resamplePix'] = "yes" if self.resamplePix else "no"

            self._params['slowSearch'] = "yes" if self.slowSearch else "no"

            self._argsCtffind4()

    def _argsCtffind3(self):
        self._program = 'export NATIVEMTZ=kk ; '
        if self.numberOfThreads.get() > 1:
            self._program += 'export NCPUS=%d ;' % self.numberOfThreads.get() + CTFFINDMP_PATH
        else:
            self._program += CTFFIND_PATH
        self._args = """   << eof > %(ctffindOut)s
%(micFn)s
%(ctffindPSD)s
%(sphericalAberration)f,%(voltage)f,%(ampContrast)f,%(magnification)f,%(scannedPixelSize)f
%(windowSize)d,%(lowRes)f,%(highRes)f,%(minDefocus)f,%(maxDefocus)f,%(step_focus)f
eof
"""

    def _argsCtffind4(self):

        # Avoid threads multiplication
        # self._program = 'export OMP_NUM_THREADS=%d; ' % self.numberOfThreads.get()
        self._program = 'export OMP_NUM_THREADS=1; '
        self._program += CTFFIND4_PATH
        self._args = """ << eof
%(micFn)s
%(ctffindPSD)s
%(sampling)f
%(voltage)f
%(sphericalAberration)f
%(ampContrast)f
%(windowSize)d
%(lowRes)f
%(highRes)f
%(minDefocus)f
%(maxDefocus)f
%(step_focus)f"""

        if getVersion('CTFFIND4') in ['4.1.5', '4.1.8']:
            if self.findPhaseShift:
                self._args += """
no
%(slowSearch)s
yes
%(astigmatism)f
%(phaseShift)s
%(minPhaseShift)f
%(maxPhaseShift)f
%(stepPhaseShift)f
yes
%(resamplePix)s
eof
"""
            else:
                self._args += """
no
%(slowSearch)s
yes
%(astigmatism)f
%(phaseShift)s
yes
%(resamplePix)s
eof
"""
        elif getVersion('CTFFIND4') == '4.0.15':
            if self.findPhaseShift:
                self._args += """
%(astigmatism)f
%(phaseShift)s
%(minPhaseShift)f
%(maxPhaseShift)f
%(stepPhaseShift)f
eof
"""
            else:
                self._args += """
%(astigmatism)f
%(phaseShift)s
eof
"""

    def _getPsdPath(self, micDir):
        return os.path.join(micDir, 'ctfEstimation.mrc')

    def _getCtfOutPath(self, micDir):
        return os.path.join(micDir, 'ctfEstimation.txt')

    def _parseOutput(self, filename):
        """ Try to find the output estimation parameters
        from filename. It search for a line containing: Final Values.
        """
        if not self.useCtffind4:
            return parseCtffindOutput(filename)
        else:
            return parseCtffind4Output(filename)

    def _getCTFModel(self, defocusU, defocusV, defocusAngle, psdFile):
        ctf = em.CTFModel()
        ctf.setStandardDefocus(defocusU, defocusV, defocusAngle)
        ctf.setPsdFile(psdFile)

        return ctf

    def _summary(self):
        summary = em.ProtCTFMicrographs._summary(self)
        if self.useCtffind4 and getVersion('CTFFIND4') == '4.1.5':
            summary.append("NOTE: ctffind4.1.5 finishes correctly (all output is generated properly),"
                           " but returns an error code. Disregard error messages until this is fixed."
                           "http://grigoriefflab.janelia.org/node/5421")
        return summary

