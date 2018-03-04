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
import sys
import pyworkflow.utils as pwutils
import pyworkflow.em as em
import pyworkflow.protocol.params as params
from pyworkflow import VERSION_1_2
from grigoriefflab import CTFTILT_PATH, CTFTILTMP_PATH, CTFFIND_HOME
from convert import readCtfModel, parseCtftiltOutput


class ProtCTFTilt(em.ProtCTFMicrographs):
    """
    Estimates CTF on a set of tilted micrographs
    using ctftilt program.
    
    """
    _label = 'ctftilt'
    _lastUpdateVersion = VERSION_1_2

    @classmethod
    def validateInstallation(cls):
        """ Check if the installation of this protocol is correct.
        Can't rely on package function since this is a "multi package" package
        Returning an empty list means that the installation is correct
        and there are no errors. If some errors are found, a list with
        the error messages will be returned.
        """
        missingPaths = []

        if not os.path.exists(CTFTILT_PATH):
            missingPaths.append("%s : ctffind3/ctftilt installation not found"
                                " - %s" % (CTFFIND_HOME, CTFTILT_PATH))
        return missingPaths

    def _defineProcessParams(self, form):
        form.addParam('astigmatism', params.FloatParam, default=100.0,
                      label='Expected astigmatism (A)',
                      expertLevel=params.LEVEL_ADVANCED,
                      help='Expected amount of astigmatism in Angstrom. ')

        line = form.addLine('Tilt angle',
                            help='Expected tilt angle value and its uncertainty in degrees.')
        line.addParam('tiltA', params.FloatParam, default=0.,
                      label='Expected value')
        line.addParam('tiltR', params.FloatParam, default=5.,
                      label='Uncertainty')

    #--------------------------- STEPS functions ------------------------------
    def _estimateCTF(self, micFn, micDir, micName):
        """ Run ctftilt with required parameters """

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
            self._params['ctftiltOut'] = self._getCtfOutPath(micDir)
            self._params['ctftiltPSD'] = self._getPsdPath(micDir)

        except Exception, ex:
            print >> sys.stderr, "Some error happened: %s" % ex
            import traceback
            traceback.print_exc()

        try:
            self.runJob(self._program, self._args % self._params)
        except Exception, ex:
            print >> sys.stderr, "ctftilt has failed with micrograph %s" % micFnMrc

        # Let's notify that this micrograph have been processed
        # just creating an empty file at the end (after success or failure)
        open(doneFile, 'w')
        # Let's clean the temporary mrc micrographs
        pwutils.cleanPath(micFnMrc)

    def _restimateCTF(self, ctfId):
        """ Run ctftilt with required parameters """

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
        self._params['ctftiltOut'] = out
        self._params['ctftiltPSD'] = psdFile

        pwutils.cleanPath(psdFile)
        try:
            self.runJob(self._program, self._args % self._params)
        except Exception, ex:
            print >> sys.stderr, "ctftilt has failed with micrograph %s" % micFnMrc
        pwutils.cleanPattern(micFnMrc)

    def _createCtfModel(self, mic, updateSampling=True):
        #  When downsample option is used, we need to update the
        # sampling rate of the micrograph associated with the CTF
        # since it could be downsampled
        if updateSampling:
            newSampling = mic.getSamplingRate() * self.ctfDownFactor.get()
            mic.setSamplingRate(newSampling)

        micDir = self._getMicrographDir(mic)
        out = self._getCtfOutPath(micDir)
        psdFile = self._getPsdPath(micDir)

        ctfModel = em.CTFModel()
        readCtfModel(ctfModel, out, ctf4=False, ctfTilt=True)
        ctfModel.setPsdFile(psdFile)
        ctfModel.setMicrograph(mic)

        return ctfModel

    def _createOutputStep(self):
        pass

    #--------------------------- INFO functions -------------------------------
    def _validate(self):
        errors = []
        thr = self.numberOfThreads.get()
        ctftilt = CTFTILT_PATH if thr > 1 else CTFTILTMP_PATH
        if not os.path.exists(ctftilt):
            errors.append('Missing %s' % ctftilt)

        return errors

    def _citations(self):
        return ['Mindell2003']

    def _methods(self):
        if self.inputMicrographs.get() is None:
            return ['Input micrographs not available yet.']
        methods = "We calculated the CTF of %s using CTFTilt. " % self.getObjectTag('inputMicrographs')
        methods += self.methodsVar.get('')
        methods += 'Output CTFs: %s' % self.getObjectTag('outputCTF')

        return [methods]

    #--------------------------- UTILS functions ------------------------------
    def _prepareCommand(self):
        sampling = self.inputMics.getSamplingRate() * self.ctfDownFactor.get()
        # Convert digital frequencies to spatial frequencies
        self._params['sampling'] = sampling
        self._params['lowRes'] = sampling / self._params['lowRes']
        if self._params['lowRes'] > 50:
            self._params['lowRes'] = 50
        self._params['highRes'] = sampling / self._params['highRes']
        self._params['astigmatism'] = self.astigmatism.get()
        self._params['step_focus'] = 500.0
        self._params['pixelAvg'] = 1  # set to 1 since we have our own downsampling
        self._params['tiltAngle'] = self.tiltA.get()
        self._params['tiltR'] = self.tiltR.get()
        self._argsCtftilt()

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
        self._params['astigmatism'] = self.astigmatism.get()
        self._params['windowSize'] = size
        self._params['pixelAvg'] = 1  # set to 1 since we have our own downsampling
        self._params['tiltAngle'] = self.tiltA.get()
        self._params['tiltR'] = self.tiltR.get()
        self._argsCtftilt()
       
    def _argsCtftilt(self):
        self._program = 'export NATIVEMTZ=kk ; '
        if self.numberOfThreads.get() > 1:
            self._program += 'export NCPUS=%d ;' % self.numberOfThreads.get() + CTFTILTMP_PATH
        else:
            self._program += CTFTILT_PATH
        self._args = """   << eof > %(ctftiltOut)s
%(micFn)s
%(ctftiltPSD)s
%(sphericalAberration)f,%(voltage)f,%(ampContrast)f,%(magnification)f,%(scannedPixelSize)f,%(pixelAvg)d
%(windowSize)d,%(lowRes)f,%(highRes)f,%(minDefocus)f,%(maxDefocus)f,%(step_focus)f,%(astigmatism)f,%(tiltAngle)f,%(tiltR)f
eof
"""

    def _getPsdPath(self, micDir):
        return os.path.join(micDir, 'ctfEstimation.mrc')

    def _getCtfOutPath(self, micDir):
        return os.path.join(micDir, 'ctfEstimation.txt')

    def _parseOutput(self, filename):
        """ Try to find the output estimation parameters
        from filename. It searches for a line containing: Final Values.
        """
        return parseCtftiltOutput(filename)

    def _getCTFModel(self, defocusU, defocusV, defocusAngle, psdFile):
        ctf = em.CTFModel()
        ctf.setStandardDefocus(defocusU, defocusV, defocusAngle)
        ctf.setPsdFile(psdFile)

        return ctf

    def _summary(self):
        summary = em.ProtCTFMicrographs._summary(self)
        if hasattr(self, 'outputCTF'):
            ctfs = self.outputCTF
            for ctf in ctfs:
                angle = float(ctf._ctftilt_tiltAngle)
                axis = float(ctf._ctftilt_tiltAxis)
                summary.append('Estimated tilt parameters:\n - tilt angle _%0.2f_\n'
                               ' - tilt axis _%0.2f_' %
                               (angle, axis))
                summary.append('If you think that tilt angle should have an '
                               'opposite sign than reported, use the following '
                               'values:\n - tilt angle _%0.2f_\n'
                               ' - tilt axis _%0.2f_' %
                               (-angle, axis + 180.0) )

        return summary
