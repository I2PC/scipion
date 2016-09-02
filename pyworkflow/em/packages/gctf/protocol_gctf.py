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
"""
This module contains the protocol for CTF estimation with gctf
"""

import os
from os.path import join, exists, basename
import sys
import pyworkflow.utils as pwutils
import pyworkflow.em as em
import pyworkflow.protocol.params as params
from pyworkflow.utils.properties import Message
from convert import (readCtfModel, parseGctfOutput)


class ProtGctf(em.ProtCTFMicrographs):
    """
    Estimates CTF on a set of micrographs
    using GPU-accelerated Gctf program.
    
    To find more information about Gctf go to:
    http://www.mrc-lmb.cam.ac.uk/kzhang
    """
    _label = 'CTF estimation on GPU'


    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_CTF_ESTI)
        form.addParam('recalculate', params.BooleanParam, default=False,
              condition='recalculate',
              label="Do recalculate ctf?")
        
        form.addParam('continueRun', params.PointerParam, allowsNull=True,
              condition='recalculate', label="Input previous run",
              pointerClass=self.getClassName())
        form.addHidden('sqliteFile', params.FileParam, condition='recalculate',
              allowsNull=True)
        
        form.addParam('inputMicrographs', params.PointerParam, important=True,
              condition='not recalculate', label=Message.LABEL_INPUT_MIC,
              pointerClass='SetOfMicrographs')
        form.addParam('ctfDownFactor', params.FloatParam, default=1.,
              label='CTF Downsampling factor',
              condition='not recalculate',
              help='Set to 1 for no downsampling. Non-integer downsample '
                   'factors are possible. This downsampling is only used for '
                   'estimating the CTF and it does not affect any further '
                   'calculation. Ideally the estimation of the CTF is optimal '
                   'when the Thon rings are not too concentrated at the origin '
                   '(too small to be seen) and not occupying the whole power '
                   'spectrum (since this downsampling might entail aliasing).')
        
        line = form.addLine('Resolution', condition='not recalculate',
              help='Give a value in digital frequency (i.e. between 0.0 and 0.5). '
                   'These cut-offs prevent the typical peak at the center of the '
                   'PSD and high-resolution terms where only noise exists, to '
                   'interfere with CTF estimation. The default lowest value is '
                   '0.05 but for micrographs with a very fine sampling this '
                   'may be lowered towards 0. The default highest value is '
                   '0.35, but it should be increased for micrographs with '
                   'signals extending beyond this value. However, if your '
                   'micrographs extend further than 0.35, you should consider '
                   'sampling them at a finer rate.')
        line.addParam('lowRes', params.FloatParam, default=0.05,
                      label='Lowest' )
        line.addParam('highRes', params.FloatParam, default=0.35,
                      label='Highest')
        # Switched (microns) by 'in microns' by fail in the identifier with jquery
        line = form.addLine('Defocus search range (microns)',
              expertLevel=params.LEVEL_ADVANCED,
              condition='not recalculate',
              help='Select _minimum_ and _maximum_ values for defocus search '
                   'range (in microns). Underfocus is represented by a '
                   'positive number.')
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
              help='The PSD is estimated from small patches of this size. '
                   'Bigger patches allow identifying more details. However, '
                   'since there are fewer windows, estimations are noisier.')
    
        form.addParam('plotResRing', params.BooleanParam, default=True,
              label='Plot a resolution ring on a PSD file',
              help='Whether to plot an estimated resolution ring on the '
                   'power spectrum',
              expertLevel=params.LEVEL_ADVANCED)
        form.addParam('GPUCore', params.IntParam, default=0,
              expertLevel=params.LEVEL_ADVANCED,
              label="Choose GPU core",
              help='GPU may have several cores. Set it to zero if you do '
                   'not know what we are talking about. First core index '
                   'is 0, second 1 and so on.')

        form.addSection(label='Advanced')
        form.addParam('bfactor', params.IntParam, default=150,
              expertLevel=params.LEVEL_ADVANCED,
              label="B-factor",
              help='B-factors used to decrease high resolution amplitude, A^2; '
              'suggested range 50~300 except using REBS method')
        form.addParam('doBasicRotave', params.BooleanParam, default=False,
              expertLevel=params.LEVEL_ADVANCED,
              label="Do rotational average",
              help='Do rotational average used for output CTF file. '
              'Only for nice output, will NOT be used for CTF determination.')
        form.addParam('doEPA', params.BooleanParam, default=False,
              expertLevel=params.LEVEL_ADVANCED,
              label="Do EPA",
              help='Do Equiphase average used for output CTF file. '
              'Only for nice output, will NOT be used for CTF determination.')
        form.addParam('overlap', params.FloatParam, default=0.5,
              expertLevel=params.LEVEL_ADVANCED,
              label="Overlap factor",
              help='Overlapping factor for grid boxes sampling, '
              'for windowsize=512, 0.5 means 256 pixels overlapping.')
        form.addParam('convsize', params.IntParam, default=85,
              expertLevel=params.LEVEL_ADVANCED,
              label="Boxsize for smoothing",
              help='Boxsize to be used for smoothing, '
                   'suggested 1/5 ~ 1/20 of boxsize in pixel, '
                   'e.g. 99 for 512 boxsize')

        group = form.addGroup('High-res refinement')
        group.addParam('doHighRes', params.BooleanParam, default=False,
              expertLevel=params.LEVEL_ADVANCED,
              label="Do high-resolution refinement",
              help='Whether to do High-resolution refinement or not, '
              'very useful for selecting high quality micrographs. '
              'Especially useful when your data has strong low-resolution bias')
        group.addParam('HighResL', params.FloatParam, default=15.0,
              expertLevel=params.LEVEL_ADVANCED,
              condition='doHighRes',
              label="Lowest resolution",
              help='Lowest resolution  to be used for High-resolution '
                   'refinement, in Angstroms')
        group.addParam('HighResH', params.FloatParam, default=4.0,
              expertLevel=params.LEVEL_ADVANCED,
              condition='doHighRes',
              label="Highest resolution",
              help='Highest resolution  to be used for High-resolution '
                   'refinement, in Angstroms')
        group.addParam('HighResBf', params.IntParam, default=50,
              expertLevel=params.LEVEL_ADVANCED,
              condition='doHighRes',
              label="B-factor",
              help='B-factor to be used for High-resolution '
                   'refinement, in Angstroms')

        form.addParam('doValidate', params.BooleanParam, default=False,
              expertLevel=params.LEVEL_ADVANCED,
              label="Do validation",
              help='Whether to validate the CTF determination.')

#        form.addParallelSection(threads=1, mpi=1)


    #--------------------------- STEPS functions -------------------------------
    def _estimateCTF(self, micFn, micDir, micName):
        """ Run Gctf with required parameters """
        # Create micrograph dir 
        pwutils.makePath(micDir)
        downFactor = self.ctfDownFactor.get()
        micFnMrc = self._getTmpPath(pwutils.replaceBaseExt(micFn, 'mrc'))
        micFnCtf = self._getTmpPath(pwutils.replaceBaseExt(micFn, 'ctf'))
        micFnCtfFit = self._getTmpPath(pwutils.removeBaseExt(micFn) + '_EPA.log')

        if downFactor != 1:
            # Replace extension by 'mrc' cause there are some formats
            # that cannot be written (such as dm3)
            import pyworkflow.em.packages.xmipp3 as xmipp3
            self.runJob("xmipp_transform_downsample",
                        "-i %s -o %s --step %f --method fourier"
                        % (micFn, micFnMrc, downFactor),
                        env=xmipp3.getEnviron())
            sps = self.inputMicrographs.get().getScannedPixelSize() * downFactor
            self._params['scannedPixelSize'] = sps
        else:
            em.ImageHandler().convert(micFn, micFnMrc, em.DT_FLOAT)
        
        # Update _params dictionary
        self._params['micFn'] = micFnMrc
        self._params['micDir'] = micDir
        self._params['gctfOut'] = self._getCtfOutPath(micDir)

        try:
            self.runJob(self._getProgram(), self._args % self._params)
        except Exception, ex:
            print >> sys.stderr, "Gctf has failed with micrograph %s" % micFnMrc

        psdFile = self._getPsdPath(micDir)
        ctffitFile = self._getCtfFitOutPath(micDir)
        pwutils.moveFile(micFnCtf, psdFile)
        pwutils.moveFile(micFnCtfFit, ctffitFile)
        pwutils.cleanPath(micFnMrc)
 
    def _restimateCTF(self, ctfId):
        """ Run Gctf with required parameters """

        ctfModel = self.recalculateSet[ctfId]
        mic = ctfModel.getMicrograph()
        micFn = mic.getFileName()
        micDir = self._getMicrographDir(mic)
        micFnCtf = self._getTmpPath(pwutils.replaceBaseExt(micFn, 'ctf'))
        micFnCtfFit = self._getTmpPath(pwutils.removeBaseExt(micFn) + '_EPA.log')

        out = self._getCtfOutPath(micDir)
        psdFile = self._getPsdPath(micDir)
        ctffitFile = self._getCtfFitOutPath(micDir)

        pwutils.cleanPath(out)
        micFnMrc = self._getTmpPath(pwutils.replaceBaseExt(micFn, 'mrc'))
        em.ImageHandler().convert(micFn, micFnMrc, em.DT_FLOAT)

        # Update _params dictionary
        self._prepareRecalCommand(ctfModel)
        self._params['micFn'] = micFnMrc
        self._params['micDir'] = micDir
        self._params['gctfOut'] = out
        
        pwutils.cleanPath(psdFile)
        try:
            self.runJob(self._getProgram(), self._args % self._params)
        except Exception, ex:
            print >> sys.stderr, "Gctf has failed with micrograph %s" % micFnMrc
        pwutils.moveFile(micFnCtf, psdFile)
        pwutils.moveFile(micFnCtfFit, ctffitFile)
        pwutils.cleanPattern(micFnMrc)
    
    def _createNewCtfModel(self, mic):
        micDir = self._getMicrographDir(mic)
        out = self._getCtfOutPath(micDir)
        psdFile = self._getPsdPath(micDir)
        ctfModel2 = em.CTFModel()
        readCtfModel(ctfModel2, out)
        ctfModel2.setPsdFile(psdFile)
        ctfModel2.setMicrograph(mic)
        return ctfModel2

    def _updateOutput(self, ctfSet):
        firstTime = not self.hasAttribute('outputCTF')
        ctfSet.setMicrographs(self.inputMics)
        self._computeDefocusRange(ctfSet)
        self._defineOutputs(outputCTF=ctfSet)
        if firstTime:  # define relation just once
            self._defineCtfRelation(self.inputMics, ctfSet)

    def _createOutputStep(self):
        pass
        
    #--------------------------- INFO functions --------------------------------
    def _validate(self):
        errors = []
        # Check that the program exists
        if not exists(self._getProgram()):
            errors.append("Binary '%s' does not exits.\n"
                          "Check configuration file: \n"
                          "~/.config/scipion/scipion.conf\n"
                          "and set GCTF variables properly."
                          % self._getProgram())
        return errors
    
    def _citations(self):
        return ['Zhang2016']

    def _methods(self):
        if self.inputMicrographs.get() is None:
            return ['Input micrographs not available yet.']
        methods = "We calculated the CTF of "
        methods += self.getObjectTag('inputMicrographs')
        methods += " using Gctf [Zhang2016]. "
        methods += self.methodsVar.get('')

        if self.hasAttribute('outputCTF'):
            methods += 'Output CTFs: %s' % self.getObjectTag('outputCTF')
        
        return [methods]
    
    #--------------------------- UTILS functions -------------------------------
    def _prepareCommand(self):
        sampling = self.inputMics.getSamplingRate() * self.ctfDownFactor.get()
        # Convert digital frequencies to spatial frequencies
        self._params['sampling'] = sampling
        self._params['lowRes'] = sampling / self._params['lowRes']
        if self._params['lowRes'] > 50:
            self._params['lowRes'] = 50
        self._params['highRes'] = sampling / self._params['highRes']
        self._params['step_focus'] = 500.0
        self._argsGctf()
    
    def _prepareRecalCommand(self, ctfModel):
        line = ctfModel.getObjComment().split()
        self._defineRecalValues(ctfModel)
        # get the size and the image of psd
    
        imgPsd = ctfModel.getPsdFile()
        imgh = em.ImageHandler()
        size, _, _, _ = imgh.getDimensions(imgPsd)
        
        mic = ctfModel.getMicrograph()
        micDir = self._getMicrographDir(mic)
        
        # Convert digital frequencies to spatial frequencies
        sampling = mic.getSamplingRate()
        self._params['step_focus'] = 1000.0
        self._params['sampling'] = sampling
        self._params['lowRes'] = sampling / float(line[3])
        self._params['highRes'] = sampling / float(line[4])
        self._params['minDefocus'] = min([float(line[0]), float(line[1])])
        self._params['maxDefocus'] = max([float(line[0]), float(line[1])])
        self._params['windowSize'] = size
        
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
        self._args += "--do_basic_rotave %d " % ( 1 if self.doBasicRotave else 0)
        self._args += "--overlap %f " % self.overlap.get()
        self._args += "--convsize %d " % self.convsize.get()
        self._args += "--do_Hres_ref %d " % (1 if self.doHighRes else 0)

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
    
    def _parseOutput(self, filename):
        """ Try to find the output estimation parameters
        from filename. It search for  lines containing: Final Values
        and Resolution limit.
        """
        return parseGctfOutput(filename)
    
    def _getCTFModel(self, defocusU, defocusV, defocusAngle, psdFile):
        ctf = em.CTFModel()
        ctf.setStandardDefocus(defocusU, defocusV, defocusAngle)
        ctf.setPsdFile(psdFile)
         
        return ctf

    def _getProgram(self):
        """ Return the program binary that will be used. """
        binary = os.environ['GCTF']
        program = join(os.environ['GCTF_HOME'], 'bin', basename(binary))

        return program

