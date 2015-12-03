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
from convert import (readCtfModel, parseGctfOutput)


class ProtGctf(em.ProtCTFMicrographs):
    """
    Estimates CTF on a set of micrographs
    using Gctf program.
    
    To find more information about Gctf go to:
    http://www.mrc-lmb.cam.ac.uk/kzhang
    """
    _label = 'CTF estimation on GPU'
    
    
    def _defineProcessParams(self, form):
        form.addParam('astigmatism', params.FloatParam, default=100.0,
              label='Expected (tolerated) astigmatism',
              help='Estimated astigmatism in Angstroms',
              expertLevel=params.LEVEL_ADVANCED)
        form.addParam('plotResRing', params.BooleanParam, default=True,
              label='Plot a resolution ring on a PSD file',
              help='Whether to plot an estimated resolution ring on the power spectrum',
              expertLevel=params.LEVEL_ADVANCED)
        form.addParam('GPUCore', params.IntParam, default=0,
              expertLevel=params.LEVEL_ADVANCED,
              label="Choose GPU core",
              help="GPU may have several cores. Set it to zero if you do not know what we are talking about. First core index is 0, second 1 and so on.")
    
    #--------------------------- STEPS functions ---------------------------------------------------
    def _estimateCTF(self, micFn, micDir, micName):
        """ Run Gctf with required parameters """
        # Create micrograph dir 
        pwutils.makePath(micDir)
        downFactor = self.ctfDownFactor.get()
        
        micFnMrc = self._getTmpPath(pwutils.replaceBaseExt(micFn, 'mrc'))
        micFnCtf = self._getTmpPath(pwutils.replaceBaseExt(micFn, "ctf"))
        if downFactor != 1:
            #Replace extension by 'mrc' cause there are some formats that cannot be written (such as dm3)
            import pyworkflow.em.packages.xmipp3 as xmipp3
            self.runJob("xmipp_transform_downsample","-i %s -o %s --step %f --method fourier" % (micFn, micFnMrc, downFactor), env=xmipp3.getEnviron())
            self._params['scannedPixelSize'] = self.inputMicrographs.get().getScannedPixelSize() * downFactor
        else:
            micFnMrc = self._getTmpPath(pwutils.replaceBaseExt(micFn, "mrc"))
            em.ImageHandler().convert(micFn, micFnMrc, em.DT_FLOAT)
        
        # Update _params dictionary
        self._params['micFn'] = micFnMrc
        self._params['micDir'] = micDir
        self._params['gctfOut'] = self._getCtfOutPath(micDir)
        self._params['gctfPSD'] = self._getPsdPath(micDir)
        try:
            self.runJob(self._getProgram(), self._args % self._params)
        except Exception, ex:
            print >> sys.stderr, "Gctf has failed with micrograph %s" % micFnMrc
        psdFile = self._getPsdPath(micDir)
        pwutils.moveFile(micFnCtf, psdFile)
        pwutils.cleanPath(micFnMrc)
 
    def _restimateCTF(self, ctfId):
        """ Run Gctf with required parameters """

        ctfModel = self.recalculateSet[ctfId]
        mic = ctfModel.getMicrograph()
        micFn = mic.getFileName()
        micDir = self._getMicrographDir(mic)
        micFnCtf = self._getTmpPath(pwutils.replaceBaseExt(micFn, "ctf"))

        out = self._getCtfOutPath(micDir)
        psdFile = self._getPsdPath(micDir)

        pwutils.cleanPath(out)
        micFnMrc = self._getTmpPath(pwutils.replaceBaseExt(micFn, "mrc"))
        em.ImageHandler().convert(micFn, micFnMrc, em.DT_FLOAT)

        # Update _params dictionary
        self._prepareRecalCommand(ctfModel)
        self._params['micFn'] = micFnMrc
        self._params['micDir'] = micDir
        self._params['gctfOut'] = out
        self._params['gctfPSD'] = psdFile
        
        pwutils.cleanPath(psdFile)
        try:
            self.runJob(self._getProgram(), self._args % self._params)
        except Exception, ex:
            print >> sys.stderr, "Gctf has failed with micrograph %s" % micFnMrc
        pwutils.moveFile(micFnCtf, psdFile)
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
    
    def _createOutputStep(self):
        ctfSet = self._createSetOfCTF()
        ctfSet.setMicrographs(self.inputMics)
        defocusList = []
        
        for _, micDir, mic in self._iterMicrographs():
            samplingRate = mic.getSamplingRate() * self.ctfDownFactor.get()
            mic.setSamplingRate(samplingRate)
            micFn = mic.getFileName()
            psdFile = self._getPsdPath(micDir)
            out = self._getCtfOutPath(micDir)
            
            ctfModel = em.CTFModel()
            readCtfModel(ctfModel, out)
            ctfModel.setPsdFile(psdFile)
            ctfModel.setMicrograph(mic)
            
            defocusList.append(ctfModel.getDefocusU())
            defocusList.append(ctfModel.getDefocusV())
            ctfSet.append(ctfModel)
        
        self._defocusMaxMin(defocusList)
        self._defineOutputs(outputCTF=ctfSet)
        self._defineCtfRelation(self.inputMics, ctfSet)
        
    #--------------------------- INFO functions ----------------------------------------------------
    def _validate(self):
        errors = []
        # Check that the program exists
        if not exists(self._getProgram()):
            errors.append("Binary '%s' does not exits. \n"
                          "Check configuration file: ~/.config/scipion/scipion.conf\n"
                          "and set GCTF variables properly." % self._getProgram())
            print "os.environ['GCTF_HOME']", os.environ['GCTF_HOME']
            print "os.environ['GCTF']", os.environ['GCTF']
        return errors
    
    def _citations(self):
        return ['Zhang2015']

    def _methods(self):
        if self.inputMicrographs.get() is None:
            return ['Input micrographs not available yet.']
        methods = "We calculated the CTF of %s using Gctf [Zhang2015]. " % self.getObjectTag('inputMicrographs')
        methods += self.methodsVar.get('')
        methods += 'Output CTFs: %s' % self.getObjectTag('outputCTF')
        
        return [methods]
    
    #--------------------------- UTILS functions ---------------------------------------------------
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
        self._args += "--defL %f " % self._params['minDefocus']
        self._args += "--defH %f " % self._params['maxDefocus']
        self._args += "--defS %f " % self._params['step_focus']
        self._args += "--astm %f " % self.astigmatism.get()
        self._args += "--resL %f " % self._params['lowRes']
        self._args += "--resH %f " % self._params['highRes']
        self._args += "--do_EPA 1 "
        self._args += "--boxsize %d " % self._params['windowSize']
        self._args += "--plot_res_ring %d " % (1 if self.plotResRing else 0)
        self._args += "--gid %d " % self.GPUCore.get()
        self._args += "%(micFn)s "
        self._args += "> %(gctfOut)s"
 
    def _getPsdPath(self, micDir):
        return os.path.join(micDir, 'ctfEstimation.mrc')
    
    def _getCtfOutPath(self, micDir):
        return os.path.join(micDir, 'ctfEstimation.txt')
    
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
        program = join(os.environ['GCTF_HOME'], basename(binary))

        return program

