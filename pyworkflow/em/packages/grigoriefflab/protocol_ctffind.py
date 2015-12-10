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
This module contains the protocol for CTF estimation with ctffind3
"""

import os
import sys
import pyworkflow.utils as pwutils
import pyworkflow.em as em
import pyworkflow.protocol.params as params
from grigoriefflab import CTFFIND_PATH, CTFFIND4_PATH
from convert import (readCtfModel, parseCtffindOutput, parseCtffind4Output)



class ProtCTFFind(em.ProtCTFMicrographs):
    """
    Estimates CTF on a set of micrographs
    using either ctffind3 or ctffind4 program.
    
    To find more information about ctffind4 go to:
    http://grigoriefflab.janelia.org/ctffind4
    """
    _label = 'ctffind'
    
    
    def _defineProcessParams(self, form):
        form.addParam('useCftfind4', params.BooleanParam, default=True,
              label="Use ctffind4 to estimate the CTF?",
              help='If is true, the protocol will use ctffind4 instead of ctffind3')
        form.addParam('astigmatism', params.FloatParam, default=100.0,
              label='Expected (tolerated) astigmatism', expertLevel=params.LEVEL_ADVANCED,
              condition='useCftfind4', )
        form.addParam('findPhaseShift', params.BooleanParam, default=False,
              label="Find additional phase shift?", condition='useCftfind4',
              expertLevel=params.LEVEL_ADVANCED,)
    
    #--------------------------- STEPS functions ---------------------------------------------------
    def _estimateCTF(self, micFn, micDir, micName):
        """ Run ctffind, 3 or 4, with required parameters """
        # Create micrograph dir 
        pwutils.makePath(micDir)
        downFactor = self.ctfDownFactor.get()
        
        micFnMrc = self._getTmpPath(pwutils.replaceBaseExt(micFn, 'mrc'))
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
        self._params['ctffindOut'] = self._getCtfOutPath(micDir)
        self._params['ctffindPSD'] = self._getPsdPath(micDir)
        try:
            self.runJob(self._program, self._args % self._params)
        except Exception, ex:
            print >> sys.stderr, "ctffind has failed with micrograph %s" % micFnMrc
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
    
    def _createNewCtfModel(self, mic):
        micDir = self._getMicrographDir(mic)
        out = self._getCtfOutPath(micDir)
        psdFile = self._getPsdPath(micDir)
        ctfModel2 = em.CTFModel()
        
        if not self.useCftfind4:
            readCtfModel(ctfModel2, out)
        else:
            readCtfModel(ctfModel2, out, True)
        
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
            psdFile = self._getPsdPath(micDir)
            out = self._getCtfOutPath(micDir)
            
            ctfModel = em.CTFModel()
            readCtfModel(ctfModel, out, ctf4=self.useCftfind4.get())
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
        if self.useCftfind4:
            ctffind = CTFFIND4_PATH
        else:
            ctffind = CTFFIND_PATH
        if not os.path.exists(ctffind):
            errors.append('Missing %s' % ctffind)
        return errors
    
    def _citations(self):
        return ['Mindell2003']

    def _methods(self):
        if self.inputMicrographs.get() is None:
            return ['Input micrographs not available yet.']
        methods = "We calculated the CTF of %s using CtfFind [Midell2003]. " % self.getObjectTag('inputMicrographs')
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
        if not self.useCftfind4:
            self._argsCtffind3()
        else:
            self._params['astigmatism'] = self.astigmatism.get()
            if self.findPhaseShift:
                self._params['phaseShift'] = "yes"
            else:
                self._params['phaseShift'] = "no"
            self._argsCtffind4()
    
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
        
        if not self.useCftfind4:
            self._argsCtffind3()
        else:
            self._params['astigmatism'] = self.astigmatism.get()
            if self.findPhaseShift:
                self._params['phaseShift'] = "yes"
            else:
                self._params['phaseShift'] = "no"
            self._argsCtffind4()
    
    def _argsCtffind3(self):
        self._program = 'export NATIVEMTZ=kk ; ' + CTFFIND_PATH
        self._args = """   << eof > %(ctffindOut)s
%(micFn)s
%(ctffindPSD)s
%(sphericalAberration)f,%(voltage)f,%(ampContrast)f,%(magnification)f,%(scannedPixelSize)f
%(windowSize)d,%(lowRes)f,%(highRes)f,%(minDefocus)f,%(maxDefocus)f,%(step_focus)f
eof
"""
    
    def _argsCtffind4(self):
        self._program = 'export OMP_NUM_THREADS=1; ' + CTFFIND4_PATH
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
%(step_focus)f
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
        if not self.useCftfind4:
            return parseCtffindOutput(filename)
        else:
            return parseCtffind4Output(filename)
    
    def _getCTFModel(self, defocusU, defocusV, defocusAngle, psdFile):
        ctf = em.CTFModel()
        ctf.setStandardDefocus(defocusU, defocusV, defocusAngle)
        ctf.setPsdFile(psdFile)
         
        return ctf
