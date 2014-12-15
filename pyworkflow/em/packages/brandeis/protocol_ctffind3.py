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


from pyworkflow.utils.path import makePath, replaceBaseExt, join, basename, cleanPath, getExt, cleanPattern
from pyworkflow.em import *
from brandeis import *
from convert import parseCtffindOutput
# from pyworkflow.utils.which import which


class ProtBaseCTFFind():
    """ This class cointains the common functionalities for all protocols to
estimate CTF on a set of micrographs using ctffind """
    
    #--------------------------- INFO functions ----------------------------------------------------
    def _validate(self):
        errors = []
        ctffind = join(os.environ['CTFFIND_HOME'], 'ctffind3.exe')
        if not exists(ctffind):
            errors.append('Missing ctffind3.exe')
        return errors
    
    def _citations(self):
        return ['Mindell2003']

    def _methods(self):
        str="We calculated the CTF using CtfFind [Midell2003]."
        if self.methodsInfo.hasValue():
            str+=" "+self.methodsInfo.get()
        return [str]
    
    #--------------------------- UTILS functions ---------------------------------------------------
    def _getPsdPath(self, micDir):
        return join(micDir, 'ctffind_psd.mrc')
    
    def _getCtfOutPath(self, micDir):
        return join(micDir, 'ctffind.out')
    
    def _parseOutput(self, filename):
        """ Try to find the output estimation parameters
        from filename. It search for a line containing: Final Values.
        """
        return parseCtffindOutput(filename)
            
    def _getCTFModel(self, defocusU, defocusV, defocusAngle, psdFile):
        ctf = CTFModel()
        ctf.setStandardDefocus(defocusU, defocusV, defocusAngle)
        ctf.setPsdFile(psdFile)
        
        return ctf

class ProtCTFFind(ProtBaseCTFFind, ProtCTFMicrographs):
    """Estimates CTF on a set of micrographs
    using the ctffind3 program"""
    _label = 'ctffind3'
    
    #--------------------------- STEPS functions ---------------------------------------------------
    def _estimateCTF(self, micFn, micDir):
        """ Run ctffind3 with required parameters """
        # Create micrograph dir 
        makePath(micDir)
        downFactor = self.ctfDownFactor.get()
        
        micFnMrc = self._getTmpPath(replaceBaseExt(micFn, 'mrc'))
        if downFactor != 1:
            #Replace extension by 'mrc' cause there are some formats that cannot be written (such as dm3)
            self.runJob("xmipp_transform_downsample","-i %s -o %s --step %f --method fourier" % (micFn, micFnMrc, downFactor))
            
            self._params['samplingRate'] = self._params['samplingRate'] * downFactor
            self._params['scannedPixelSize'] = self._params['scannedPixelSize'] * downFactor
        else:
            micFnMrc = self._getTmpPath(replaceBaseExt(micFn, "mrc"))
            ImageHandler().convert(micFn, micFnMrc, DT_FLOAT)
        
        # Update _params dictionary
        self._params['micFn'] = micFnMrc
        self._params['micDir'] = micDir
        self._params['ctffindOut'] = self._getCtfOutPath(micDir)
        self._params['ctffindPSD'] = self._getPsdPath(micDir)
        try:
            self.runJob(self._program, self._args % self._params)
        except Exception:
            raise
        cleanPath(micFnMrc)
    
    def createOutputStep(self):
        ctfSet = self._createSetOfCTF()
        ctfSet.setMicrographs(self.inputMics)
        defocusList = []
        
        for fn, micDir, mic in self._iterMicrographs():
            
            mic.setSamplingRate(self._params['samplingRate'])
            
            out = self._getCtfOutPath(micDir)
            psdFile = self._getPsdPath(micDir)
            result = self._parseOutput(out)
            defocusU, defocusV, defocusAngle = result
            # save the values of defocus for each micrograph in a list
            ctfModel = self._getCTFModel(defocusU, defocusV, defocusAngle, psdFile)
            
            ctfModel.setMicrograph(mic)
            
            defocusList.append(ctfModel.getDefocusU())
            defocusList.append(ctfModel.getDefocusV())
            ctfSet.append(ctfModel)
        
        self._defocusMaxMin(defocusList)
        self._defineOutputs(outputCTF=ctfSet)
        self._defineCtfRelation(self.inputMics, ctfSet)
    
    #--------------------------- UTILS functions ---------------------------------------------------
    def _prepareCommand(self):
        self._params['step_focus'] = 1000.0
        # Convert digital frequencies to spatial frequencies
        sampling = self.inputMics.getSamplingRate()
        self._params['lowRes'] = sampling / self._params['lowRes']
        self._params['highRes'] = sampling / self._params['highRes']        
        self._program = 'export NATIVEMTZ=kk ; ' + CTFFIND_PATH
        self._args = """   << eof > %(ctffindOut)s
%(micFn)s
%(ctffindPSD)s
%(sphericalAberration)f,%(voltage)f,%(ampContrast)f,%(magnification)f,%(scannedPixelSize)f
%(windowSize)d,%(lowRes)f,%(highRes)f,%(minDefocus)f,%(maxDefocus)f,%(step_focus)f
eof
"""


class ProtRecalculateCTFFind(ProtBaseCTFFind, ProtRecalculateCTF):
    """Re-estimate CTF on a set of micrographs
    using the ctffind3 program"""
    _label = 'ctffind re-estimation'
    
    def __init__(self, **args):
        ProtRecalculateCTF.__init__(self, **args)
    
    #--------------------------- STEPS functions ---------------------------------------------------
    def _estimateCTF(self, id):
        """ Run ctffind3 with required parameters """

        ctfModel = self.recalculateSet[id]
        mic = ctfModel.getMicrograph()
        micFn = mic.getFileName()
        micDir = self._getMicrographDir(mic)
        
        out = self._getCtfOutPath(micDir)
        psdFile = self._getPsdPath(micDir)
        
        cleanPath(out)
        cleanPath(psdFile)
        micFnMrc = self._getTmpPath(replaceBaseExt(micFn, "mrc"))
        ImageHandler().convert(micFn, micFnMrc, DT_FLOAT)

        # Update _params dictionary
        self._prepareCommand(ctfModel)
        self._params['micFn'] = micFnMrc
        self._params['micDir'] = micDir
        self._params['ctffindOut'] = out
        self._params['ctffindPSD'] = psdFile
        
        try:
            self.runJob(self._program, self._args % self._params)
        except Exception:
            raise
        cleanPattern(micFnMrc)
    
    def _createNewCtfModel(self, mic):
        micDir = self._getMicrographDir(mic)                    
        out = self._getCtfOutPath(micDir)
        psdFile = self._getPsdPath(micDir)
        result = self._parseOutput(out)
        defocusU, defocusV, defocusAngle = result
        # save the values of defocus for each micrograph in a list
        ctfModel2 = self._getCTFModel(defocusU, defocusV, defocusAngle, psdFile)
        ctfModel2.setMicrograph(mic)
        return ctfModel2
    
    #--------------------------- UTILS functions ---------------------------------------------------
    def _prepareCommand(self, ctfModel):
        line = ctfModel.getObjComment().split()
        self._defineValues(ctfModel)
        # get the size and the image of psd

        imgPsd = ctfModel.getPsdFile()
        imgh = ImageHandler()
        size, _, _, _ = imgh.getDimensions(imgPsd)
        
        mic = ctfModel.getMicrograph()
        micDir = self._getMicrographDir(mic)
        
        # Convert digital frequencies to spatial frequencies
        sampling = mic.getSamplingRate()
        self._params['step_focus'] = 1000.0
        self._params['lowRes'] = sampling / float(line[3])
        self._params['highRes'] = sampling / float(line[4])
        self._params['minDefocus'] = min([float(line[0]), float(line[1])])
        self._params['maxDefocus'] = max([float(line[0]), float(line[1])])
        self._params['windowSize'] = size
        
        self._program = 'export NATIVEMTZ=kk ; ' + CTFFIND_PATH
        self._args = """   << eof > %(ctffindOut)s
%(micFn)s
%(ctffindPSD)s
%(sphericalAberration)f,%(voltage)f,%(ampContrast)f,%(magnification)f,%(scannedPixelSize)f
%(windowSize)d,%(lowRes)f,%(highRes)f,%(minDefocus)f,%(maxDefocus)f,%(step_focus)f
eof
"""
    
    