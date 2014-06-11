# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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


from pyworkflow.utils.path import makePath, replaceBaseExt, join, basename
from pyworkflow.em import *
from brandeis import *
# from pyworkflow.utils.which import which


class ProtCTFFind(ProtCTFMicrographs):
    """Estimates CTF on a set of micrographs
    using the ctffind3 program"""
    _label = 'ctffind3'
    
    #--------------------------- STEPS functions ---------------------------------------------------
    def _estimateCTF(self, micFn, micDir):
        """ Run ctffind3 with required parameters """
        # Create micrograph dir 
        makePath(micDir)
        # Update _params dictionary
        self._params['micFn'] = micFn
        self._params['micDir'] = micDir
        self._params['ctffindOut'] = join(micDir, 'ctffind.out')
        self._params['ctffindPSD'] = self._getPsdPath(micDir)
        self.runJob(self._program, self._args % self._params)
    
    def createOutputStep(self):
        ctfSet = self._createSetOfCTF()
        defocusList = []
        
        for fn, micDir, mic in self._iterMicrographs():
            out = join(micDir, 'ctffind.out')
            psdfile = join(micDir, 'ctffind_psd.mrc')
            result = self._parseOutput(out)
            defocusU, defocusV, defocusAngle = result
            # save the values of defocus for each micrograph in a list
            ctfModel = self._getCTFModel(defocusU, defocusV, defocusAngle, psdfile)
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
        ctffind = join(os.environ['CTFFIND_HOME'], 'ctffind3.exe')
        if not exists(ctffind):
            errors.append('Missing ctffind3.exe')
        return errors
    
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
    
    def _getPsdPath(self, micDir):
        return join(micDir, 'ctffind_psd.mrc')

    def _parseOutput(self, filename):
        """ Try to find the output estimation parameters
        from filename. It search for a line containing: Final Values.
        """
        f = open(filename)
        result = None
        for line in f:
            if 'Final Values' in line:
                # Take DefocusU, DefocusV and Angle as a tuple
                # that are the first three values in the line
                result = tuple(map(float, line.split()[:3]))
                break
        f.close()
        return result
            
    def _getCTFModel(self, defocusU, defocusV, defocusAngle, psdFile):
        ctf = CTFModel()
        ctf.setDefocusU(defocusU)
        ctf.setDefocusV(defocusV)
        ctf.setDefocusAngle(defocusAngle)
        ctf.setPsdFile(psdFile)
        
        return ctf


    def _citations(self):
        return ['Mindell2003']