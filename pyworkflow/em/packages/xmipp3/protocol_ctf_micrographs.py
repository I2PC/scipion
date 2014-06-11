# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              Laura del Cano (ldelcano@cnb.csic.es)
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
This sub-package contains the XmippCtfMicrographs protocol
"""


from pyworkflow.em import *  
from pyworkflow.utils.path import makePath, moveFile

from convert import *
from xmipp3 import XmippMdRow
#import xmipp


class XmippProtCTFMicrographs(ProtCTFMicrographs):
    """Protocol to estimate CTF on a set of micrographs using xmipp3"""
    _label = 'ctf estimation'
    __prefix = join('%(micDir)s','xmipp_ctf')
    _templateDict = {
        # This templates are relative to a micDir
        'micrographs': 'micrographs.xmd',
        'prefix': __prefix,
        'ctfparam': __prefix +  '.ctfparam',
        'psd': __prefix + '.psd',
        'enhanced_psd': __prefix + '_enhanced_psd.xmp',
        'ctfmodel_quadrant': __prefix + '_ctfmodel_quadrant.xmp',
        'ctfmodel_halfplane': __prefix + '_ctfmodel_halfplane.xmp'
#        'ctffind_ctfparam': join('%(micDir)s', 'ctffind.ctfparam'),
#        'ctffind_spectrum': join('%(micDir)s', 'ctffind_spectrum.mrc')
        }
    
    #--------------------------- STEPS functions ---------------------------------------------------
    def _estimateCTF(self, micFn, micDir):
        """ Run the estimate CTF program """        
        # Create micrograph dir under extra directory
        print "creating path micDir=", micDir
        makePath(micDir)
        if not exists(micDir):
            raise Exception("No created dir: %s " % micDir)
        # Update _params dictionary with mic and micDir
        self._params['micFn'] = micFn
        self._params['micDir'] = self._getFilename('prefix', micDir=micDir)
        # CTF estimation with Xmipp                
        self.runJob(self._program, self._args % self._params)    
    
    def createOutputStep(self):
        ctfSet = self._createSetOfCTF()
        defocusList = []
        
        for _, micDir, mic in self._iterMicrographs():
            ctfparam = self._getFilename('ctfparam', micDir=micDir)
            
            ctfModel = readCTFModel(ctfparam, mic)
            ctfModel.setObjId(mic.getObjId())
            ctfModel.setMicrograph(mic)
            ctfSet.append(self._setPsdFiles(ctfModel, micDir))
            
            # save the values of defocus for each micrograph in a list
            defocusList.append(ctfModel.getDefocusU())
            defocusList.append(ctfModel.getDefocusV())
        
        self._defineOutputs(outputCTF=ctfSet)
        self._defineCtfRelation(self.inputMics, ctfSet)
        self._defocusMaxMin(defocusList)
        
#         # Write as a Xmipp metadata
#         mdFn = self._getPath('micrographs_ctf.xmd')
#         writeSetOfMicrographs(micSet, mdFn, self.setupMicRow)
#                    
#         tmpFn = self._getPath('micrographs_backup.xmd')
#         writeSetOfMicrographs(micSet, tmpFn, self.setupMicRow)
#   
#         # Evaluate the PSD and add some criterias
#         auxMdFn = self._getTmpPath('micrographs.xmd')
#         self.runJob("xmipp_ctf_sort_psds","-i %s -o %s" % (mdFn, auxMdFn))
#         # Copy result to output metadata
#         moveFile(auxMdFn, mdFn)
        

    def _citations(self):
        return ['Vargas2013']
    
    #--------------------------- UTILS functions ---------------------------------------------------
    def _prepareCommand(self):
        self._program = 'xmipp_ctf_estimate_from_micrograph'       
        self._args = "--micrograph %(micFn)s --oroot %(micDir)s --fastDefocus"
        
        # Mapping between base protocol parameters and the package specific command options
        self.__params = {'kV': self._params['voltage'],
                'Cs': self._params['sphericalAberration'],
                'sampling_rate': self._params['samplingRate'], 
                'ctfmodelSize': self._params['windowSize'],
                'Q0': self._params['ampContrast'],
                'min_freq': self._params['lowRes'],
                'max_freq': self._params['highRes'],
                'pieceDim': self._params['windowSize'],
                'defocus_range': (self._params['maxDefocus']-self._params['minDefocus'])/2,
                'defocusU': (self._params['maxDefocus']+self._params['minDefocus'])/2
                }
        
        for par, val in self.__params.iteritems():
            self._args += " --%s %s" % (par, str(val))

    def _setPsdFiles(self, ctfModel, micDir):
        
        ctfModel._psdFile = String(self._getFilename('psd', micDir=micDir))
        ctfModel._xmipp_enhanced_psd = String(self._getFilename('enhanced_psd', micDir=micDir))
        ctfModel._xmipp_ctfmodel_quadrant = String(self._getFilename('ctfmodel_quadrant', micDir=micDir))
        ctfModel._xmipp_ctfmodel_halfplane = String(self._getFilename('ctfmodel_halfplane', micDir=micDir))
        
        return ctfModel
    