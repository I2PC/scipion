# **************************************************************************
# *
# * Authors:     Javier Vargas (jvargas@cnb.csic.es)
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


from glob import glob
from os.path import exists, basename

import pyworkflow.em.metadata as md
from pyworkflow.object import String, Float
from pyworkflow.protocol.constants import STEPS_PARALLEL, LEVEL_ADVANCED, STATUS_FINISHED
from pyworkflow.protocol.params import (PointerParam, EnumParam, FloatParam, IntParam, 
                                        BooleanParam, RelationParam, Positive)
from pyworkflow.em.data import SetOfMicrographs
from pyworkflow.em.constants import RELATION_CTF
from pyworkflow.utils.path import removeBaseExt, replaceBaseExt, removeExt
from pyworkflow.em.protocol import ProtMicrographs 
from convert import micrographToCTFParam
from xmipp3 import XmippProtocol
from pyworkflow.em import metadata
from xmipp import MetaData, MD_APPEND, MD_OVERWRITE
from pyworkflow.em.metadata.constants import MDL_MICROGRAPH, MDL_CTF_DEFOCUSU, MDL_CTF_DEFOCUSV, MDL_IMAGE, \
    MDL_CTF_DEFOCUS_ANGLE

class XmippProtHilbertPhasePlate(ProtMicrographs):
    """Protocol to extract particles from a set of coordinates"""
    _label = 'hilbert_phase_plate'
            
       
        #--------------------------- DEFINE param functions --------------------------------------------   
    def _defineParams(self, form):
        form.addSection(label='Input')
        
        form.addParam('inputMicrographs', PointerParam, label="Micrographs", 
                      pointerClass='SetOfMicrographs',
                      help='Select the original SetOfMicrographs')
       
        form.addParam('ctfRelations', RelationParam, allowsNull=True,
                      relationName=RELATION_CTF, attributeName='getInputMicrographs',
                      label='CTF estimation',
                      help='Choose some CTF estimation related to input micrographs. \n'
                           'CTF estimation is needed if you want to do phase flipping or \n'
                           'you want to associate CTF information to the particles.')

        
        form.addParallelSection(threads=4, mpi=1)
                
    #--------------------------- INSERT steps functions --------------------------------------------  
    def _insertAllSteps(self):
        """for each micrograph do CTF Wiener Filtering and then Obtain the Spiral Hilbert Transform
        """
        self._createTempCTFFile()        
        self._insertFunctionStep('wienerStep')
       
        
               
    def wienerStep(self):
        
        # Generate projections from this reconstruction        
        nproc = self.numberOfMpi.get()
        nT=self.numberOfThreads.get() 
        
        mdCTF = MetaData()
        #self._getTmpPath('ctf_file.xmd')
        
        self.ctfDict = {}
        for mic in self.getInputMicrographs():            
            self.ctfDict[removeBaseExt(mic.getFileName())] = mic.getFileName()

        objIdCorr = mdCTF.addObject()
        
        prefix=1
                
        for ctf in self.ctfRelations.get():

            mdCTF.setValue(MDL_IMAGE,  self.ctfDict[self._getMicrographNameFromCTF(ctf)],objIdCorr)
            mdCTF.setValue(MDL_CTF_DEFOCUSU, ctf.getDefocusU(),objIdCorr)
            mdCTF.setValue(MDL_CTF_DEFOCUSV, ctf.getDefocusV(),objIdCorr)
            mdCTF.setValue(MDL_CTF_DEFOCUS_ANGLE, ctf.getDefocusAngle(),objIdCorr)
            
            mdCTF.write(self._getExtraPath('ctf_file.xmd'), MD_OVERWRITE)
                       
            params =  '  -i %s' % self._getExtraPath('ctf_file.xmd')
            params +=  '  -o %s' % self._getExtraPath(('vol_%03d_' % prefix)+'wiener_corrected.tif')
#            params +=  '  --sampling_rate %s' % self.getInputMicrographs().getSamplingRate()
            params +=  '  --isIsotropic'
            params +=  '  --pad 1'
            
            self.runJob("xmipp_ctf_correct_wiener2d", params, numberOfMpi=nproc, numberOfThreads=nT)
            prefix = prefix+1         
                        
        #--------------------------- UTILS functions --------------------------------------------
    def _createTempCTFFile(self):

        mdCTF = MetaData()
        #self._getTmpPath('ctf_file.xmd')
        
        self.ctfDict = {}
        for mic in self.getInputMicrographs():            
            self.ctfDict[removeBaseExt(mic.getFileName())] = mic.getFileName()
        
        for ctf in self.ctfRelations.get():
            objIdCorr = mdCTF.addObject()
            mdCTF.setValue(MDL_IMAGE,  self.ctfDict[self._getMicrographNameFromCTF(ctf)],objIdCorr)
            mdCTF.setValue(MDL_CTF_DEFOCUSU, ctf.getDefocusU(),objIdCorr)
            mdCTF.setValue(MDL_CTF_DEFOCUSV, ctf.getDefocusV(),objIdCorr)
            mdCTF.setValue(MDL_CTF_DEFOCUS_ANGLE, ctf.getDefocusAngle(),objIdCorr)
            
        mdCTF.write(self._getExtraPath('ctf_file.xmd'))
                
                
    def getInputMicrographs(self):
        """ Return the micrographs """
        return self.inputMicrographs.get()
    
    def _getMicrographDir(self, mic):
        """ Return an unique dir name for results of the micrograph. """
        return self._getExtraPath(removeBaseExt(mic.getFileName()))  
    
    def _getMicrographNameFromCTF(self, ctf):
        """ Return an unique  name for results of the micrograph. """
        return removeExt(ctf.getMicrograph().getMicName())
        

    
