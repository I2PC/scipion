# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              Laura del Cano (ldelcano@cnb.csic.es)
# *              Josue Gomez Blanco (jgomez@cnb.csic.es)
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
from pyworkflow.utils.path import makePath, moveFile, removeBaseExt
from convert import *
from xmipp3 import XmippMdRow
from pyworkflow.protocol.constants import LEVEL_EXPERT, STEPS_PARALLEL
import xmipp


class XmippProtCTFDiscrepancy(ProtCTFMicrographs):
    """Protocol to estimate CTF on a set of micrographs using xmipp3"""
    _label = 'ctf discrepancy'
    
    def __init__(self, **args):
        ProtCTFMicrographs.__init__(self, **args)
        #uncomment if you do not inherit from ProtCTFMicrographs
        #self.methodsInfo = String()
        #self.stepsExecutionMode = STEPS_PARALLEL
        self._freqDict = {}

    def _defineParams(self, form):
        form.addSection(label='Input')
        # Read N ctfs estimations
        form.addParam('inputCTFs', MultiPointerParam, label="setOfCTF1", 
                      pointerClass='SetOfCTF',
                      help='Select the first set of CTFs to compare')        
        
        #read two CTFs sets...
#        form.addParam('inputCTF1', PointerParam, label="setOfCTF1", 
#                      pointerClass='SetOfCTF',
#                      help='Select the first set of CTFs to compare')
#        form.addParam('inputCTF2', PointerParam, label="setOfCTF2", 
#                      pointerClass='SetOfCTF',
#                      help='Select the second set of CTFs to compare')
#        form.addParam('phaseMaxDifference', FloatParam, default=90., 
#                      expertLevel=LEVEL_EXPERT,
#                      label='phaseMaxDifference',
#                      help='Two CTFs are identical at frequency f'
#                           ' if the phase (wave aberration function) difference '
#                           ' is smaller than this parameter in degrees')
        form.addParallelSection(threads=4, mpi=0)       
#--------------------------- INSERT steps functions --------------------------------------------  
    def _insertAllSteps(self):
        """for each ctf insert the steps to compare it
        """
        self.setOfCTF1 = self.inputCTFs[0].get()
        self.setOfCTF2 = self.inputCTFs[1].get()
#        self.setOfCTF1 = self.inputCTF1.get()
#        self.setOfCTF2 = self.inputCTF2.get()
        
        deps = [] # Store all steps ids, final step createOutput depends on all of them
        # For each ctf pair insert the steps to process it
        # check same size, same micrographs
        for ctf in self.setOfCTF1:
            # CTF estimation with Xmipp
            stepId = self._insertFunctionStep('_computeCTFDiscrepancy' 
                                              , ctf.getObjId()
                                              ,prerequisites=[]) # Make estimation steps indepent between them
            deps.append(stepId)
        # Insert step to create output objects       
        self._insertFunctionStep('createOutputStep', prerequisites=deps)
    
    def _computeCTFDiscrepancy(self, ctfId):
        #TODO must me same micrographs
        #convert to md
        mdList = [xmipp.MetaData(), xmipp.MetaData()]
        ctfList = [self.setOfCTF1[ctfId], self.setOfCTF2[ctfId]]
        ctfRow = XmippMdRow()
        
        for md, ctf in izip(mdList, ctfList):
            objId = md.addObject()
            ctfModelToRow(ctf, ctfRow)
            micrographToRow(ctf.getMicrograph(), ctfRow)
            ctfRow.writeToMd(md, objId)

        self._freqDict[ctfId] = xmipp.errorMaxFreqCTFs2D(*mdList)
        
    def createOutputStep(self):
        ctfSet = self._createSetOfCTF()
        ctfSet.setMicrographs(self.setOfCTF1.getMicrographs())
        
        ctf = CTFModel()
        print self._freqDict
        for ctf1, ctf2 in izip(self.setOfCTF1, self.setOfCTF2):
            ctfId = ctf1.getObjId()
            ctf.setDefocusU( (ctf1.getDefocusU() + ctf2.getDefocusU())/2.0 )
            ctf.setDefocusV( (ctf1.getDefocusV() + ctf2.getDefocusV())/2.0 )
            ctf.setDefocusAngle( (ctf1.getDefocusAngle() + ctf2.getDefocusAngle())/2.0 )
            ctf.setMicrograph(ctf1.getMicrograph())
            ctf.setObjId(ctfId)
            ctf.discrepancy = Float(self._freqDict[ctfId])
            # save the values of defocus for each micrograph in a list
            ctfSet.append(ctf)
        
        self._defineOutputs(outputCTF=ctfSet)
        
    def _citations(self):
        return ['Marabini2014a']
    
    def _summary(self):
        message = []
        #TODO size de la cosa calculada
        ####message.append("Comparered <%d> micrograph" % (size,'micrographs'))
        return message    
    
    def _methods(self):
        pass#nothing here
    
    def _validate(self):
        """ The function of this hook is to add some validation before the protocol
        is launched to be executed. It should return a list of errors. If the list is
        empty the protocol can be executed.
        """
        #same micrographs in both CTF??
        errors = [ ] 
        # Add some errors if input is not valid
        return errors
