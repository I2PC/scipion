# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              Roberto Marabini (roberto@cnb.csic.es)
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

import collections
from itertools import izip

from pyworkflow.utils.path import cleanPath, removeBaseExt
from pyworkflow.object import Set, Float, String, Object
import pyworkflow.protocol.params as params
import pyworkflow.em as em
from pyworkflow.em.metadata import Row, MetaData

import convert 
import xmipp



class XmippProtCTFDiscrepancy(em.ProtCTFMicrographs):
    """
    Protocol to estimate the agreement between different estimation of the CTF
    for the same set of micrographs. The algorithm assumes that two CTF are consistent
    if the phase (wave aberration function) of the two CTFs are closer than 90 degrees.
    The reported resolution is the resolution at which the two CTF phases differ in 90 degrees.
    """
    _label = 'ctf discrepancy'
    
    def __init__(self, **args):
        em.ProtCTFMicrographs.__init__(self, **args)
        self._freqResol = {}
        self.stepsExecutionMode = params.STEPS_SERIAL

    def _defineParams(self, form):
        form.addSection(label='Input')
        # Read N ctfs estimations
        form.addParam('inputCTF1', params.PointerParam, pointerClass='SetOfCTF',
                      label="Reference CTF",
                      help='Reference CTF in comparison')
        form.addParam('inputCTF2', params.PointerParam, pointerClass='SetOfCTF',
                      label="target CTF",
                      help='CTF to be compared with reference CTF')
        form.addParallelSection(threads=0, mpi=0)
        
#--------------------------- INSERT steps functions --------------------------------------------  
                                
    def _insertAllSteps(self):
        """for each ctf insert the steps to compare it
        """
        deps = [] # Store all steps ids, final step createOutput depends on all of them
        # For each ctf pair insert the steps to process it
        # check same size, same micrographs
        stepId = self._insertFunctionStep('_computeCTFDiscrepancyStep'
                                              ,prerequisites=[]) # Make estimation steps independent between them
        deps.append(stepId)
        # Insert step to create output objects       
        self._insertFunctionStep('createAnalyzeFilesStep', prerequisites=deps)
    
    def _computeCTFDiscrepancyStep(self):
        self._computeCTFDiscrepancy(self.inputCTF1.get(), self.inputCTF2.get())

    def _ctfToMd(self, ctf, ctfMd):
        """ Write the proper metadata for Xmipp from a given CTF """
        ctfMd.clear()
        ctfRow = Row()
        convert.ctfModelToRow(ctf, ctfRow)
        convert.micrographToRow(ctf.getMicrograph(), ctfRow, alignType=convert.ALIGN_NONE)
        ctfRow.addToMd(ctfMd)

    def _computeCTFDiscrepancy(self, method1, method2):
        #TODO must be same micrographs
        #move to a single step, each step takes 5 sec while the function takes 0.03 sec
        #convert to md
        md1 = MetaData()
        md2 = MetaData()

        for ctf1 in method1:#reference CTF
            ctfId = ctf1.getObjId()
            self._ctfToMd(ctf1, md1)
            ctf2 = method2[ctfId]#.target
            self._ctfToMd(ctf2, md2)
            key = ctfId
            self._freqResol[key] = xmipp.errorMaxFreqCTFs2D(md1, md2)

    def createAnalyzeFilesStep(self):
        """ This method will add a column with a relative consistence table
        """

        ctfs = self._createSetOfCTF()
        inputCTF1= self.inputCTF1.get()
        ctfs.copyInfo(inputCTF1)
        ctfs.setMicrographs(inputCTF1.getMicrographs())
        for ctf in inputCTF1:
            
            ctfAux = ctf.clone()
            ctfId = ctf.getObjId()
            key = ctfId
            resolution = self._freqResol[key]
            ctfAux._discrepancy_resolution  = Float(resolution)
            ctfAux._discrepancy_astigmatism = Float(ctf.getDefocusU() - ctf.getDefocusV())
            ctfs.append(ctfAux)
        self._defineOutputs(outputCTF=ctfs)
        self._defineTransformRelation(self.inputCTF1.get(), ctfs)

    def _citations(self):
        return ['Marabini2014a']
    
    def _summary(self):
        message = []
        for i, ctfs in enumerate([self.inputCTF1, self.inputCTF2]):
            protocol = self.getMapper().getParent(ctfs.get())
            message.append("Method %d: %s" % (i+1, protocol.getClassLabel()))

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

    def _stepsCheck(self):
        # Just to avoid the stream checking inherited from ProtCTFMicrographs
        pass
    
        
