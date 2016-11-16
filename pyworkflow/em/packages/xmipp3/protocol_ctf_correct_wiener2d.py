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

import pyworkflow.em.metadata as md
import pyworkflow.protocol.params as params

from pyworkflow.em.protocol import ProtProcessParticles
from pyworkflow.em.packages.xmipp3.convert import (writeSetOfParticles, xmippToLocation)


class XmippProtCTFCorrectWiener2D(ProtProcessParticles):
    """    
    Perform CTF correction by Wiener filtering.
    """
    _label = 'ctf_correct_wiener2d'
    
    def __init__(self, *args, **kwargs):
        ProtProcessParticles.__init__(self, *args, **kwargs)
        #self.stepsExecutionMode = STEPS_PARALLEL
        
    #--------------------------- DEFINE param functions --------------------------------------------   
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputParticles', params.PointerParam, pointerClass='SetOfParticles', 
                      label="Input particles",  
                      help='Select the input projection images .') 
        form.addParam('isIsotropic', params.BooleanParam, default='True',
                      label="Isotropic Correction", 
                      help='If true, Consider that there is not astigmatism and then it is performed an isotropic correction.') 
        form.addParam('padding_factor', params.IntParam, default=2,expertLevel=params.LEVEL_ADVANCED,
                      label="Padding factor",  
                      help='Padding factor for Wiener correction ')
        form.addParam('wiener_constant', params.FloatParam, default=-1,expertLevel=params.LEVEL_ADVANCED,
                      label="Wiener constant",  
                      help=' Wiener-filter constant (if < 0: use FREALIGN default)')
        form.addParam('correctEnvelope', params.BooleanParam, default='False',expertLevel=params.LEVEL_ADVANCED,
                      label="Correct for CTF envelope",  
                      help=' Only in cases where the envelope is well estimated correct for it')                       
        form.addParallelSection(threads=1, mpi=1)

    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('convertInputStep',self.inputParticles.get().getObjId())        
        self._insertFunctionStep('wienerStep')
        self._insertFunctionStep('createOutputStep')
        
    def convertInputStep(self, particlesId):
        """ Write the input images as a Xmipp metadata file. 
        particlesId: is only need to detect changes in
        input particles and cause restart from here.
        """
        writeSetOfParticles(self.inputParticles.get(), 
                            self._getPath('input_particles.xmd'))
    
    def wienerStep(self):
        params =  '  -i %s' % self._getPath('input_particles.xmd')
        params +=  '  -o %s' % self._getPath('corrected_ctf_particles.stk')
        params +=  '  --save_metadata_stack %s' % self._getPath('corrected_ctf_particles.xmd')
        params +=  '  --pad %s' % self.padding_factor.get()
        params +=  '  --wc %s' % self.wiener_constant.get()
        params +=  '  --sampling_rate %s' % self.inputParticles.get().getSamplingRate()

        if (self.inputParticles.get().isPhaseFlipped()):
            params +=  '  --phase_flipped '
            
        if (self.correctEnvelope):
            params +=  '  --correct_envelope '

        print params 
        nproc = self.numberOfMpi.get()
        nT=self.numberOfThreads.get() 

        self.runJob('xmipp_ctf_correct_wiener2d', 
                    params, numberOfMpi=nproc,numberOfThreads=nT)
    
    def createOutputStep(self):
        imgSet = self.inputParticles.get()
        partSet = self._createSetOfParticles()
        imgFn = self._getPath('corrected_ctf_particles.xmd')
        
        partSet.copyInfo(imgSet)
        partSet.setIsPhaseFlipped(True)
        partSet.copyItems(imgSet,
                            updateItemCallback=self._updateLocation,
                            itemDataIterator=md.iterRows(imgFn, sortByLabel=md.MDL_ITEM_ID))
        
        self._defineOutputs(outputParticles=partSet)
        self._defineSourceRelation(imgSet, partSet)
    #--------------------------- INFO functions -------------------------------------------- 
    def _validate(self):
        pass
    
    def _summary(self):
        pass
    
    def _methods(self):
        messages = []
        return messages
    
    def _citations(self):
        return ['Vargas2014a']
    
    #--------------------------- UTILS functions -------------------------------------------- 
    def _updateLocation(self, item, row):
        index, filename = xmippToLocation(row.getValue(md.MDL_IMAGE))
        item.setLocation(index, filename)

