# **************************************************************************
# *
# * Authors:  J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es), Nov 2014
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


from pyworkflow.protocol.params import PointerParam, FileParam
from pyworkflow.em.protocol import BatchProtocol
from pyworkflow.em.data import SetOfParticles, Volume
from pyworkflow.em.packages.xmipp3.convert import writeSetOfParticles



class BatchProtNMACluster(BatchProtocol):
    """ Protocol executed when a cluster is created
    from NMA images and theirs deformations.
    """
    _label = 'nma cluster'
    
    def _defineParams(self, form):
        form.addHidden('inputNmaDimred', PointerParam, pointerClass='EMObject')
        form.addHidden('sqliteFile', FileParam)
        
    #--------------------------- INSERT steps functions --------------------------------------------
        
    def _insertAllSteps(self):
        imagesMd = self._getExtraPath('images.xmd')
        outputVol = self._getExtraPath('reconstruction.vol')
        
        self._insertFunctionStep('convertInputStep', imagesMd)
        params = '-i %(imagesMd)s -o %(outputVol)s ' % locals()
        self._insertFunctionStep('reconstructStep', params)
        self._insertFunctionStep('createOutputStep', outputVol)
        
    #--------------------------- STEPS functions --------------------------------------------   
        
    def convertInputStep(self, imagesMd):
        # It is unusual to create the output in the convertInputStep,
        # but just to avoid reading twice the sqlite with particles
        inputSet = self.inputNmaDimred.get().getInputParticles()
        partSet = self._createSetOfParticles()
        partSet.copyInfo(inputSet)
        
        tmpSet = SetOfParticles(filename=self.sqliteFile.get())        
        partSet.appendFromImages(tmpSet)
        # Register outputs
        partSet.setAlignmentProj()
        self._defineOutputs(outputParticles=partSet)
        self._defineTransformRelation(inputSet, partSet)
        
        writeSetOfParticles(partSet, imagesMd)

    def reconstructStep(self, params):
        self.runJob('xmipp_reconstruct_fourier', params)
    
    def createOutputStep(self, outputVol):
        vol = Volume()
        vol.setFileName(outputVol)
        vol.setSamplingRate(self.outputParticles.getSamplingRate())
        self._defineOutputs(outputVol=vol)
        
    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        summary = []
        return summary
    
    def _validate(self):
        errors = []
        return errors
    
    def _citations(self):
        return []
    
    def _methods(self):
        return []
    
