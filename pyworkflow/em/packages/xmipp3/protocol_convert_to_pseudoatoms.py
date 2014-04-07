# **************************************************************************
# *
# * Authors:  Carlos Oscar Sanchez Sorzano (coss@cnb.csic.es), May 2013
# *           Slavica Jonic                (jonic@impmc.upmc.fr)
# * Ported to Scipion:
# *           J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es), Jan 2014
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

from protocol_convert_to_pseudoatoms_base import *

class XmippProtConvertToPseudoAtoms(XmippProtConvertToPseudoAtomsBase):
    """ Converts an EM volume into pseudoatoms """
    _label = 'convert to pseudoatoms'
    
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputStructure', PointerParam, label="Input structure", important=True, 
                      pointerClass='Volume')
        XmippProtConvertToPseudoAtomsBase._defineParams(self,form)
        form.addParallelSection(threads=4, mpi=1)    
             
    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        inputStructure = self.inputStructure.get()
        self.sampling = inputStructure.getSamplingRate()
        self.fnIn=getImageLocation(inputStructure)

        fnMask = self._insertMaskStep(self.fnIn)
        self._insertFunctionStep('convertToPseudoAtomsStep', self.fnIn, fnMask)
        self._insertFunctionStep('createChimeraScriptStep',inputStructure,self.fnIn)
        self._insertFunctionStep('createOutputStep')
        
        
    #--------------------------- STEPS functions --------------------------------------------
    def createOutputStep(self):
        pdb = PdbFile(self._getPath('pseudoatoms.pdb'), pseudoatoms=True)
        self._defineOutputs(outputPdb=pdb)
        self._defineSourceRelation(self.inputStructure.get(), self.outputPdb)


    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        summary = []
        summary.append('Pseudoatom radius (voxels): %f'%self.pseudoAtomRadius.get())
        summary.append('Approximation target error (%%): %f'%self.pseudoAtomTarget.get())
        return summary

    def _methods(self):
        summary = []
        summary.append('We converted the volume %s into a pseudoatomic representation with Gaussian atoms (sigma=%f A and a target error'\
                       ' of %f%%) [Nogales2013].'%(self.inputStructure.get().getNameId(),
                                     self.pseudoAtomRadius.get()*self.inputStructure.get().getSamplingRate(),
                                     self.pseudoAtomTarget.get()));
        if self.hasAttribute('outputPdb'):
            summary.append('We refer to the pseudoatomic model as %s.'%self.outputPdb.getNameId())
        return summary

    def _citations(self):
        return ['Nogales2013']
        