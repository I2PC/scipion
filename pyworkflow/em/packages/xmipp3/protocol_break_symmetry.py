# **************************************************************************
# *
# * Authors:     Roberto Marabini (roberto@cnb.csic.es)
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

from pyworkflow.object import String
from pyworkflow.protocol.params import StringParam
from pyworkflow.em.protocol import ProtProcessParticles
from convert import writeSetOfParticles
import pyworkflow.em.metadata as md

 
class XmippProtAngBreakSymmetry(ProtProcessParticles):
    """
    Given an input set of particles with angular assignment, find an
    equivalent angular assignment for a given symmetry.

    Be aware that input symmetry values follows Xmipp conventions as described in:
    http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Symmetry
    """
    _label = 'break symmetry'

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineProcessParams(self, form):
        
        form.addParam('symmetryGroup', StringParam, default="c1",
                      label='Symmetry group',
                      help="See http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Symmetry"
                           " for a description of the symmetry groups format in Xmipp.\n"
                           "If no symmetry is present, use _c1_.")
    
    def _getDefaultParallel(self):
        """This protocol doesn't have mpi version"""
        return (0, 0)
     
    #--------------------------- INSERT steps functions --------------------------------------------            
    def _insertAllSteps(self):
        """ Mainly prepare the command line for call brak symmetry program"""
        # Create a metadata with the geometrical information
        # as expected by Xmipp
        imgsFn = self._getPath('input_particles.xmd')
        self._insertFunctionStep('convertInputStep', imgsFn)
        self._insertFunctionStep('breakSymmetryStep', imgsFn)
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions --------------------------------------------

    def convertInputStep(self, outputFn):
        """ Create a metadata with the images and geometrical information. """
        writeSetOfParticles(self.inputParticles.get(), outputFn)

    #--------------------------- STEPS functions --------------------------------------------
    def breakSymmetryStep(self, imgsFn):
        outImagesMd = self._getPath('images.xmd')
        args = "-i Particles@%s --sym %s -o %s" % (imgsFn,
                                                 self.symmetryGroup.get(),
                                                 outImagesMd )
        self.runJob("xmipp_angular_break_symmetry", args)
        self.outputMd = String(outImagesMd)

    def createOutputStep(self):
        imgSet = self.inputParticles.get()
        partSet = self._createSetOfParticles()
        partSet.copyInfo(imgSet)
        partSet.copyItems(imgSet,
                          updateItemCallback=self._createItemMatrix,
                          itemDataIterator=md.iterRows(self.outputMd.get(), sortByLabel=md.MDL_ITEM_ID))
        
        self._defineOutputs(outputParticles=partSet)
        self._defineSourceRelation(imgSet, partSet)

    #--------------------------- INFO functions --------------------------------------------                
    def _summary(self):
        import os
        summary = []
        if not hasattr(self, 'outputParticles'):
            summary.append("Output particles not ready yet.")
        else:
            summary.append("Symmetry: %s"% self.symmetryGroup.get())
        return summary
    
    def _validate(self):
        pass
        
    def _citations(self):
        return []#['Vargas2013b']
    
    def _methods(self):
        methods = []
#        if hasattr(self, 'outputParticles'):
#            outParticles = len(self.outputParticles) if self.outputParticles is not None else None
#            particlesRejected = len(self.inputParticles.get())-outParticles if outParticles is not None else None
#            particlesRejectedText = ' ('+str(particlesRejected)+')' if particlesRejected is not None else ''
#            rejectionText = [
#                             '',# REJ_NONE
#                             ' and removing those not reaching %s%s' % (str(self.maxZscore.get()), particlesRejectedText),# REJ_MAXZSCORE
#                             ' and removing worst %s percent%s' % (str(self.percentage.get()), particlesRejectedText)# REJ_PERCENTAGE
#                             ]
#            methods.append('Input dataset %s of %s particles was sorted by'
#                           ' its ZScore using xmipp_image_sort_by_statistics'
#                           ' program%s. ' % (self.getObjectTag('inputParticles'), len(self.inputParticles.get()), rejectionText[self.autoParRejection.get()]))
#            methods.append('Output set is %s.'%self.getObjectTag('outputParticles'))
        return methods

    #--------------------------- Utils functions --------------------------------------------                
    def _createItemMatrix(self, item, row):
        from pyworkflow.em.packages.xmipp3.convert import createItemMatrix
        import pyworkflow.em as em
        
        createItemMatrix(item, row, align=em.ALIGN_PROJ)

