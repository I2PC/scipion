# **************************************************************************
# *
# * Authors:     I. Foche (ifoche@cnb.csic.es)
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

from pyworkflow.object import String
from pyworkflow.protocol.params import IntParam, PointerParam, LEVEL_ADVANCED
from pyworkflow.em.protocol import ProtProcessParticles

from convert import writeSetOfParticles, writeSetOfClasses2D, xmippToLocation


        
class XmippProtDenoiseParticles(ProtProcessParticles):
    """ Remove particles noise by filtering them. 
    This filtering process is based on a projection over a basis created
    from some averages (extracted from classes). This filtering is not 
    intended for processing particles. The huge filtering they will be 
    passed through is known to remove part of the signal with the noise. 
    However this is a good method for clearly see which particle are we 
    going to process before it's done.
    """
    _label = 'denoise particles'

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineProcessParams(self, form):
        # First we customize the inputParticles param to fit our needs in this protocol
        form.getParam('inputParticles').pointerCondition = String('hasAlignment')
        form.getParam('inputParticles').help = String('Input images you want to filter. It is important that the images have alignment information with '
                                                      'respect to the chosen set of classes. This is the standard situation '
                                                      'after CL2D or ML2D.')
        form.addParam('inputClasses', PointerParam, label='Input Classes', important=True,
                      pointerClass='SetOfClasses', 
                      help='Select the input classes for the basis construction against images will be projected to.')
        
        form.addSection(label='Basis construction')
        form.addParam('maxClasses', IntParam, default=128,
                      label='Max. number of classes', expertLevel=LEVEL_ADVANCED,
                      help='Maximum number of classes.')
        form.addParam('maxPCABases', IntParam, default=200,
                      label='Number of PCA bases', expertLevel=LEVEL_ADVANCED,
                      help='Number of PCA bases.')
        form.addSection(label='Denoising')
        form.addParam('PCABases2Project', IntParam, default=200,
                      label='Number of PCA bases on which to project', expertLevel=LEVEL_ADVANCED,
                      help='Number of PCA bases on which to project.')
        
    def _getDefaultParallel(self):
        """ Return the default value for thread and MPI
        for the parallel section definition.
        """
        return (2, 4)
     
    #--------------------------- INSERT steps functions --------------------------------------------            
    def _insertAllSteps(self):
        """ Insert every step of the protocol"""
        
        # Convert input images if necessary
        self._insertFunctionStep('denoiseImages', self.inputParticles.getObjId(), self.inputClasses.getObjId()) 
        
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions --------------------------------------------
    def denoiseImages(self, inputId, inputClassesId):
        # We start preparing writing those elements we're using as input to keep them untouched
        imagesMd = self._getPath('images.xmd')
        writeSetOfParticles(self.inputParticles.get(), imagesMd)
        classesMd = self._getPath('classes.xmd')
        writeSetOfClasses2D(self.inputClasses.get(), classesMd)

        fnRoot = self._getExtraPath('pca')
        fnRootDenoised = self._getExtraPath('imagesDenoised')

        args = '-i Particles@%s --oroot %s --eigenvectors %d --maxImages %d' % (imagesMd, fnRoot, self.maxPCABases.get(), self.maxClasses.get())
        self.runJob("xmipp_image_rotational_pca", args)

        N=min(self.maxPCABases.get(), self.PCABases2Project.get())
        args='-i %s -o %s.stk --save_metadata_stack %s.xmd --basis %s.stk %d'\
             % (imagesMd, fnRootDenoised, fnRootDenoised, fnRoot, N)

        self.runJob("xmipp_transform_filter", args)

        self.outputMd = String('%s.stk' % fnRootDenoised)

    def createOutputStep(self):
        imgSet = self.inputParticles.get()
        partSet = self._createSetOfParticles()
        
        partSet.copyInfo(imgSet)
        partSet.copyItems(imgSet,
                            updateItemCallback=self._updateLocation,
                            itemDataIterator=md.iterRows(self.outputMd.get(), sortByLabel=md.MDL_ITEM_ID))
        
        self._defineOutputs(outputParticles=partSet)
        self._defineSourceRelation(imgSet, partSet)
    
    #--------------------------- INFO functions --------------------------------------------                
    def _summary(self):
        summary = []
        if not hasattr(self, 'outputParticles'):
            summary.append("Output particles not ready yet.")
        else:
            summary.append('PCA basis created by using %d classes' % len(self.inputClasses.get()))
            summary.append('Max. number of classes defined for PCA basis creation: %d' % self.maxClasses.get())
            summary.append('Max. number of PCA bases defined for PCA basis creation: %d' % self.maxPCABases.get())
            summary.append('PCA basis on which to project for denoising: %d' % self.PCABases2Project.get())
        return summary
    
    def _validate(self):
        pass
        
    def _citations(self):
        return ['zhao2013', 'ponce2011']
    
    def _methods(self):
        methods = []
        if not hasattr(self, 'outputParticles'):
            methods.append("Output particles not ready yet.")
        else:
            methods.append('An input dataset of %d particles was filtered creating a PCA basis (%d components) with '
                           'xmipp_image_rotational_pca and projecting the dataset into that base with xmipp_transform_filter.'\
                           % (len(self.inputParticles.get()), len(self.inputClasses.get())))
        return methods
    
    #--------------------------- UTILS functions --------------------------------------------
    def _updateLocation(self, item, row):
        index, filename = xmippToLocation(row.getValue(md.MDL_IMAGE))
        item.setLocation(index, filename)

