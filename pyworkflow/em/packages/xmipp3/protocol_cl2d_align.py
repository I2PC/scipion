# **************************************************************************
# *
# * Authors:     Josue Gomez Blanco (jgomez@cnb.csic.es)
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
# *  e-mail address 'jgomez@cnb.csic.es'
# *
# **************************************************************************
from pyworkflow.utils.process import runJob
"""
This sub-package contains wrapper around align2d Xmipp program
"""

from os.path import join, dirname, exists
from pyworkflow.em import *  
import xmipp

from convert import createXmippInputImages, readSetOfParticles
from glob import glob
       
        
class XmippProtCL2DAlign(ProtAlign):
    """ Protocol to align a set of particles. """
    _label = 'cl2d align'
    _references = ['[[http://www.ncbi.nlm.nih.gov/pubmed/20362059][Sorzano, et.al,  JSB (2010)]]']

    def _defineAlignParams(self, form):
        form.addParam('useReferenceImage', BooleanParam, default=True,
                      label='Use a Reference Image ?', 
                      help='If you set to *Yes*, you should provide a reference image'
                           'If *No*, the default generation is done by averaging'
                           'subsets of the input images.')
        form.addParam('ReferenceImage', StringParam , condition='useReferenceImage',
                      label="Reference image(s)", 
                      help='Image that will serve as class reference')
        
        form.addParam('maximumShift', IntParam, default=10,
                      label='Maximum shift:',
                      help='in pixels')
        form.addParam('numberOfIterations', IntParam, default=10, expertLevel=LEVEL_ADVANCED,
                      label='Number of iterations:',
                      help='Maximum number of iterations')
        form.addParallelSection(threads=1, mpi=3)
        
        
    def _defineSteps(self):
        """ Mainly prepare the command line for call cl2d align program"""
        
        # Convert input images if necessary
        imgsFn = createXmippInputImages(self, self.inputParticles.get())
        
        # Prepare arguments to call program: xmipp_classify_CL2D
        self._params = {'imgsFn': imgsFn, 
                        'extraDir': self._getExtraPath(),
                        'maxshift': self.maximumShift.get(),
                        'iter': self.numberOfIterations.get(),
                        }
        args = '-i %(imgsFn)s --odir %(extraDir)s --nref 1 --iter %(iter)d --maxShift %(maxshift)d'
        if self.useReferenceImage:
            args += " --ref0 " + self.ReferenceImage.get()
        else:
            args += " --nref0 1"
            
        self._defineClassifySteps("xmipp_classify_CL2D", args)
              
    def _defineClassifySteps(self, program, args, subset=''):
        self._insertRunJobStep(program, args % self._params)
        self._insertFunctionStep('createOutput')        
        
    def createOutput(self):
        """ Store the setOfParticles object 
        as result of the protocol. 
        """
        #TODO: Generate a set of Transforms
        if self.writeAlignedParticles:
            lastMdFn = self._getExtraPath("images.xmd")
            alignedMd = self._getPath("particles_aligned.xmd")
            alignedStk = self._getPath("particles_aligned.stk")
            params = "%(lastMdFn)s -o %(alignedStk)s --save_metadata_stack %(alignedMd)s"
            params += " --apply_transform --keep_input_columns" 
            
            # Apply transformation
            self.runJob(None, "xmipp_transform_geometry", params % locals())
            md = xmipp.MetaData(alignedMd)
            for label in [xmipp.MDL_ANGLE_PSI, xmipp.MDL_SHIFT_X, xmipp.MDL_SHIFT_Y, xmipp.MDL_FLIP]:
                md.removeLabel(label)
            if md.containsLabel(xmipp.MDL_ZSCORE):
                md.sort(xmipp.MDL_ZSCORE)
            md.write(alignedMd)
            
            particles = self.inputParticles.get()
            # Define the output average
            avgFile = self._getExtraPath('level_00', 'class_classes.stk')
            avg = Particle()
            avg.setLocation(1, avgFile)
            avg.copyInfo(particles)
            self._defineOutputs(outputAverage=avg)
            # Define the output of aligned particles
            imgSet = self._createSetOfParticles()
            imgSet.copyInfo(particles)
            imgSet.setHasAlignment(True)
            readSetOfParticles(alignedMd, imgSet, imgSet.hasCTF())
            imgSet.write()
            self._defineOutputs(outputParticles=imgSet)
            self._defineTransformRelation(particles, imgSet)

    def validate(self):
        errors = []
        refImage = self.ReferenceImage.get()
        if self.useReferenceImage:
            if len(refImage) > 0:
                if not exists(refImage):
                    errors.append("Cannot find the file " + refImage)
            else:
                errors.append("Please, enter an Image file")
        return errors
        
    def _summary(self):
        summary = []
        if not hasattr(self, 'outputParticles'):
            summary.append("Output alignment not ready yet.")
        else:
            summary.append("Input Images: %s" % self.inputParticles.get().getNameId())
            summary.append("Output Aligned Images: %s" % self.outputParticles.get())
        return summary
