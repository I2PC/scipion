# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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
This sub-package contains protocol for particles filters operations
"""

from pyworkflow.em import ProtCreateMask2D, Mask
from pyworkflow.protocol.params import PointerParam, FloatParam
from pyworkflow.em.convert import ImageHandler 

from ..spider import runCustomMaskScript
from protocol_base import SpiderProtocol


class SpiderProtCustomMask(ProtCreateMask2D, SpiderProtocol):
    """ 
    This step creates a custom mask. 
    
    In the step following this one, dimension-reduction, the covariance of 
    the pixels in all images will be computed. Only pixels under a given 
    mask will be analyzed. If this step is performed, a mask that follows 
    closely the contour the particle of interest will be used. Absent a 
    custom-made mask, a circular mask will be used.  For non-globular structures, 
    this customized mask will reduce computational demand and the likelihood 
    of numerical inaccuracy in the next dimension-reduction step. On the other 
    hand, given the power of modern computers, this step may be unnecessary.
    """
    _label = 'custom mask'
    
    def __init__(self, **kwargs):
        ProtCreateMask2D.__init__(self, **kwargs)
        SpiderProtocol.__init__(self, **kwargs)
        
        self._params = {'ext': 'stk',
                        'inputImage': 'input_image',
                        'outputMask': 'output_mask'
                        }
    
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputImage', PointerParam, label="Input image", important=True, 
                      pointerClass='Particle',
                      help='Select the input image to create the mask. \n'
                           'It is recommended to used an average image.')        
        form.addParam('filterRadius1', FloatParam, default=0.1,
                      label='Fourier radius for input image (range: 0 - 0.5 px^-1)',
                      help='The input image will be low-pass filtered to smooth any jagged edges.')        
        form.addParam('sdFactor', FloatParam, default=0.6,
                      label='First threshold (units of standard deviations)',
                      help='The filtered image will be thresholded at the average plus this number * st.dev.')
        form.addParam('filterRadius2', FloatParam, default=0.1,
                      label='Fourier radius for intermediate mask (range: 0 - 0.5)',
                      help='The intermediate thresholded mask will be again filtered for further smoothing.')       
        form.addParam('maskThreshold', FloatParam, default=0.01,
                      label='Mask threshold (range: approx. 0 - 1)',
                      help='The filtered intermediate mask will be thresholded to generate the final mask.')
        
    def _insertAllSteps(self):
        # Store references of input converted image filename
        # and the input image object
        self.outFn = self._getFileName('inputImage')
        self.inputImg = self.inputImage.get()
        # Convert the input image to Spider format
        self._insertFunctionStep('convertInputStep', self.inputImg.getLocation(), 
                                 self._getFileName('inputImage'))
        # Run Spider script to generate the custom mask
        self._insertFunctionStep('createMaskStep', 
                                 self.filterRadius1.get(), self.sdFactor.get(),
                                 self.filterRadius2.get(), self.maskThreshold.get())
        # Create the output Mask object
        self._insertFunctionStep('createOutputStep')
        
    def convertInputStep(self, inputLoc, outputFn):
        """ Convert the input image to a Spider (with stk extension). """
        ImageHandler().convert(inputLoc, (1, outputFn))
        
    def createMaskStep(self, filterRadius1, sdFactor, filterRadius2, maskThreshold):
        """ Apply the selected filter to particles. 
        Create the set of particles.
        """
        runCustomMaskScript(filterRadius1, sdFactor,
                            filterRadius2, maskThreshold,
                            workingDir=self._getPath(), ext=self.getExt(),
                            inputImage=self._params['inputImage']+'@1',
                            outputMask=self._params['outputMask'])
                            
    def createOutputStep(self):
        maskFn = self._getFileName('outputMask')
        mask = Mask()
        mask.copyInfo(self.inputImg)
        mask.setLocation(4, maskFn)
        self._defineOutputs(outputMask=mask)
        self._defineSourceRelation(self.inputImg, mask)
            
    def _summary(self):
        pixelSize = self.inputImage.get().getSamplingRate()
        filter1Angstroms = pixelSize/self.filterRadius1.get()
        filter2Angstroms = pixelSize/self.filterRadius2.get()
        
        summary = []
        summary.append('Radius for initial Fourier filter: *%s px^-1*' % self.filterRadius1.get() )
        summary.append('              On the object scale: *%s Angstroms*' % filter1Angstroms )
        summary.append('Threshold for filtered image: *avg + %s s.d.*' % self.sdFactor.get() )
        summary.append('Radius for Fourier filter for intermediate mask: *%s px^-1*' % self.filterRadius2.get() )
        summary.append('                            On the object scale: *%s Angstroms*' % filter2Angstroms )
        summary.append('Threshold for filtered mask: *%s*' % self.maskThreshold.get() )
        
        return summary
    
    def _methods(self):
        
        pixelSize = self.inputImage.get().getSamplingRate()
        filter1Angstroms = pixelSize/self.filterRadius1.get()
        filter2Angstroms = pixelSize/self.filterRadius2.get()
        
        methods = []
        
        msg  = '\nFor multivariate data analysis, a custom mask was generated by low-pass filtering the average image ' 
        msg += 'to 1/%s Angstroms^-1, thresholding it at its average plus %s * s.d., low-pass filtering this intermediate mask ' % \
                (filter1Angstroms, self.sdFactor.get() )
        msg += 'to 1/%s Angstroms^-1, and finally thresholding this filtered mask at a value of %s.' % \
                (filter2Angstroms, self.maskThreshold.get() )
        
        methods.append(msg)
        return methods

