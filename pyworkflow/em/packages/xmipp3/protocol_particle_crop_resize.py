# **************************************************************************
# *
# * Authors:     Airen Zaldivar (azaldivar@cnb.csic.es)
# *              Joaquin Oton   (oton@cnb.csic.es)
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

from pyworkflow.em import *  
from pyworkflow.utils import *  
import xmipp
import xmipp3
from protocol_filters import XmippProcess
from convert import createXmippInputImages, readSetOfParticles

RESIZE_OP_AUTO    = 0
RESIZE_OP_FOURIER = 1
RESIZE_OP_SPLINES = 2
RESIZE_OP_PYRAMID = 3


class XmippProtResize(XmippProcess):
    """ This class implement a protocol to change dimensions of the particles with Xmipp.
    """
    _label = "resize particles"
    
    def __init__(self):
        XmippProcess.__init__(self)
#         self._program = "xmipp_transform_window"
        self._programWindow = "xmipp_transform_window"
        self._programResize = "xmipp_image_resize"
        self._resizeFactor = Float(1)
        self._resizeOperation = Integer(0)
        
    def _defineProcessParams(self, form):        
        form.addParam('useWindowOperation', BooleanParam, default=False,
                      label='Apply a window operation?',
                      help='If you set to *Yes*, you should provide a window option.')
        form.addParam('windowOperation', EnumParam,
                      choices=['crop', 'window'],
                      condition='useWindowOperation',
                      default=0,
                      label="Window operation", display=EnumParam.DISPLAY_COMBO,
                      help='Select how do you want to change the size of the particles. \n '
                      '_resize_: you will provide the new size (in pixels) for your particles. \n '
                      '_crop_: you choose how many pixels you want to crop from each border. \n ')
        
        form.addParam('cropSize', IntParam, default=0,
                      condition='useWindowOperation and windowOperation == 0',
                      label='Crop size (px)',
                      help='This is the amount of pixels cropped. Half of the pixels are cropped from \n '
                      ' each side of the image.') 

        form.addParam('windowSize', IntParam, default=0,
                      condition='useWindowOperation and windowOperation == 1',
                      label='Window size (px)',
                      help='This is the size in pixels of the particle images.')
        
        form.addParam('useResizeOperation', BooleanParam, default=False,
                      label='Resize particles?',
                      help='If you set to *Yes*, you should provide a resize option.')
                
        form.addParam('resizeOperation', EnumParam,
                      choices=['Auto', 'Fourier', 'Splines', 'Pyramid'],
                      condition='useResizeOperation',
                      default=RESIZE_OP_AUTO,
                      label="Resize operation", display=EnumParam.DISPLAY_COMBO, expertLevel=LEVEL_ADVANCED,
                      help='Select the method to resize the images: \n '
                      '_Auto_: The method is automatically chosen according to the resize factor. \n '
                      '_Fourier_: Resize is done in Fourier space. Not available for resize factor>1. \n'
                      '_Splines_: Use splines interpolation. \n'
                      '_Pyramid_: Use positive level value to expand and negative to reduce. \n')
               
        form.addParam('resizeFactor', FloatParam, default=0.5,
                      condition='useResizeOperation and resizeOperation==%d' % (RESIZE_OP_AUTO),
                      label='Resize factor',
                      help='New image size is the old one x resize factor.')
        
        form.addParam('resizeSize', IntParam, default=0,
                      condition='useResizeOperation and resizeOperation!=%d and resizeOperation!=%d' % (RESIZE_OP_PYRAMID, RESIZE_OP_AUTO),
                      label='New image size (px)',
                      help='Size in pixels of the particle images <x> <y=x> <z=x>.')

        form.addParam('resizeLevel', IntParam, default=0,
                      condition='useResizeOperation and resizeOperation==%d' % (RESIZE_OP_PYRAMID),
                      label='Pyramid level',
                      help='Use positive value to expand and negative to reduce.')
                

    def _insertProcessStep(self, inputFn, outputFn, outputMd):
        
        windowSize = 0
        inPutFnOrig = inputFn
        outputFnOrig = outputFn
        outputMdOrig = outputMd
        outputFnTmp = self._getTmpPath(basename(outputFn))
        outputMdTmp = self._getTmpPath(basename(outputMd))
        
        applyBoth = self.useWindowOperation.get() and self.useResizeOperation
        
        imgSet = self.inputParticles.get()
        x, y, z, n = imgSet.getDimensions()
        
        if self.useWindowOperation.get():
            if self.getEnumText('windowOperation') == "crop":
                cropSize = self.cropSize.get()
                windowSize = x - 2*cropSize
                args = self._args + " --crop %(cropSize)s "
            else:
                windowSize = self.windowSize.get()
                args = self._args + " --size %(windowSize)s"
            
            if applyBoth:
                outputFn = outputFnTmp
                outputMd = outputMdTmp
            
            args += " -o " + outputFn
                
            self._insertRunJobStep(self._programWindow, args % locals())
        
        if self.useResizeOperation:
            
            if applyBoth:
                inputFn = outputMdTmp
                outputFn = outputFnOrig
                outputMd = outputMdOrig
                x = windowSize
            
            resizeSize = 0
            resizeFactor = 1.0
            resizeLevel = 0
            resizeOperation = 0
            
            if self.resizeOperation.get() == RESIZE_OP_AUTO:
                resizeFactor = self.resizeFactor.get()
                if resizeFactor == 0.5:
                    resizeOperation = RESIZE_OP_PYRAMID
                    resizeLevel = -1
                else:
                    resizeSize = round(x*resizeFactor)
                    if resizeFactor < 1:
                        resizeOperation = RESIZE_OP_FOURIER
                    else:
                        resizeOperation = RESIZE_OP_SPLINES
            else:
                resizeOperation = self.resizeOperation.get()
                if resizeOperation == RESIZE_OP_PYRAMID:
                    resizeLevel = self.resizeLevel.get()
                    resizeFactor = float(2**resizeLevel)
                else:
                    resizeSize = self.resizeSize.get()
                    resizeFactor = resizeSize/float(x)
                    
                    
            if resizeOperation == RESIZE_OP_FOURIER:
                args = self._args + " --fourier %s" % (resizeSize)
            elif resizeOperation == RESIZE_OP_SPLINES:
                args = self._args + " --dim %s" % (resizeSize)
            elif resizeOperation == RESIZE_OP_PYRAMID:
                args = self._args + " --pyramid %s" % (resizeLevel)
            
            
                
            args += " -o " + outputFn
            
            self._insertRunJobStep(self._programResize, args % locals())
            
            self._resizeFactor.set(resizeFactor)
            self._resizeOperation.set(resizeOperation)
       
    def _processOutput(self, imgSet):
        imgSet.setDownsample(1/self._resizeFactor.get())
            
    def _validate(self):
        errors = []
        
        if self.useResizeOperation and self.resizeOperation.get()==RESIZE_OP_FOURIER:
            imgSet = self.inputParticles.get()
            x, y, z, n = imgSet.getDimensions()
            if self.resizeSize.get()/float(x) > 1:
                errors.append('Fourier resize method cannot be used to increase the resolution')

        return errors
    
    def _summary(self):
        summary = []
        
        resizeTypeText = {RESIZE_OP_FOURIER:'Fourier',
                          RESIZE_OP_SPLINES:'Splines',
                          RESIZE_OP_PYRAMID:'Pyramid'}
        
        if not hasattr(self, 'outputParticles'):
            summary.append("Output images not ready yet.") 
        else:
            if self.useWindowOperation.get():
                if self.getEnumText('windowOperation') == "crop":
                    summary.append("*Crop operation*: New size %s" % (self.cropSize.get()))
                else:
                    summary.append("*Window operation*: New size %s" % (self.windowSize.get()))
            
            if self.useResizeOperation:
                summary.append("*Resize operation*: method _%s_" % (resizeTypeText.get(self._resizeOperation.get())) )
                summary.append("_Resize factor: %0.2f_" % self._resizeFactor.get() ) 
        
        return summary