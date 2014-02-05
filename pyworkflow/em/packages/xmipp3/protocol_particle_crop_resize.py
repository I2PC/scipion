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


class XmippProtResize(XmippProcess):
    """ This class implement a protocol to change dimensions of the particles with Xmipp.
    """
    _label = "resize particles"
    
    def __init__(self):
        XmippProcess.__init__(self)
#         self._program = "xmipp_transform_window"
        self._programWindow = "xmipp_transform_window"
        self._programResize = "xmipp_image_resize"
        
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
        
        form.addParam('resizeSize', IntParam, default=0,
                      condition='useResizeOperation',
                      label='New image size resize (px)',
                      help='This is the size in pixels of the particle images.')
               
#         form.addParam('resizeOperation', EnumParam, 
#                       choices=['resize', 'crop'], 
#                       default=0,
#                       label="Resize operation", display=EnumParam.DISPLAY_COMBO,
#                       help='Select how do you want to change the size of the particles. \n '
#                       '_resize_: you will provide the new size (in pixels) for your particles. \n '
#                       '_crop_: you choose how many pixels you want to crop from each border. \n ')
#         

    def _insertProcessStep(self, inputFn, outputFn, outputMd):
        
        cropSize = self.cropSize.get()
        windowSize = self.windowSize.get()
        
        print self.inputFn
        if self.useWindowOperation.get():
            print 'crop: ',self.getEnumText('windowOperation') 
#             args = ' -i %s -o %s --save_metadata_stack %s --keep_input_colums' % (self.inputFn, self.outputStk, self.outputMd)
            if self.getEnumText('windowOperation') == "crop":
                args = self._args + " --crop %(cropSize)s "
            else:
                args = self._args + " --size %(windowSize)s"
                
            if outputFn != inputFn:
                args += " -o " + outputFn
            
            self._insertRunJobStep(self._programWindow, args % locals())
        
#         if self.useResizeOperation:
#             args = ""
            
