# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              Laura del Cano (ldelcano@cnb.csic.es)
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
This sub-package contains the XmippPreprocessMicrographs protocol
"""


from pyworkflow.em import *  
from pyworkflow.utils import *  
from xmipp import MetaData, MDL_MICROGRAPH, MDL_MICROGRAPH_ORIGINAL, MDL_MICROGRAPH_TILTED, MDL_MICROGRAPH_TILTED_ORIGINAL
from pyworkflow.em.packages.xmipp3.data import *


class XmippDefPreprocessMicrograph(Form):
    """Create the definition of parameters for
    the XmippPreprocessMicrographs protocol"""
    def __init__(self):
        Form.__init__(self)
    
        self.addSection(label='Input')
        self.addParam('inputMicrographs', PointerParam, label="Micrographs", 
                      pointerClass='SetOfMicrographs',
                      help='Select the SetOfMicrograph ')
        
        self.addSection(label='Preprocess')
        self.addParam('doCrop', BooleanParam, default=False, important=True,
                      label='Crop borders?', 
                      help='Crop a given amount of pixels from each border.')
        self.addParam('cropPixels', IntParam, default=10, condition='doCrop',
                      label='Pixels to crop',
                      help='Amount of pixels you want to crop from borders.')
        self.addParam('doLog', BooleanParam, default=False, important=True,
                      label='Take logarithm?', 
                      help='Depending on your acquisition system you may need to take the logarithm '
                           'of the pixel values in order to have a linear relationship between '
                           'the gray values in the image and those in the volume. a - b ln(x+c) '
                           'by default 4.431-0.4018*LN((P1+336.6)) is applied (right one for nikon coolscan 9000)')
        self.addParam('logA', FloatParam, default=4.431, condition='doLog',
                      label='a',
                      help='Parameter a in a - b ln(x+c).')
        self.addParam('logB', FloatParam, default=0.4018, condition='doLog',
                      label='b',
                      help='Parameter b in a - b ln(x+c).')
        self.addParam('logC', FloatParam, default=336.6, condition='doLog',
                      label='c',
                      help='Parameter c in a - b ln(x+c).')
        self.addParam('doRemoveBadPix', BooleanParam, default=False, important=True,
                      label='Remove bad pixels?',
                      help='Values will be thresholded to this multiple of standard deviations. '
                           'Typical values are about 5, i.e., pixel values beyond 5 times the '
                           'standard deviation will be substituted by the local median. '
                           'Set this option to -1 for not applying it.')
        self.addParam('mulStddev', IntParam, default=5, condition='doRemoveBadPix',
                      label='Multiple of Stddev',
                      help='Multiple of standard devitation.')    
        self.addParam('doDownsample', BooleanParam, default=False, important=True,
                      label='Downsample micrographs?',
                      help='Downsample micrographs by a given factor.')
        self.addParam('downFactor', FloatParam, default=2., condition='doDownsample',
                      label='Downsampling factor',
                      help='Non-integer downsample factors are possible. Must be larger than 1.')
    

class XmippProtPreprocessMicrographs(ProtPreprocessMicrographs):
    """Protocol to preprocess a set of micrographs in the project"""
    _definition = XmippDefPreprocessMicrograph()
    
    def __init__(self, **args):
        
        Protocol.__init__(self, **args)
        
    def _defineSteps(self):
        '''for each micrograph insert the steps to preprocess it
        '''
        print "En defineSteps"
        # Get pointer to input micrographs 
        self.inputMics = self.inputMicrographs.get() 
        
        IOTable = {}
        
        # Parameters needed to preprocess the micrographs
        self.params = {'downFactor': self.downFactor.get(),
                       'cropPixels': self.cropPixels.get(),
                       'logA': self.logA.get(),
                       'logB': self.logB.get(),
                       'logC': self.logC.get(),
                       'stddev': self.mulStddev.get()}
        
        # For each micrograph insert the steps to preprocess it
        for mic in self.inputMics:
            fn = mic.getFileName()
            print "found mic ", mic
            fnOut = self._getPath(os.path.basename(fn))
            self.__insertStepsForMicrograph(fn, fnOut)
            IOTable[fn] = fnOut
        
        # Insert step to create output objects       
        self._insertFunctionStep('createOutput', IOTable)
        
    def __insertOneStep(self, condition, program, arguments):
        """Insert operation if the condition is met.
        Possible conditions are: doDownsample, doCrop...etc"""
        if condition.get():
            # If the input micrograph and output micrograph differss, 
            # add the -o option
            if self.params['inputMic'] != self.params['outputMic']:
                arguments += " -o %(outputMic)s"
            # Insert the command with the formatted parameters
            self._insertRunJobStep(program, arguments % self.params)
            # Update inputMic for next step as outputMic
            self.params['inputMic'] = self.params['outputMic']
            
    def __insertStepsForMicrograph(self, inputMic, outputMic):
        self.params['inputMic'] = inputMic
        self.params['outputMic'] = outputMic
        # Downsample
        self.__insertOneStep(self.doDownsample, "xmipp_transform_downsample",
                            "-i %(inputMic)s --step %(downFactor)f --method fourier")
                    
        # Crop
        self.__insertOneStep(self.doCrop, "xmipp_transform_window",
                            " -i %(inputMic)s --crop %(cropPixels)d -v 0")
        
            
        # Take logarithm
        self.__insertOneStep(self.doLog, "xmipp_transform_filter",
                            " -i %(inputMic)s --log --fa %(logA)f --fb %(logB)f --fc %(logC)f")
                    
        # Remove bad pixels
        self.__insertOneStep(self.doRemoveBadPix, "xmipp_transform_filter",
                            " -i %(inputMic)s --bad_pixels outliers %(stddev)f -v 0")
                
    def createOutput(self, IOTable):
        
        mdOut = self._getPath("micrographs.xmd")    
            
        # Create the xmipp metadata micrographs.xmd           
        md = MetaData()      
        for i, v in IOTable.iteritems():
            objId = md.addObject()
            md.setValue(MDL_MICROGRAPH,v,objId)
            md.setValue(MDL_MICROGRAPH_ORIGINAL,i,objId)
            #TODO: Handle Tilted micrographs
#            if tiltPairs:
#                MD.setValue(xmipp.MDL_MICROGRAPH_TILTED,IOTable[fnMicrographTilted],objId)
#                MD.setValue(xmipp.MDL_MICROGRAPH_TILTED_ORIGINAL,fnMicrographTilted,objId)
        md.write("micrographs"+"@"+mdOut)
                
        # Create the SetOfMicrographs object on the database
        self.outputMicrographs = XmippSetOfMicrographs(filename=mdOut)     
        self.outputMicrographs.copyInfo(self.inputMics)

        self._defineOutputs(micrograph=self.outputMicrographs)
