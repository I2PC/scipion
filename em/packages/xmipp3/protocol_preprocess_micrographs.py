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
        self.addParam('inputMicrographs', PointerParam, label="Micrographs", pointerClass='SetOfMicrographs')
        
        self.addSection(label='Preprocess')
        self.addParam('doCrop', BooleanParam, default=False,
                   label='Crop borders?')
        self.addParam('cropPixels', IntParam, default=10,
                   label='Pixels to crop')
        self.addParam('doLog', BooleanParam, default=False,
                   label='Take logarithm?')
        self.addParam('logA', FloatParam, default=4.431,
                   label='a')
        self.addParam('logB', FloatParam, default=0.4018,
                   label='b')
        self.addParam('logC', FloatParam, default=336.6,
                   label='c')
        self.addParam('doRemoveBadPix', BooleanParam, default=False,
                   label='Remove bad pixels?')
        self.addParam('mulStddev', IntParam, default=5,
                   label='Multiple of Stddev')    
        self.addParam('doDownsample', BooleanParam, default=False,
                   label='Downsample micrographs?')
        self.addParam('downFactor', FloatParam, default=2.,
                   label='Downsampling factor')
    

class XmippProtPreprocessMicrographs(ProtPreprocessMicrographs):
    """Protocol to preprocess a set of micrographs in the project"""
    _definition = XmippDefPreprocessMicrograph()
    
    def __init__(self, **args):
        
        Protocol.__init__(self, **args)
        
    def defineSteps(self):
        '''for each micrograph insert the steps to preprocess it
        '''
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
            fnOut = self.getPath(os.path.basename(fn))
            self.insertStepsForMicrograph(fn, fnOut)
            IOTable[fn] = fnOut
        
        # Insert step to create output objects       
        self.insertFunctionStep('createOutput', IOTable)
        


        
    def _insertOneStep(self, condition, program, arguments):
        """Insert operation if the condition is met.
        Possible conditions are: doDownsample, doCrop...etc"""
        if condition.get():
            # If the input micrograph and output micrograph differss, 
            # add the -o option
            if self.params['inputMic'] != self.params['outputMic']:
                arguments += " -o %(outputMic)s"
            # Insert the command with the formatted parameters
            self.insertRunJobStep(program, arguments % self.params)
            # Update inputMic for next step as outputMic
            self.params['inputMic'] = self.params['outputMic']
            
        
    def insertStepsForMicrograph(self, inputMic, outputMic):
        self.params['inputMic'] = inputMic
        self.params['outputMic'] = outputMic
        # Downsample
        self._insertOneStep(self.doDownsample, "xmipp_transform_downsample",
                            "-i %(inputMic)s --step %(downFactor)f --method fourier")
                    
        # Crop
        self._insertOneStep(self.doCrop, "xmipp_transform_window",
                            " -i %(inputMic)s --crop %(cropPixels)d -v 0")
        
            
        # Take logarithm
        self._insertOneStep(self.doLog, "xmipp_transform_filter",
                            " -i %(inputMic)s --log --fa %(logA)f --fb %(logB)f --fc %(logC)f")
                    
        # Remove bad pixels
        self._insertOneStep(self.doRemoveBadPix, "xmipp_transform_filter",
                            " -i %(inputMic)s --bad_pixels outliers %(stddev)f -v 0")
                
    def createOutput(self, IOTable):

        mdOut = self.getPath("micrographs.xmd")
                
        self.outputMicrographs = XmippSetOfMicrographs(value=mdOut)     
        self.outputMicrographs.microscope.voltage.set(self.inputMics.microscope.voltage.get())
        self.outputMicrographs.microscope.sphericalAberration.set(self.inputMics.microscope.sphericalAberration.get())
        self.outputMicrographs.samplingRate.set(self.inputMics.samplingRate.get())
        
        #Create the xmipp metadata micrographs.xmd  
         
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
        
        self.outputMicrographs.setFileName(mdOut)

        self.defineOutputs(micrograph=self.outputMicrographs)
