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

from pyworkflow.em import *  
from pyworkflow.em.constants import NO_INDEX
from pyworkflow.utils import removeExt, removeBaseExt, makePath, getLastFile

      
class ProtResMap(ProtValidate3D):
    """ 
    """
        
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputVolume', PointerParam, label="Input volume", important=True, 
                      pointerClass='Volume',
                      help='Select the input volume.')
        form.addParam('useSplitVolume', BooleanParam, default=False,
                      label="Use split volume?", 
                      help='Use to split volumes for gold-standard FSC.')                
        form.addParam('splitVolume', PointerParam, label="Split volume", important=True, 
                      pointerClass='Volume', condition="useSplitVolume",
                      help='Select the second split volume.')  
        form.addParam('applyMask', BooleanParam, default=False,
                      label="Mask input volume?", 
                      help='If set to <No> ResMap will automatically compute a mask.')                
        form.addParam('maskVolume', PointerParam, label="Mask volume",  
                      pointerClass='Volume', condition="applyMask",
                      help='Select a volume to apply as a mask.')               
        
        form.addSection(label='Optional')
        form.addParam('stepRes', FloatParam, default=1,
                      label='Step size (Ang):',
                      help='in Angstroms (min 0.25, default 1.0)')
        form.addParam('minRes', FloatParam, default=0,
                      label='Min resolution (Ang):',
                      help='Default: algorithm will start a just above 2*voxelSize')
        form.addParam('maxRes', FloatParam, default=0,
                      label='Max resolutoin (Ang):',
                      help='Default: algorithm will start a just above 2*voxelSize')
        form.addParam('pVal', FloatParam, default=0.05,
                      label='Confidence level:',
                      help='usually between [0.01, 0.05]')   
                     
    def _insertAllSteps(self):
        # Insert processing steps
        inputs = [self.inputVolume.get().getLocation()]
        if self.useSplitVolume:
            inputs.append(self.splitVolume.get().getLocation())
            
        self._insertFunctionStep('convertInput', *inputs)
        self._insertFunctionStep('estimateResolution', self.pVal.get(), 
                                 self.minRes.get(), self.maxRes.get(), self.stepRes.get())
        self._insertFunctionStep('createOutput')

    def convertInput(self, volLocation1, volLocation2=None):
        """ Convert input volume to .mrc as expected by ResMap. 
        Params:
            volLocation1: a tuple containing index and filename of the input volume.
            volLocation2: if not None, a tuple like volLocation1 for the split volume.
        """
        ih = ImageHandler()
        outputLocation = (NO_INDEX, self._getPath('volume1.map'))
        ih.convert(volLocation1, outputLocation)
        if volLocation2 is not None:
            outputLocation = (NO_INDEX, self._getPath('volume2.map'))
            ih.convert(volLocation2, outputLocation) 
            

    def estimateResolution(self, pVal, minRes, maxRes, stepRes):
        """ Call ResMap.py with the appropiate parameters. """
        self._enterWorkingDir()
        
        inputVol = 'volume1.map'
        vxSize = self.inputVolume.get().getSamplingRate()
        
        path = os.environ['RESMAP_HOME']
        os.environ['PATH'] = path + os.pathsep + os.environ['PATH']
        os.environ['PYTHONPATH'] = path + os.pathsep + os.environ['PYTHONPATH']
        args = join(path, 'ResMap.py')
        if self.useSplitVolume:
            args += ' --noguiSplit volume1.map volume2.map '
        else:
            args += ' --noguiSingle volume1.map '
        args += '--vxSize=%(vxSize)f --pVal=%(pVal)f '
        args += '--minRes=%(minRes)f --maxRes=%(maxRes)f --stepRes=%(stepRes)f '
        args += '--vis2D' # TODO: move this to viewer
        self.runJob('xmipp_python', args % locals())
                    
        self._leaveWorkingDir()
        
    def createOutput(self):
        """ Register the output (the alignment and the aligned particles.)
        """
        pass
        
    def _summary(self):
        summary = []
        return summary
    
    def _validate(self):
        errors = []
                
        return errors
    

