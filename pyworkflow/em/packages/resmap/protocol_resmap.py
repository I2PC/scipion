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
Protocol wrapper around the ResMap tool for local resolutio
"""

import os
import sys

from pyworkflow.object import String
from pyworkflow.protocol.params import PointerParam, BooleanParam, FloatParam
from pyworkflow.em.protocol import ProtAnalysis3D
from pyworkflow.em.convert import ImageHandler



class ProtResMap(ProtAnalysis3D):
    """
    ResMap is software tool for computing the local resolution of 3D
    density maps studied in structural biology, primarily by cryo-electron
    microscopy (cryo-EM).
     
    Please find the manual at http://resmap.sourceforge.net 
    """
    _label = 'local resolution'
    
    def __init__(self, **kwargs):
        ProtAnalysis3D.__init__(self, **kwargs)
        self.histogramData = String()
             
    #--------------------------- DEFINE param functions --------------------------------------------   
    
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputVolume', PointerParam, pointerClass='Volume',  
                      label="Input volume", important=True,
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
                      pointerClass='VolumeMask', condition="applyMask",
                      help='Select a volume to apply as a mask.')
        
        line = form.addLine('Pre-whitening')
        line.addParam('prewhitenAng', FloatParam, label="Angstroms")
        line.addParam('prewhitenRamp', FloatParam, label='Ramp')                      
        
        group = form.addGroup('Extra parameters')
        #form.addSection(label='Optional')
        group.addParam('stepRes', FloatParam, default=1,
                      label='Step size (Ang):',
                      help='in Angstroms (min 0.25, default 1.0)')
        line = group.addLine('Resolution Range (A)', 
                            help='Default: algorithm will start a just above 2*voxelSize')
        line.addParam('minRes', FloatParam, default=0, label='Min')
        line.addParam('maxRes', FloatParam, default=0, label='Max')
        group.addParam('pVal', FloatParam, default=0.05,
                      label='Confidence level:',
                      help='usually between [0.01, 0.05]')   
                     
    #--------------------------- INSERT steps functions --------------------------------------------  
    
    def _insertAllSteps(self):
        # Insert processing steps
        inputs = [self.inputVolume.get().getLocation()]
        if self.useSplitVolume:
            inputs.append(self.splitVolume.get().getLocation())
            
        self._insertFunctionStep('convertInputStep', *inputs)
        self._insertFunctionStep('estimateResolutionStep', 
                                 self.pVal.get(), 
                                 self.minRes.get(), 
                                 self.maxRes.get(), 
                                 self.stepRes.get(),
                                 self.prewhitenAng.get(),
                                 self.prewhitenRamp.get())

    #--------------------------- STEPS functions --------------------------------------------       
    
    def convertInputStep(self, volLocation1, volLocation2=None):
        """ Convert input volume to .mrc as expected by ResMap. 
        Params:
            volLocation1: a tuple containing index and filename of the input volume.
            volLocation2: if not None, a tuple like volLocation1 for the split volume.
        """
        ih = ImageHandler()
        ih.convert(volLocation1, self._getPath('volume1.map'))
        if volLocation2 is not None:
            ih.convert(volLocation2, self._getPath('volume2.map')) 
            

    def estimateResolutionStep(self, pValue, minRes, maxRes, stepRes, ang, rampWeight):
        """ Call ResMap.py with the appropiate parameters. """
        results = self.runResmap(self._getPath())
        
        from cPickle import dumps
        self.histogramData.set(dumps(results['resHisto']))
        self._store(self.histogramData)
        
    #--------------------------- INFO functions -------------------------------------------- 
    
    def _summary(self):
        summary = []
        return summary
    
    def _validate(self):
        errors = []
                
        return errors
    
    #--------------------------- UTILS functions --------------------------------------------
 
    def runResmap(self, workingDir, wizardMode=False):
        """ Prepare the args dictionary to be used
        and call the ResMap algorithm.
        Params:
            workingDir: where to run ResMap
            wizardMode: some custom params to be used by the wizard
                to display the pre-whitening GUI and only that.
        with the  """
        self._enterDir(workingDir)
        
        volumes = ['volume1.map', 'volume2.map']
        
        # Add resmap libraries to the path
        sys.path.append(os.environ['RESMAP_HOME'])
        from ResMap_algorithm import ResMap_algorithm
        from ResMap_fileIO import MRC_Data
        
        # Always read the first volume as mrc data
        data1 = MRC_Data(volumes[0],'ccp4')
        
        prewhitenArgs = {'display': wizardMode,
                         'force-stop': wizardMode                         
                         }
        if (self.prewhitenAng.hasValue() and 
            self.prewhitenRamp.hasValue()):
            prewhitenArgs['newElbowAngstrom'] = self.prewhitenAng.get()
            prewhitenArgs['newRampWeight'] = self.prewhitenRamp.get()
            
        args = {'vxSize': self.inputVolume.get().getSamplingRate(),
                'pValue': self.pVal.get(),
                'minRes': self.minRes.get(), 
                'maxRes': self.maxRes.get(),
                'stepRes': self.stepRes.get(),
                'chimeraLaunch': False, # prevent ResMap to launch some graphical analysis
                'graphicalOutput': False, 
                'scipionPrewhitenParams': prewhitenArgs
                }
        
        if self.useSplitVolume:
            # Read the second splitted volume
            data2 = MRC_Data(volumes[0],'ccp4')
            args.update({'inputFileName1': 'volume1.map',
                         'inputFileName2': 'volume2.map',
                         'data1': data1,
                         'data2': data2,
                         })
        else:
            args.update({'inputFileName': 'volume1.map',
                         'data': data1,
                         })   
            
        results = ResMap_algorithm(**args)  
        self._leaveDir()
        
        return results   
    
     

