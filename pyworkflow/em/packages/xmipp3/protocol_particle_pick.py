# **************************************************************************
# *
# * Authors:     Jose Gutierrez Tabuenca (jose.gutierrez@cnb.csic.es)
# *              J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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
This sub-package contains the XmippParticlePicking protocol
"""

from pyworkflow.em import *  
from pyworkflow.utils.path import *  
import xmipp
from pyworkflow.em.packages.xmipp3.data import *
from xmipp3 import XmippProtocol
from data import XmippSetOfMicrographs


class XmippDefParticlePicking(Form):
    """Create the definition of parameters for
    the XmippParticlePicking protocol"""
    def __init__(self):
        Form.__init__(self)
    
        self.addSection(label='Input')
        self.addParam('inputMicrographs', PointerParam, label="Micrographs",
                      pointerClass='SetOfMicrographs',
                      help='Select the SetOfMicrograph ')
        self.addParam('memory', FloatParam, default=2,
                   label='Memory to use (In Gb)', expertLevel=2)        


class XmippProtParticlePicking(ProtParticlePicking, XmippProtocol):
    """Protocol to pick particles manually of a set of micrographs in the project"""
    _definition = XmippDefParticlePicking()
    
    def __init__(self, **args):        
        ProtParticlePicking.__init__(self, **args)
        # The following attribute is only for testing
        self.importFolder = String(args.get('importFolder', None))
        
    def _defineSteps(self):
        """The Particle Picking proccess is realized for a set of micrographs"""
        
        # Get pointer to input micrographs 
        self.inputMics = self.inputMicrographs.get()
        # Parameters needed 
        self._params = {'memory': self.memory.get(),
                        'pickingMode': 'manual',
                        'extraDir': self._getExtraPath()
                        }
        
        # Convert input SetOfMicrographs to Xmipp if needed
        self._insertConvertStep('inputMics', XmippSetOfMicrographs, 
                                 self._getPath('micrographs.xmd'))
        # Launch Particle Picking GUI
        if not self.importFolder.hasValue():
            self._insertFunctionStep('launchParticlePickGUI', isInteractive=True)
        else: # This is only used for test purposes
            self._insertFunctionStep('_importFromFolder')       
        # Insert step to create output objects       
        self._insertFunctionStep('createOutput')
        
        
    def launchParticlePickGUI(self):
        # Get the converted input micrographs in Xmipp format
        # if not exists, means the input was already in Xmipp
        inputMicsXmipp = self.getConvertedInput('inputMics')
        self._params['inputMicsXmipp'] = inputMicsXmipp.getFileName()
        # Launch the particle picking GUI
        program = "xmipp_micrograph_particle_picking"
        arguments = "-i %(inputMicsXmipp)s -o %(extraDir)s --mode %(pickingMode)s --memory %(memory)dg"
        # TiltPairs
        if inputMicsXmipp.hasTiltPairs():
            self._params['inputMicsXmipp'] = "TiltedPairs@" + inputMicsXmipp.getFileName()
            program = "xmipp_micrograph_tiltpair_picking"
        # Run the command with formatted parameters
        self.runJob(None, program, arguments % self._params)
        
    def _createSetOfCoordinates(self, size):
        inputMicsXmipp = self.getConvertedInput('inputMics')
        # Create a md with the coordinates files for each micrograph
        coordId = 0L
        posMd = xmipp.MetaData()
        for mic in inputMicsXmipp:
            micPosFn = self._getExtraPath(replaceBaseExt(mic.getFileName(), 'pos'))
            micPosBlockFn = 'particles@' + micPosFn
            micPosMd = xmipp.MetaData(micPosBlockFn)
            #TODO:  micPosMd.fillLinear
            for objId in micPosMd:
                coordId += 1
                micPosMd.setValue(xmipp.MDL_ITEM_ID, coordId, objId)
            micPosMd.write(micPosBlockFn, xmipp.MD_APPEND)
            posId = posMd.addObject()
            posMd.setValue(xmipp.MDL_ITEM_ID, mic.getId(), posId)
            posMd.setValue(xmipp.MDL_MICROGRAPH_PARTICLES, micPosFn, posId)
        coordsFn = self._getExtraPath('scipion_micrographs_coordinates.xmd')
        posMd.write(coordsFn)  
        coords = XmippSetOfCoordinates(filename=coordsFn)
        coords.setMicrographs(inputMicsXmipp)
        coords.boxSize.set(size)
        
        return coords                    
        
    def _importFromFolder(self):
        """ This function will copy Xmipp .pos files for
        simulating an particle picking run...this is only
        for testing purposes.
        """
        from pyworkflow.utils.path import getFiles
        import shutil
        
        for f in getFiles(self.importFolder.get()):
            shutil.copy(f, self._getExtraPath())
        
    def createOutput(self):
        fn = self._getExtraPath('config.xmd')
        md = xmipp.MetaData('properties@%s' % fn)
        size = md.getValue(xmipp.MDL_PICKING_PARTICLE_SIZE, md.firstObject())
        coords = self._createSetOfCoordinates(size)
        self._defineOutputs(outputCoordinates=coords)

