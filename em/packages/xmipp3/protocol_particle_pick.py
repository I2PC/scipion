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
from pyworkflow.utils.process import runJob
from xmipp import MetaData, MDL_MICROGRAPH, MDL_MICROGRAPH_ORIGINAL, MDL_MICROGRAPH_TILTED, MDL_MICROGRAPH_TILTED_ORIGINAL, MDL_PICKING_FAMILY, MDL_PICKING_PARTICLE_SIZE
from pyworkflow.em.packages.xmipp3.data import *
from convert import convertSetOfMicrographs
from xmipp3 import XmippProtocol


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
        #self._insertFunctionStep('convertToXmippSetOfMicrograph')
        self._insertFunctionStep('convertInputToXmipp', 'inputMics', 
                                 convertSetOfMicrographs, self._getPath('micrographs.xmd'))
        # Launch Particle Picking GUI
        self._insertFunctionStep('launchParticlePickGUI', isInteractive=True)       
        # Insert step to create output objects       
        self._insertFunctionStep('createOutput')
        
    def convertToXmippSetOfMicrograph(self):
        """ We need to ensure the micrograph.xmd metadata is available
        to Xmipp picking program before launching it.
        """
        # Convert from SetOfMicrographs to XmippSetOfMicrographs
        micFn = self._getPath('micrographs.xmd')
        inputMicsXmipp = convertSetOfMicrographs(self.inputMics, micFn)

        if inputMicsXmipp != self.inputMics:  # If distintic, micrographs.xmd should be produced
            self._insertChild('inputMicsXmipp', inputMicsXmipp)
            return [micFn]
        
    def launchParticlePickGUI(self):
        # Get the converted input micrographs in Xmipp format
        # if not exists, means the input was already in Xmipp
        inputMicsXmipp = getattr(self, 'inputMicsXmipp', self.inputMics)
        self._params['inputMicsXmipp'] = self.getConvertedInput('inputMics').getFileName()
        # Launch the particle picking GUI
        program = "xmipp_micrograph_particle_picking"
        arguments = "-i %(inputMicsXmipp)s -o %(extraDir)s --mode %(pickingMode)s --memory %(memory)dg"
        # TiltPairs
        if inputMicsXmipp.hasTiltPairs():
            program = "xmipp_micrograph_tiltpair_picking"
        # Run the command with formatted parameters
        runJob(None, program, arguments % self._params)
        
    def _createSetOfCoordinates(self, family, size):
        inputMicsXmipp = getattr(self, 'inputMicsXmipp', self.inputMics)
        print "createSetOfCoordinates for family: ", family, size
        coords = XmippSetOfCoordinates(filename=self._getExtraPath())
        coords.setMicrographs(inputMicsXmipp)
        coords.family.set(family)
        coords.boxSize.set(size)
        return coords                    
        
    def createOutput(self):
        fn = self._getExtraPath('families.xmd')
        md = MetaData(fn)
        for objId in md:
            family = md.getValue(MDL_PICKING_FAMILY, objId)
            size = md.getValue(MDL_PICKING_PARTICLE_SIZE, objId)
            coords = self._createSetOfCoordinates(family, size)
            self._defineOutputs(outputCoordinates=coords)

