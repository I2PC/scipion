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
from convert import readSetOfCoordinates




class BsoftProtParticlePicking(ProtParticlePicking):

    """Bsoft protocol to pick particles in a set of micrographs of the project"""
    _label = 'particle picking'

    
    def __init__(self, **args):        
        ProtParticlePicking.__init__(self, **args)
        # The following attribute is only for testing

        
    def _defineParams(self, form):
    
        form.addSection(label='Input')
        form.addParam('inputMicrographs', PointerParam, label="Micrographs",
                      pointerClass='SetOfMicrographs',
                      help='Select the SetOfMicrograph ')
        form.addParam('memory', FloatParam, default=2,
                   label='Memory to use (In Gb)', expertLevel=2)    
        

        
    def _insertAllSteps(self):
        """The Particle Picking proccess is realized for a set of micrographs"""
        
        # Get pointer to input micrographs 
        self.inputMics = self.inputMicrographs.get()
        # Launch Particle Picking GUI
        self._insertFunctionStep('launchParticlePickGUIStep', isInteractive=True)
        # Insert step to create output objects       
        self._insertFunctionStep('createOutputStep')
        
        
    def launchParticlePickGUIStep(self):
        self._enterWorkingDir()
        # Launch the particle picking GUI
        program = "bshow"  
       
        for mic in self.inputMics:
            micfile = os.path.relpath(mic.getFileName(), self.workingDir.get())
            print micfile
            self.runJob(program, micfile)
        self._leaveWorkingDir()
   
        
    def createOutputStep(self):
        outputDir = self._getExtraPath()
        coordSet = self._createSetOfCoordinates()
        coordSet.setMicrographs(self.inputMics)
        readSetOfCoordinates(outputDir, self.inputMics, coordSet)
        self._defineOutputs(outputCoordinates=coordSet)        
        self._defineSourceRelation(self.inputMics, coordSet)
        

    def __str__(self):
        """ String representation of a Supervised Picking run """
    
        if not hasattr(self, 'outputCoordinates'):
            picked = 0
        else:
            picked = self.outputCoordinates.getSize()
        return  "Number of particles picked: %d (from %d micrographs)" % (picked, self.inputMicrographs.get().getSize())
    
  

