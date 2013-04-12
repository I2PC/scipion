# **************************************************************************
# *
# * Authors:     José Gutiérrez Tabuenca (jose.gutierrez@cnb.csic.es)
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
from pyworkflow.utils import *  
from xmipp import MetaData, MDL_MICROGRAPH, MDL_MICROGRAPH_ORIGINAL, MDL_MICROGRAPH_TILTED, MDL_MICROGRAPH_TILTED_ORIGINAL
from pyworkflow.em.packages.xmipp3.data import *


class XmippDefParticlePicking(Form):
    """Create the definition of parameters for
    the XmippParticlePicking protocol"""
    def __init__(self):
        Form.__init__(self)
    
        self.addSection(label='Input')
        self.addParam('inputMicrographs', PointerParam, label="Micrographs", pointerClass='SetOfMicrographs')
        self.addParam('memory', FloatParam, default=2,
                   label='Memory to use (In Gb)', expertLevel=2)        

class XmippProtParticlePicking(ProtParticlePicking):
    """Protocol to pick particles manually of a set of micrographs in the project"""
    _definition = XmippDefParticlePicking()
    
    def __init__(self, **args):
        
        Protocol.__init__(self, **args)
        
    def defineSteps(self):
        '''The Particle Picking proccess is realized for a set of micrographs'''
        
        # Get pointer to input micrographs 
        self.inputMics = self.inputMicrographs.get()
        
        # Convert from SetOfMicrographs to SetOfMicrographsXmipp
        self.inputMicsXmipp = self.inputMics
        pass
        
        # Parameters needed 
        self.params = {'inputMicsXmipp': self.inputMicsXmipp, #pillarlo del micrxmipp.metadata.filename
                       'memory': self.memory.get(),
                       'pickingMode': 'PM_MANUAL',
                       }
        
        # For a list of micrograph is launched the Particle Picking GUI
        self.launchParticlePickGUI()

        # Insert step to create output objects       
        self.insertFunctionStep('createOutput')
        
    def launchParticlePickGUI(self):
        # (log, InputMicrographs, ExtraDir, PickingMode=PM_MANUAL,TiltPairs=False, Memory=2, Family=""):
        
        arguments = "-i %(inputMicsXmipp)s -o %(ExtraDir)s --mode %(PickingMode)s --memory %(Memory)dg"
        program = "xmipp_micrograph_particle_picking"
                     
        # TiltPairs
        if self.params['tiltPairs'] != "":
            program = "xmipp_micrograph_tiltpair_picking"
               
        # Insert the command with the formatted parameters
        self.insertRunJobStep(program, arguments % self.params)
        
    def createOutput(self):

        #Para cada familia crear un SetOfCoordinate.

        mdOut = self.getPath("micrographs.xmd")
                
        self.outputMicrographs = XmippSetOfMicrographs(value=mdOut)     
        self.outputMicrographs.microscope.voltage.set(self.inputMics.microscope.voltage.get())
        self.outputMicrographs.microscope.sphericalAberration.set(self.inputMics.microscope.sphericalAberration.get())
        self.outputMicrographs.samplingRate.set(self.inputMics.samplingRate.get())
        
        # Create the xmipp metadata micrographs.xmd  
         
        md = MetaData()      
        for i, v in IOTable.iteritems():
            objId = md.addObject()
            md.setValue(MDL_MICROGRAPH, v, objId)
            md.setValue(MDL_MICROGRAPH_ORIGINAL, i, objId)
            # TODO: Handle Tilted micrographs
#            if tiltPairs:
#                MD.setValue(xmipp.MDL_MICROGRAPH_TILTED,IOTable[fnMicrographTilted],objId)
#                MD.setValue(xmipp.MDL_MICROGRAPH_TILTED_ORIGINAL,fnMicrographTilted,objId)
        md.write("micrographs" + "@" + mdOut)
        
        self.outputMicrographs.setFileName(mdOut)

        self.defineOutputs(micrograph=self.outputMicrographs)
