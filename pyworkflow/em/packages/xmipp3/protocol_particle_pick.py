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
from xmipp3 import XmippProtocol
from pyworkflow.em.showj import runJavaIJapp

from convert import createXmippInputMicrographs, readSetOfCoordinates


class XmippProtParticlePicking(ProtParticlePicking, XmippProtocol):

    """Xmipp protocol to pick particles in a set of micrographs 
    either manually or using supervised picking support  """
    _label = 'supervised picking'

    
    def __init__(self, **args):        
        ProtParticlePicking.__init__(self, **args)
        # The following attribute is only for testing
        self.importFolder = String(args.get('importFolder', None))
    
    #--------------------------- DEFINE param functions --------------------------------------------    
    def _defineParams(self, form):
    
        form.addSection(label='Input')
        form.addParam('inputMicrographs', PointerParam, label="Micrographs",
                      pointerClass='SetOfMicrographs',
                      help='Select the SetOfMicrograph ')
        form.addParam('memory', FloatParam, default=2,
                   label='Memory to use (In Gb)', expertLevel=2)  
              
    #--------------------------- INSERT steps functions --------------------------------------------    
    def _insertAllSteps(self):
        """The Particle Picking proccess is realized for a set of micrographs"""
        
        # Get pointer to input micrographs 
        self.inputMics = self.inputMicrographs.get()
        # Parameters needed 
        
        
        # Convert input SetOfMicrographs to Xmipp if needed
#        self._insertConvertStep('inputMics', XmippSetOfMicrographs, 
#                                 self._getPath('micrographs.xmd'))
        # Launch Particle Picking GUI
        if not self.importFolder.hasValue():
            self._insertFunctionStep('launchParticlePickGUIStep', interactive=True)
        else: # This is only used for test purposes
            self._insertFunctionStep('_importFromFolderStep')       
        # Insert step to create output objects       
        self._insertFunctionStep('createOutputStep')
        
    
    #--------------------------- STEPS functions --------------------------------------------
    def launchParticlePickGUIStep(self):
        
        # Get the converted input micrographs in Xmipp format
        # if not exists, means the input was already in Xmipp
        micFn = createXmippInputMicrographs(self, self.inputMics)
        
        # Launch the particle picking GUI
        extraDir = self._getExtraPath()
        scipion =  "%s %s \"%s\" %s" % ( pw.PYTHON, pw.join('apps', 'pw_create_coords.py'), self.getDbPath(), self.strId())
        app = "xmipp.viewer.particlepicker.training.SupervisedPickerRunner"
        args = "--input %(micFn)s --output %(extraDir)s --mode manual  --scipion %(scipion)s"%locals()
        # TiltPairs
#        if self.inputMics.hasTiltPairs():
#            self._params['inputMicsXmipp'] = "TiltedPairs@" + fn
#            program = "xmipp_micrograph_tiltpair_picking"
        # Run the command with formatted parameters
        
        #self.runJob(program, arguments % self._params)
        
        runJavaIJapp("%dg" % self.memory.get(), app, args)
        
    def _importFromFolderStep(self):
        """ This function will copy Xmipp .pos files for
        simulating a particle picking run...this is only
        for testing purposes.
        """
        for f in getFiles(self.importFolder.get()):
            copyFile(f, self._getExtraPath())
        
    def createOutputStep(self):
        posDir = self._getExtraPath()
        coordSet = self._createSetOfCoordinates(self.inputMics)
        readSetOfCoordinates(posDir, self.inputMics, coordSet)
        self._defineOutputs(outputCoordinates=coordSet)
        
        self._defineSourceRelation(self.inputMics, coordSet)
        

    #--------------------------- INFO functions --------------------------------------------
    def _citations(self):
        return ['Abrishami2013']


    #--------------------------- UTILS functions --------------------------------------------
    def __str__(self):
        """ String representation of a Supervised Picking run """
        if not hasattr(self, 'outputCoordinates'):
            msg = "No particles picked yet."
        else:
            picked = self.outputCoordinates.getSize()
            msg = "Number of particles picked: %d (from %d micrographs)" % (picked, self.inputMicrographs.get().getSize())
    
        return msg
    
    def _methods(self):
        methodsMsgs = self.summary()
        #TODO: Provide summary with more details
        return methodsMsgs
    
    def _summary(self):
        summary = []
        if not hasattr(self, 'outputCoordinates'):
            summary.append("Output coordinates not ready yet.") 
        else:
            for key, output in self.iterOutputAttributes(EMObject):
                summary.append("Particles picked: %d (from %d micrographs)" % (output.getSize(), self.inputMicrographs.get().getSize()))
            # Read the picking state from the config.xmd metadata
            configfile = join(self._getExtraPath(), 'config.xmd')
            if exists(configfile):
                
                md = xmipp.MetaData('properties@' + configfile)
                configobj = md.firstObject()
                size = md.getValue(xmipp.MDL_PICKING_PARTICLE_SIZE, configobj)
                state = md.getValue(xmipp.MDL_PICKING_STATE, configobj)
                if(state is "Manual"):
                    autopick = "No"
                else:
                    autopick = "Yes"
                activemic = md.getValue(xmipp.MDL_MICROGRAPH, configobj)
                summary.append("Particle size:%d" % size)
                summary.append("Autopick: " + autopick)
                summary.append("Last micrograph: " + activemic)
        return summary
    
    def getInputMicrographs(self):
        return self.inputMicrographs.get()
    
    
    

