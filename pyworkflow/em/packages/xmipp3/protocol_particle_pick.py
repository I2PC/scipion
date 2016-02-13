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

from pyworkflow.em import *
from pyworkflow.protocol.launch import launch
from pyworkflow.utils.path import *  
import xmipp
from xmipp3 import XmippProtocol
from pyworkflow.em.showj import launchSupervisedPickerGUI
from convert import writeSetOfMicrographs, readSetOfCoordinates


class XmippProtParticlePicking(ProtParticlePicking, XmippProtocol):
    """ Picks particles in a set of micrographs
    either manually or in a supervised mode.
    """
    _label = 'manual-picking (step 1)'

    def __init__(self, **args):        
        ProtParticlePicking.__init__(self, **args)
        # The following attribute is only for testing
        self.importFolder = String(args.get('importFolder', None))

    #--------------------------- DEFINE param functions --------------------------------------------    
    def _defineParams(self, form):
        ProtParticlePicking._defineParams(self, form)
        form.addParam('memory', FloatParam, default=2,
                   label='Memory to use (In Gb)', expertLevel=2)  
              
    #--------------------------- INSERT steps functions --------------------------------------------    
    def _insertAllSteps(self):
        """The Particle Picking process is realized for a set of micrographs"""
        # Get pointer to input micrographs
        self.inputMics = self.inputMicrographs.get()
        micFn = self.inputMics.getFileName()

        # Launch Particle Picking GUI
        if not self.importFolder.hasValue():
            self._insertFunctionStep('launchParticlePickGUIStep', micFn, interactive=True)
        else: # This is only used for test purposes
            self._insertFunctionStep('_importFromFolderStep')
            # Insert step to create output objects
            self._insertFunctionStep('createOutputStep')

    def launchParticlePickGUIStep(self, micFn):
        # Launch the particle picking GUI
        extraDir = self._getExtraPath()
        memory = '%dg'%self.memory.get()
        process = launchSupervisedPickerGUI(micFn, extraDir, self, memory=memory)
        process.wait()

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
        self._defineSourceRelation(self.inputMicrographs, coordSet)
        
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
        if self.getOutputsSize() > 0:
            return ProtParticlePicking._methods(self)
        else:
            return [self.getMethods(None)]
    
    def getMethods(self, output):#output is not used but to overwrite getMethods it is used
        configfile = join(self._getExtraPath(), 'config.xmd')
        existsConfig = exists(configfile)
        msg = ''
        
        if existsConfig:
            md = xmipp.MetaData('properties@' + configfile)
            configobj = md.firstObject()
            pickingState = md.getValue(xmipp.MDL_PICKING_STATE, configobj)
            particleSize = md.getValue(xmipp.MDL_PICKING_PARTICLE_SIZE, configobj)
            isAutopick = pickingState != "Manual"
            manualParticlesSize = md.getValue(xmipp.MDL_PICKING_MANUALPARTICLES_SIZE, configobj)
            autoParticlesSize = md.getValue(xmipp.MDL_PICKING_AUTOPARTICLES_SIZE, configobj)

            if manualParticlesSize is None:
                manualParticlesSize = 0

            if autoParticlesSize is None:
                autoParticlesSize = 0
            msg = 'User picked %d particles with a particle size of %d.' % (autoParticlesSize + manualParticlesSize, particleSize)

            if isAutopick:
                msg += "Automatic picking was used ([Abrishami2013]). %d particles were picked automatically and %d  manually."%(autoParticlesSize, manualParticlesSize)
        return msg

    def getSummary(self, coordSet):
        summary = []
        configfile = join(self._getExtraPath(), 'config.xmd')
        existsConfig = exists(configfile)
        if existsConfig:
            md = xmipp.MetaData('properties@' + configfile)
            configobj = md.firstObject()
            pickingState = md.getValue(xmipp.MDL_PICKING_STATE, configobj)
            particleSize = md.getValue(xmipp.MDL_PICKING_PARTICLE_SIZE, configobj)
            activeMic = md.getValue(xmipp.MDL_MICROGRAPH, configobj)
            isAutopick = pickingState != "Manual"
            manualParticlesSize = md.getValue(xmipp.MDL_PICKING_MANUALPARTICLES_SIZE, configobj)
            autoParticlesSize = md.getValue(xmipp.MDL_PICKING_AUTOPARTICLES_SIZE, configobj)

            summary.append("Manual particles picked: %d"%manualParticlesSize)
            summary.append("Particle size:%d" %(particleSize))
            autopick = "Yes" if isAutopick else "No"
            summary.append("Autopick: " + autopick)
            if isAutopick:
                summary.append("Automatic particles picked: %d"%autoParticlesSize)
            summary.append("Last micrograph: " + activeMic)
        return "\n".join(summary)

    
    def getCoordsDir(self):
        return self._getExtraPath()
    
    





