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

from pyworkflow.em import ProtParticlePicking, PointerParam, FloatParam, String, CoordinatesTiltPair, EMObject
from pyworkflow.utils.path import pw, getFiles, copyFile, join, exists
from xmipp3 import XmippProtocol
from pyworkflow.em.showj import runJavaIJapp
from pyworkflow.em.packages.xmipp3 import readSetOfCoordinates, readAnglesFromMicrographs
from convert import writeSetOfMicrographsPairs

import xmipp


class XmippProtParticlePickingPairs(ProtParticlePicking, XmippProtocol):

    """Xmipp protocol to pick particles in a set of micrographs 
    either manually or using supervised picking support  """
    _label = 'tilt pairs particle picking'

    
    def __init__(self, **args):        
        ProtParticlePicking.__init__(self, **args)
        # The following attribute is only for testing
        self.importFolder = String(args.get('importFolder', None))
               
    
    #--------------------------- DEFINE param functions --------------------------------------------    
    def _defineParams(self, form):
    
        form.addSection(label='Input')
        form.addParam('inputMicrographsTiltedPair', PointerParam, label="Micrographs tilt pair",
                      pointerClass='MicrographsTiltPair',
                      help='Select the MicrographsTiltPair ')
        form.addParam('memory', FloatParam, default=2,
                   label='Memory to use (In Gb)', expertLevel=2)  
              
    #--------------------------- INSERT steps functions --------------------------------------------    
    def _insertAllSteps(self):
        """The Particle Picking proccess is realized for a pair of set of micrographs"""
        
        self.micsFn = self._getPath('input_micrographs.xmd')
        # Convert input into xmipp Metadata format
        self._insertFunctionStep('convertInputStep')
        
        # Launch Particle Picking GUI
        if not self.importFolder.hasValue():
            self._insertFunctionStep('launchParticlePickGUIStep', interactive=True)
        else: # This is only used for test purposes
            self._insertFunctionStep('_importFromFolderStep')     
        
    
    #--------------------------- STEPS functions --------------------------------------------
    def convertInputStep(self):
        
        # Get the converted input micrographs in Xmipp format
        writeSetOfMicrographsPairs(self.inputMicrographsTiltedPair.get().getUntilted(), 
                                   self.inputMicrographsTiltedPair.get().getTilted(), 
                                   self.micsFn)       
        
    def launchParticlePickGUIStep(self):
        
        # Launch the particle picking GUI
        micFn = self.micsFn
        extraDir = self._getExtraPath()
        scipion =  "%s %s \"%s\" %s" % ( pw.PYTHON, pw.join('apps', 'pw_create_coords.py'), self.getDbPath(), self.strId())
        app = "xmipp.viewer.particlepicker.tiltpair.TiltPairPickerRunner"
        args = " --input %(micFn)s --output %(extraDir)s --mode manual --scipion %(scipion)s"%locals()
        
        runJavaIJapp("%dg" % self.memory.get(), app, args)
        
    def _importFromFolderStep(self):
        """ This function will copy Xmipp .pos files for
        simulating a particle picking run...this is only
        for testing purposes.
        """
        for f in getFiles(self.importFolder.get()):
            copyFile(f, self._getExtraPath())

        extradir = self._getExtraPath()  
        
        inputset = self.inputMicrographsTiltedPair.get()
        uSet = inputset.getUntilted()
        tSet = inputset.getTilted()

        # Create Untilted and Tilted SetOfCoordinates
        uCoordSet = self._createSetOfCoordinates(uSet, suffix='Untilted')
        readSetOfCoordinates(extradir, uSet, uCoordSet)
        uCoordSet.write()
        tCoordSet = self._createSetOfCoordinates(tSet, suffix='Tilted')
        readSetOfCoordinates(extradir, tSet, tCoordSet)
        tCoordSet.write()
        
        # Read Angles from faked input micrographs
        micsFn = self._getExtraPath('input_micrographs.xmd')
        setAngles = self._createSetOfAngles()
        readAnglesFromMicrographs(micsFn, setAngles)
        setAngles.write()
        # Create CoordinatesTiltPair object
        outputset = CoordinatesTiltPair()
        outputset.setTilted(tCoordSet)
        outputset.setUntilted(uCoordSet)
        outputset.setAngles(setAngles)
        outputset.setMicsPair(inputset)
        
        self._defineOutputs(outputCoordinatesTiltPair=outputset)
        self._defineSourceRelation(inputset, outputset)
        
    #--------------------------- INFO functions --------------------------------------------
    def _citations(self):
        return ['Abrishami2013']


    #--------------------------- UTILS functions --------------------------------------------
    def __str__(self):
        """ String representation of a Particle Picking Tilt run """
        if not hasattr(self, 'outputCoordinatesTiltPair'):
            msg = "No particles picked yet."
        else:
            picked = self.outputCoordinatesTiltPair.getTilted().getSize()
            msg = "Number of particles picked: %d (from %d micrographs)" % (picked, self.inputMicrographsTiltedPair.get().getTilted().getSize())
        return msg

    


    def getMethods(self, output):#output is not used but to overwrite getMethods it is used

        configfile = join(self._getExtraPath(), 'config.xmd')
        existsConfig = exists(configfile)
        if existsConfig:
            md = xmipp.MetaData('properties@' + configfile)
            configobj = md.firstObject()
            particleSize = md.getValue(xmipp.MDL_PICKING_PARTICLE_SIZE, configobj)
            manualParticlesSize = md.getValue(xmipp.MDL_PICKING_MANUALPARTICLES_SIZE, configobj)
            msg = 'User picked %d particles with a particle size of %d.' % (manualParticlesSize, particleSize)

        return msg

    def _summary(self):
        if self.getOutputsSize() > 0:
            return ProtParticlePicking._summary(self)
        else:
            return self.getSummary(self.getCoords())

    def _methods(self):
        if self.getOutputsSize() > 0:
            return ProtParticlePicking._methods(self)
        else:
            return [self.getMethods(None)]

    def getInputMicrographs(self):
        return self.inputMicrographsTiltedPair.get().getTilted()

    def getCoords(self):
        count = self.getOutputsSize()
        suffix = str(count) if count > 1 else ''
        outputName = 'outputCoordinatesTiltPair' + suffix
        return getattr(self, outputName)
    
    def getSummary(self):
        summary = []
        configfile = join(self._getExtraPath(), 'config.xmd')
        if exists(configfile):
            md = xmipp.MetaData('properties@' + configfile)
            configobj = md.firstObject()
            activemic = md.getValue(xmipp.MDL_MICROGRAPH, configobj)
            manualParticlesSize = md.getValue(xmipp.MDL_PICKING_MANUALPARTICLES_SIZE, configobj)
            summary.append("Particles picked: %d"%manualParticlesSize)
            summary.append("Last micrograph: " + activemic)
        return "\n".join(summary)
       
    

