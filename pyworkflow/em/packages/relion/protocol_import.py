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
This module contains the protocol for 3d classification with relion.
"""

from pyworkflow.protocol.params import FileParam
from pyworkflow.em.protocol import ProtImport

from protocol_base import *



class ProtRelionImport(ProtImport, ProtRelionBase):
    """    
    Protocol to import existing Relion runs.
    """
    _label = 'import relion'
    
    CHANGE_LABELS = [xmipp.MDL_AVG_CHANGES_ORIENTATIONS, 
                     xmipp.MDL_AVG_CHANGES_OFFSETS]
    
    def __init__(self, **args):        
        ProtRelionBase.__init__(self, **args)
        
    def _initialize(self):
        """ This function is mean to be called after the 
        working dir for the protocol have been set. (maybe after recovery from mapper)
        """
        ProtRelionBase._initialize(self)

    #--------------------------- DEFINE param functions --------------------------------------------   
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputStar', FileParam, 
                      label="Input data STAR file",  
                      help='Select the input data STAR file from a Relion run.'
                           'Also the *optimiser.star and *sampling.star files '
                           'should be present.')  
            
    #--------------------------- INSERT steps functions --------------------------------------------  
    def _insertAllSteps(self): 
        self._insertFunctionStep('createOutputStep', self.inputStar.get())
    
    #--------------------------- STEPS functions --------------------------------------------
    def createOutputStep(self, dataFile):
        
        # TODO: Register set of Micrographs
        partSet = self._createParticles(dataFile)
        self._defineOutputs(outputParticles=partSet)
        
        classes = self._createClasses(dataFile, partSet)
        self._defineOutputs(outputClasses=classes)
        self._defineSourceRelation(partSet, classes)
        # TODO: Register input volume and also clases if necesary
        
    def _createParticles(self, dataFile):
        self.info('Creating the set of particles...')
        from convert import readSetOfParticles
        # Create the set of particles
        partSet = self._createSetOfParticles()
        readSetOfParticles(dataFile, partSet)
        particle = partSet.getFirstItem() 
        # Copy acquisition from first element
        partSet.setAcquisition(particle.getAcquisition())
        
        # Grap the sampling rate from the --angpix option
        optimiserFile = dataFile.replace('_data', '_optimiser')
        opts = self._parseCommand(optimiserFile)
        partSet.setSamplingRate(float(opts['--angpix']))
        
        return partSet   
    
    def _createClasses(self, dataFile, partSet):     
        self.info('Creating the set of classes...')
        from convert import readSetOfClasses3D
        # Create the set of classes 2D or 3D  
        classesSqlite = self._getTmpPath('classes.sqlite')
        classTemplate = dataFile.replace('_data.star', '_class%(ref)03d.mrc:mrc')
        createClassesFromImages(partSet, dataFile, classesSqlite, 
                                self.OUTPUT_TYPE, self.CLASS_LABEL, classTemplate, 0)      
        # FIXME: Check whether create classes 2D or 3D
        classes = self._createSetOfClasses3D(partSet)
        readSetOfClasses3D(classes, classesSqlite)
        
        return classes
    
    #--------------------------- INFO functions -------------------------------------------- 
    def _validate(self):
        """ Should be overriden in subclasses to 
        return summary message for NORMAL EXECUTION. 
        """
        return []
    
    def _summary(self):
        """ Should be overriden in subclasses to 
        return summary message for NORMAL EXECUTION. 
        """
        return []
    
    #--------------------------- UTILS functions --------------------------------------------
    
    def _parseCommand(self, optimiserFile):
        """ Read from the optimiser.star the second line which should 
        contain the exact command line used to launch relion and 
        grap the parameters in a dictionary way. """
        opts = {}
        self.info('Parsing parameters from optimiser file: %s' % optimiserFile)
        f = open(optimiserFile)
        for line in f:
            if '--angpix' in line:
                parts = line.split()
                for i, p in enumerate(parts):
                    if p.startswith('--'): # it an option
                        opts[p] = parts[i+1] # take what follows the option
                break
        f.close()
        return opts