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

import os

from pyworkflow.protocol.params import FileParam, FloatParam, BooleanParam, IntParam
from pyworkflow.em.protocol import ProtImport
from pyworkflow.utils.properties import Message

from protocol_base import ProtRelionBase



class ProtRelionImport(ProtImport, ProtRelionBase):
    """    
    Protocol to import existing Relion runs.
    """
    _label = 'import'
    
    #--------------------------- DEFINE param functions --------------------------------------------   
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputStar', FileParam, 
                      label="Input data STAR file",  
                      help='Select the input data STAR file from a Relion run.'
                           'Also the *optimiser.star and *sampling.star files '
                           'should be present.')  
        form.addParam('samplingRate', FloatParam, 
                      label=Message.LABEL_SAMP_RATE,
                      help='Provide the sampling rate of your particles. (in Angstroms per pixel)')
        form.addParam('magnification', IntParam, default=50000,
           label=Message.LABEL_MAGNI_RATE)
        form.addParam('isPhaseFlipped', BooleanParam, default=False,
                      label='The particles are phase flipped?')
        form.addParam('hasAlignment', BooleanParam, default=False,
                      label='The particles has alignment?')
        
        form.addParam('importClasses', BooleanParam, default=True,
                      label='Import also classes?')
        
        #TODO: Add an option to allow the user to decide if copy binary files or not        
            
    #--------------------------- INSERT steps functions --------------------------------------------  
    def _insertAllSteps(self): 
        self._insertFunctionStep('createOutputStep', self.inputStar.get())
    
    #--------------------------- STEPS functions --------------------------------------------
    def createOutputStep(self, dataFile):
        partSet = self._createParticles(dataFile)
        self._defineOutputs(outputParticles=partSet)
        firstParticle = partSet.getFirstItem()
        if firstParticle.getMicId() is None:
            if firstParticle.hasAttribute("_micrograph"):
                #create micID: aggregate function, argument function, group by
                micIdList = partSet.aggregate(['count'],'_micrograph',['_micrograph'])
                micIdMap={}
                counter = 0;
                for mic in micIdList:
                    micIdMap[mic['_micrograph']]=counter
                    counter = counter +1
                for par in partSet:
                    par.setMicID(micIdMap[mic['_micrograph']])
    #               MDL_MICROGRAPH:        'rlnMicrographName'

            self.warning("Micrograph ID was not set for particles!!!")
        if not firstParticle.hasAlignment():
            self.warning("Alignment was not read from particles!!!")
        
        if self.importClasses:
            classes = self._createClasses(dataFile, partSet)
            self._defineOutputs(outputClasses=classes)
            self._defineSourceRelation(partSet, classes)
        # TODO: Register input volume and also clases if necesary
        
    def _createParticles(self, dataFile):
        self.info('Creating the set of particles...')
        from convert import readSetOfParticles
        # Create the set of particles
        partSet = self._createSetOfParticles()
        self._findImagesPath(dataFile)
        # Copy acquisition from first element
        partSet.setSamplingRate(self.samplingRate.get())
        readSetOfParticles(dataFile, partSet, 
                           preprocessImageRow=self._preprocessImageRow, 
                           magnification=self.magnification.get())
        particle = partSet.getFirstItem()
        partSet.setAcquisition(particle.getAcquisition())
        
        
        if self.isPhaseFlipped:
            partSet.setIsPhaseFlipped(True)

        if self.hasAlignment:
            #FIXME: Detect 2D or 3D alignment
            partSet.setAlignment3D()
        return partSet   
    
    def _createClasses(self, dataFile, partSet):     
        self.info('Creating the set of classes...')
        from convert import readSetOfClasses3D, createClassesFromImages
        # Create the set of classes 2D or 3D  
        classesSqlite = self._getTmpPath('classes.sqlite')
        relDataFile = os.path.relpath(dataFile)
        classTemplate = relDataFile.replace('_data.star', '_class%(ref)03d.mrc:mrc')
        self.info('  Using classes template: %s' % classTemplate)
        createClassesFromImages(partSet, dataFile, classesSqlite, 
                                self.OUTPUT_TYPE, self.CLASS_LABEL, classTemplate, 
                                0, preprocessImageRow=self._preprocessImageRow)      
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
    
    def copyOfLinkReferences(self, refPattern):
        foundRef = True
        ref = 0
        while foundRef:
            refPath = refPattern % ref
            foundRef = os.path.exists(refPath)
        
    def _findImagesPath(self, starFile):
        from convert import findImagesPath
        self._imagesPath = findImagesPath(starFile)
        
        if self._imagesPath:
            self.info('Images path found in: %s' % self._imagesPath)
        else:
            self.warning('Images binaries not found!!!')
            

    def _preprocessImageRow(self, img, imgRow):
        from convert import setupCTF, copyOrLinkFileName
        if self._imagesPath is not None:
            copyOrLinkFileName(imgRow, self._imagesPath, self._getExtraPath())
        setupCTF(imgRow, self.samplingRate.get())
#     def _preprocessRow(self, img, imgRow):
#         from convert import setupCTF, prependToFileName
#         if self._imagesPath is not None:
#             prependToFileName(imgRow, self._imagesPath)
#         
#         setupCTF(imgRow, self.samplingRate.get())
        
