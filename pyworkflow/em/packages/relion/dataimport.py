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
from pyworkflow.utils.path import findRootFrom

import xmipp


class RelionImport():
    """    
    Protocol to import existing Relion runs.
    """
    def __init__(self, protocol, starFile):
        self.protocol = protocol
        self._starFile = starFile
        self.copyOrLink = protocol.getCopyOrLink()
            
    def importParticles(self):
        """ Import particles from a metadata 'images.xmd' """
        self._imgDict = {} # store which images stack have been linked/copied and the new path
        self._findImagesPath(label=xmipp.RLN_IMAGE_NAME)
        if self._micIdOrName:
            # If MDL_MICROGRAPH_ID or MDL_MICROGRAPH then
            # create a set to link from particles
            self.micSet = self.protocol._createSetOfMicrographs()
            self.protocol.setSamplingRate(self.micSet)
            self.protocol.fillAcquisition(self.micSet.getAcquisition())

        partSet = self.protocol._createSetOfParticles()
        partSet.setObjComment('Particles imported from Relion star file:\n%s' % self._starFile)
        
        # Update both samplingRate and acquisition with parameters
        # selected in the protocol form
        self.protocol.setSamplingRate(partSet)
        self.protocol.fillAcquisition(partSet.getAcquisition())
        # Read the micrographs from the 'self._starFile' metadata
        # but fixing the filenames with new ones (linked or copy to extraDir)
        from convert import readSetOfParticles
        readSetOfParticles(self._starFile, partSet, 
                           preprocessImageRow=self._preprocessImageRow, 
                           readAcquisition=False)
        if self._micIdOrName:
            self.protocol._defineOutputs(outputMicrographs=self.micSet)
        self.protocol._defineOutputs(outputParticles=partSet)
        
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
    def validateParticles(self):
        """ Should be overriden in subclasses to 
        return summary message for NORMAL EXECUTION. 
        """
        return []
    
    def summaryParticles(self):
        """ Should be overriden in subclasses to 
        return summary message for NORMAL EXECUTION. 
        """
        return []
        
    def _findImagesPath(self, label, warnings=True):
        from pyworkflow.em.packages.xmipp3.utils import getMdFirstRow

        row = getMdFirstRow(self._starFile)
        
        if row is None:
            raise Exception("Can not import from an empty metadata: %s" % self._starFile)
        
        if not row.containsLabel(label):
            raise Exception("Label *%s* is missing in metadata: %s" % (xmipp.label2Str(label), 
                                                                         self._starFile))
            
        self._imgPath = findRootFrom(self._starFile, row.getValue(label))
        
        if warnings and self._imgPath is None:
            self.protocol.warning("Binary data was not found from metadata: %s" % self._starFile)

        # Check if the MetaData contains either MDL_MICROGRAPH_ID
        # or MDL_MICROGRAPH, this will be used when imported
        # particles to keep track of the particle's micrograph
        self._micIdOrName = False#(row.containsLabel(xmipp.MDL_MICROGRAPH_ID) or
                                 #row.containsLabel(xmipp.MDL_MICROGRAPH))
        #init dictionary. It will be used in the preprocessing
        self.micDict = {}


    def _preprocessImageRow(self, img, imgRow):
        from convert import setupCTF, copyOrLinkFileName
        if self._imgPath is not None:
            copyOrLinkFileName(imgRow, self._imagesPath, self.protocol._getExtraPath())
        setupCTF(imgRow, self.protocol.samplingRate.get())
        
