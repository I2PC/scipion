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

from os.path import join, basename, exists
from collections import OrderedDict

from pyworkflow.utils.path import findRootFrom
from pyworkflow.em.data import SetOfParticles, SetOfMicrographs
from pyworkflow.em.constants import ALIGN_NONE



class ScipionImport():
    """ Import 
    """
    def __init__(self, protocol, sqliteFile):
        self.protocol = protocol
        self._sqliteFile = sqliteFile
        self.copyOrLink = protocol.getCopyOrLink()
    
    def importMicrographs(self):
        """ Import a SetOfMicrographs from a given sqlite file. """
        inputSet = SetOfMicrographs(filename=self._sqliteFile)
        self._findImagesPath(inputSet)
        
        micSet = self.protocol._createSetOfMicrographs()
        micSet.setObjComment('Micrographs imported from sqlite file:\n%s'
                             % self._sqliteFile)

        # Update both samplingRate and acquisition with parameters
        # selected in the protocol form.
        # Acquisition should be updated before, just to ensure that
        # scannedPixedSize will be computed properly when calling
        # setSamplingRate
        self.protocol.fillAcquisition(micSet.getAcquisition())
        self.protocol.setSamplingRate(micSet)
        micSet.setIsPhaseFlipped(self.protocol.haveDataBeenPhaseFlipped.get())
        # Read the micrographs from the 'self._sqliteFile' metadata
        # but fixing the filenames with new ones (linked or copy to extraDir)
        micSet.copyItems(inputSet, updateItemCallback=self._updateParticle)
        self.protocol._defineOutputs(outputMicrographs=micSet)

    def importParticles(self):
        """ Import particles from a metadata 'images.xmd' """
        inputSet = SetOfParticles(filename=self._sqliteFile)
        inputSet.loadProperty('_alignment', ALIGN_NONE)
        inputSet.loadProperty('_hasCtf', False)
        self._findImagesPath(inputSet)

        partSet = self.protocol._createSetOfParticles()
        partSet.copyInfo(inputSet)
        
        partSet.setObjComment('Particles imported from Scipion sqlite file:\n%s'
                              % self._sqliteFile)
        partSet.copyItems(inputSet, updateItemCallback=self._updateParticle)
        
        # Update both samplingRate and acquisition with parameters
        # selected in the protocol form
        self.protocol.setSamplingRate(partSet)
        self.protocol.fillAcquisition(partSet.getAcquisition())
        partSet.setIsPhaseFlipped(self.protocol.haveDataBeenPhaseFlipped.get())
        self.protocol._defineOutputs(outputParticles=partSet)
        
    def _findImagesPath(self, inputSet):
        """ Find the relative path from which the images exists
        repect to the sqlite location. 
        """
        # Store which images stack have been linked/copied and the new path
        self._imgDict = {}
        img = inputSet.getFirstItem()
        self._imgPath = findRootFrom(self._sqliteFile, img.getFileName())
        
        if self._imgPath is None:
            self.protocol.warning("Binary data was not found from sqlite: %s"
                                  % self._sqliteFile)
    
    def validate(self):
        """ Try to find errors on import. """
        errors = []
        return errors
        
    def validateMicrographs(self):
        return self.validate()
    
    def validateParticles(self):
        return self.validate()
        
    def _updateParticle(self, particle, partRow):
        if self._imgPath:
            # Create a link or copy files to extraPath
            # and update the Row properly
            fn = particle.getFileName()
            imgBase = basename(fn)
            imgDst = self.protocol._getExtraPath(imgBase)
            if not exists(imgDst):
                self.copyOrLink(join(self._imgPath, fn), imgDst)
            particle.setFileName(imgDst)

    def loadAcquisitionInfo(self):
        """ Return a dictionary with acquisition values and 
        the sampling rate information.
        In the case of Scipion, these values will be read from the
        'Properties' table of the particles.sqlite file.
        """
        acq = OrderedDict()
        inputSet = SetOfParticles(filename=self._sqliteFile)
        def _get(key):
            return inputSet.getProperty(key)
        acq['samplingRate'] = _get('_samplingRate')
        acq['voltage'] = _get('_acquisition._voltage')
        acq['amplitudeContrast'] = _get('_acquisition._amplitudeContrast')
        acq['sphericalAberration'] = float(_get('_acquisition._sphericalAberration'))

        return acq
