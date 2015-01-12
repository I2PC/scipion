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
from pyworkflow.em.data import SetOfParticles
from pyworkflow.em.constants import ALIGN_NONE



class ScipionImport():
    """ Import 
    """
    def __init__(self, protocol, sqliteFile):
        self.protocol = protocol
        self._sqliteFile = sqliteFile
        self.copyOrLink = protocol.getCopyOrLink()
    
    def importMicrographs(self):
        """ Import a SetOfMicrographs from a given micrograph metadata.
        (usually the result "micrographs.xmd" from Xmipp protocols)
        If the CTF is found, a SetOfCTF will be also created.
        """
        pass
            
    def importParticles(self):
        """ Import particles from a metadata 'images.xmd' """
        self._imgDict = {} # store which images stack have been linked/copied and the new path

        
        inputSet = self._findPathAndCtf()
        
        partSet = self.protocol._createSetOfParticles()
        partSet.copyInfo(inputSet)
        
        partSet.setObjComment('Particles imported from Scipion sqlite file:\n%s' % self._sqliteFile)
        partSet.copyItems(inputSet, updateItemCallback=self._updateParticle)
        
        # Update both samplingRate and acquisition with parameters
        # selected in the protocol form
        self.protocol.setSamplingRate(partSet)
        self.protocol.fillAcquisition(partSet.getAcquisition())

        self.protocol._defineOutputs(outputParticles=partSet)
        
    def _findPathAndCtf(self, warnings=True):
        """ Find the relative path from which the micrographs exists
        repect to the metadata location. Also check if it contains
        CTF information and their relative root.
        """
        inputSet = SetOfParticles(filename=self._sqliteFile)
        inputSet.loadProperty('_alignment', ALIGN_NONE)
        inputSet.loadProperty('_hasCtf', False)
        
        particle = inputSet.getFirstItem()
        self._imgPath = findRootFrom(self._sqliteFile, particle.getFileName())
        
        if warnings and self._imgPath is None:
            self.protocol.warning("Binary data was not found from metadata: %s" % self._sqliteFile)

        return inputSet
    
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
        In the case of Xmipp, they are stored in files:
        acquisition_info.xmd and microscope.xmd 
        """
        acquisitionDict = OrderedDict()
            
        return acquisitionDict
          
                