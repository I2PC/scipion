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

from os.path import join, basename

from pyworkflow.utils.path import findRootFrom
from pyworkflow.em.packages.xmipp3.convert import readSetOfMicrographs
from pyworkflow.em.packages.xmipp3.utils import getMdFirstRow

import xmipp


class XmippImport():
    """ Class used to import different kind of objects
    from Xmipp projects into Scipion.
    """
    def __init__(self, protocol):
        self.protocol = protocol
        self.copyOrLink = protocol.getCopyOrLink()
    
    def importMicrographs(self, micsMd):
        """ Import a SetOfMicrographs from a given micrograph metadata.
        (usually the result "micrographs.xmd" from Xmipp protocols)
        If the CTF is found, a SetOfCTF will be also created.
        Params:
            protocol: the protocol that will be used to register
                the outputs.
            micsMd: the metadata with micrographs.
        """
        self._findPathAndCtf(micsMd)
        micSet = self.protocol._createSetOfMicrographs()
        micSet.setObjComment('Micrographs imported from Xmipp metadata:\n%s' % micsMd)
        # Update both samplingRate and acquisition with parameters
        # selected in the protocol form
        self.protocol.setSamplingRate(micSet)
        self.protocol.fillAcquisition(micSet.getAcquisition())
        # Read the micrographs from the 'micsMd' metadata
        # but fixing the filenames with new ones (linked or copy to extraDir)
        readSetOfMicrographs(micsMd, micSet, 
                           preprocessImageRow=self._preprocessImageRow, 
                           readAcquisition=False)
        self.protocol._defineOutputs(outputMicrographs=micSet)
        
    def _findPathAndCtf(self, micsMd):
        """ Find the relative path from which the micrographs exists
        repect to the metadata location. Also check if it contains
        CTF information and their relative root.
        """
        row = getMdFirstRow(micsMd)
        
        if row is None:
            raise Exception("Can not import Micrographs from an empty metadata: %s" % micsMd)
        
        if not row.containsLabel(xmipp.MDL_MICROGRAPH):
            raise Exception("Label 'micrograph' is missing from metadata: %s" % micsMd)
            
        self._micsPath = findRootFrom(micsMd, row.getValue(xmipp.MDL_MICROGRAPH))
        if self._micsPath is None:
            self.protocol.warning("Micrographs binary data was not found from metadata: %s" % micsMd)
        self._hasCtf = row.containsLabel(xmipp.MDL_CTF_MODEL)
    
        
    def _preprocessImageRow(self, img, imgRow):
        if self._micsPath:
            # Create a link or copy files to extraPath
            # and update the Row properly
            micFile = imgRow.getValue(xmipp.MDL_MICROGRAPH)
            micBase = basename(micFile)
            micDst = self.protocol._getExtraPath(micBase)
            self.copyOrLink(join(self._micsPath, micFile), micDst)
            imgRow.setValue(xmipp.MDL_MICROGRAPH, micDst)
    
    def _postprocessImageRow(self, img, imgRow):
        pass  
                