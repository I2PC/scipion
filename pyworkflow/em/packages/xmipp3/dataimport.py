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

from os.path import join, basename, dirname, exists
from collections import OrderedDict

from pyworkflow.utils.path import findRootFrom, copyTree, createLink
from pyworkflow.em.data import Micrograph
from pyworkflow.em.packages.xmipp3.convert import (readSetOfMicrographs, readSetOfParticles,
                                                   xmippToLocation, locationToXmipp, 
                                                   CTF_PSD_DICT)
from pyworkflow.em.packages.xmipp3.utils import getMdFirstRow
from pyworkflow.em.packages.xmipp3 import XmippMdRow

import xmipp


class XmippImport():
    """ Class used to import different kind of objects
    from Xmipp projects into Scipion.
    """
    def __init__(self, protocol, mdFile):
        self.protocol = protocol
        self._mdFile = mdFile
        self.copyOrLink = protocol.getCopyOrLink()
    
    def importMicrographs(self):
        """ Import a SetOfMicrographs from a given micrograph metadata.
        (usually the result "micrographs.xmd" from Xmipp protocols)
        If the CTF is found, a SetOfCTF will be also created.
        """
        self._findPathAndCtf(label=xmipp.MDL_MICROGRAPH)
        micSet = self.protocol._createSetOfMicrographs()
        micSet.setObjComment('Micrographs imported from Xmipp metadata:\n%s' % self._mdFile)
        
        # Update both samplingRate and acquisition with parameters
        # selected in the protocol form
        self.protocol.setSamplingRate(micSet)
        self.protocol.fillAcquisition(micSet.getAcquisition())
        # Read the micrographs from the 'self._mdFile' metadata
        # but fixing the filenames with new ones (linked or copy to extraDir)
        readSetOfMicrographs(self._mdFile, micSet, 
                           preprocessImageRow=self._preprocessMicrographRow, 
                           readAcquisition=False)
        self.protocol._defineOutputs(outputMicrographs=micSet)
        
        # Also create a SetOfCTF if the present
        if self._ctfPath:
            ctfSet = self.protocol._createSetOfCTF()
            for mic in micSet:
                ctf = mic.getCTF()
                ctf.copyObjId(mic)
                ctfSet.append(ctf)
                
            self.protocol._defineOutputs(outputCTF=ctfSet)
            self.protocol._defineCtfRelation(micSet, ctfSet)
            
    def importParticles(self):
        """ Import particles from a metadata 'images.xmd' """
        self._imgDict = {} # store which images stack have been linked/copied and the new path
        self._findPathAndCtf(label=xmipp.MDL_IMAGE)
        if self._micIdOrName:
            # If MDL_MICROGRAPH_ID or MDL_MICROGRAPH then
            # create a set to link from particles
            self.micSet = self.protocol._createSetOfMicrographs()
            self.protocol.setSamplingRate(self.micSet)
            self.protocol.fillAcquisition(self.micSet.getAcquisition())

        partSet = self.protocol._createSetOfParticles()
        partSet.setObjComment('Particles imported from Xmipp metadata:\n%s' % self._mdFile)
        
        # Update both samplingRate and acquisition with parameters
        # selected in the protocol form
        self.protocol.setSamplingRate(partSet)
        self.protocol.fillAcquisition(partSet.getAcquisition())
        # Read the micrographs from the 'self._mdFile' metadata
        # but fixing the filenames with new ones (linked or copy to extraDir)
        readSetOfParticles(self._mdFile, partSet, 
                           preprocessImageRow=self._preprocessParticleRow, 
                           readAcquisition=False)
        if self._micIdOrName:
            self.protocol._defineOutputs(outputMicrographs=self.micSet)
        self.protocol._defineOutputs(outputParticles=partSet)
        
    def _findPathAndCtf(self, label, warnings=True):
        """ Find the relative path from which the micrographs exists
        repect to the metadata location. Also check if it contains
        CTF information and their relative root.
        """
        row = getMdFirstRow(self._mdFile)
        
        if row is None:
            raise Exception("Can not import from an empty metadata: %s" % self._mdFile)
        
        if not row.containsLabel(label):
            raise Exception("Label *%s* is missing in metadata: %s" % (xmipp.label2Str(label), 
                                                                         self._mdFile))

        # take only the filename part after the @
        index, fn = xmippToLocation(row.getValue(label))
        self._imgPath = findRootFrom(self._mdFile, fn)
        
        if warnings and self._imgPath is None:
            self.protocol.warning("Binary data was not found from metadata: %s" % self._mdFile)
        
        if row.containsLabel(xmipp.MDL_CTF_MODEL):
            self._ctfPath = findRootFrom(self._mdFile, row.getValue(xmipp.MDL_CTF_MODEL))
        else:
            self._ctfPath = None # means no CTF info from micrographs metadata

        # Check if the MetaData contains either MDL_MICROGRAPH_ID
        # or MDL_MICROGRAPH, this will be used when imported
        # particles to keep track of the particle's micrograph
        self._micIdOrName = (row.containsLabel(xmipp.MDL_MICROGRAPH_ID) or
                             row.containsLabel(xmipp.MDL_MICROGRAPH))
        #init dictionary. It will be used in the preprocessing
        self.micDict = {}
    
    def validate(self, label):
        """ Try to find errors on import. """
        errors = []
        try:
            self._findPathAndCtf(label, warnings=False)
        except Exception, ex:
            errors.append(str(ex))
            
        return errors
        
    def validateMicrographs(self):
        return self.validate(xmipp.MDL_MICROGRAPH)
    
    def validateParticles(self):
        return self.validate(xmipp.MDL_IMAGE)
        
    def _preprocessMicrographRow(self, img, imgRow):
        if self._imgPath:
            # Create a link or copy files to extraPath
            # and update the Row properly
            micFile = imgRow.getValue(xmipp.MDL_MICROGRAPH)
            micBase = basename(micFile)
            micDst = self.protocol._getExtraPath(micBase)
            self.copyOrLink(join(self._imgPath, micFile), micDst)
            imgRow.setValue(xmipp.MDL_MICROGRAPH, micDst)
            
        if self._ctfPath:
            # Read Xmipp ctfModel parameters and add
            # to the original micrograph row
            ctfFile = imgRow.getValue(xmipp.MDL_CTF_MODEL) 
            ctfPath = join(self._imgPath, ctfFile)
            ctfRow = XmippMdRow()
            ctfRow.readFromFile(ctfPath)
            imgRow.copyFromRow(ctfRow)
            # Also copy or link to the result micrograph 
            # folder output by Xmipp containing the PSD and other images
            ctfSrcDir = dirname(ctfPath)
            ctfBaseDir = basename(ctfSrcDir)
            ctfDstDir = self.protocol._getExtraPath(ctfBaseDir)
            
            if self.copyOrLink == createLink:
                createLink(ctfSrcDir, ctfDstDir)
            else: # use copyTree instead of copyFile
                copyTree(ctfSrcDir, ctfDstDir)
            # Fix the path to psd files
            for label in CTF_PSD_DICT.values():
                filePath = imgRow.getValue(label)
                # Take the last part of the path including
                # the filename and the folder up to that
                fileName = basename(filePath)
                newFilePath = join(ctfDstDir, fileName)
                imgRow.setValue(label, newFilePath)
        
    def _preprocessParticleRow(self, img, imgRow):
        if self._imgPath:
            # Create a link or copy files to extraPath
            # and update the Row properly
            imgFn = imgRow.getValue(xmipp.xmipp.MDL_IMAGE)
            index, fn = xmippToLocation(imgFn)
            imgBase = basename(fn)
            imgDst = self.protocol._getExtraPath(imgBase)
            if not exists(imgDst):
                self.copyOrLink(join(self._imgPath, fn), imgDst)
            imgRow.setValue(xmipp.MDL_IMAGE, locationToXmipp(index, imgDst))

        if self._micIdOrName:
            micId = imgRow.getValue(xmipp.MDL_MICROGRAPH_ID, None)
            micName = imgRow.getValue(xmipp.MDL_MICROGRAPH, None)

            # Check which is the key to identify micrographs (id or name)
            if micId is not None:
                micKey = micId
            else:
                micKey = micName

            mic = self.micDict.get(micKey, None)

            # First time I found this micrograph (either by id or name)
            if mic is None:
                mic = Micrograph()
                mic.setObjId(micId)
                if micName is None:
                    micName = self.protocol._getExtraPath('fake_micrograph%6d' % micId)
                mic.setFileName(micName)
                self.micSet.append(mic)
                # Update dict with new Micrograph
                self.micDict[micKey] = mic

            # Update the row to set a MDL_MICROGRAPH_ID
            imgRow.setValue(xmipp.MDL_MICROGRAPH_ID, long(mic.getObjId()))

    def loadAcquisitionInfo(self):
        """ Return a dictionary with acquisition values and 
        the sampling rate information.
        In the case of Xmipp, they are stored in files:
        acquisition_info.xmd and microscope.xmd 
        """
        acquisitionDict = OrderedDict()
        
        if exists(self._mdFile):
            dirName = dirname(self._mdFile)
            acquisitionFile = join(dirName, 'acquisition_info.xmd')
            microscopeFile = join(dirName, 'microscope.xmd')
                
            if exists(microscopeFile):
                row = getMdFirstRow(microscopeFile)
                acquisitionDict['voltage'] = row.getValue(xmipp.MDL_CTF_VOLTAGE)
                acquisitionDict['sphericalAberration'] = row.getValue(xmipp.MDL_CTF_CS)
            
            if exists(acquisitionFile):
                row = getMdFirstRow(acquisitionFile)
                acquisitionDict['samplingRate'] = row.getValue(xmipp.MDL_SAMPLINGRATE)
            
        return acquisitionDict
          
                