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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

from os.path import join, basename, dirname, exists
from collections import OrderedDict

from pyworkflow.utils.path import findRootFrom, copyTree, createLink, replaceExt
from pyworkflow.em.data import Micrograph
import pyworkflow.em.metadata as md
from pyworkflow.em.packages.xmipp3.convert import *



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
        self._findPathAndCtf(label=md.MDL_MICROGRAPH)
        micSet = self.protocol._createSetOfMicrographs()
        micSet.setObjComment('Micrographs imported from Xmipp metadata:\n'
                             '%s' % self._mdFile)

        # Update both samplingRate and acquisition with parameters
        # selected in the protocol form
        self.protocol.setSamplingRate(micSet)
        micSet.setIsPhaseFlipped(self.protocol.haveDataBeenPhaseFlipped.get())
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
        # Store which images stack have been linked/copied and the new path
        self._imgDict = {}
        # Keep a dictionary of which ctfparams have been found
        self._ctfDict = {}

        self._findPathAndCtf(label=md.MDL_IMAGE)

        if self._micIdOrName:
            # If MDL_MICROGRAPH_ID or MDL_MICROGRAPH then
            # create a set to link from particles
            self.micSet = self.protocol._createSetOfMicrographs()
            self.protocol.setSamplingRate(self.micSet)
            self.protocol.fillAcquisition(self.micSet.getAcquisition())

        partSet = self.protocol._createSetOfParticles()
        partSet.setObjComment('Particles imported from Xmipp metadata:\n'
                              '%s' % self._mdFile)

        # Update both samplingRate and acquisition with parameters
        # selected in the protocol form
        self.protocol.setSamplingRate(partSet)
        partSet.setIsPhaseFlipped(self.protocol.haveDataBeenPhaseFlipped.get())
        self.protocol.fillAcquisition(partSet.getAcquisition())
        # Read the micrographs from the 'self._mdFile' metadata
        # but fixing the filenames with new ones (linked or copy to extraDir)
        readSetOfParticles(self._mdFile, partSet,
                           preprocessImageRow=self._preprocessParticleRow,
                           readAcquisition=False)
        if self._micIdOrName:
            self.protocol._defineOutputs(outputMicrographs=self.micSet)

        self.protocol._defineOutputs(outputParticles=partSet)

        # Also create classes if MDL_REF or MDL_REF3D was found
        if self._classFunc is not None:
            clsSet = self._classFunc(partSet)
            fillClasses(clsSet)
            self.protocol._defineOutputs(outputClasses=clsSet)
            self.protocol._defineSourceRelation(partSet, clsSet)

    def _findPathAndCtf(self, label, warnings=True):
        """ Find the relative path from which the micrographs exists
        repect to the metadata location. Also check if it contains
        CTF information and their relative root.
        """
        row = md.getFirstRow(self._mdFile)

        if row is None:
            raise Exception("Can not import from an empty metadata: "
                            "%s" % self._mdFile)

        if not row.containsLabel(label):
            raise Exception("Label *%s* is missing in metadata: "
                            "%s" % (md.label2Str(label), self._mdFile))

        # take only the filename part after the @
        index, fn = xmippToLocation(row.getValue(label))
        self._imgPath = findRootFrom(self._mdFile, fn)

        if warnings and self._imgPath is None:
            self.protocol.warning("Binary data was not found from metadata: "
                                  "%s" % self._mdFile)

        if row.containsLabel(md.MDL_CTF_MODEL):
            self._ctfPath = findRootFrom(self._mdFile,
                                         row.getValue(md.MDL_CTF_MODEL))
        else:
            self._ctfPath = None # means no CTF info from micrographs metadata

        if row.containsLabel(md.MDL_REF):
            self._classFunc = self.protocol._createSetOfClasses2D
        elif row.containsLabel(md.MDL_REF3D):
            self._classFunc = self.protocol._createSetOfClasses3D
        else:
            self._classLabel = None
            self._classFunc = None

        # Check if the MetaData contains either MDL_MICROGRAPH_ID
        # or MDL_MICROGRAPH, this will be used when imported
        # particles to keep track of the particle's micrograph
        self._micIdOrName = (row.containsLabel(md.MDL_MICROGRAPH_ID) or
                             row.containsLabel(md.MDL_MICROGRAPH))
        #init dictionary. It will be used in the preprocessing
        self.micDict = {}

        return row

    def validate(self, label):
        """ Try to find errors on import. """
        errors = []
        try:
            self._findPathAndCtf(label, warnings=False)
        except Exception, ex:
            errors.append(str(ex))

        return errors

    def validateMicrographs(self):
        return self.validate(md.MDL_MICROGRAPH)

    def validateParticles(self):
        return self.validate(md.MDL_IMAGE)

    def _preprocessMicrographRow(self, img, imgRow):
        if self._imgPath:
            # Create a link or copy files to extraPath
            # and update the Row properly
            micFile = imgRow.getValue(md.MDL_MICROGRAPH)
            micBase = basename(micFile)
            micDst = self.protocol._getExtraPath(micBase)
            self.copyOrLink(join(self._imgPath, micFile), micDst)
            
            imgRow.setValue(md.MDL_MICROGRAPH, micDst)
            self._fillMicName(img, micBase)

        if self._ctfPath:
            # Read Xmipp ctfModel parameters and add
            # to the original micrograph row
            ctfFile = imgRow.getValue(md.MDL_CTF_MODEL)
            ctfPath = join(self._imgPath, ctfFile)
            ctfRow = md.Row()
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
            index, fn = xmippToLocation(imgRow.getValue(md.MDL_IMAGE))
            imgBase = basename(fn)
            imgDst = self.protocol._getExtraPath(imgBase)
            if not exists(imgDst):
                self.copyOrLink(join(self._imgPath, fn), imgDst)
            imgRow.setValue(md.MDL_IMAGE, locationToXmipp(index, imgDst))

        if self._micIdOrName:
            micId = imgRow.getValue(md.MDL_MICROGRAPH_ID, None)
            micName = imgRow.getValue(md.MDL_MICROGRAPH, None)

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
                    micName = self.protocol._getExtraPath('fake_micrograph%6d'
                                                          % micId)
                mic.setFileName(micName)
                self.micSet.append(mic)
                # Update dict with new Micrograph
                self.micDict[micKey] = mic

            # Update the row to set a MDL_MICROGRAPH_ID
            imgRow.setValue(md.MDL_MICROGRAPH_ID, long(mic.getObjId()))

        # JMRT: This means that the metadata contains MDL_CTF_MODEL
        # and the files path were found from some root
        # In Xmipp 3.1 the ctfparam metadata in particles
        # was replaced with directly seeting the CTF values
        # so we need to fill those in the particle row
        if self._ctfPath:
            ctfModel = imgRow.getValue(md.MDL_CTF_MODEL)
            if ctfModel in self._ctfDict:
                ctfRow = self._ctfDict[ctfModel]
            else:
                ctfRow = md.Row()
                ctfRow.readFromFile(join(self._ctfPath, ctfModel))
                self._ctfDict[ctfModel] = ctfRow
            imgRow.copyFromRow(ctfRow)

    def loadAcquisitionInfo(self):
        """ Return a dictionary with acquisition values and
        the sampling rate information.
        In the case of Xmipp, they are stored in files:
        acquisition_info.xmd and microscope.xmd
        """
        acqDict = OrderedDict()

        if exists(self._mdFile):
            dirName = dirname(self._mdFile)
            acquisitionFile = join(dirName, 'acquisition_info.xmd')
            microscopeFile = join(dirName, 'microscope.xmd')

            if exists(microscopeFile):
                row = md.getFirstRow(microscopeFile)
                acqDict['voltage'] = row.getValue(md.MDL_CTF_VOLTAGE)
                acqDict['sphericalAberration'] = row.getValue(md.MDL_CTF_CS)

            if exists(acquisitionFile):
                row = md.getFirstRow(acquisitionFile)
                acqDict['samplingRate'] = row.getValue(md.MDL_SAMPLINGRATE)
            
        return acqDict


    def importCoordinates(self, fileName, addCoordinate):
        posMd = readPosCoordinates(fileName)
        for objId in posMd:
            coord = rowToCoordinate(rowFromMd(posMd, objId))
            addCoordinate(coord)
            
    def getBoxSize(self, coordFile):
        """ Try to infer the box size from the given coordinate file.
        """
        configFile = join(dirname(coordFile), 'config.xmd')
        if exists(configFile):
            firstRow = md.getFirstRow('properties@' + configFile)
            return firstRow.getValue(md.MDL_PICKING_PARTICLE_SIZE)
        
        return None

    def importCTF(self, mic, fileName):
        ctf = readCTFModel(fileName, mic)
        ctf.setPsdFile(replaceExt(fileName, 'psd'))
        return ctf

    def _fillMicName(self, mic, filename):
        micName = filename.replace("/", "_")
        mic.setMicName(micName)
