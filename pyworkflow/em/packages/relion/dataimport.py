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

import os
from os.path import exists
from collections import OrderedDict

from pyworkflow.object import Float
from pyworkflow.em.constants import ALIGN_PROJ, ALIGN_2D, ALIGN_NONE
from pyworkflow.em.data import Micrograph
import pyworkflow.em.metadata as md
from pyworkflow.em.packages.relion.convert import relionToLocation
from pyworkflow.utils.path import findRootFrom

# import xmipp



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
        self.ignoreIds = self.protocol.ignoreIdColumn.get()
        self._imgDict = {} # store which images stack have been linked/copied and the new path
        self._findImagesPath(label=md.RLN_IMAGE_NAME)
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
                           postprocessImageRow=self._postprocessImageRow,
                           readAcquisition=False, alignType=self.alignType)
        if self._micIdOrName:
            self.protocol._defineOutputs(outputMicrographs=self.micSet)
        self.protocol._defineOutputs(outputParticles=partSet)

        if self._classesFunc is not None:
            self._createClasses(partSet)


    def _updateClass(self, item):
        classId = item.getObjId()
        if  classId in self._classesDict:
            index, fn, row = self._classesDict[classId]
            if fn.endswith('.mrc'):
                fn += ':mrc' # Specify that are volumes to read them properly in xmipp
            item.getRepresentative().setLocation(index, fn)
            item._rlnclassDistribution = Float(row.getValue('rlnClassDistribution'))
            item._rlnAccuracyRotations = Float(row.getValue('rlnAccuracyRotations'))
            item._rlnAccuracyTranslations = Float(row.getValue('rlnAccuracyTranslations'))

    def _createClasses(self, partSet):
        self._classesDict = {} # store classes info, indexed by class id
        pathDict = {}

        self.protocol.info('Loading classes info from: %s' % self._modelStarFile)
        modelMd = md.MetaData('model_classes@' + self._modelStarFile)
        for classNumber, objId in enumerate(modelMd):
            row = md.Row()
            row.readFromMd(modelMd, objId)
            index, fn = relionToLocation(row.getValue('rlnReferenceImage'))

            if fn in pathDict:
                newFn = pathDict.get(fn)
            else:
                clsPath = findRootFrom(self._modelStarFile, fn)
                if clsPath is None:
                    newFn = fn
                else:
                    newFn = self.protocol._getExtraPath(os.path.basename(fn))
                    self.copyOrLink(os.path.join(clsPath, fn), newFn)
                pathDict[fn] = newFn

            self._classesDict[classNumber+1] = (index, newFn, row)

        clsSet = self._classesFunc(partSet)
        clsSet.classifyItems(updateClassCallback=self._updateClass)

        self.protocol._defineOutputs(outputClasses=clsSet)
        self.protocol._defineSourceRelation(partSet, clsSet)

    #--------------------------- INFO functions -------------------------------------------- 
    def validateParticles(self):
        """ Should be overriden in subclasses to 
        return summary message for NORMAL EXECUTION. 
        """
        errors = []
        try:
            self._findImagesPath(label="rlnImageName", warnings=False)
        except Exception, ex:
            errors.append(str(ex))

        return errors

    def summaryParticles(self):
        """ Should be overriden in subclasses to 
        return summary message for NORMAL EXECUTION. 
        """
        return []

    def _findImagesPath(self, label, warnings=True):

        row = md.getFirstRow(self._starFile)

        if row is None:
            raise Exception("Can not import from an empty metadata: %s" % self._starFile)

        if not row.containsLabel(label):
            raise Exception("Label *%s* is missing in metadata: %s" % (md.label2Str(label),
                                                                       self._starFile))

        index, fn = relionToLocation(row.getValue(label))
        self._imgPath = findRootFrom(self._starFile, fn)

        if warnings and self._imgPath is None:
            self.protocol.warning("Binary data was not found from metadata: %s" % self._starFile)


        if self._starFile.endswith('_data.star'):

            modelStarFile = self._starFile.replace('_data.star', '_model.star')

            if exists(modelStarFile):
                self._modelStarFile = modelStarFile
            else:
                modelHalfStarFile = self._starFile.replace('_data.star', '_half1_model.star')
                if exists(modelHalfStarFile):
                    self._modelStarFile = modelHalfStarFile
                else:
                    raise Exception("Missing required model star file, search for\n%s\nor\n%s" % (modelStarFile,
                                                                                                  modelHalfStarFile))

            modelRow = md.getFirstRow(self._modelStarFile)
            classDimensionality = modelRow.getValue('rlnReferenceDimensionality')

            self._optimiserFile = self._starFile.replace('_data.star', '_optimiser.star')
            if not exists(self._optimiserFile):
                raise Exception("Missing required optimiser star file: %s" % self._optimiserFile)
            optimiserRow = md.getFirstRow(self._optimiserFile)
            autoRefine = optimiserRow.containsLabel('rlnModelStarFile2')


            self.alignType = ALIGN_PROJ

            if not autoRefine:
                if classDimensionality == 3:
                    self._classesFunc = self.protocol._createSetOfClasses3D
                else:
                    self._classesFunc = self.protocol._createSetOfClasses2D
                    self.alignType = ALIGN_2D
            else:
                self._classesFunc = None
        else:
            self.alignType = ALIGN_NONE
            self._classesFunc = None
            self._modelStarFile = None
        # Check if the MetaData contains either MDL_MICROGRAPH_ID
        # or MDL_MICROGRAPH, this will be used when imported
        # particles to keep track of the particle's micrograph
        self._micIdOrName = (row.containsLabel('rlnMicrographName') or
                             row.containsLabel('rlnMicrographId'))
        #init dictionary. It will be used in the preprocessing
        self.micDict = {}

        return row, modelRow


    def _preprocessImageRow(self, img, imgRow):
        from convert import setupCTF, copyOrLinkFileName
        if self._imgPath is not None:
            copyOrLinkFileName(imgRow, self._imgPath, self.protocol._getExtraPath())
        setupCTF(imgRow, self.protocol.samplingRate.get())

        if self._micIdOrName:
            micId = imgRow.getValue('rlnMicrographId', None)
            micName = imgRow.getValue('rlnMicrographName', None)

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
            imgRow.setValue('rlnMicrographId', long(mic.getObjId()))

    def _postprocessImageRow(self, img, imgRow):
        if self.ignoreIds:
            img.setObjId(None) # Force to generate a new id in Set

    def loadAcquisitionInfo(self):
        """ Return a dictionary with acquisition values and 
        the sampling rate information.
        In the case of Xmipp, they are stored in files:
        acquisition_info.xmd and microscope.xmd 
        """
        acquisitionDict = OrderedDict()

        try:
            row, modelRow = self._findImagesPath(label=md.RLN_IMAGE_NAME, warnings=False)

            if row.containsLabel(md.RLN_CTF_VOLTAGE):
                acquisitionDict['voltage'] = row.getValue(md.RLN_CTF_VOLTAGE)

            if row.containsLabel('rlnAmplitudeContrast'):
                acquisitionDict['amplitudeContrast'] = row.getValue('rlnAmplitudeContrast')

            if row.containsLabel('rlnSphericalAberration'):
                acquisitionDict['sphericalAberration'] = row.getValue('rlnSphericalAberration')

            if modelRow.containsLabel('rlnPixelSize'):
                acquisitionDict['samplingRate'] = modelRow.getValue('rlnPixelSize')

        except Exception, ex:
            print "Error loading acquisition: ", str(ex)

        return acquisitionDict

    def importCoordinates(self, mic, fileName, coordSet):
        print 'importFromRelion'  + fileName

            