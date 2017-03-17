# -*- coding: utf-8 -*-
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

import sys
import os
from os.path import basename, exists, isdir, join, dirname
import time
from datetime import timedelta, datetime

import pyworkflow.utils as pwutils
from pyworkflow.utils.properties import Message
import pyworkflow.protocol.params as params
from pyworkflow.em.convert import ImageHandler
from pyworkflow.em.data import Acquisition

from base import ProtImportFiles



class ProtImportImages(ProtImportFiles):
    """ Common protocol to import a set of images into the project. """
    # The following class property should be set in each import subclass
    # for example, if set to SetOfParticles, this will the output classes
    # It is also assumed that a function with the name _createSetOfParticles
    # exists in the EMProtocol base class
    _outputClassName = 'None'  
    # If set to True, each binary file will be inspected to
    # see if it is a binary stack containing more items
    _checkStacks = True
        
    #--------------------------- DEFINE param functions ------------------------

    def _defineAcquisitionParams(self, form):
        """ Define acquisition parameters, it can be overriden
        by subclasses to change what parameters to include.
        """
        group = form.addGroup('Acquisition info')
        group.addParam('haveDataBeenPhaseFlipped', params.BooleanParam,
                       default=False,
                      label='Have data been phase-flipped?',
                      help='Set this to Yes if the images have been ctf-phase '
                           'corrected.')
        group.addParam('acquisitionWizard', params.LabelParam, important=True,
                       condition=self._acquisitionWizardCondition(),
                       label='Use the wizard button to import acquisition.',
                       help='Depending on the import Format, the wizard\n'
                            'will try to import the acquisition values.\n'
                            'If not found, required ones should be provided.')
        group.addParam('voltage', params.FloatParam, default=200,
                   label=Message.LABEL_VOLTAGE, 
                   help=Message.TEXT_VOLTAGE)
        group.addParam('sphericalAberration', params.FloatParam, default=2,
                   label=Message.LABEL_SPH_ABERRATION, 
                   help=Message.TEXT_SPH_ABERRATION)
        group.addParam('amplitudeContrast', params.FloatParam, default=0.1,
                      label=Message.LABEL_AMPLITUDE,
                      help=Message.TEXT_AMPLITUDE)
        group.addParam('magnification', params.IntParam, default=50000,
                   label=Message.LABEL_MAGNI_RATE, 
                   help=Message.TEXT_MAGNI_RATE)
        return group

    def _acquisitionWizardCondition(self):
        """ By default this wizard will appears only when we import from
        a format that is not from files.
        But movie-import also can have a wizard to read from FEI xml files. """
        return 'importFrom != %d' % self.IMPORT_FROM_FILES
    
    #--------------------------- INSERT functions ------------------------------
    def _insertAllSteps(self):
        
        if self.dataStreaming:
            funcName = 'importImagesStreamStep' 
        else:
            funcName = 'importImagesStep'
        
        self._insertFunctionStep(funcName, self.getPattern(),
                                 self.voltage.get(),
                                 self.sphericalAberration.get(),
                                 self.amplitudeContrast.get(),
                                 self.magnification.get())
        
    #--------------------------- STEPS functions -------------------------------
    def importImagesStep(self, pattern, voltage, sphericalAberration, 
                         amplitudeContrast, magnification):
        """ Copy images matching the filename pattern
        Register other parameters.
        """
        self.info("Using pattern: '%s'" % pattern)
        
        createSetFunc = getattr(self, '_create' + self._outputClassName)
        imgSet = createSetFunc()
        imgSet.setIsPhaseFlipped(self.haveDataBeenPhaseFlipped.get())
        acquisition = imgSet.getAcquisition()
        
        self.fillAcquisition(acquisition)
        
        # Call a function that should be implemented by each subclass
        self.setSamplingRate(imgSet)
        
        outFiles = [imgSet.getFileName()]
        imgh = ImageHandler()
        img = imgSet.ITEM_TYPE()
        img.setAcquisition(acquisition)
        n = 1
        copyOrLink = self.getCopyOrLink()

        for i, (fileName, fileId) in enumerate(self.iterFiles()):
            dst = self._getExtraPath(basename(fileName))
            copyOrLink(fileName, dst)
            # Handle special case of Imagic images, copying also .img or .hed
            self.handleImgHed(copyOrLink, fileName, dst)
            
            if self._checkStacks:
                _, _, _, n = imgh.getDimensions(dst)
                
            if n > 1:
                for index in range(1, n+1):
                    img.cleanObjId()
                    img.setMicId(fileId)
                    img.setFileName(dst)
                    img.setIndex(index)
                    self._addImageToSet(img, imgSet)
            else:
                img.setObjId(fileId)
                img.setFileName(dst)
                # Fill the micName if img is a Micrograph.
                self._fillMicName(img, fileName)
                self._addImageToSet(img, imgSet)

            outFiles.append(dst)
            
            sys.stdout.write("\rImported %d/%d" % (i+1, self.numberOfFiles))
            sys.stdout.flush()
            
        print "\n"
        
        args = {}
        outputSet = self._getOutputName()
        args[outputSet] = imgSet
        self._defineOutputs(**args)
        
        return outFiles

    def _addImageToSet(self, img, imgSet):
        """ Add an image to a set, check if the first image to handle
         some special condition to read dimensions such as .txt or compressed
         movie files.
        """
        if imgSet.isEmpty():
            self._setupFirstImage(img, imgSet)
        imgSet.append(img)

    def _setupFirstImage(self, img, imgSet):
        pass

    def iterNewInputFiles(self):
        """ Iterate over input files that have not been imported.
        This function uses the self.importedFiles dict.
        """
        for fileName, fileId in self.iterFiles():
            # If file already imported, skip it
            if fileName not in self.importedFiles:
                yield fileName, fileId

    def importImagesStreamStep(self, pattern, voltage, sphericalAberration, 
                         amplitudeContrast, magnification):
        """ Copy images matching the filename pattern
        Register other parameters.
        """
        self.info("Using pattern: '%s'" % pattern)
        createSetFunc = getattr(self, '_create' + self._outputClassName)
        imgSet = createSetFunc()
        imgSet.setIsPhaseFlipped(self.haveDataBeenPhaseFlipped.get())
        acquisition = imgSet.getAcquisition()
        self.fillAcquisition(acquisition)
        # Call a function that should be implemented by each subclass
        self.setSamplingRate(imgSet)
        outFiles = [imgSet.getFileName()]
        imgh = ImageHandler()
        img = imgSet.ITEM_TYPE()
        img.setAcquisition(acquisition)
        n = 1
        copyOrLink = self.getCopyOrLink()
        outputName = self._getOutputName()

        finished = False
        self.importedFiles = set()
        # this is only used when creating stacks from frame files
        self.createdStacks = set()

        i = 0
        lastDetectedChange = datetime.now()
        timeout = timedelta(seconds=self.timeout.get())
        fileTimeout = timedelta(seconds=self.fileTimeout.get())

        while not finished:
            time.sleep(3) # wait 3 seconds before check for new files
            someNew = False
            someAdded = False

            for fileName, fileId in self.iterNewInputFiles():
                someNew = True

                if self.fileModified(fileName, fileTimeout):
                    continue

                self.info('Importing file: %s' % fileName)
                self.importedFiles.add(fileName)
                dst = self._getExtraPath(basename(fileName))
                copyOrLink(fileName, dst)

                if self._checkStacks:
                    _, _, _, n = imgh.getDimensions(dst)

                someAdded = True
                self.debug('Appending file to DB...')
                if self.importedFiles: # enable append after first append
                    imgSet.enableAppend()

                if n > 1:
                    for index in range(1, n+1):
                        img.cleanObjId()
                        img.setMicId(fileId)
                        img.setFileName(dst)
                        img.setIndex(index)
                        self._addImageToSet(img, imgSet)
                else:
                    img.setObjId(fileId)
                    img.setFileName(dst)
                    # Fill the micName if img is a Micrograph.
                    self._fillMicName(img, fileName)
                    self._addImageToSet(img, imgSet)

                outFiles.append(dst)
                self.debug('After append. Files: %d' % len(outFiles))

            if someAdded:
                self.debug('Updating output...')
                self._updateOutputSet(outputName, imgSet,
                                      state=imgSet.STREAM_OPEN)
                self.debug('Update Done.')

            self.debug('Checking if finished...someNew: %s' % someNew)

            now = datetime.now()

            if not someNew:
                # If there are no new detected files, we should check the
                # inactivity time elapsed (from last event to now) and
                # if it is greater than the defined timeout, we conclude
                # the import and close the output set
                # Another option is to check if the protocol have some
                # special stop condition, this can be used to manually stop
                # some protocols such as import movies
                finished = (now - lastDetectedChange > timeout or
                            self.streamingHasFinished())
                self.debug("Checking if finished:")
                self.debug("   Now - Last Change: %s"
                           % pwutils.prettyDelta(now - lastDetectedChange))
                self.debug("Finished: %s" % finished)
            else:
                # If we have detected some files, we should update
                # the timestamp of the last event
                lastDetectedChange = now

        self._updateOutputSet(outputName, imgSet,
                              state=imgSet.STREAM_CLOSED)
        return outFiles
    
    #--------------------------- INFO functions --------------------------------
    def _validateImages(self):
        errors = []
        ih = ImageHandler()
        
        for imgFn, _ in self.iterFiles():
            
            if isdir(imgFn):
                errors.append("Folders can not be selected.")
                errors.append('  %s' % imgFn)
            else:
                # Check if images are correct by reading the header of the
                # imported files:
                # Exceptions: 
                #  - Compressed movies (bz2 or tbz extensions)
                #  - Importing in streaming, since files may be incomplete
                if (not self.dataStreaming and 
                    not (imgFn.endswith('bz2') or 
                         imgFn.endswith('tbz') or 
                         ih.isImageFile(imgFn))): 
                    if not errors: # if empty add the first line
                        errors.append("Error reading the following images:")
                    errors.append('  %s' % imgFn)
        
        return errors
        
    def _validate(self):
        errors = ProtImportFiles._validate(self)
        # Check that files are proper EM images, only when importing from
        # files and not using streamming. In the later case we could
        # have partial files not completed.
        if (self.importFrom == self.IMPORT_FROM_FILES and
            not self.dataStreaming):
            errors += self._validateImages()
        
        return errors
        
    def _summary(self):
        summary = []
        outputSet = self._getOutputSet()
        if outputSet is None:
            summary.append("Output " + self._outputClassName + " not ready yet.") 
            if self.copyFiles:
                summary.append("*Warning*: You select to copy files into your project.\n"
                               "This will make another copy of your data and may take \n"
                               "more time to import. ")
        else:
            summary.append("*%d* %s imported from %s" % (outputSet.getSize(),
                                                         self._getOutputItemName(),
                                                         self.getPattern()))
            summary.append("Is the data phase flipped : %s" % outputSet.isPhaseFlipped())
            summary.append("Sampling rate : *%0.2f* Å/px" % outputSet.getSamplingRate())
        
        return summary
    
    def _methods(self):
        methods = []
        outputSet = self._getOutputSet()
        if outputSet is not None:
            methods.append("*%d* %s were imported with a sampling rate of "
                           "*%0.2f* Å/px (microscope voltage %d kV, "
                           "magnification %dx). Output set is %s."
                           % (outputSet.getSize(), self._getOutputItemName(),
                              outputSet.getSamplingRate(),
                              round(self.voltage.get()),
                              round(self.magnification.get()),
                              self.getObjectTag(self._getOutputName())))
            
        return methods
    
    #--------------------------- UTILS functions -------------------------------
    def getFiles(self):
        outputSet = self._getOutputSet()
        if outputSet is not None:
            return self._getOutputSet().getFiles()
        else:
            return []
    
    def _getOutputName(self):
        # We assume that the import output is always a 'SetOfSomething'
        return self._outputClassName.replace('SetOf', 'output')
    
    def _getOutputItemName(self):
        return self._outputClassName.replace('SetOf', '')
    
    def _getOutputSet(self):
        return getattr(self, self._getOutputName(), None)
    
    def fillAcquisition(self, acquisition):
        """ Fill the acquition object with protocol params. """
        acquisition.setVoltage(self.voltage.get())
        acquisition.setSphericalAberration(self.sphericalAberration.get())
        acquisition.setAmplitudeContrast(self.amplitudeContrast.get())
        acquisition.setMagnification(self.magnification.get())
        
    def getAcquisition(self):
        """ Build and fill an acquisition object. """
        acquisition = Acquisition()
        self.fillAcquisition(acquisition)
        
        return acquisition   
    
    def loadAcquisitionInfo(self):
        """ Return a proper acquistionInfo (dict)
        or an error message (str).
        """
        result = "Could not import acquisition info."
        # Check if the user selected to import from a given package
        # the property self.importFilePath should be set
        ci = self.getImportClass()

        if self.hasAttribute('importFilePath'):
            if not self.importFilePath:
                result = "Select import file first"
            elif not exists(self.importFilePath):
                result = ("*Import file doesn't exist.*\n"
                         "File:\n%s" % self.importFilePath)
            else:

                acq = ci.loadAcquisitionInfo()
                if acq is not None:
                    result = acq

        return result

    def _fillMicName(self, img, filename):
        from pyworkflow.em import Micrograph
        if isinstance(img, Micrograph):
            filePaths = self.getMatchFiles()
            commPath = pwutils.commonPath(filePaths)
            micName = filename.replace(commPath + "/", "").replace("/", "_")
            img.setMicName(micName)
            
    def handleImgHed(self, copyOrLink, src, dst):
        """ Check the special case of Imagic files format
        composed by two files: .hed and .img.
        When copying or linking we need to take care
        of the other one.
        """
        if src.endswith('.hed'):
            src2 = src.replace('.hed', '.img')
            dst2 = dst.replace('.hed', '.img')
        elif src.endswith('.img'):
            src2 = src.replace('.img', '.hed')
            dst2 = dst.replace('.img', '.hed')
        else:
            src2 = None
            dst2 = None
        
        if src2 is not None and exists(src2):
            copyOrLink(src2, dst2)
    
    def _onNewFile(self, newFile):
        """ This method will be called we a new files is found in
        streaming mode. """
        pass

    def _createOutputSet(self):
        """ Create the output set that will be populated as more data is
        imported. """

    # --------------- Streaming special functions -----------------------
    def _getStopStreamingFilename(self):
        return self._getExtraPath("STOP_STREAMING.TXT")

    def getActions(self):
        """ This method will allow that the 'Stop import' action to appears
        in the GUI when the user right-click in the protocol import box.
        It will allow a user to manually stop the streaming.
        """
        # Only allow to stop if running and in streaming mode
        if self.dataStreaming and self.isRunning():
            return [('STOP STREAMING', self.stopImport)]
        else:
            return []

    def stopImport(self):
        """ Since the actual protocol that is running is in a different
        process that the one that this method will be invoked from the GUI,
        we will use a simple mechanism to place an special file to stop
        the streaming.
        """
        # Just place an special file into the run folder
        f = open(self._getStopStreamingFilename(), 'w')
        f.close()

    def streamingHasFinished(self):
        return os.path.exists(self._getStopStreamingFilename())
