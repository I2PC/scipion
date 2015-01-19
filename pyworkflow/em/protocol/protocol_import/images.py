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
In this module are protocol base classes related to EM imports of Micrographs, Particles, Volumes...
"""

import sys
from os.path import basename

from pyworkflow.utils.properties import Message
from pyworkflow.protocol.params import FloatParam, IntParam, LabelParam
from pyworkflow.em.convert import ImageHandler
from pyworkflow.em.data import Acquisition

from base import ProtImportFiles



class ProtImportImages(ProtImportFiles):
    """Common protocol to import a set of images into the project"""
    # This label should be set in subclasses
    _label = 'None'
    # The following class property should be set in each import subclass
    # for example, if set to SetOfParticles, this will the output classes
    # It is also assumed that a function with the name _createSetOfParticles
    # exists in the EMProtocol base class
    _outputClassName = 'None'  
    # If set to True, each binary file will be inspected to
    # see if it is a binary stack containing more items
    _checkStacks = True
        
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        ProtImportFiles._defineParams(self, form)
        self._defineAcquisitionParams(form)
        
    def _defineAcquisitionParams(self, form):
        """ Define acquisition parameters, it can be overriden
        by subclasses to change what parameters to include.
        """
        group = form.addGroup('Acquisition info')
        group.addParam('acquisitionWizard', LabelParam, important=True,
                       condition='importFrom != %d' % self.IMPORT_FROM_FILES,
                       label='Use the wizard button to import acquisition.',
                       help='Depending on the import Format, the wizard\n'
                            'will try to import the acquisition values.\n'
                            'If not found, required ones should be provided.')
        group.addParam('voltage', FloatParam, default=200,
                   label=Message.LABEL_VOLTAGE, 
                   help=Message.TEXT_VOLTAGE)
        group.addParam('sphericalAberration', FloatParam, default=2,
                   label=Message.LABEL_SPH_ABERRATION, 
                   help=Message.TEXT_SPH_ABERRATION)
        group.addParam('amplitudeContrast', FloatParam, default=0.1,
                      label=Message.LABEL_AMPLITUDE,
                      help=Message.TEXT_AMPLITUDE)
        group.addParam('magnification', IntParam, default=50000,
                   label=Message.LABEL_MAGNI_RATE, 
                   help=Message.TEXT_MAGNI_RATE)
        return group
        
        
    #--------------------------- INSERT functions ---------------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('importImagesStep', self.getPattern(), 
                                 self.voltage.get(), self.sphericalAberration.get(), 
                                 self.amplitudeContrast.get(), self.magnification.get()) #, self.samplingRate.get(),
        
    #--------------------------- STEPS functions ---------------------------------------------------
    def importImagesStep(self, pattern, voltage, sphericalAberration, amplitudeContrast, magnification):
        """ Copy images matching the filename pattern
        Register other parameters.
        """
        self.info("Using pattern: '%s'" % pattern)
        
        createSetFunc = getattr(self, '_create' + self._outputClassName)
        imgSet = createSetFunc()
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
            
            if self._checkStacks:
                _, _, _, n = imgh.getDimensions(dst)
                
            if n > 1:
                for index in range(1, n+1):
                    img.cleanObjId()
                    img.setMicId(fileId)
                    img.setFileName(dst)
                    img.setIndex(index)
                    imgSet.append(img)
            else:
                img.setObjId(fileId)
                img.setFileName(dst)
                imgSet.append(img)
            outFiles.append(dst)
            
            sys.stdout.write("\rImported %d/%d" % (i+1, self.numberOfFiles))
            sys.stdout.flush()
            
        print "\n"
        
        args = {}
        outputSet = self._getOutputName()
        args[outputSet] = imgSet
        self._defineOutputs(**args)
        
        return outFiles
    
    #--------------------------- INFO functions ----------------------------------------------------
    
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
            summary.append("*%d* %s imported from %s" % (outputSet.getSize(), self._getOutputItemName(), self.getPattern()))
            summary.append("Sampling rate : *%0.2f* A/px" % outputSet.getSamplingRate())
        
        return summary
    
    def _methods(self):
        methods = []
        outputSet = self._getOutputSet()
        if outputSet is not None:
            methods.append("*%d* %s were imported" % (outputSet.getSize(), self._getOutputItemName())+\
                           " with a sampling rate of *%0.2f* A/px (microscope voltage %d kV, magnification %dx). Output set is %s."%
                            (outputSet.getSamplingRate(),round(self.voltage.get()),round(self.magnification.get()), self.getObjectTag(outputSet)))
            
        return methods
    
    #--------------------------- UTILS functions ---------------------------------------------------
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
        """ Override to import acquisition from the specified format. """
        return None 


