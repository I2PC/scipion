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
In this module are two protocols to Import/Export data from/to EMX.
"""
import emxlib

from pyworkflow.em.protocol import *
from pyworkflow.protocol.params import *
from pyworkflow.em.data import Acquisition



class ProtEmxImport(ProtImport):
    """
    Import micrographs, coordinates or particles from EMX file.
    
    EMX is a joint initiative for data exchange format between different 
    EM software packages. See more about [[http://i2pc.cnb.csic.es/emx][EMX format]]
    """
    _label = 'emx import'
    
    # Mapping betwween form parameters and EMX tags
    PARAM_DICT = {'voltage': 'acceleratingVoltage',
                  'sphericalAberration': 'cs',
                  'amplitudeContrast': 'amplitudeContrast',
                  'samplingRate': 'pixelSpacing__X',                  
                  }
        
    #--------------------------- DEFINE param functions --------------------------------------------   
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputEMX', FileParam, 
                      label="Input EMX file",
                      help='Provide the path to a valid EMX file.')
        form.addParam('doCopyFiles', BooleanParam, default=False,
                      label="Copy files?",
                      help='Copy files into project folder. '
                           'VERY useful if you plan to move '
                           'your project to another computer '
                           'and the data is not inside the working directory')
        
        
        form.addParam('provideAcquisition', BooleanParam, default=False,
                      label='Provide acquision info?',
                      help='You can provide acquisition infomation such as:\n'
                           'Voltage, spherical aberration, sampling rate, etc...\n'
                           'and override the values from the EMX file.\n\n'
                           'If such parameters are not present in the EMX file\n'
                           'their values are mandatory and should be provided\n'
                           'in the form entries.')
        group = form.addGroup('Acquisition info', condition='provideAcquisition')
        group.addParam('voltage', FloatParam, allowsNull=True,
                   label=Message.LABEL_VOLTAGE)
        group.addParam('sphericalAberration', FloatParam, allowsNull=True,
                   label=Message.LABEL_SPH_ABERRATION)
        group.addParam('amplitudeContrast', FloatParam, allowsNull=True,
                      label=Message.LABEL_AMPLITUDE,
                      help=Message.TEXT_AMPLITUDE)
        group.addParam('samplingRate', FloatParam, allowsNull=True,  
                   label=Message.LABEL_SAMP_RATE)
        group.addParam('magnification', FloatParam, default=10000,  
                   label=Message.LABEL_MAGNI_RATE)
        
    def _loadEmxInfo(self):
        """ Load the EMX file and get some information about the type
        of objects contained and the binary data.
        """
        self.emxFn = self.inputEMX.get()
        emxData = emxlib.EmxData()
        
        emxDir = dirname(self.emxFn)
        self.classElement = None
        self.binaryFile = None
        self.object = None
        self.objDict = {}
        
        for classElement in emxlib.CLASSLIST:
            obj = emxData.readFirstObject(classElement, self.emxFn)
            if obj is not None:
                self.objDict[classElement] = obj
                #is the binary file of this type
                binaryFile = join(emxDir, obj.get(emxlib.FILENAME))
                
                if exists(binaryFile):
                    self.object = obj
                    self.binaryFile = binaryFile
                    self.classElement = classElement
                    
                for k, v in self.PARAM_DICT.iteritems():
                    # Read from the EMX files the values not set in the Form
                    if not getattr(self, k).hasValue():
                        if obj.get(v) is not None:
                            self.setAttributeValue(k, obj.get(v))

    #--------------------------- INSERT steps functions --------------------------------------------  
    def _insertAllSteps(self):
        self._insertFunctionStep('importDataStep', self.inputEMX.get())       

    #--------------------------- STEPS functions --------------------------------------------       
    def importDataStep(self, emxFile):
        """ Export micrographs to EMX file.
        micsId is only passed to force redone of this step if micrographs change.
        """
        acquisition = Acquisition()
        # Setting Acquisition properties
        acquisition.setMagnification(self.magnification.get())
        acquisition.setVoltage(self.voltage.get())
        acquisition.setSphericalAberration(self.sphericalAberration.get())
        acquisition.setAmplitudeContrast(self.amplitudeContrast.get())
        
        from convert import importData
        importData(self, emxFile, self._getPath(), 
                   acquisition, self.samplingRate.get(),self.doCopyFiles.get())
    
    #--------------------------- INFO functions -------------------------------------------- 
    def _validate(self):
        errors = []
        self.emxFn = self.inputEMX.get()
        # Check that input EMX file exist
        if not exists(self.emxFn):
                errors.append("Input EMX file *%s* doesn't exists" % self.emxFn)
        else:
            self._loadEmxInfo()   
            
            if self.object is None:
                errors.append('Cannot find any object in EMX file *%s*' % self.emxFn)
            else:
                if self.binaryFile is None:
                    errors.append('Cannot find binary data *%s* associated with EMX metadata file.\n' % self.binaryFile)
                
                for k, v in self.PARAM_DICT.iteritems():
                    if not getattr(self, k).hasValue():
                        errors.append('*%s* param was left empty and *%s* does not have attribute *%s*.\n' % (k, self.classElement, v))
        
        return errors
    
    def _citations(self):
        cites = []
        return cites
    
    def _summary(self):
        summary = []
        return summary
    
    def _methods(self):
        return self._summary()  # summary is quite explicit and serve as methods
    
    
class ProtEmxExport(EMProtocol):
    """
    Export micrographs, coordinates or particles to EMX format.
    
    EMX is a joint initiative for data exchange format between different 
    EM software packages. See more about [[http://i2pc.cnb.csic.es/emx][EMX format]]
    """
    _label = 'emx export'
            
    #--------------------------- DEFINE param functions --------------------------------------------   
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputSet', PointerParam, 
                      pointerClass='SetOfMicrographs,SetOfCoordinates,SetOfParticles', 
                      label="Set to export",
                      help='Select the microgrpahs, coordinates or particles set to be exported to EMX.')
        form.addParam('ctfEstimation', RelationParam, 
                      allowsNull=True, relationName=RELATION_CTF, attributeName='getInputSet', 
                      label='Include CTF from', 
                      help='You can select a CTF estimation associated with these\n'
                           'micrographs to be included in the EMX file')
        
        form.addParam('outputPrefix', StringParam, default='data',
                      label='EMX files prefix',
                      help='Select how do you want to name the EMX files.'
                           'For example, if you use "data" as prefix, two'
                           'files will be generated:\n'
                           '_data.emx_ and _data.mrc_')
                 
    def getInputSet(self):
        return self.inputSet.get()
    
    #--------------------------- INSERT steps functions --------------------------------------------  
    def _insertAllSteps(self):
        self._insertFunctionStep('exportDataStep', self.inputSet.get().getObjId())       

    #--------------------------- STEPS functions --------------------------------------------       
    def exportDataStep(self, micsId):
        """ Export micrographs to EMX file.
        micsId is only passed to force redone of this step if micrographs change.
        """
        from convert import exportData
        emxDir = self._getPath('emxData')
        xmlFile = self.outputPrefix.get() + '.emx'
        binaryFile = self.outputPrefix.get() + '.mrc'
        exportData(emxDir, self.inputSet.get(), ctfSet=self.ctfEstimation.get(), 
                   xmlFile=xmlFile, binaryFile=binaryFile)
        
        self._defineOutputs(emxOutput=EMXObject(join(emxDir, xmlFile), 
                                                join(emxDir, binaryFile)))
    
    #--------------------------- INFO functions -------------------------------------------- 
    def _validate(self):
        errors = []
        return errors
    
    def _citations(self):
        cites = []
        return cites
    
    def _summary(self):
        summary = []
        return summary
    
    def _methods(self):
        return self._summary()  # summary is quite explicit and serve as methods
            
        
