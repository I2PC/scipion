# **************************************************************************
# *
# * Authors:     Mohsen Kazemi  (mkazemi@cnb.csic.es)
# *              
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

import sys
from pyworkflow.em.protocol.protocol_import import ProtImportImages
import pyworkflow.protocol.params as params
from os.path import basename
from h5py import File
from pyworkflow.utils.path import removeExt, cleanPattern, cleanPath, expandPattern, copyFile
from pyworkflow.em import ImageHandler
import numpy as np
import xmipp
from pyworkflow.em.packages.sxt.data import TiltSeries, SetOfTiltSeries
from pyworkflow.utils import getFloatListFromValues
import pyworkflow.em.metadata as md
from pyworkflow.mapper.sqlite_db import SqliteDb



class ProtImportTiltSeries(ProtImportImages):
    """    
    This prtocol is to import tilt seies and related info included tilt angles.
    """
    _label = 'import tilt series'     
           
    #--------------------------- DEFINE param functions --------------------------------------------    
    def _defineParams(self, form):        
        form.addSection(label='Import')

        form.addParam('importFocalSeries', params.BooleanParam, default=False,
                      label='Import focal tilt series?',
                      help="This is used to add more parameters "
                           "related to the focal tilt series")
        form.addParam('filesPath', params.PathParam,                      
                      label="Files directory",
                      help="Directory with the .hdf5 file you want to import.\n\n"
                           "To import tilt series, select the related hdf5 file.\n"
                           "To import focal tilt series, select the directory "
                           "includes hdf5 files, then define the file patterns "
                           "in the related input parameters field")        
        form.addParam('filesPattern', params.StringParam, 
                      condition="importFocalSeries",
                      label='Pattern', 
                      help="Pattern of the files to be imported.\n\n"
                           "The pattern can contain standard wildcards such as\n"
                           "*, ?, etc, or special ones like ### to mark some\n"
                           "digits in the filename as ID.")
        form.addParam('copyFiles', params.BooleanParam, default = False, 
                      expertLevel = params.LEVEL_ADVANCED, 
                      label = "Copy files?",
                      help = "By default the files are not copied into the\n"
                             "project to avoid data duplication and to save\n"
                             "disk space. Instead of copying, symbolic links are\n"
                             "created pointing to original files. This approach\n"
                             "has the drawback that if the project is moved to\n"
                             "another computer, the links need to be restored.\n")        
        form.addParam('doCrop', params.BooleanParam, default=False,
                      label='Crop input tilt series?',
                      help="This is used to avoid some black pixels of some cameras")
        form.addParam('cropSize', params.IntParam, default=20, condition="doCrop",
                      label='Crop size (px)',
                      help="Number of pixels to crop from each side")
        
        
        """
        
        group = form.addGroup('Focal tilt series info', 
                              condition="importFocalSeries")
        
        group.addParam('defocusValue', params.NumericListParam,
                      default="0 -2 +2",
                      label="Defocus value (micron)",
                      help="Defocus value related to each tilt series "
                           "in the focal series.\n"
                           "Note: the orders of these values should be the same "
                           "as the selected files are importing.")
        group.addParam('reference', params.IntParam, default=1,
                      label='Reference tilt series',
                      help="tilt series (file) with defocus value of 0")        
        """
        
        
        
        group = form.addGroup('Acquisition info')
        group.addParam('lensLabel', params.StringParam, default = 'lens1',
                   label = "Lens label", 
                   help = "Type of the lens has been used for imaging")
        group.addParam('energy', params.FloatParam, default = 00.00,
                   label = "Energy (ev)")        
        group.addParam('date', params.StringParam, default = '15072016',
                      label = "Date",
                      help = "Date of imaging in ddmmyyyy format")
        group.addParam('samplingRate', params.FloatParam,
                      label = "Sampling rate (nm/px)")         
                      
        form.addParallelSection(threads=0, mpi=0)
        
    #--------------------------- INSERT steps functions --------------------------------------------    
    def _insertAllSteps(self):        
        fnOutMd = self._defineOutputMdName()
        if not self.importFocalSeries.get():              
            fnIn = self.filesPath.get()
            copyOrLink = self.getCopyOrLink()
            dst = self._getExtraPath(basename(fnIn))        
            copyOrLink(fnIn, dst)            
            fnStack = removeExt(dst) + '.mrc'            
            self._insertFunctionStep('createOutputStepTiltSeries', fnStack, dst, fnOutMd)
        else:
            self._insertFunctionStep('createOutputStepSetOfTiltSeries', self.getPattern(), fnOutMd)    
          
    #--------------------------- STEPS functions --------------------------------------------    
    def createOutputStepTiltSeries(self, fnStack, fnIn, fnAngles):          
        self._importHdf5Step(fnIn, removeExt(fnIn), fnAngles, fnStack)
        
        pattern = expandPattern(fnIn)
        self.info("Using input file: '%s'" % pattern)
                        
        tiltSeries = TiltSeries()
        tiltSeries.setFileName(fnStack)        
        acquisition = tiltSeries.getXrayAcquisition()        
        self._fillXrayAcquisition(acquisition)    
        tiltSeries.setSamplingRate(self.samplingRate.get() * 10)
        
        # Read angles from metadata
        mdObj = xmipp.MetaData(fnAngles)
        angles = mdObj.getColumnValues(xmipp.MDL_ANGLE_TILT)
        tiltSeries.setAngles(angles)
        tiltSeries.setXrayAcquisition(acquisition)
        tiltSeries.setSize(np.shape(angles)[0])
        
        sys.stdout.write("\rImported %d/%d" % (1, 1))
        sys.stdout.flush()            
        print "\n"
                
        self._defineOutputs(outputTiltSeries=tiltSeries)
        
    def createOutputStepSetOfTiltSeries(self, pattern, fnOutMd):
        self.info("Using pattern: '%s'" % pattern)
        
        #fnStack = self._getExtraPath('output_focalSeries.mrc')
        #focalSeries = SetOfTiltSeries(fnStack)
        #focalSeries.setFileName(fnStack)
        focalSeries = self._createSetOfTiltSeries()
        focalSeries.setSamplingRate(self.samplingRate.get()* 10)
        
        mdOut = xmipp.MetaData()
        tiltIndex = 0
        refTilt = 0
        
        #outFiles = [focalSeries.getFileName()]
        imgh = ImageHandler()
        
        copyOrLink = self.getCopyOrLink()
        for i, (fileName, fileId) in enumerate(self.iterFiles()):
            dst = self._getExtraPath(basename(fileName))
            copyOrLink(fileName, dst)
            fnRoot = removeExt(dst)
            fnStackTilt = fnRoot + '.mrc'
            fnOutMdTilt = fnRoot + '_imgs_angles.xmd'
            self._importHdf5Step(dst, fnRoot, fnOutMdTilt, fnStackTilt)
            _, _, _, n = imgh.getDimensions(fnStackTilt)
                                    
            tiltSeries = TiltSeries()
            tiltSeries.setFileName(fnStackTilt)
            tiltSeries.setSize(n)
            tiltSeries.setSamplingRate(self.samplingRate.get() * 10)
            mdObj = xmipp.MetaData(fnOutMdTilt)
            angles = mdObj.getColumnValues(xmipp.MDL_ANGLE_TILT)
            tiltSeries.setAngles(angles)
            acquisition = tiltSeries.getXrayAcquisition()
            self._fillXrayAcquisition(acquisition)
            tiltSeries.setXrayAcquisition(acquisition)
            
            
            fileBaseName = basename(fileName)
            for focalId in range(100):
                if '_tomo_%02d' % focalId in fileBaseName:
                    tiltIndex += 1
                    if tiltIndex == 4:
                        tiltIndex = 1
                    break
            if '_0.hdf5' in fileBaseName:
                defocusValue = 0.0
                refTilt = tiltIndex
            elif '_m2.hdf5' in fileBaseName:
                defocusValue = -2.0
                if tiltIndex == 1:
                    refTilt = 2                    
            else:
                defocusValue = +2.0
                if tiltIndex == 1:
                    refTilt = 2
                
                   
            
            
            
            focalInfo = tiltSeries.getFocalSeries()            
            focalInfo.settiltSeriesGroup(focalId)
            focalInfo.setIndex(tiltIndex)
            focalInfo.setDefocus(defocusValue)
            focalInfo.setReference(refTilt)
            tiltSeries.setFocalSeries(focalInfo)
            
            
            objId = mdOut.addObject()
            mdOut.setValue(xmipp.MDL_TOMOGRAMMD, fnOutMdTilt, objId)
            mdOut.setValue(xmipp.MDL_XRAY_FOCAL_IDX, focalId, objId)
            mdOut.setValue(xmipp.MDL_XRAY_TILT_IDX, tiltIndex, objId)
            mdOut.setValue(xmipp.MDL_XRAY_DEFOCUS, defocusValue, objId)
            mdOut.setValue(xmipp.MDL_XRAY_REF_IDX, refTilt, objId)
            
            
            
            focalSeries.append(tiltSeries)
            sys.stdout.write("\rImported %d/%d\n\n" % (i+1, self.numberOfFiles))
            sys.stdout.flush()
            #outFiles.append(fnStackTilt)
        
        mdOut.write(fnOutMd)
        self._defineOutputs(outputFocalSeries=focalSeries)    

            
        
        
        
        
        
        
        
        
        
        
        
        """
        angles, sampling, index, size, reference, acquisition, defocus ... moshabehe acquisition, baraye focal series
        focalSeries = FocalSeries()
        focalSeries.setFileName(fnStack)        
        acquisition = tiltSeries.getXrayAcquisition()        
        self._fillXrayAcquisition(acquisition)      
        tiltSeries.setSamplingRate(self.samplingRate.get() * 10)
        """    
    
    #--------------------------- INFO functions --------------------------------        
    def _validate(self):
        errors = []
        if self.filesPath.get():
            if not self.filesPath.get().endswith('hdf5') and not self.filesPattern.get():
                errors.append ("Expected hdf5 files for importing or indicate the files pattern!!!") 
        else:
            errors.append("The path can not be empty!!!")   
        if self.getPattern():
            # Just check the number of files matching the pattern
            self.getMatchFiles()
            if self.numberOfFiles == 0:
                errors.append("There are no files matching the pattern %s"
                                % self.getPattern())     
        return errors
        
    def _summary(self):
        summary = []
        outputSet1 = self._getOutputSetTiltSeries()
        outputSet2 = self._getOutputFocalSeries()
        if not self.importFocalSeries.get() and outputSet1 is None:
            summary.append("Output TiltSeries is not ready yet.") 
            if self.copyFiles:
                summary.append("*Warning*: You select to copy files into your project.\n"
                               "This will make another copy of your data and may take \n"
                               "more time to import.")
        elif self.importFocalSeries.get() and outputSet2 is None:
            summary.append("Output FocalSeries is not ready yet.") 
            if self.copyFiles:
                summary.append("*Warning*: You select to copy files into your project.\n"
                               "This will make another copy of your data and may take \n"
                               "more time to import.")     
        if  outputSet1 is not None:  
            summary.append("*%d* %s imported from:  %s" % (
                                outputSet1.getSize(),
                                'images related to the selected TiltSeries',
                                basename(self.filesPath.get())))            
            summary.append("Sampling rate : *%0.2f* A/px" % (
                                outputSet1.getSamplingRate()*10))
            summary.append("*Imaging info*:\n" +
                           "Lens label: %s \n" % self.lensLabel.get() +
                           "Energy (ev): %f \n" % self.energy.get() +
                           "Date of imaging (ddmmyyy): %s" % self.date.get()) 
        elif outputSet2 is not None:
            summary.append("*%d* %s imported from:  %s" % (
                                outputSet2.getSize(),
                                'tiltSeries related to the selected FocalSeries',
                                self.filesPath.get()))            
            summary.append("Sampling rate : *%0.2f* A/px" % (
                                outputSet2.getSamplingRate()*10))
            summary.append("*Imaging info*:\n" +
                           "Lens label: %s \n" % self.lensLabel.get() +
                           "Energy (ev): %f \n" % self.energy.get() +
                           "Date of imaging (ddmmyyy): %s" % self.date.get())                  
        return summary
    
    def _methods(self):
        methods = []
        outputSet1 = self._getOutputSetTiltSeries()
        outputSet2 = self._getOutputFocalSeries()
        if outputSet1 is not None:
            methods.append("*%d* %s were imported with a sampling rate of "
                           "*%0.2f* A/px (Lens label: %s, "
                           "Energy (ev): %f, Date of imaging (ddmmyyy): %s)."
                           " Output set is %s."
                           % (outputSet1.getSize(), 'images related to the selected TiltSeries',
                              outputSet1.getSamplingRate()*10,
                              self.lensLabel.get(),
                              self.energy.get(), self.date.get(),
                              self.getObjectTag('outputTiltSeries'))) 
        if outputSet2 is not None:
            methods.append("*%d* %s were imported with a sampling rate of "
                           "*%0.2f* A/px (Lens label: %s, "
                           "Energy (ev): %f, Date of imaging (ddmmyyy): %s)."
                           " Output set is %s."
                           % (outputSet2.getSize(), 'tiltSeries related to the selected FocalSeries',
                              outputSet2.getSamplingRate()*10,
                              self.lensLabel.get(),
                              self.energy.get(), self.date.get(),
                              self.getObjectTag('outputFocalSeries')))           
        return methods 
       
    #--------------------------- UTILS functions -------------------------------    
    def _importHdf5Step(self, fnIn, fnRoot, fnOutMd, fnStack):
        fhHdf5 = File(fnIn, 'r')
        pixels = self.cropSize.get()
        if not "TomoNormalized" in fhHdf5:
            print "Input data need to be normalized ..."
            if self.doCrop.get():                
                self.runJob("xmipp_xray_import", 
                            "--mistral %s --oroot %s --crop %d" % (
                                                                   fnIn, 
                                                                   fnRoot,
                                                                   pixels))
            else:
                self.runJob("xmipp_xray_import", 
                            "--mistral %s --oroot %s" % (fnIn, fnRoot))                                    
            
            anglesArray = fhHdf5 ["NXtomo/data/rotation_angle"][()]            
            mdOut = xmipp.MetaData()            
            for k in range(np.shape(anglesArray)[0]):
                objId = mdOut.addObject()
                mdOut.setValue(xmipp.MDL_IMAGE, "%d@%s"%(k+1,fnStack), objId)
                mdOut.setValue(xmipp.MDL_ANGLE_TILT, anglesArray[k], objId)
            mdOut.write(fnOutMd)
            
            cleanPattern(fnRoot + '.tlt')
            cleanPattern(fnRoot + '*.xmp')
            cleanPattern(fnRoot + '.xmd')
                         
        else:
            print "Input data were already normalized ..."
            imgArray = fhHdf5 ["TomoNormalized/TomoNormalized"][()]            
            ih = ImageHandler()
            outputImg = ih.createImage()         
            
            i = 0
            for j in range(np.shape(imgArray)[0]):
                outputImg.setData(imgArray[j, :, :])
                i += 1
                outputImg.write((i,fnStack))
                
            if self.doCrop.get():
                fnStacknew = fnRoot + '_new.mrc'
                self.runJob("xmipp_transform_window", 
                            "-i %s -o %s --crop %d" % (fnStack, fnStacknew, pixels))
                copyFile(fnStacknew, fnStack)
                cleanPattern(fnStacknew)          
                
            anglesArray = fhHdf5 ["TomoNormalized/rotation_angle"][()]            
            mdOut = xmipp.MetaData()            
            for k in range(np.shape(anglesArray)[0]):
                objId = mdOut.addObject()
                mdOut.setValue(xmipp.MDL_IMAGE, "%d@%s"%(k+1,fnStack), objId)
                mdOut.setValue(xmipp.MDL_ANGLE_TILT, anglesArray[k], objId)            
            mdOut.write(fnOutMd)
    
    def _fillXrayAcquisition(self, acquisition):
        """ Fill the xrayAcquition object with protocol params. """
        acquisition.setLensLabel(self.lensLabel.get())
        acquisition.setEnergy(self.energy.get())
        acquisition.setDate(self.date.get())
    
    def _getOutputSetTiltSeries(self):
        return getattr(self, 'outputTiltSeries', None)
    
    def _getOutputFocalSeries(self):
        return getattr(self, 'outputFocalSeries', None)
    
    def _defineOutputMdName(self):
        return self._getExtraPath('normalized_tiltSeries_plus_related_info.xmd')
        
    def _createSetOfTiltSeries(self):
        """ Create a set and set the filename. 
        If the file exists, it will be delete. """
        setFn = self._getPath('focalSeries.sqlite')
        # Close the connection to the database if
        # it is open before deleting the file
        cleanPath(setFn)
        
        SqliteDb.closeConnection(setFn)        
        setObj = SetOfTiltSeries(filename=setFn)
        return setObj
    
    