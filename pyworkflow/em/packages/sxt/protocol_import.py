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
from pyworkflow.em.packages.sxt.data import TiltSeries




class ProtImportTiltSeries(ProtImportImages):
    """    
    This prtocol is to import tilt seies and related info included tilt angles.
    """
    _label = 'import tilt series'     
           
    #--------------------------- DEFINE param functions --------------------------------------------    
    def _defineParams(self, form):        
        form.addSection(label='Import')

        form.addParam('filesPath', params.PathParam,                      
                      label="Files directory",
                      help="Directory with the .hdf5 file you want to import.")        
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
        fnIn = self.filesPath.get()
        copyOrLink = self.getCopyOrLink()
        dst = self._getExtraPath(basename(fnIn))        
        copyOrLink(fnIn, dst)
        fnOutMd = self._defineOutputMdName()
        fnRoot = removeExt(dst)
        fnStack = fnRoot + '.mrc'
                      
        self._insertFunctionStep('importHdf5Step', dst, fnRoot, fnOutMd, fnStack)
        self._insertFunctionStep('createOutputStep', fnStack, dst, fnOutMd)     
          
    #--------------------------- STEPS functions --------------------------------------------    
    def importHdf5Step(self, fnIn, fnRoot, fnOutMd, fnStack):
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
              
    def createOutputStep(self, fnStack, fnIn, fnAngles):        
        pattern = expandPattern(fnIn)
        self.info("Using input file: '%s'" % pattern)
                        
        tiltSeries = TiltSeries()
        tiltSeries.setFileName(fnStack)        
        acquisition = tiltSeries.getXrayAcquisition()        
        self.fillXrayAcquisition(acquisition) 
        if tiltSeries.hasXrayAcquisition():
            tiltSeries.setXrayAcquisition(acquisition)      
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
    
    #--------------------------- INFO functions --------------------------------        
    def _validate(self):
        errors = []
        if self.filesPath.get():
            if not self.filesPath.get().endswith('hdf5'):
                errors.append ("Expected hdf5 files for importing!!!") 
        else:
            errors.append("The path can not be empty!!!")     
        return errors
        
    def _summary(self):
        summary = []
        outputSet = self._getOutputSet()
        if outputSet is None:
            summary.append("Output TiltSeries is not ready yet.") 
            if self.copyFiles:
                summary.append("*Warning*: You select to copy files into your project.\n"
                               "This will make another copy of your data and may take \n"
                               "more time to import.")
        else:
            summary.append("*%d* %s imported from:  %s" % (
                                outputSet.getSize(),
                                'images related to the selected TiltSeries',
                                basename(self.filesPath.get())))            
            summary.append("Sampling rate : *%0.2f* A/px" % (
                                outputSet.getSamplingRate()*10))
            summary.append("*Imaging info*:\n" +
                           "Lens label: %s \n" % self.lensLabel.get() +
                           "Energy (ev): %f \n" % self.energy.get() +
                           "Date of imaging (ddmmyyy): %s" % self.date.get())       
        return summary
    
    def _methods(self):
        methods = []
        outputSet = self._getOutputSet()
        if outputSet is not None:
            methods.append("*%d* %s were imported with a sampling rate of "
                           "*%0.2f* A/px (Lens label: %s, "
                           "Energy (ev): %f, Date of imaging (ddmmyyy): %s)."
                           " Output set is %s."
                           % (outputSet.getSize(), 'images related to the selected TiltSeries',
                              outputSet.getSamplingRate()*10,
                              self.lensLabel.get(),
                              self.energy.get(), self.date.get(),
                              self.getObjectTag('outputTiltSeries')))            
        return methods 
       
    #--------------------------- UTILS functions -------------------------------    
    def fillXrayAcquisition(self, acquisition):
        """ Fill the acquition object with protocol params. """
        acquisition.setLensLabel(self.lensLabel.get())
        acquisition.setEnergy(self.energy.get())
        acquisition.setDate(self.date.get())
    
    def _getOutputSet(self):
        return getattr(self, 'outputTiltSeries', None)
    
    def _defineOutputMdName(self):
        return self._getExtraPath('normalized_image_plus_related_Angles.xmd')
    