# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:  Erney Ramírez Aportela (eramirez@cnb.csic.es), Feb 2018
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

import math
from glob import glob
import numpy

import pyworkflow.em as em 
from pyworkflow.em.protocol import ProtPreprocessVolumes 
from pyworkflow.utils import * 
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pyworkflow.protocol.params import (PointerParam, StringParam, BooleanParam, FloatParam, LabelParam, EnumParam, FileParam)
from phenix import  *
from pyworkflow.utils.path import createLink
from pyworkflow.em.convert import ImageHandler
import pyworkflow.em.packages.xmipp3 as xmipp3
from pyworkflow.em.convert_header.CCP4.convert import (
    adaptFileToCCP4, runCCP4Program, START, getProgram)

OUTPUT_SHARP = 'outputvolume'
OUTPUT_HALF1 = 'outputhalf1'
OUTPUT_HALF2 = 'outputhalf2'
#COEF_SHARP = 'coefficients'
#SHIFTED_INPUT = 'a'
#SHIFTED_SHARP = 'b'

class PhenixProtAutomatedSharpening(ProtPreprocessVolumes):
    """ Protocol for make sharpening to a input cryoEM map
    
    This is actually a program phenix.auto_sharpen from Phenix package.
    See documentation at:
       https://www.phenix-online.org/documentation/reference/auto_sharpen.html
    """
    _label = 'automated sharpening'
    
    #IMPORT PDB
    IMPORT_FROM_ID = 0
    IMPORT_OBJ = 1
    IMPORT_FROM_FILES = 2 
    
    #For Sharpening
    SHARP_ISOCUT = 0
    SHARP_ISO = 1
    SHARP_TARGET = 2
    SHARP_RESDEPEND = 3
    SHARP_SPLITV = 4
    SHARP_PDB = 5
    SHARP_NOSHARP = 6  
    
    #EVALUATION
    METHOD_ADJUSTED = 0
    METHOD_KURTOSIS = 1
        
    #--------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputMap', PointerParam, pointerClass='Volume',
                      label="Input map", important=True,  
                      help='File with CCP4-style map')  
        form.addParam('sharpSplitVolumes', BooleanParam, default=False,
                      label="Sharpening also for half volumes?",
                      help='In addition to the main volume it makes sharpening for half volumes.')        
        form.addParam('useSplitVolumes', BooleanParam, default=False,
                      label="Use half volumes for FSC?",
                      help='Use split volumes for FSC calculation.' 
                      ' Then you must select halfMapSharp in sharpenig methods')                          
        form.addParam('usePDB', BooleanParam, default=False,
                      label="Introduce atomic model?",
                      help='If a model is supplied, the map will be adjusted' 
                      ' to maximize map-model correlation.'
                      ' Then you must select modelSharp in sharpenig methods')
        #if usePDB      
        form.addParam('inputPdbData', EnumParam, choices=['id', 'object', 'file'],
                      label="Retrieve PDB from", default=self.IMPORT_FROM_ID,
                      display=EnumParam.DISPLAY_HLIST, condition='usePDB',
                      help='Retrieve PDB data from server,' 
                      ' use a pdb Object, or a local file')
        form.addParam('pdbId', StringParam, 
                      label="Pdb Id ", allowsNull=True,
                      condition='inputPdbData == IMPORT_FROM_ID and usePDB', 
                      help='Type a pdb Id (four alphanumeric characters).')
        form.addParam('pdbObj', PointerParam, pointerClass='PdbFile',
                      label="Input pdb ",  allowsNull=True,
                      condition='inputPdbData == IMPORT_OBJ and usePDB',
                      help='Specify a pdb object.')
        form.addParam('pdbFile', FileParam,  
                      label="File path", allowsNull=True, 
                      condition='inputPdbData == IMPORT_FROM_FILES and usePDB',
                      help='Specify a path to desired PDB structure.')       
  
        form.addParam('resolution', FloatParam, default=-1, 
                      label='resolution', 
                      help="Optimal nominal resolution of the map") 
               
        form.addParam('sharpening', EnumParam, 
                      choices=['b_iso_to_d_cut', 'b_iso', 'target_b_iso_to_d_cut', 
                               'resolDependent',  'halfMapSharp', 'modelSharp', 'no_sharpening'],
                      label='sharpening methods', default=self.SHARP_ISOCUT,                    
                      help="Methods to use in sharpening. b_iso searches for b_iso" 
                      " to maximize sharpening target (kurtosis or adjusted_sa)."
                      " b_iso_to_d_cut applies b_iso only up to resolution specified,"
                      " with fall-over of k_sharpen. Resolution dependent adjusts 3 parameters"
                      " to sharpen variably over resolution range. Default is b_iso_to_d_cut" )     
    

        form.addParam('sharpTarget', EnumParam, choices=['adjusted_sa', 'kurtosis'],
                      label='sharpening_target', default=self.METHOD_ADJUSTED,                    
                       help='Overall target for sharpening.' 
                      ' Can be "kurtosis" or "adjusted" (adjusted surface area).' 
                      ' Used to decide which sharpening approach is used.') 
             
        form.addParam('bIso',  FloatParam, default=-1, 
                      expertLevel=LEVEL_ADVANCED, 
                      label='target b_iso',
                      help='Target B-value for map' 
                      ' (sharpening will be applied to yield this value of b_iso)' ) 
        form.addParam('bSharp',  FloatParam, default=-1, 
                      expertLevel=LEVEL_ADVANCED,                       
                      label='b_sharpen',
                      help='Sharpen with this b-value.' 
                      ' Contrast with b_iso that yield a targeted value of b_iso' )             
                                
           
    #--------------------------- INSERT steps functions ------------------------
    def _createFilenameTemplates(self):
        """ Centralize how files are called """
        myDict = {
                 OUTPUT_SHARP: self._getExtraPath('sharpened_map_file.mrc'),
                 OUTPUT_HALF1: self._getExtraPath('sharpened_half1_file.mrc'),
                 OUTPUT_HALF2: self._getExtraPath('sharpened_half2_file.mrc'),                                                  
                 }
        self._updateFilenamesDict(myDict)    
    
    def _insertAllSteps(self):
        self._createFilenameTemplates()

                               
        if self.inputPdbData == self.IMPORT_FROM_ID:
             self._insertFunctionStep('pdbDownloadStep')    
        self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('autoSharpeningStep')
        self._insertFunctionStep('createOutputStep')    
        
    #--------------------------- STEPS functions -------------------------------
    
    def pdbDownloadStep(self):
        """Download all pdb files in file_list and unzip them."""
        if self.usePDB: 
             em.downloadPdb(self.pdbId.get(), self._getPdbFileName(), self._log)
       
    def convertInputStep(self):
        fnVol = self.inputMap.get()
        inFileName = fnVol.getFileName()
        if inFileName.endswith(":mrc"):
            inFileName.replace(":mrc", "")
        self.vol0Fn = self._getTmpPath("inputMap.mrc")
        origin = fnVol.getOrigin(force=True).getShifts()
        sampling = fnVol.getSamplingRate()
        adaptFileToCCP4(inFileName, self.vol0Fn, origin, sampling,
                        START)
        
        if self.useSplitVolumes.get() is True or self.sharpSplitVolumes.get() is True:          
             self.vol1Fn, self.vol2Fn = self.inputMap.get().getHalfMaps()
             halfFilename1=self.vol1Fn
             if halfFilename1.endswith(":mrc"):
                 halfFilename1.replace(":mrc", "")
             halfFilename2=self.vol2Fn
             if halfFilename2.endswith(":mrc"):
                 halfFilename2.replace(":mrc", "")                        
             self.vol1FnM = self._getTmpPath("halfMap1.mrc")       
             self.vol2FnM = self._getTmpPath("halfMap2.mrc")       
             adaptFileToCCP4(halfFilename1, self.vol1FnM, origin, sampling,
                        START)   
             adaptFileToCCP4(halfFilename2, self.vol2FnM, origin, sampling,
                        START)                    
        
#         if self.useSplitVolumes.get() is True or self.sharpSplitVolumes.get() is True:  
#              self.vol1Fn, self.vol2Fn = self.inputMap.get().getHalfMaps()
#              self.vol1FnM = self._getTmpPath("halfMap1.mrc")
#              img.convert(self.vol1Fn, self.vol1FnM)            
#              self.vol2FnM = self._getTmpPath("halfMap2.mrc")  
#              img.convert(self.vol2Fn, self.vol2FnM) 
#              
#              xmipp3.runXmippProgram("xmipp_image_header","-i %s --sampling_rate %f"
#                                    %(self.vol1FnM+":mrc",self.inputMap.get().getSamplingRate()))           
#              xmipp3.runXmippProgram("xmipp_image_header","-i %s --sampling_rate %f"
#                                    %(self.vol2FnM+":mrc",self.inputMap.get().getSamplingRate()))                    
                
        if self.usePDB:  
             self.pdbFn = self._getPdbFileName()                           
    
    def autoSharpeningStep(self):
        params=''           
        if self.useSplitVolumes.get() is True:    
             params += ' half_map_file=%s' % self.vol1FnM      
             params += ' half_map_file=%s' % self.vol2FnM  
        if self.usePDB.get() is True:                                            
             params += ' pdb_file=%s' % self.pdbFn  
        if self.resolution.get() != -1:              
             params += ' resolution=%f' % self.resolution.get()  
d={}
d[self.SHARP_ISOCUT]='b_iso_to_d_cut'

d[self.sharpening]
          #auto_sharpen_methods        
        if self.sharpening == self.SHARP_ISOCUT:                                            
             params += ' auto_sharpen_methods=b_iso_to_d_cut'                  
        elif self.sharpening == self.SHARP_ISO:                                            
             params += ' auto_sharpen_methods=b_iso'  
        elif self.sharpening == self.SHARP_TARGET:                                             
             params += ' auto_sharpen_methods=target_b_iso_to_d_cut'                    
        elif self.sharpening == self.SHARP_RESDEPEND:                                             
             params += ' auto_sharpen_methods=resolution_dependent'  
        elif self.sharpening == self.SHARP_SPLITV:                                            
             params += ' auto_sharpen_methods=half_map_sharpening'              
        elif self.sharpening == self.SHARP_PDB:                                          
             params += ' auto_sharpen_methods=model_sharpening'             
        elif self.sharpening == self.SHARP_NOSHARP:                                             
             params += ' auto_sharpen_methods=no_sharpening'         
              
            
        if self.sharpTarget == self.METHOD_ADJUSTED:              
             params += ' sharpening_target=adjusted_sa  residual_target=adjusted_sa'    
        else:      
             params += ' sharpening_target=kurtosis  residual_target=kurtosis'     
                   
        if self.bIso.get() != -1:               
             params += ' b_iso=%f' % self.bIso.get()   
        if self.bSharp.get() != -1:                                       
             params += ' b_sharpen=%f' % self.bSharp.get()   
            
        self.runJob("phenix.auto_sharpen map_file=%s sharpened_map_file=%s" 
                    % (self.vol0Fn, self._getFileName(OUTPUT_SHARP)) , params)
        xmipp3.runXmippProgram("xmipp_transform_mirror","-i %s --flipZ"
                               %self._getFileName(OUTPUT_SHARP))
        xmipp3.runXmippProgram("xmipp_transform_geometry","-i %s --rotate_volume euler 180 90 0"
                               %self._getFileName(OUTPUT_SHARP))
        
        if self.sharpSplitVolumes:
             self.runJob("phenix.auto_sharpen map_file=%s sharpened_map_file=%s" 
                         % (self.vol1FnM, self._getFileName(OUTPUT_HALF1)) , params)
             xmipp3.runXmippProgram("xmipp_transform_mirror","-i %s --flipZ"
                                    %self._getFileName(OUTPUT_HALF1))
             xmipp3.runXmippProgram("xmipp_transform_geometry","-i %s --rotate_volume euler 180 90 0"
                                    %self._getFileName(OUTPUT_HALF1))
        
             self.runJob("phenix.auto_sharpen map_file=%s sharpened_map_file=%s" 
                         % (self.vol2FnM, self._getFileName(OUTPUT_HALF2)) , params)
             xmipp3.runXmippProgram("xmipp_transform_mirror","-i %s --flipZ"
                                    %self._getFileName(OUTPUT_HALF2))
             xmipp3.runXmippProgram("xmipp_transform_geometry","-i %s --rotate_volume euler 180 90 0"
                                    %self._getFileName(OUTPUT_HALF2))
 
       
    def createOutputStep(self):
        inputMap = self.inputMap.get()
               
        vol = em.Volume()          
        vol.setLocation(self._getFileName(OUTPUT_SHARP))        
        vol.setSamplingRate(inputMap.getSamplingRate())     
        
        fnVol = self.inputMap.get()
        origin = fnVol.getOrigin(force=True)
        vol.setOrigin(origin)
        self._defineOutputs(outputVol=vol)       
        self._defineSourceRelation(self.inputMap, vol)    
        
        if self.sharpSplitVolumes:            
            half1 = self._getFileName(OUTPUT_HALF1)
            half2 = self._getFileName(OUTPUT_HALF2)
            vol.setHalfMaps([half1, half2])   
            self._defineOutputs(outputVol=vol)
            self._defineSourceRelation(self.inputMap, vol)                      
                  
    #--------------------------- INFO functions --------------------------------
    def _summary(self):
        summary = []
        summary.append('Input Map for sharpening: %s' % self.inputMap.get().getFileName())
        return summary

    def _methods(self):
        summary = []
        summary.append('We computed the sharpening for initial cryoEM map %s.'
                       % self.inputMap.get().getFileName())
        summary.append('The method used is described '
                       'in [Terwilliger2018].')
        return summary
    
    def _validate(self):
        errors = []
        if which('phenix.auto_sharpen') is '':
            errors.append('You should have the directory'
            ' phenix-xxx/buils/bin in the PATH')  
            
        if self.inputMap.get():                  
             if self.resolution.get() == -1.0:  
                 errors.append('If you introduce input map,'
                          ' you should introduce nominal resolution of the map') 
         
        if self.inputMap.get():
             if self.sharpSplitVolumes.get() is True and not self.inputMap.get().hasHalfMaps():                     
                 errors.append('The input Map you use does not have any associated half volumes ')   
                 
        if self.inputMap.get():
             if self.useSplitVolumes.get() is True and not self.inputMap.get().hasHalfMaps():                     
                 errors.append('The input Map you use does not have any associated half volumes ')                   
             if self.useSplitVolumes.get() is True and self.inputMap.get().hasHalfMaps(): 
                 if  self.sharpening != self.SHARP_SPLITV:   
                     errors.append('If you use half Maps for FSC calculation,'
                                ' you can use halfMapSharp in sharpening methods.')    
                     
        if self.usePDB.get() is True and  self.sharpening != self.SHARP_PDB:    
            errors.append('If you use atomic model,' 
                          ' you can use modelSharp in sharpening methods.')     
                                                      
        return errors

    def _citations(self):
        return ['Terwilliger2018']
    
    #--------------------------- UTLIS functions --------------------------------------------
    def _getPdbFileName(self):
        if self.inputPdbData == self.IMPORT_FROM_ID:
            return self._getExtraPath('%s.pdb' % self.pdbId.get())
        elif self.inputPdbData == self.IMPORT_OBJ:
            return self.pdbObj.get().getFileName()
        else:
            return self.pdbFile.get()
    

        