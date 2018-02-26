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

OUTPUT_SHARP = 'outputvolume'
COEF_SHARP = 'coefficients'
SHIFTED_INPUT = 'a'
SHIFTED_SHARP = 'b'

class PhenixProtAutomatedSharpening(ProtPreprocessVolumes):
    """ Protocol for make sharpening to a input cryoEM map
    
    This is actually a program phenix.auto_sharpen from Phenix package.
    See documentation at:
       https://www.phenix-online.org/documentation/reference/auto_sharpen.html
    """
    _label = 'automated sharpening'
    IMPORT_FROM_ID = 0
    IMPORT_OBJ = 1
    IMPORT_FROM_FILES = 2 
        
    #--------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputMap', PointerParam, pointerClass='Volume',
                      label="Input map", important=True,  
                      help='File with CCP4-style map')  
        form.addParam('useSplitVolumes', BooleanParam, default=False,
                      label="Use half volumes?",
                      help='Use to split volumes for FSC calculation.' 
                      ' Must have grid identical to map file')                          
        form.addParam('volumeHalf1', PointerParam,
                      label="Volume half 1", important=True,
                      pointerClass='Volume', condition="useSplitVolumes")
        form.addParam('volumeHalf2', PointerParam,
                      pointerClass='Volume', condition="useSplitVolumes",
                      label="Volume half 2", important=True)
        form.addParam('usePDB', BooleanParam, default=False,
                      label="Introduce atomic model?",
                      help='If a model is supplied, the map will be adjusted' 
                      ' to maximize map-model correlation')
        #if usePDB      
        form.addParam('inputPdbData', EnumParam, choices=['id', 'object', 'file'],
                      label="Retrieve PDB from", default=self.IMPORT_FROM_ID,
                      display=EnumParam.DISPLAY_HLIST, condition='usePDB',
                      help='Retrieve PDB data from server, use a pdb Object, or a local file')
        form.addParam('pdbId', StringParam, condition='inputPdbData == IMPORT_FROM_ID and usePDB', 
                      label="Pdb Id ", allowsNull=True,
                      help='Type a pdb Id (four alphanumeric characters).')
        form.addParam('pdbObj', PointerParam, pointerClass='PdbFile',
                      label="Input pdb ", condition='inputPdbData == IMPORT_OBJ and usePDB', allowsNull=True,
                      help='Specify a pdb object.')
        form.addParam('pdbFile', FileParam,
                      label="File path", condition='inputPdbData == IMPORT_FROM_FILES and usePDB', allowsNull=True,
                      help='Specify a path to desired PDB structure.')       
  
        form.addParam('resolution', FloatParam, default=-1, 
                      label='resolution', 
                      help="Optimal nominal resolution of the map") 
          
               
        form.addSection(label='Sharpening Methods') 
        form.addParam('sharpeningLabel', LabelParam, important=True,
                      label="Sharpening methods ",
                      help="Methods to use in sharpening. b_iso searches for b_iso" 
                      " to maximize sharpening target (kurtosis or adjusted_sa)."
                      " b_iso_to_d_cut applies b_iso only up to resolution specified,"
                      " with fall-over of k_sharpen. Resolution dependent adjusts 3 parameters"
                      " to sharpen variably over resolution range. Default is b_iso_to_d_cut" )
        form.addParam('halfMapSharp', BooleanParam, default=True,
                      label='Use half_map_sharpening', condition = 'useSplitVolumes')                          
        form.addParam('modelSharp', BooleanParam, default=True,
                      label='Use model_sharpening', condition = 'usePDB')              
        form.addParam('bIsoCut', BooleanParam, default=False,
                      label='Use b_iso_to_d_cut')  
        form.addParam('bIsoSharp', BooleanParam, default=False,
                      label='Use b_iso') 
        form.addParam('resolDependent', BooleanParam, default=False,
                      label='Use resolution_dependent')                               
        form.addParam('targetIsoCut', BooleanParam, default=False,
                      label='Use target_b_iso_to_d_cut')         
        form.addParam('noSharp', BooleanParam, default=False,
                      label='Use no_sharpening')   
        form.addParam('sharpTarget', StringParam, default="adjusted_sa",
                      label="sharpening_target", allowsNull=True,
                      help='Overall target for sharpening.' 
                      ' Can be "kurtosis" or "adjusted_sa" (adjusted surface area).' 
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
                 COEF_SHARP: self._getExtraPath('sharpened_map_file.mtz'),
                 SHIFTED_INPUT: self._getExtraPath('shifted_map_file.ccp4'),
                 SHIFTED_SHARP: self._getExtraPath('shifted_sharpened_map_file.ccp4')                                  
                 }
        self._updateFilenamesDict(myDict)    
    
    def _insertAllSteps(self):
        self._createFilenameTemplates()
        
        #self.micsFn = self._getPath()

        if self.useSplitVolumes:        
            self.vol1Fn = self.volumeHalf1.get().getFileName()
            self.vol2Fn = self.volumeHalf2.get().getFileName() 
        else:
            self.volumeHalf1.set(None)
            self.volumeHalf2.set(None)
            
                   
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
        
        self.vol0Fn = self.inputMap.get().getFileName()
        extVol0 = getExt(self.vol0Fn)
        if (extVol0 == '.mrc:mrc'):
            filenamebase = len(self.vol0Fn)
            self.vol0Fn = self.vol0Fn[0:(filenamebase-4)]          

        if self.useSplitVolumes.get() is True:
            extVol1 = getExt(self.vol1Fn)
            extVol2 = getExt(self.vol2Fn)
        if (extVol1 == '.mrc:mrc'):
            filenamebase1 = len(self.vol1Fn)
            self.vol1Fn = self.vol1Fn[0:(filenamebase1-4)]             
        if (extVol2 == '.mrc:mrc'):
            filenamebase2 = len(self.vol2Fn)
            self.vol2Fn = self.vol2Fn[0:(filenamebase1-4)]           
            
#             if (extVol1 == '.mrc') or (extVol1 == '.ccp4'):
#                 self.vol1Fn = self.vol1Fn + ':ccp4'
#             if (extVol2 == '.mrc') or (extVol2 == '.ccp4'):
#                 self.vol2Fn = self.vol2Fn + ':ccp4'
 
                
        if self.usePDB:  
            self.pdbFn = self._getPdbFileName()                           
    
    def autoSharpeningStep(self):
        params = ' map_file=%s' % self.vol0Fn            
        if self.useSplitVolumes.get() is True:    
            params += ' half_map_file=%s' % self.vol1Fn      
            params += ' half_map_file=%s' % self.vol2Fn  
        if self.usePDB.get() is True:                                            
            params += ' pdb_file=%s' % self.pdbFn  
        if self.resolution.get() != -1:              
            params += ' resolution=%f' % self.resolution.get()  
        #auto_sharpen_methods         
        if self.bIsoCut.get() is True:                                            
            params += ' auto_sharpen_methods=b_iso_to_d_cut'                  
        elif self.bIsoSharp.get() is True:                                            
            params += ' auto_sharpen_methods=b_iso'  
        elif self.targetIsoCut.get() is True:                                            
            params += ' auto_sharpen_methods=target_b_iso_to_d_cut'                    
        elif self.resolDependent.get() is True:                                            
            params += ' auto_sharpen_methods=resolution_dependent'  
        elif self.useSplitVolumes.get() is True:                                            
            params += ' auto_sharpen_methods=half_map_sharpening'              
        elif self.usePDB.get() is True:                                            
            params += ' auto_sharpen_methods=model_sharpening'             
        elif self.noSharp.get() is True:                                            
            params += ' auto_sharpen_methods=no_sharpening' 
        else: 
            params += ' auto_sharpen_methods=b_iso_to_d_cut' 
                         
        params += ' sharpening_target=%s  residual_target=%s' % (self.sharpTarget.get(), self.sharpTarget.get())     
        if self.bIso.get() != -1:               
            params += ' b_iso=%f' % self.bIso.get()   
        if self.bSharp.get() != -1:                                       
            params += ' b_sharpen=%f' % self.bSharp.get()   

        params += ' sharpened_map_file=%s' % self._getFileName(OUTPUT_SHARP) 
        params += ' sharpened_map_coeffs_file=%s' % self._getFileName(COEF_SHARP)         
        params += ' shifted_map_file=%s' % self._getFileName(SHIFTED_INPUT)  
        params += ' shifted_sharpened_map_file=%s' % self._getFileName(SHIFTED_SHARP)       
          

        self.runJob("phenix.auto_sharpen", params)      
        
    def createOutputStep(self):
        inputMap = self.inputMap.get()
               
        vol = em.Volume()          
        vol.setFileName(self._getFileName(OUTPUT_SHARP))        
        vol.setSamplingRate(inputMap.getSamplingRate())       
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
        if self.resolution.get() == -1.0:  
            errors.append('If you introduce input map, you should introduce nominal resolution '
                          'of the map')                           
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
    

        