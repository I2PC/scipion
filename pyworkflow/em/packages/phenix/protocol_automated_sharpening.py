# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:  Erney Ramírez Aportela (eramirez@cnb.csic.es), May 2018
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

from pyworkflow.em import *  
from pyworkflow.utils import * 
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pyworkflow.protocol.params import (PointerParam, StringParam, BooleanParam, FloatParam, LabelParam)
import phenix
from pyworkflow.utils.path import createLink
from pyworkflow.em.convert import ImageHandler

# TODO: Move to 3D Tools
class PhenixProtAutomatedSharpening(ProtPreprocessVolumes):
    """ Protocol for make sharpening to a input cryoEM map
    
    This is actually a program phenix.auto_sharpen from Phenix package.
    See documentation at:
       https://www.phenix-online.org/documentation/reference/auto_sharpen.html
    """
    _label = 'automated sharpening'
    
    #--------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputMap', PointerParam, pointerClass='Volume',
                      label="Input map", important=True,  
                      help='File with CCP4-style map')  
        form.addParam('useSplitVolumes', BooleanParam, default=False,
                      label="Use half volumes?",
                      help='Use to split volumes for FSC calculation. Must have grid identical to map file')                          
        form.addParam('volumeHalf1', PointerParam,
                      label="Volume half 1", important=True,
                      pointerClass='Volume', condition="useSplitVolumes")
        form.addParam('volumeHalf2', PointerParam,
                      pointerClass='Volume', condition="useSplitVolumes",
                      label="Volume half 2", important=True)
        form.addParam('usePDB', BooleanParam, default=False,
                      label="Introduce atomic model?",
                      help='If a model is supplied, the map will be adjusted to maximize map-model correlation')           
        form.addParam('inputModel', PointerParam, pointerClass='PdbFile',
                      label="Input atomic model", important=True, condition="usePDB")  
        form.addParam('resolution', FloatParam, default=-1, 
                      label='resolution', 
                      help="Optimal nominal resolution of the map") 
          
               
        form.addSection(label='Sharpening Methods') 
        form.addParam('sharpeningLabel', LabelParam, important=True,
                      label="Sharpening methods ",
                      help="Methods to use in sharpening. b_iso searches for b_iso" 
                      "to maximize sharpening target (kurtosis or adjusted_sa)."
                      "b_iso_to_d_cut applies b_iso only up to resolution specified,"
                      "with fall-over of k_sharpen. Resolution dependent adjusts 3 parameters"
                      "to sharpen variably over resolution range. Default is b_iso_to_d_cut" )      
        form.addParam('b_iso_to_d_cut', BooleanParam, default=True,
                      label='Use b_iso_to_d_cut')  
        form.addParam('b_iso_sharp', BooleanParam, default=False,
                      label='Use b_iso')                       
        form.addParam('target_b_iso_to_d_cut', BooleanParam, default=False,
                      label='Use target_b_iso_to_d_cut')         
        form.addParam('half_map_sharpening', BooleanParam, default=False,
                      label='Use half_map_sharpening')                  
        form.addParam('resolution_dependent', BooleanParam, default=False,
                      label='Use resolution_dependent')         
        form.addParam('model_sharpening', BooleanParam, default=False,
                      label='Use model_sharpening') 
        form.addParam('no_sharpening', BooleanParam, default=False,
                      label='Use no_sharpening')         
        form.addParam('none', BooleanParam, default=False,
                      label='none') 
        form.addParam('b_iso',  FloatParam, default=-1, 
                      label='target b_iso',
                      help="Target B-value for map (sharpening will be applied to yield this value of b_iso)" ) 
        form.addParam('b_sharpen',  FloatParam, default=-1, 
                      label='b_sharpen',
                      help="Sharpen with this b-value. Contrast with b_iso that yield a targeted value of b_iso" )         
        form.addParam('b_blur_hires',  FloatParam, default=200, 
                      label='b_blur_hires',
                      help="Blur high_resolution data (higher than d_cut) with this b-value." 
                      "Contrast with b_sharpen applied to data up to d_cut" )       
                                
 
            
    #--------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):
        
        self.micsFn = self._getPath()

        if self.useSplitVolumes:          
            self.vol1Fn = self.volumeHalf1.get().getFileName()
            self.vol2Fn = self.volumeHalf2.get().getFileName()

        else:
            self.volumeHalf1.set(None)
            self.volumeHalf2.set(None)
            
        self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('autoSharpeningStep')
        self._insertFunctionStep('createOutputStep')    
        
    #--------------------------- STEPS functions -------------------------------
    def convertInputStep(self):
        " Convert input volume to .ccp4 as expected by phenix.auto_sharpen."
        
        self.vol0Fn = self.inputMap.get().getFileName()   
        extVol0 = getExt(self.vol0Fn)
        if (extVol0 == '.mrc') or (extVol0 == '.ccp4'):
            self.vol0Fn = self.vol0Fn + ':ccp4'             

        if self.useSplitVolumes.get() is True:
            extVol1 = getExt(self.vol1Fn)
            extVol2 = getExt(self.vol2Fn)
            if (extVol1 == '.mrc') or (extVol1 == '.ccp4'):
                self.vol1Fn = self.vol1Fn + ':ccp4'
            if (extVol2 == '.mrc') or (extVol2 == '.ccp4'):
                self.vol2Fn = self.vol2Fn + ':ccp4'
                
        if self.usePDB:            
            inputModel=os.path.abspath(self.inputModel.get().getFileName())       
    
    def autoSharpeningStep(self):
        params = ' map_file=%s' % self.vol0Fn            
        if self.useSplitVolumes.get() is True:    
            params += ' half_map_file=%s' % self.vol1Fn      
            params += ' half_map_file=%s' % self.vol2Fn #I do not Know  
        if self.usePDB.get() is True:                                            
            params += ' pdb_file=%s' % inputModel  
        if self.resolution.get() != -1:              
            params += ' resolution=%f' % self.resolution.get()  
        #auto_sharpen_methods         
        if self.b_iso_to_d_cut.get() is True:                                            
            params += ' auto_sharpen_methods=b_iso_to_d_cut'                  
        if self.b_iso_sharp.get() is True:                                            
            params += ' auto_sharpen_methods=b_iso'  
        if self.target_b_iso_to_d_cut.get() is True:                                            
            params += ' auto_sharpen_methods=target_b_iso_to_d_cut'    
        if self.half_map_sharpening.get() is True:                                            
            params += ' auto_sharpen_methods=half_map_sharpening'                  
        if self.resolution_dependent.get() is True:                                            
            params += ' auto_sharpen_methods=resolution_dependent'  
        if self.model_sharpening.get() is True:                                            
            params += ' auto_sharpen_methods=model_sharpening'             
        if self.no_sharpening.get() is True:                                            
            params += ' auto_sharpen_methods=no_sharpening'  
        if self.none.get() is True:                                            
            params += ' auto_sharpen_methods=none'  
        if self.b_iso.get() != -1:               
            params += ' b_iso=%f' % self.b_iso.get()   
        if self.b_sharpen.get() != -1:                                       
            params += ' b_sharpen=%f' % self.b_sharpen.get()   
        if self.b_blur_hires.get() != -1:                    
            params += ' b_blur_hires=%f' % self.b_blur_hires.get()         
        params += 'sharpened_map_file=%s' % self._getExtraPath("sharpened_map_file.ccp4") 
        params += 'sharpened_map_coeffs_file=%s' % self._getExtraPath("sharpened_map_file.mtz")         
        params += 'shifted_map_file=%s' % self._getExtraPath("shifted_map_file.ccp4")  
        params += 'shifted_sharpened_map_file=%s' % self._getExtraPath("shifted_sharpened_map_file.ccp4")  
 
        print("HOLAAA111")
        self.runJob("source /home/erney/git/scipion/software/em/phenix-1.13-2998/phenix_env.sh")       
        print("HOLAAA222")                            
        self.runJob("phenix.auto_sharpen", params)
        
        
    def createOutputStep(self):
        sharpVol = self._getExtraPath("sharpened_map.ccp4")
        vol = em.Volume()        
        vol.setLocation(SharpVol)
        vol.setSamplingRate(inputMap.getSamplingRate())       
        self._defineOutputs(outputVol=vol)       
        self._defineSourceRelation(self.inputMap, vol)    
          
    #--------------------------- INFO functions --------------------------------
    def _summary(self):
        summary = []
        summary.append('Input Map for sharpening: %s' % self.vol0Fn)
        return summary

    def _methods(self):
        summary = []
        summary.append('We computed the sharpening for initial cryoEM map %s.'
                       % self.vol0Fn)
        summary.append('We make the automated map sharpening using the method described '
                       'in [Terwilliger2018].')
        return summary

    def _citations(self):
        return ['Terwilliger2018']
    
    #def _validate(self):
        #errors = []
        #if which('phenix.auto_sharpen') is '':
        #    errors.append('You should have the program phenix.auto_sharpen in the PATH')
        #return errors
        