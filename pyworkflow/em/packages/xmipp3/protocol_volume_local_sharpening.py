# -*- coding: utf-8 -*-
# **************************************************************************
# *
#* Authors:    Erney Ramirez-Aportela,                                    eramirez@cnb.csic.es
# *             Jose Luis Vilas,                                           jlvilas@cnb.csic.es
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
from pyworkflow import VERSION_1_1
from pyworkflow.protocol.params import (PointerParam, StringParam, 
                                        BooleanParam, FloatParam, IntParam, LEVEL_ADVANCED)
from pyworkflow.em.protocol.protocol_3d import ProtAnalysis3D
from pyworkflow.object import Float
from pyworkflow.em import ImageHandler
from pyworkflow.utils import getExt
from pyworkflow.em.data import Volume
import numpy as np
import os
import pyworkflow.em.metadata as md
from pyworkflow.em.metadata.constants import (MDL_COST, MDL_ITER, MDL_SCALE)
from ntpath import dirname

LOCALDEBLUR_METHOD_URL= 'http://github.com/I2PC/scipion/wiki/XmippProtLocSharp'

CHIMERA_RESOLUTION_VOL = 'MG_Chimera_resolution.vol'
BINARY_MASK = 'binaryMask.vol'
OUTPUT_RESOLUTION_FILE = 'resolutionMonoRes.vol'
OUTPUT_RESOLUTION_FILE_CHIMERA = 'MonoResChimera.vol'
METADATA_MASK_FILE = 'mask_data.xmd'
METADATA_PARAMS_SHARPENING = 'params.xmd'


class XmippProtLocSharp(ProtAnalysis3D):
    """    
    Given a resolution map the protocol calculate the sharpened map.
    """
    _label = 'localdeblur sharpening'
    _lastUpdateVersion = VERSION_1_1
    
    def __init__(self, **args):
        ProtAnalysis3D.__init__(self, **args)
    
    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')      
        
        form.addParam('inputVolume', PointerParam, pointerClass='Volume',
                      label="Input Map", important=True,
                      help='Select a volume for sharpening.')

        form.addParam('resolutionVolume', PointerParam, pointerClass='Volume',
                      label="Resolution Map", important=True,
                      help='Select a local resolution map.')
                
        form.addParam('const', FloatParam, default=1, 
                      expertLevel=LEVEL_ADVANCED,
                      label="lambda",
                      help='Regularization Param.' 
                        'The software determines this parameter automatically.'
                         'However, if the user wishes it can be modified')
        
        form.addParallelSection(threads = 4, mpi = 0)
  
    # --------------------------- INSERT steps functions --------------------------------------------


    def _insertAllSteps(self):
            # Convert input into xmipp Metadata format
        self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('createMaskStep')      
        self._insertFunctionStep('sharpeningAndMonoResStep')       
        self._insertFunctionStep('createOutputStep')


    def convertInputStep(self):
        """ Read the input volume.
        """

        self.volFn = self.inputVolume.get().getFileName()
        self.resFn = self.resolutionVolume.get().getFileName()      
        extVol = getExt(self.volFn)
        extRes = getExt(self.resFn)        
        if (extVol == '.mrc') or (extVol == '.map'):
            self.volFn = self.volFn + ':mrc'
        if (extRes == '.mrc') or (extRes == '.map'):
            self.resFn = self.resFn + ':mrc'     
    
    
    def createMaskStep(self):
        
        params = ' -i %s' % self.resFn
        params += ' -o %s' % self._getTmpPath(BINARY_MASK)
        params += ' --select below %f' % 0.1
        params += ' --substitute binarize'
             
        self.runJob('xmipp_transform_threshold', params)
        
    
    def MonoResStep(self, iter):
        sampling = self.inputVolume.get().getSamplingRate()
        
        if (iter == 1):
            imageFile = self.resolutionVolume.get().getFileName()
        else:
            imageFile = self._getTmpPath(OUTPUT_RESOLUTION_FILE)
            
        pathres=dirname(imageFile)
               
        img = ImageHandler().read(imageFile)
        imgData = img.getData()
        max_res = np.amax(imgData)
        
        significance = 0.95

        mtd = md.MetaData()
        mtd.read(pathres+"/"+(METADATA_MASK_FILE))  
        radius = mtd.getValue(MDL_SCALE,1)
        
        params = ' --vol %s' % self._getExtraPath('sharpenedMap_'+str(iter)+'.mrc')
        params += ' --mask %s' % self._getTmpPath(BINARY_MASK)
        params += ' --sampling_rate %f' % sampling
        params += ' --minRes %f' % (2*sampling)
        params += ' --maxRes %f' % max_res           
        params += ' --step %f' % 0.25
        params += ' --mask_out %s' % self._getTmpPath('refined_mask.vol')
        params += ' -o %s' % self._getTmpPath(OUTPUT_RESOLUTION_FILE)
        params += ' --volumeRadius %f' % radius
        params += ' --exact'
        params += ' --chimera_volume %s' % self._getTmpPath(
                                                OUTPUT_RESOLUTION_FILE_CHIMERA)
        params += ' --sym %s' % 'c1'
        params += ' --significance %f' % significance
        params += ' --md_outputdata %s' % self._getTmpPath(METADATA_MASK_FILE)
        params += ' --threads %i' % self.numberOfThreads.get()

        self.runJob('xmipp_resolution_monogenic_signal', params)
        
        
    def sharpenStep(self, iter):   
        sampling = self.inputVolume.get().getSamplingRate()   
        #params = ' --vol %s' % self.volFn
        if (iter == 1):
            params = ' --resolution_map %s' % self.resFn
        else:
            params = ' --resolution_map %s' % self._getTmpPath(
                                                   OUTPUT_RESOLUTION_FILE)
            
        params += ' --sampling %f' % sampling
        params += ' -n %i' %  self.numberOfThreads.get()   
        params += ' --md %s' % self._getTmpPath(METADATA_PARAMS_SHARPENING)  
        if (iter == 1 and self.const!=1):
             params += ' -l %f' % self.const
        #params += ' -o %s' % self._getExtraPath('sharpenedMap.vol')
         
        if (iter == 1):
            invol = ' --vol %s' % self.volFn
        else:
            invol = ' --vol %s' % self._getExtraPath('sharpenedMap_'+str(iter-1)+'.mrc')
            
        params +=invol    
            
        self.runJob("xmipp_volume_local_sharpening  -o %s"
                    %(self._getExtraPath('sharpenedMap_'+str(iter)+'.mrc')), params)
        self.runJob("xmipp_image_header -i %s -s %f"
                    %(self._getExtraPath('sharpenedMap_'+str(iter)+'.mrc'), sampling), "")


    def sharpeningAndMonoResStep(self):
        last_Niters = -1
        last_lambda_sharpening = 1e38
        nextIter = True
        iteration = 0
        while nextIter is True:
            iteration = iteration + 1
            #print iteration
            print ("Iteration  %s"  % (iteration))
            self.sharpenStep(iteration)           
            mtd = md.MetaData()
            mtd.read(self._getTmpPath(METADATA_PARAMS_SHARPENING))
            
            lambda_sharpening = mtd.getValue(MDL_COST,1)
            Niters = mtd.getValue(MDL_ITER,1)
            
#             if (Niters == last_Niters):
#                 nextIter = False
#                 break
               
            if (abs(lambda_sharpening - last_lambda_sharpening)<= 0.2):
                nextIter = False   
                break

            last_Niters = Niters
            last_lambda_sharpening = lambda_sharpening
            
            self.MonoResStep(iteration)
            
            imageFile = self._getTmpPath(OUTPUT_RESOLUTION_FILE)

            img = ImageHandler().read(imageFile)
            imgData = img.getData()
            max_res = np.amax(imgData)
            min_res = 2*self.inputVolume.get().getSamplingRate()
            
            #print ("minres %s  y maxres  %s"  % (min_res, max_res))          
        
            if (max_res-min_res<0.75):
                nextIter = False
                break
                
        os.system('cp '  +self._getExtraPath('sharpenedMap_'+str(iteration)+'.mrc')+
                   ' '  +self._getExtraPath('sharpenedMap.mrc'))   
  
             
    def createOutputStep(self):
        volume=Volume()
        volume.setFileName(self._getExtraPath('sharpenedMap.mrc'))
        volume.setSamplingRate(self.inputVolume.get().getSamplingRate())

        self._defineOutputs(sharpened_map=volume)
        self._defineSourceRelation(self.inputVolume, volume)          
            
    # --------------------------- INFO functions ------------------------------

    def _methods(self):
        messages = []
        if hasattr(self, 'sharpened_map'):
            messages.append(
                'Information about the method/article in ' + LOCALDEBLUR_METHOD_URL)
        return messages
    
    def _summary(self):
        summary = []
        summary.append("LocalDeblur Map")
        return summary

    def _citations(self):
        return ['Ramirez-Aportela2018']

