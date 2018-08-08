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
from os.path import exists, lexists

LOCALDEBLUR_METHOD_URL='http://github.com/I2PC/scipion/wiki/XmippProtLocSharp' 

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
                      help='Select a local resolution map.'
                      ' LocalDeblur has been specially designed to work with'
                      ' resolution maps obtained with MonoRes,' 
                      ' however resolution map from ResMap and BlocRes are also accepted.')
                
        form.addParam('const', FloatParam, default=1, 
                      expertLevel=LEVEL_ADVANCED,
                      label="lambda",
                      help='Regularization Param.' 
                        'The method determines this parameter automatically.'
                        ' This parameter is directly related to the convergence.'
                        ' Increasing it would accelerate the convergence,' 
                        ' however it presents the risk of falling into local minima.')
        form.addParam('K', FloatParam, default=0.025, 
                      expertLevel=LEVEL_ADVANCED,
                      label="K",
                      help='K = 0.025 works well for all tested cases.'
                      ' K should be in the 0.01-0.05 range.'
                      ' For maps with FSC resolution lower than 6Å,'
                      ' K = 0.01 can be a good alternative.')        
        
        form.addParallelSection(threads = 4, mpi = 0)
        
    # --------------------------- INSERT steps functions --------------------------------------------

    def _createFilenameTemplates(self):
        """ Centralize how files are called """
        myDict= {           
                'BINARY_MASK': self._getTmpPath('binaryMask.vol'),
                'OUTPUT_RESOLUTION_FILE': self._getTmpPath('resolutionMonoRes.vol'),                
                'OUTPUT_RESOLUTION_FILE_CHIMERA': self._getTmpPath('MonoResChimera.vol'),
                'METADATA_PARAMS_SHARPENING': self._getTmpPath('params.xmd'),                                                                                 
                }
        self._updateFilenamesDict(myDict) 
    

    def _insertAllSteps(self):
        self._createFilenameTemplates() 
        self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('checkBackgroundStep')
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

    def checkBackgroundStep(self): 
        
        initRes=self.resolutionVolume.get().getFileName() 

        img = ImageHandler().read(initRes)
        imgData = img.getData()
        max_value = np.amax(imgData)
        min_value = np.amin(imgData) 
        #print ("minvalue %s  y maxvalue  %s"  % (min_value, max_value)) 
        
        if (min_value > 0.01):
            params = ' -i %s' % self.resFn  
            params += ' -o %s' % self.resFn    
            params += ' --select above %f' % (max_value-1)   
            params += ' --substitute value 0'
            
            self.runJob('xmipp_transform_threshold', params)  
                             
    def createMaskStep(self):
        
        params = ' -i %s' % self.resFn
        params += ' -o %s' % self._getFileName('BINARY_MASK')
        params += ' --select below %f' % 0.1
        params += ' --substitute binarize'
             
        self.runJob('xmipp_transform_threshold', params)
        
    
    def MonoResStep(self, iter):
        sampling = self.inputVolume.get().getSamplingRate()
        
        if (iter == 1):
            resFile = self.resolutionVolume.get().getFileName()
        else:
            resFile = self._getFileName('OUTPUT_RESOLUTION_FILE')
            
        pathres=dirname(resFile)
                       
        img = ImageHandler().read(resFile)
        imgData = img.getData()
        max_res = np.amax(imgData)
        
        significance = 0.95

        mtd = md.MetaData()
        if exists(os.path.join(pathres,'mask_data.xmd')):    
            mtd.read(os.path.join(pathres,'mask_data.xmd')) 
            radius = mtd.getValue(MDL_SCALE,1)
        else:
            xdim, _ydim, _zdim = self.inputVolume.get().getDim()
            radius = xdim*0.5 
        
        params = ' --vol %s' % self._getExtraPath('sharpenedMap_'+str(iter)+'.mrc')
        params += ' --mask %s' % self._getFileName('BINARY_MASK')
        params += ' --sampling_rate %f' % sampling
        params += ' --minRes %f' % (2*sampling)
        params += ' --maxRes %f' % max_res           
        params += ' --step %f' % 0.25
        params += ' --mask_out %s' % self._getTmpPath('refined_mask.vol')
        params += ' -o %s' % self._getFileName('OUTPUT_RESOLUTION_FILE')
        params += ' --volumeRadius %f' % radius
        params += ' --exact'
        params += ' --chimera_volume %s' % self._getFileName(
                                                'OUTPUT_RESOLUTION_FILE_CHIMERA')
        params += ' --sym %s' % 'c1'
        params += ' --significance %f' % significance
        params += ' --md_outputdata %s' % self._getTmpPath('mask_data.xmd')
        params += ' --threads %i' % self.numberOfThreads.get()

        self.runJob('xmipp_resolution_monogenic_signal', params)
        
        
    def sharpenStep(self, iter):   
        sampling = self.inputVolume.get().getSamplingRate()   
        #params = ' --vol %s' % self.volFn
        if (iter == 1):
            params = ' --resolution_map %s' % self.resFn
        else:
            params = ' --resolution_map %s' % self._getFileName(
                                                   'OUTPUT_RESOLUTION_FILE')
            
        params += ' --sampling %f' % sampling
        params += ' -n %i' %  self.numberOfThreads.get()   
        params += ' -k %f' %  self.K
        params += ' --md %s' % self._getFileName('METADATA_PARAMS_SHARPENING')  
        if (iter == 1 and self.const!=1):
             params += ' -l %f' % self.const
         
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
            print ('\n====================\n'
                'Iteration  %s'  % (iteration))
            self.sharpenStep(iteration)           
            mtd = md.MetaData()
            mtd.read(self._getFileName('METADATA_PARAMS_SHARPENING'))
            
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
            
            imageFile = self._getFileName('OUTPUT_RESOLUTION_FILE')

            img = ImageHandler().read(imageFile)
            imgData = img.getData()
            max_res = np.amax(imgData)
            min_res = 2*self.inputVolume.get().getSamplingRate()
            
            #print ("minres %s  y maxres  %s"  % (min_res, max_res))          
        
            if (max_res-min_res<0.75):
                nextIter = False
                break
                
        os.system('cp '  +self._getExtraPath('sharpenedMap_'+str(iteration)+'.mrc')+
                   ' '  +self._getExtraPath('sharpenedMap_last.mrc'))
        
        resFile = self.resolutionVolume.get().getFileName()        
        pathres=dirname(resFile)
        if  not exists(os.path.join(pathres,'mask_data.xmd')):     

            print ('\n====================\n' 
                   ' WARNING---This is not the ideal case because resolution map has been imported.'
                   ' The ideal case is to calculate it previously' 
                   ' in the same project using MonoRes.'
                   '\n====================\n')           
  
             
    def createOutputStep(self):
        volume=Volume()
        volume.setFileName(self._getExtraPath('sharpenedMap_last.mrc'))
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
        return ['Ramirez-Aportela 2018']

