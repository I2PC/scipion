# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              Jose Gutierrez (jose.gutierrez@cnb.csic.es)
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
This module implement some wizards
"""

import os
import Tkinter as tk
import ttk

import xmipp

from pyworkflow.em.constants import *
from constants import *

from pyworkflow.em import SetOfImages, SetOfMicrographs, Volume, ProtCTFMicrographs
from protocol_projmatch import XmippProtProjMatch 
from protocol_preprocess_micrographs import XmippProtPreprocessMicrographs
from protocol_filter import XmippProtFilterParticles, XmippProtFilterVolumes
from protocol_mask import XmippProtMaskParticles, XmippProtMaskVolumes

from pyworkflow.em.wizard import * 


class XmippDownsampleWizard(downsampleWizard):
    _targets = [(XmippProtPreprocessMicrographs, ['downFactor'])]
    
    def _getParameters(self, protocol):
        protParams = {}
        protParams['input']= protocol.inputMicrographs
        protParams['label']= "downFactor"
        protParams['value']= protocol.downFactor.get()
        
        return protParams

    def _getProvider(self, protocol):
        _objs = self._getParameters(protocol)['input']
        return downsampleWizard._getListProvider(self, _objs)
    
    def show(self, form):
        params = self._getParameters(form.protocol)
        _value = params['value']
        _label = params['label']
        downsampleWizard.show(self, form, _value, _label, UNIT_PIXEL)
        
    @classmethod
    def getView(self):
        return "wiz_xmipp_downsampling"
        
    
class XmippCTFWizard(ctfWizard):
    _targets = [(ProtCTFMicrographs, ['lowRes', 'highRes'])]
    
    def _getParameters(self, protocol):
        protParams = {}
        protParams['input']= protocol.inputMicrographs
        protParams['label']= ["lowRes", "highRes"]
        protParams['value']= [protocol.lowRes.get(), protocol.highRes.get()]
        return protParams
    
    def _getProvider(self, protocol):
        _objs = self._getParameters(protocol)['input']
        return ctfWizard._getListProvider(self, _objs)
    
    def show(self, form):
        params = self._getParameters(form.protocol)
        _value = params['value']
        _label = params['label']
        ctfWizard.show(self, form, _value, _label, UNIT_PIXEL)
    
    @classmethod
    def getView(self):
        return "wiz_xmipp_ctf"
        

class XmippParticleMaskRadiusWizard(particleMaskRadiusWizard):
    _targets = [(XmippProtMaskParticles, ['radius'])]
    
    def _getProvider(self, protocol):
        _objs = protocol.inputParticles
        return particleMaskRadiusWizard._getListProvider(self, _objs)
    
    def show(self, form):
        _value = form.protocol.radius.get()
        _label = "radius"
        particleMaskRadiusWizard.show(self, form, _value, _label, UNIT_PIXEL)
    
    @classmethod
    def getView(self):
        return "wiz_xmipp_particle_mask_radius"
    
class XmippVolumeMaskRadiusWizard(volumeMaskRadiusWizard):
    _targets = [(XmippProtMaskVolumes, ['radius'])]
      
    def _getProvider(self, protocol):
        _objs = protocol.inputVolumes
        return volumeMaskRadiusWizard._getListProvider(self, _objs)
    
    def show(self, form):
        _label = "radius"
        _value = form.protocol.radius.get()
        volumeMaskRadiusWizard.show(self, form, _value, _label, UNIT_PIXEL)
 
    @classmethod
    def getView(self):
        return "wiz_xmipp_volume_mask_radius"
 
class XmippParticleMaskRadiiWizard(particlesMaskRadiiWizard):
    _targets = [(XmippProtMaskParticles, ['innerRadius', 'outerRadius'])]
    
    def _getProvider(self, protocol):
        _objs = protocol.inputParticles
        return particlesMaskRadiiWizard._getListProvider(self, _objs)
    
    def show(self, form):
        _value = [form.protocol.innerRadius.get(), form.protocol.outerRadius.get()]
        _label = ["innerRadius", "outerRadius"]
        particlesMaskRadiiWizard.show(self, form, _value, _label, UNIT_PIXEL)
    
    @classmethod
    def getView(self):
        return "wiz_xmipp_particle_mask_radii"

class XmippVolumeRadiiWizard(volumeMaskRadiiWizard):
    _targets = [(XmippProtMaskVolumes, ['innerRadius', 'outerRadius'])]
 
    def _getProvider(self, protocol):
        _objs = protocol.inputVolumes  
        return volumeMaskRadiiWizard._getListProvider(self, _objs)
    
    def show(self, form):
        _value = [form.protocol.innerRadius.get(), form.protocol.outerRadius.get()]
        _label = ["innerRadius", "outerRadius"]
        volumeMaskRadiiWizard.show(self, form, _value, _label, UNIT_PIXEL)
    
    @classmethod
    def getView(self):
        return "wiz_xmipp_volume_mask_radii"
        
class XmippFilterParticlesWizard(filterParticlesWizard):   
    _targets = [(XmippProtFilterParticles, ['lowFreq', 'highFreq', 'freqDecay'])]
    
    def _getProvider(self, protocol):
        _objs = protocol.inputParticles
        return filterParticlesWizard._getListProvider(self, _objs)
    
    def show(self, form):
        _value = [form.protocol.lowFreq.get(), 
                  form.protocol.highFreq.get(), 
                  form.protocol.freqDecay.get()]
        _label = ["lowFreq", 
                  "highFreq", 
                  "freqDecay"]
        _mode = form.protocol.fourierMode.get()
        
        filterParticlesWizard.show(self, form, _value, _label, _mode, UNIT_PIXEL_FOURIER)
    
    @classmethod
    def getView(self):
        return "wiz_xmipp_filter_particle"

class XmippGaussianParticlesWizard(gaussianParticlesWizard):
    _targets = [(XmippProtFilterParticles, ['freqSigma'])]
    
    def _getProvider(self, protocol):
        _objs = protocol.inputParticles
        return gaussianParticlesWizard._getListProvider(self, _objs)
    
    def show(self, form):
        _value = form.protocol.freqSigma.get()
        _label = "freqSigma"
        gaussianParticlesWizard.show(self, form, _value, _label, UNIT_PIXEL_FOURIER)
    
    @classmethod
    def getView(self):
        return "wiz_xmipp_gaussian_particle"
    
class XmippFilterVolumesWizard(filterVolumesWizard):   
    _targets = [(XmippProtFilterVolumes, ['lowFreq', 'highFreq', 'freqDecay'])]
    
    def _getProvider(self, protocol):
        _objs = protocol.inputVolumes    
        return filterVolumesWizard._getListProvider(self, _objs)
    
    def show(self, form):
        _value = [form.protocol.lowFreq.get(), 
                  form.protocol.highFreq.get(), 
                  form.protocol.freqDecay.get()]
        _label = ["lowFreq", 
                  "highFreq", 
                  "freqDecay"]
        _mode = form.protocol.fourierMode.get()
        
        filterVolumesWizard.show(self, form, _value, _label, _mode, UNIT_PIXEL_FOURIER)
     
    @classmethod
    def getView(self):
        return "wiz_xmipp_filter_volume"
     
class XmippGaussianVolumesWizard(gaussianVolumesWizard):
    _targets = [(XmippProtFilterVolumes, ['freqSigma'])]
    
    def _getProvider(self, protocol):
        _objs = protocol.inputVolumes
        return gaussianVolumesWizard._getListProvider(self, _objs)
    
    def show(self, form):
        _value = form.protocol.freqSigma.get()
        _label = "freqSigma"
        gaussianVolumesWizard.show(self, form, _value, _label, UNIT_PIXEL_FOURIER)

    @classmethod
    def getView(self):
        return "wiz_xmipp_gaussian_volume"
  
