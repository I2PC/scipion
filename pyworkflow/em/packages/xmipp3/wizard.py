# **************************************************************************
# *
# * Authors:     Jose Gutierrez (jose.gutierrez@cnb.csic.es)
# *
# * Unidad de Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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

#===============================================================================
# DOWNSAMPLING
#===============================================================================

class XmippDownsampleWizard(DownsampleWizard):
    _targets = [(XmippProtPreprocessMicrographs, ['downFactor'])]
    
    def _getParameters(self, protocol):
        protParams = {}
        protParams['input']= protocol.inputMicrographs
        protParams['label']= "downFactor"
        protParams['value']= protocol.downFactor.get()
        
        return protParams

    def _getProvider(self, protocol):
        _objs = self._getParameters(protocol)['input']
        return DownsampleWizard._getListProvider(self, _objs)
    
    def show(self, form):
        params = self._getParameters(form.protocol)
        _value = params['value']
        _label = params['label']
        DownsampleWizard.show(self, form, _value, _label, UNIT_PIXEL)

#===============================================================================
# CTFS
#===============================================================================
        
class XmippCTFWizard(CtfWizard):
    _targets = [(ProtCTFMicrographs, ['lowRes', 'highRes'])]
    
    def _getParameters(self, protocol):
        protParams = {}
        protParams['input']= protocol.inputMicrographs
        protParams['label']= ["lowRes", "highRes"]
        protParams['value']= [protocol.lowRes.get(), protocol.highRes.get()]
        return protParams
    
    def _getProvider(self, protocol):
        _objs = self._getParameters(protocol)['input']
        return CtfWizard._getListProvider(self, _objs)
    
    def show(self, form):
        params = self._getParameters(form.protocol)
        _value = params['value']
        _label = params['label']
        CtfWizard.show(self, form, _value, _label, UNIT_PIXEL)

#===============================================================================
# MASKS 
#===============================================================================

class XmippParticleMaskRadiusWizard(ParticleMaskRadiusWizard):
    _targets = [(XmippProtMaskParticles, ['radius'])]
    
    def _getParameters(self, protocol):
        protParams = {}
        protParams['input']= protocol.inputParticles
        protParams['label']= "radius"
        protParams['value']= protocol.radius.get()
        return protParams
    
    def _getProvider(self, protocol):
        _objs = self._getParameters(protocol)['input']
        return ParticleMaskRadiusWizard._getListProvider(self, _objs)
    
    def show(self, form):
        params = self._getParameters(form.protocol)
        _value = params['value']
        _label = params['label']
        ParticleMaskRadiusWizard.show(self, form, _value, _label, UNIT_PIXEL)
    
    
class XmippParticleMaskRadiiWizard(ParticlesMaskRadiiWizard):
    _targets = [(XmippProtMaskParticles, ['innerRadius', 'outerRadius'])]
    
    def _getParameters(self, protocol):
        protParams = {}
        protParams['input']= protocol.inputParticles
        protParams['label']= ["innerRadius", "outerRadius"]
        protParams['value']= [protocol.innerRadius.get(), protocol.outerRadius.get()]
        return protParams  
    
    def _getProvider(self, protocol):
        _objs = self._getParameters(protocol)['input']
        return ParticlesMaskRadiiWizard._getListProvider(self, _objs)
    
    def show(self, form):
        params = self._getParameters(form.protocol)
        _value = params['value']
        _label = params['label']
        ParticlesMaskRadiiWizard.show(self, form, _value, _label, UNIT_PIXEL)
    

class XmippVolumeMaskRadiusWizard(VolumeMaskRadiusWizard):
    _targets = [(XmippProtMaskVolumes, ['radius'])]
    
    def _getParameters(self, protocol):
        protParams = {}
        protParams['input']= protocol.inputVolumes
        protParams['label']= "radius"
        protParams['value']= protocol.radius.get()
        return protParams
    
    def _getProvider(self, protocol):
        _objs = self._getParameters(protocol)['input']
        return VolumeMaskRadiusWizard._getListProvider(self, _objs)
    
    def show(self, form):
        params = self._getParameters(form.protocol)
        _value = params['value']
        _label = params['label']
        VolumeMaskRadiusWizard.show(self, form, _value, _label, UNIT_PIXEL)
    
 

class XmippVolumeRadiiWizard(VolumeMaskRadiiWizard):
    _targets = [(XmippProtMaskVolumes, ['innerRadius', 'outerRadius'])]
    
    def _getParameters(self, protocol):
        protParams = {}
        protParams['input']= protocol.inputVolumes
        protParams['label']= ["innerRadius", "outerRadius"]
        protParams['value']= [protocol.innerRadius.get(), protocol.outerRadius.get()]
        return protParams  
    
    def _getProvider(self, protocol):
        _objs = self._getParameters(protocol)['input'] 
        return VolumeMaskRadiiWizard._getListProvider(self, _objs)
    
    def show(self, form):
        params = self._getParameters(form.protocol)
        _value = params['value']
        _label = params['label']
        VolumeMaskRadiiWizard.show(self, form, _value, _label, UNIT_PIXEL)
    

#===============================================================================
#  FILTERS
#===============================================================================
        
class XmippFilterParticlesWizard(FilterParticlesWizard):   
    _targets = [(XmippProtFilterParticles, ['lowFreq', 'highFreq', 'freqDecay'])]
    
    def _getParameters(self, protocol):
        protParams = {}
        protParams['input']= protocol.inputParticles
        protParams['label']= ["lowFreq", "highFreq","freqDecay"]
        protParams['value']= [protocol.lowFreq.get(), protocol.highFreq.get(), protocol.freqDecay.get()]
        protParams['mode'] = protocol.fourierMode.get()
        return protParams  
    
    def _getProvider(self, protocol):
        _objs = self._getParameters(protocol)['input']
        return FilterParticlesWizard._getListProvider(self, _objs)
    
    def show(self, form):
        params = self._getParameters(form.protocol)
        _value = params['value']
        _label = params['label']
        _mode = params['mode']
        FilterParticlesWizard.show(self, form, _value, _label, _mode, UNIT_PIXEL_FOURIER)
    

    
class XmippFilterVolumesWizard(FilterVolumesWizard):   
    _targets = [(XmippProtFilterVolumes, ['lowFreq', 'highFreq', 'freqDecay'])]
    
    def _getParameters(self, protocol):
        protParams = {}
        protParams['input']= protocol.inputVolumes
        protParams['label']= ["lowFreq", "highFreq","freqDecay"]
        protParams['value']= [protocol.lowFreq.get(), protocol.highFreq.get(), protocol.freqDecay.get()]
        protParams['mode'] = protocol.fourierMode.get()
        return protParams  
    
    def _getProvider(self, protocol):
        _objs = self._getParameters(protocol)['input']
        return FilterVolumesWizard._getListProvider(self, _objs)
    
    def show(self, form):
        params = self._getParameters(form.protocol)
        _value = params['value']
        _label = params['label']
        _mode = params['mode']
        FilterVolumesWizard.show(self, form, _value, _label, _mode, UNIT_PIXEL_FOURIER)



class XmippGaussianParticlesWizard(GaussianParticlesWizard):
    _targets = [(XmippProtFilterParticles, ['freqSigma'])]
    
    def _getParameters(self, protocol):
        protParams = {}
        protParams['input']= protocol.inputParticles
        protParams['label']= "freqSigma"
        protParams['value']= protocol.freqSigma.get()
        return protParams  
    
    def _getProvider(self, protocol):
        _objs = self._getParameters(protocol)['input']
        return GaussianParticlesWizard._getListProvider(self, _objs)
    
    def show(self, form):
        params = self._getParameters(form.protocol)
        _value = params['value']
        _label = params['label']
        GaussianParticlesWizard.show(self, form, _value, _label, UNIT_PIXEL_FOURIER)



class XmippGaussianVolumesWizard(GaussianVolumesWizard):
    _targets = [(XmippProtFilterVolumes, ['freqSigma'])]
    
    def _getParameters(self, protocol):
        protParams = {}
        protParams['input']= protocol.inputVolumes
        protParams['label']= "freqSigma"
        protParams['value']= protocol.freqSigma.get()
        return protParams  
    
    def _getProvider(self, protocol):
        _objs = self._getParameters(protocol)['input']
        return GaussianVolumesWizard._getListProvider(self, _objs)
    
    def show(self, form):
        params = self._getParameters(form.protocol)
        _value = params['value']
        _label = params['label']
        GaussianVolumesWizard.show(self, form, _value, _label, UNIT_PIXEL_FOURIER)
    

    
