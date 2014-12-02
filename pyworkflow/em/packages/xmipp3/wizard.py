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

from protocol_ctf_micrographs import XmippProtCTFMicrographs
from protocol_projmatch import XmippProtProjMatch 
from protocol_preprocess_micrographs import XmippProtPreprocessMicrographs
from protocol_preprocess import XmippProtPreprocessParticles, XmippProtPreprocessVolumes
from protocol_extract_particles import XmippProtExtractParticles
from protocol_extract_particles_pairs import XmippProtExtractParticlesPairs
from protocol_preprocess import (XmippProtFilterParticles, XmippProtFilterVolumes,
                                 XmippProtMaskParticles, XmippProtMaskVolumes)
from protocol_align_volume import XmippProtAlignVolume
from protocol_cl2d import XmippProtCL2D

from pyworkflow.em.wizard import *

#===============================================================================
# DOWNSAMPLING
#===============================================================================

class XmippDownsampleWizard(DownsampleWizard):
    _targets = [(XmippProtPreprocessMicrographs, ['downFactor'])]
    
    def _getParameters(self, protocol):
        
        label, value = self._getInputProtocol(self._targets, protocol)
        
        protParams = {}
        protParams['input']= protocol.inputMicrographs
        protParams['label']= label
        protParams['value']= value
        
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
    _targets = [(XmippProtCTFMicrographs, ['ctfDownFactor', 'lowRes', 'highRes'])]
    
    def _getParameters(self, protocol):
        
        label, value = self._getInputProtocol(self._targets, protocol)
        
        protParams = {}
        protParams['input']= protocol.inputMicrographs
        protParams['label']= label
        protParams['value']= value
        return protParams
    
    def _getProvider(self, protocol):
        _objs = self._getParameters(protocol)['input']
        return CtfWizard._getListProvider(self, _objs)
        
    def show(self, form):
        protocol = form.protocol
        params = self._getParameters(protocol)
        _value = params['value']
        _label = params['label']
        
#        form.setParamFromVar('inputMicrographs') # update selected input micrographs
        provider = self._getProvider(protocol)
        
        if provider is not None:
            args = {'unit': UNIT_PIXEL,
                    'downsample': _value[0],
                    'lf': _value[1],
                    'hf': _value[2]
                    }
            d = CtfDownsampleDialog(form.root, provider, **args)

            if d.resultYes():
                form.setVar(_label[0], d.getDownsample())
                form.setVar(_label[1], d.getLowFreq())
                form.setVar(_label[2], d.getHighFreq())
        else:
            dialog.showWarning("Empty input", "Select elements first", form.root)    
    
    @classmethod    
    def getView(self):
        return "wiz_ctf_downsampling"  

#===============================================================================
# BOXSIZE
#===============================================================================
class XmippBoxSizeWizard(Wizard):
    _targets = [(XmippProtExtractParticles, ['boxSize'])]

    def _getBoxSize(self, protocol):

        # Get input coordinates from protocol and if they have not value exit
        inputCoords = protocol.getCoords()
        if  inputCoords is None:
            return 0

        # Get boxSize from coordinates and sampling input from associated micrographs
        boxSize = inputCoords.getBoxSize()
        samplingInput = inputCoords.getMicrographs().getSamplingRate()

        # If downsampling type same as picking sampling does not change
        if protocol.downsampleType.get() == SAME_AS_PICKING:
            samplingFinal = samplingInput
        else:
            # In case is not same as picking get input micrographs and its sampling rate
            if issubclass(protocol.getClass(), XmippProtExtractParticles):
                inputMics = protocol.inputMicrographs.get()
            if issubclass(protocol.getClass(), XmippProtExtractParticlesPairs):
                inputMics = inputCoords.getMicrographs()

            samplingMics = inputMics.getSamplingRate()

            if protocol.downsampleType.get() == ORIGINAL:
                # If 'original' get sampling rate from original micrographs
                samplingFinal = samplingMics
            else:
                # IF 'other' multiply the original sampling rate by the factor provided
                samplingFinal = samplingMics*protocol.downFactor.get()

        downFactor = samplingFinal/samplingInput

        return int(boxSize/downFactor)



        return boxSize


    def show(self, form):
        form.setVar('boxSize', self._getBoxSize(form.protocol))

#===============================================================================
# NUMBER OF CLASSES
#===============================================================================
class XmippCL2DNumberOfClassesWizard(Wizard):
    _targets = [(XmippProtCL2D, ['numberOfClasses'])]

    def _getNumberOfClasses(self, protocol):

        numberOfClasses = 64

        if protocol.inputParticles.hasValue():
            from protocol_cl2d import IMAGES_PER_CLASS
            numberOfClasses = int(protocol.inputParticles.get().getSize()/IMAGES_PER_CLASS)

        return numberOfClasses


    def show(self, form):
        form.setVar('numberOfClasses', self._getNumberOfClasses(form.protocol))

#===============================================================================
# MASKS 
#===============================================================================

class XmippParticleMaskRadiusWizard(ParticleMaskRadiusWizard):
    _targets = [(XmippProtMaskParticles, ['radius']),
                (XmippProtPreprocessParticles, ['backRadius'])]
    
    def _getParameters(self, protocol):
        
        label, value = self._getInputProtocol(self._targets, protocol)
        
        protParams = {}
        protParams['input']= protocol.inputParticles
        protParams['label']= label
        protParams['value']= value
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
        
        label, value = self._getInputProtocol(self._targets, protocol)
        
        protParams = {}
        protParams['input']= protocol.inputParticles
        protParams['label']= label
        protParams['value']= value
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
    _targets = [(XmippProtMaskVolumes, ['radius']), 
                (XmippProtAlignVolume, ['maskRadius']),
                (XmippProtPreprocessVolumes, ['backRadius'])]
    
    def _getParameters(self, protocol):
        
        label, value = self._getInputProtocol(self._targets, protocol)
        
        protParams = {}
        protParams['input']= protocol.inputVolumes
        protParams['label']= label
        protParams['value']= value
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
        
        label, value = self._getInputProtocol(self._targets, protocol)
        
        protParams = {}
        protParams['input']= protocol.inputVolumes
        protParams['label']= label
        protParams['value']= value
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
        
        label, value = self._getInputProtocol(self._targets, protocol)
        
        protParams = {}
        protParams['input']= protocol.inputParticles
        protParams['label']= label
        protParams['value']= value
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
        
        label, value = self._getInputProtocol(self._targets, protocol)
        
        protParams = {}
        protParams['input']= protocol.inputVolumes
        protParams['label']= label
        protParams['value']= value
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
        
        label, value = self._getInputProtocol(self._targets, protocol)
        
        protParams = {}
        protParams['input']= protocol.inputParticles
        protParams['label']= label
        protParams['value']= value
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
        
        label, value = self._getInputProtocol(self._targets, protocol)
        
        protParams = {}
        protParams['input']= protocol.inputVolumes
        protParams['label']= label
        protParams['value']= value
        return protParams  
    
    def _getProvider(self, protocol):
        _objs = self._getParameters(protocol)['input']
        return GaussianVolumesWizard._getListProvider(self, _objs)
    
    def show(self, form):
        params = self._getParameters(form.protocol)
        _value = params['value']
        _label = params['label']
        GaussianVolumesWizard.show(self, form, _value, _label, UNIT_PIXEL_FOURIER)
    

    
