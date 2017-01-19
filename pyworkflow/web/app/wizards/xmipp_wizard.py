# **************************************************************************
# *
# * Authors:    Jose Gutierrez (jose.gutierrez@cnb.csic.es)
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

import os
from os.path import basename
import xmipp

from pyworkflow.em.wizard import EmWizard
from django.shortcuts import render_to_response
from django.http import HttpResponse

from pyworkflow.em.packages.xmipp3.wizard import * 
from pyworkflow.web.app.em_wizard import *
from tools import *
from pyworkflow.web.app.views_base import base_wiz
from pyworkflow.wizard import WEB_DJANGO

#===============================================================================
# DOWNSAMPLING
#===============================================================================

class XmippDownsamplingWeb(XmippDownsampleWizard):
    _environments = [WEB_DJANGO]
    
    def _run(self, protocol, request):
        params = self._getParameters(protocol)     
        objs = params['input'].get()
        
        res = validateSet(objs)
        
        if res is not 1:
            return HttpResponse(res)
        else:           
            
            context = {'objects': self._getMics(objs),
                       'params': params}
    
            context = base_wiz(request, context)
            return render_to_response('wizards/wiz_downsampling.html', context)

#===============================================================================
# CTFS
#===============================================================================

class XmippCTFWeb(XmippCTFWizard):
    _environments = [WEB_DJANGO]
    
    def _run(self, protocol, request):
        params = self._getParameters(protocol)     
        objs = params['input'].get()
        
        res = validateSet(objs)
        
        if res is not 1:
            return HttpResponse(res)
        else:           

            context = {'objects': self._getMics(objs),
                       'params': params}
        
            context = base_wiz(request, context)
            return render_to_response('wizards/wiz_ctf_downsampling.html', context)
        
#===============================================================================
# MASKS
#===============================================================================

class XmippParticleMaskRadiusWeb(XmippParticleMaskRadiusWizard):
    _environments = [WEB_DJANGO]
    
    def _run(self, protocol, request):
        params = self._getParameters(protocol)
        objs = params['input'].get()
        
        res = validateParticles(objs) 
        
        if res is not 1:
            return HttpResponse(res)
        else:
            particles = self._getParticles(objs)
            
            if len(particles) == 0:
                return HttpResponse("errorIterate")
            
            xdim = getImageXdim(request, particles[0].text)
            
            params['value'] = validateMaskRadius(params['value'], xdim, radius=1)
            
            context = {'objects': particles,
                       'xdim': int(xdim/2),
                       'params': params}
        
            context = base_wiz(request, context)
            return render_to_response('wizards/wiz_particle_mask_radius.html', context)    


class XmippParticleMaskRadiiWeb(XmippParticleMaskRadiiWizard):
    _environments = [WEB_DJANGO]
    
    def _run(self, protocol, request):
        params = self._getParameters(protocol)
        objs = params['input'].get()
        
        res = validateParticles(objs) 
        
        if res is not 1:
            return HttpResponse(res)
        else:
            particles = self._getParticles(objs)
            
            if len(particles) == 0:
                return HttpResponse("errorIterate")
            
            xdim = getImageXdim(request, particles[0].text)

            params['value'] = validateMaskRadius(params['value'], xdim, radius=2)               
            
            context = {'objects': particles,
                       'xdim': int(xdim/2),
                       'params': params}
        
            context = base_wiz(request, context)
            return render_to_response('wizards/wiz_particles_mask_radii.html', context)    


class XmippVolumeMaskRadiusWeb(XmippVolumeMaskRadiusWizard):
    _environments = [WEB_DJANGO]
    
    def _run(self, protocol, request):
        params = self._getParameters(protocol)
        objs = params['input'].get()
        
        res = validateSet(objs) 
        
        if res is not 1:
            return HttpResponse(res)
        else:
            vols = self._getVols(objs)
            
            if len(vols) == 0:
                return HttpResponse("errorIterate")
            
            xdim = getImageXdim(request, vols[0].getFileName())
            
            params['value'] = validateMaskRadius(params['value'], xdim, radius=1)
               
            context = {'objects': vols,
                       'xdim': int(xdim/2),
                       'params': params}
        
            context = base_wiz(request, context)
            return render_to_response('wizards/wiz_volume_mask_radius.html', context)    

#This classes are not used in any place and cause error
# class XmippVolumeRadiusWizardWeb(XmippVolumeRadiusWizard, XmippVolumeMaskRadiusWeb):
#     pass
# 
# class XmippVolumeRadiusWizardPWeb(XmippVolumeMaskRadiusProjMWizard, XmippVolumeMaskRadiusWeb):
#     pass

class XmippVolumeMaskRadiiWeb(XmippVolumeRadiiWizard):
    _environments = [WEB_DJANGO]

    def _run(self, protocol, request):
        params = self._getParameters(protocol)
        objs = params['input'].get()

        res = validateSet(objs)

        if res is not 1:
            return HttpResponse(res)
        else:
            vols = self._getVols(objs)

            if len(vols) == 0:
                return HttpResponse("errorIterate")

            xdim = getImageXdim(request, vols[0].getFileName())

            params['value'] = validateMaskRadius(params['value'], xdim, radius=2)

            context = {'objects': vols,
                       'xdim': xdim/2,
                       'params': params}

            context = base_wiz(request, context)
            return render_to_response('wizards/wiz_volumes_mask_radii.html', context)

class XmippVolumeMaskRadiiPWeb(XmippVolumeRadiiProjMWizard, XmippVolumeMaskRadiiWeb):
    pass

#===============================================================================
# FILTERS
#===============================================================================

class XmippFilterParticlesWeb(XmippFilterParticlesWizard):
    _environments = [WEB_DJANGO]

    def _run(self, protocol, request):
        params = self._getParameters(protocol)
        objs = params['input'].get()

        res = validateParticles(objs)

        if res is not 1:
            return HttpResponse(res)
        else:
            particles = self._getParticles(objs)

            if len(particles) == 0:
                return HttpResponse("errorIterate")

            params['value'] = proccessModeFilter(params['mode'], params['value'])

            params['samplingRate'] = particles[1].getSamplingRate()

            itemDim,_,_ = particles[1].getDim()

            if params['unit'] == UNIT_ANGSTROM:
                params['min'] = 2.*params['samplingRate']
                params['max'] = 2.*itemDim*params['samplingRate']
            else:
                params['min'] = 0.01
                params['max'] = 0.5

            params['value'] = setValueOnRange(params['value'], params['min'], params['max'])

            context = {'objects': particles,
                       'params':params}

            context = base_wiz(request, context)

            return render_to_response('wizards/wiz_filter_particles.html', context)


class XmippFilterVolumesWeb(XmippFilterVolumesWizard):
    _environments = [WEB_DJANGO]

    def _run(self, protocol, request):
        params = self._getParameters(protocol)
        objs = params['input'].get()

        res = validateSet(objs)

        if res is not 1:
            return HttpResponse(res)
        else:
            volumes = self._getVols(objs)

            if len(volumes) == 0:
                return HttpResponse("errorIterate")

            params['value'] = proccessModeFilter(params['mode'], params['value'])
            params['label'] = params['label']

            params['samplingRate'] = objs.getSamplingRate()

            itemDim,_,_ = objs.getDim()

            if params['unit'] == UNIT_ANGSTROM:
                params['min'] = 2.*params['samplingRate']
                params['max'] = 2.*itemDim*params['samplingRate']
            else:
                params['min'] = 0.01
                params['max'] = 0.5

            params['value'] = setValueOnRange(params['value'], params['min'], params['max'])

            context = {'objects': volumes,
                       'params':params}

            context = base_wiz(request, context)

            return render_to_response('wizards/wiz_filter_volumes.html', context)


class XmippGaussianParticlesWeb(XmippGaussianParticlesWizard):
    _environments = [WEB_DJANGO]

    def _run(self, protocol, request):
        params = self._getParameters(protocol)
        objs = params['input'].get()

        res = validateParticles(objs)

        if res is not 1:
            return HttpResponse(res)
        else:
            particles = self._getParticles(objs)

            if len(particles) == 0:
                return HttpResponse("errorIterate")

            context = {'objects': particles,
                       'params':params}

            context = base_wiz(request, context)

            return render_to_response('wizards/wiz_gaussian_particle.html', context)


class XmippGaussianVolumesWeb(XmippGaussianVolumesWizard):
    _environments = [WEB_DJANGO]

    def _run(self, protocol, request):
        params = self._getParameters(protocol)
        objs = params['input'].get()

        res = validateSet(objs)

        if res is not 1:
            return HttpResponse(res)
        else:
            volumes = self._getVols(objs)

            if len(volumes) == 0:
                return HttpResponse("errorIterate")

            context = {'objects': volumes,
                       'params':params}

            context = base_wiz(request, context)

            return render_to_response('wizards/wiz_gaussian_vol.html', context)


class XmippBoxSizeWizardWeb(XmippBoxSizeWizard):
    _environments = [WEB_DJANGO]

    def _run(self, protocol, request):
        boxSize = protocol.getBoxSize()
        context = '{boxSize='+ str(boxSize) +'}'
        return HttpResponse('auto_wizard:'+context)


class XmippParticleConsensusRadiusWizardWeb(XmippParticleConsensusRadiusWizard):
    _environments = [WEB_DJANGO]

    def _run(self, protocol, request):
        radius = self._getRadius(protocol)
        context = '{consensusRadius='+ str(radius) +'}'
        return HttpResponse('auto_wizard:'+context)

class XmippCL2DNumberOfClassesWizardWeb(XmippCL2DNumberOfClassesWizard):
    _environments = [WEB_DJANGO]

    def _run(self, protocol, request):
        numClasses = self._getNumberOfClasses(protocol)
        context = '{numberOfClasses='+ str(numClasses) +'}'
        return HttpResponse('auto_wizard:'+context)
