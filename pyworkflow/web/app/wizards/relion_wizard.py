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

from pyworkflow.em.packages.relion.wizard import * 
from pyworkflow.web.app.em_wizard import *
from tools import *
from pyworkflow.web.app.views_base import base_wiz

#===============================================================================
# MASKS
#===============================================================================

class RelionPartMaskRadiusWeb(RelionPartMaskDiameterWizard):
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
            
            # getting dimension in Pixels
            xdim = getImageXdim(request, particles[0].text) 
            
            # Multiply by sampling rate to convert to angstroms
            xdim *= objs.getSamplingRate()
            
            params['value'] = validateMaskRadius(params['value'], xdim, radius=1)
            
            context = {'objects': particles,
                       'xdim': int(xdim/2),
                       'params': params}
        
            context = base_wiz(request, context)
            return render_to_response('wizards/wiz_particle_mask_radius.html', context)    



class RelionVolMaskRadiusWeb(RelionBackRadiusWizard):
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


#===============================================================================
# FILTERS
#===============================================================================

class RelionVolFilterWizardWeb(RelionVolFilterWizard):
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
            params['label'] = ['',params['label'],'']

            params['samplingRate'] = objs.getSamplingRate()
            params['unit'] = UNIT_ANGSTROM

            # Correct default values to angstroms
            params['value'][0] = params['samplingRate'] / params['value'][0]
            params['value'][2] = params['samplingRate'] / params['value'][2]

            itemDim,_,_ = objs.getDim()

            params['min'] = 2.*params['samplingRate']
            params['max'] = 2.*itemDim*params['samplingRate']

            params['value'] = setValueOnRange(params['value'], params['min'], params['max'])

            context = {'objects': volumes,
                       'params':params}

            context = base_wiz(request, context)

            return render_to_response('wizards/wiz_filter_volumes.html', context)


    