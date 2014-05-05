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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************

import os
from os.path import basename
import xmipp

from pyworkflow.em.wizard import EmWizard
from django.shortcuts import render_to_response
from django.http import HttpResponse
from pyworkflow.web.app.em_wizard import *



def wiz_relion_bandpass(protocol, request):
    particles = protocol.inputParticles.get()
    
    highFreq = protocol.iniLowPassFilter.get()
    
    res = validateParticles(particles)
    
    if res is not 1:
        return HttpResponse(res)
    else:
        parts = getParticleSubset(particles.clone(),100)
        
        if len(parts) == 0:
            return HttpResponse("errorIterate")
        else:
            context = {'objects': parts,
                       'lowFreq': 0,
                       'highFreq': highFreq,
                       'decayFreq': 0,
                       'unit': UNIT_PIXEL
                       }
            
            context = wiz_base(request, context)
            
            return render_to_response('wizards/wiz_relion_bandpass.html', context)
    