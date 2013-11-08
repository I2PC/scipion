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
# *  e-mail address 'jose.gutierrez@cnb.csic.es'
# *
# **************************************************************************

import os
from os.path import basename
import xmipp
from views_util import * 
from views_protocol import updateProtocolParams
from pyworkflow.manager import Manager
from pyworkflow.project import Project
from django.shortcuts import render_to_response
from pyworkflow.gui import getImage, getPILImage
from django.http import HttpResponse

def wizard(request):
    # Get the Wizard Name
    requestDict = getattr(request, "POST")
    functionName = requestDict.get("wizName")
    function = globals().get(functionName, None)
    
    print "======================= in wizard: " + functionName
    
    # Get the protocol object
    project, protocol = loadProtocolProject(request)
    updateProtocolParams(request, protocol, project)
    
    if function is None:
        pass  # redirect to error: wizard not found
    elif not callable(function):
        pass  # redirect to error: name is not a function
    else:
        return function(protocol, request)

def wiz_downsampling(protocol, request):
    micrographs = protocol.inputMicrographs.get()
    res = validateSet(micrographs)
    
    if res is not 1:
        return HttpResponse(res)
    else:
        mics = [mic for mic in micrographs]
        for m in mics:
            m.basename = basename(m.getFileName())
            
        context = {'objects': mics,
                   'downFactor': protocol.downFactor.get(),
                   }
        
        return render_to_response('wiz_downsampling.html', context)

def wiz_ctf(protocol, request):
    micrographs = protocol.inputMicrographs.get()
    
    res = validateSet(micrographs)
    
    if res is not 1:
        return HttpResponse(res)
    else:
        mics = [mic for mic in micrographs]
        for m in mics:
            m.basename = basename(m.getFileName())    
        
        context = {'objects': mics,
                   'raphael':getResourceJs('raphael'),
                   'high_res' : protocol.highRes.get(),
                   'low_res': protocol.lowRes.get()
                   }
        
        return render_to_response('wiz_ctf.html', context)

def wiz_particle_mask(protocol, request):
    particles = protocol.inputParticles.get()
    res = validateParticles(particles) 
    
    if res is not 1:
        return HttpResponse(res)
    else:
        parts = getParticleSubset(particles,100)
        
        if len(parts) == 0:
            return HttpResponse("errorIterate");
        else:
            xdim = getImageXdim(request, parts[0].text)
    
            mask_radius = protocol.maskRadius.get()
             
            if mask_radius > xdim :
                mask_radius = xdim
            elif mask_radius == -1 :
                mask_radius = xdim/2
                
            context = {'objects': parts,
                       'raphael': getResourceJs('raphael'),
                       'maskRadius': mask_radius,
                       'xdim':xdim
                       }
            
            return render_to_response('wiz_particle_mask.html', context)

def wiz_particle_mask_radii(protocol, request):
    particles = protocol.inputParticles.get()
    res = validateParticles(particles) 
    
    if res is not 1:
        return HttpResponse(res)
    else:
        parts = getParticleSubset(particles,100)
        
        if len(parts) == 0:
            return HttpResponse("errorIterate");
        else:
            xdim = getImageXdim(request, parts[0].text)
            inner_radius = protocol.innerRadius.get()
            outer_radius = protocol.outerRadius.get()
                
            context = {'objects': parts,
                       'raphael': getResourceJs('raphael'),
                       'innerRadius': inner_radius,
                       'outerRadius': outer_radius,
                       'xdim':xdim
                       }
            
            return render_to_response('wiz_particle_mask_radii.html', context)

def wiz_volume_mask(protocol, request):
    volumes = protocol.input3DReferences.get()
    res = validateSet(volumes)
    
    if res is not 1:
        return HttpResponse(res)
    else:
        vols = [vol for vol in volumes]
        
        for v in vols:
            v.basename = basename(v.getFileName())
                    
        xdim = getImageXdim(request, vols[0].getFileName())
    
        mask_radius = protocol.maskRadius.get()
         
        if mask_radius > xdim :
            mask_radius = xdim
        elif mask_radius == -1 :
            mask_radius = xdim/2
            
        context = {'objects': vols,
                   'raphael': getResourceJs('raphael'),
                   'maskRadius': mask_radius,
                   'xdim': xdim
                   }
        
        return render_to_response('wiz_volume_mask.html', context)

def wiz_volume_mask_radii(protocol, request):
    volumes = protocol.input3DReferences.get()
    res = validateSet(volumes)
    
    if res is not 1:
        return HttpResponse(res)
    else:
        vols = [vol for vol in volumes]
        
        for v in vols:
            v.basename = basename(v.getFileName())
                    
        xdim = getImageXdim(request, vols[0].getFileName())
    
        inner_radius = protocol.innerRadius.get()
        outer_radius = protocol.outerRadius.get()
            
        context = {'objects': vols,
                   'raphael': getResourceJs('raphael'),
                   'innerRadius': inner_radius,
                   'outerRadius': outer_radius,
                   'xdim': xdim
                   }
        
        return render_to_response('wiz_volume_mask_radii.html', context)

def wiz_filter_spider(protocol, request):
    particles = protocol.inputParticles.get()
    res = validateParticles(particles) 
    
    if res is not 1:
        return HttpResponse(res)
    else:
        parts = getParticleSubset(particles,100)
        
        if len(parts) == 0:
            return HttpResponse("errorIterate");
        else:
            context = {'objects': parts,
                       'raphael': getResourceJs('raphael'),
                       'protocol': protocol
                       }
            
            return render_to_response('wiz_filter_spider.html', context)

def wiz_bandpass(protocol, request):
    particles = protocol.inputParticles.get()
    res = validateParticles(particles) 
    
    if res is not 1:
        return HttpResponse(res)
    else:
        parts = getParticleSubset(particles,100)
        
        if len(parts) == 0:
            return HttpResponse("errorIterate");
        else:
            context = {'objects': parts,
                       'lowFreq': protocol.lowFreq.get(),
                       'highFreq': protocol.highFreq.get(),
                       'decayFreq': protocol.freqDecay.get()
                       }
            
            return render_to_response('wiz_bandpass.html', context)

def wiz_gaussian(protocol, request):
    particles = protocol.inputParticles.get()
    res = validateParticles(particles) 
    
    if res is not 1:
        return HttpResponse(res)
    else:
        parts = getParticleSubset(particles,100)
        
        if len(parts) == 0:
            return HttpResponse("errorIterate");
        else:
            context = {'objects': parts,
                       'freqSigma': protocol.freqSigma.get()
                       }
            
            return render_to_response('wiz_gaussian.html', context)

def getParticleSubset(particles, num):
    """
    Method to prepare the particles to use a wizard
    """
    particleList = []
    for i, particle in enumerate(particles):
        if i == num: # Only load up to NUM particles
            break
        particle.text = particle.getFileName()
        particle.basename = basename(particle.text)
        
        index = particle.getIndex()
        if index:
            particle.text = "%03d@%s" % (index, particle.text)
            particle.basename = "%03d@%s" % (index, particle.basename)
        
        particleList.append(particle)

    return particleList

def validateSet(setOf):
    """
    Validation for a set of micrographs
    """
    if setOf is None:
        res = "errorInput"
    else:
        res = 1
    return res

def validateParticles(particles):
    """
    Validation for a set of particles
    """
    if particles is None:
        res = "errorInput"
    elif particles.getSize() == 0:
        res = "errorEmpty"
    else:
        res = 1
#        else:
#            res = parts
    return res

def get_image_psd(request):
    """
    Function to get the computing psd image
    """
    imagePath = request.GET.get('image', None)
    downsample = request.GET.get('downsample', None)
    dim = request.GET.get('dim', None)
    
    # create a xmipp image empty
    imgXmipp = xmipp.Image()
    
    # compute the PSD image
    xmipp.fastEstimateEnhancedPSD(imgXmipp, str(imagePath), float(downsample), int(dim), 2)
        
    # from PIL import Image
    img = getPILImage(imgXmipp, dim)
        
    response = HttpResponse(mimetype="image/png")
    img.save(response, "PNG")
    return response

def get_image_bandpass(request):
    """
    Function to get the computing image with a fourier filter applied
    """
    imagePath = request.GET.get('image', None)
    lowFreq = request.GET.get('lowFreq', None)
    highFreq = request.GET.get('highFreq', None)
    decay = request.GET.get('decayFreq', None)
    dim = request.GET.get('dim', None)
    
    # create a xmipp image empty
    imgXmipp = xmipp.Image()
    
    # compute the Fourier Filter in the image
    xmipp.bandPassFilter(imgXmipp, str(imagePath), float(lowFreq), float(highFreq), float(decay), int(dim))
    
    # from PIL import Image
    img = getPILImage(imgXmipp, dim)
        
    response = HttpResponse(mimetype="image/png")
    img.save(response, "PNG")
    
    return response

def get_image_gaussian(request):
    """
    Function to get the computing image with a gaussian filter applied
    """
    imagePath = request.GET.get('image', None)
    freqSigma = request.GET.get('sigmaFreq', None)
    dim = request.GET.get('dim', None)
    
    # create a xmipp image empty
    imgXmipp = xmipp.Image()
    
    # compute the Gaussian Filter in the image
    xmipp.gaussianFilter(imgXmipp, str(imagePath), float(freqSigma), int(dim))
        
    # from PIL import Image
    img = getPILImage(imgXmipp, dim)
        
    response = HttpResponse(mimetype="image/png")
    img.save(response, "PNG")
    return response

def get_image_filter_spider(request):
    """
    Function to get the computing image with a spider filter applied
    """
    imagePath = request.GET.get('image', None)
    dim = request.GET.get('dim', None)
    mode = request.GET.get('mode', None)
    
    # Instance of the protocol necessary
    
    # check values
    if mode < 2:
        radius = request.GET.get('radius', None)
    else: 
        highFreq = request.GET.get('highFreq', None)
        lowFreq = request.GET.get('lowFreq', None)
        if mode == 2:
            temperature = request.GET.get('temperature', None)    

    # Get output image and update filtered image
    imgXmipp = xmipp.Image()
    
    
    # compute the Spider Filter in the image
    #
    #
    #
    
    
    
    # from PIL import Image
    img = getPILImage(imgXmipp, dim)
    
    response = HttpResponse(mimetype="image/png")
    img.save(response, "PNG")
    return response
    
    
