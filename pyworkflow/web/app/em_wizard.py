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
        index = particle.getIndex()
        text = particle.getFileName()
        particle.basename = basename(text)
        if index:
            particle.text = "%03d@%s" % (index, text)
            particle.basename = "%03d@%s" % (index, basename(text))
        particleList.append(particle)
        
    return particleList

"""
Validation for a set of micrographs
"""
def validateSet(setOf):
    if setOf is None:
        res = "errorInput"
    else:
        res = 1
    return res

"""
Validation for a set of particles
"""
def validateParticles(particles):
    if particles is None:
        res = "errorInput"
    elif particles.getSize() == 0:
        res = "errorEmpty"
    else:
        res = 1
#        else:
#            res = parts
    return res

"""
Function to get the computing psd image
"""
def get_image_psd(request):
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

"""
Function to get the computing image with a fourier filter applied
"""
def get_image_bandpass(request):
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

"""
Function to get the computing image with a gaussian filter applied
"""
def get_image_gaussian(request):
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
