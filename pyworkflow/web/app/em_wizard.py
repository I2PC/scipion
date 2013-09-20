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
    mics = [mic for mic in protocol.inputMicrographs.get()]
    for m in mics:
        m.basename = basename(m.getFileName())
        
    context = {'objects': mics,
               'downFactor': protocol.downFactor.get(),
               }
    
    return render_to_response('wiz_downsampling.html', context)

def wiz_ctf(protocol, request):
    mics = [mic for mic in protocol.inputMicrographs.get()]
    for m in mics:
        m.basename = basename(m.getFileName())    
    
    context = {'objects': mics,
               'raphael':getResourceJs('raphael'),
               'high_res' : protocol.highRes.get(),
               'low_res': protocol.lowRes.get()
               }
    
    return render_to_response('wiz_ctf.html', context)

def wiz_particle_mask(protocol, request):
    
    parts = [p for p in protocol.inputParticles.get()] 
    
    for p in parts:
        index = p.getIndex()
        text = p.getFileName()
        p.basename = basename(text)
        if index:
            p.text = "%03d@%s" % (index, text)
            p.basename = "%03d@%s" % (index, basename(text))
            
    xdim = getImageXdim(request, parts[0].text)

    mask_radius = protocol.maskRadius.get()
     
    if mask_radius > xdim :
        mask_radius = xdim
    elif mask_radius == -1 :
        mask_radius = xdim/2
        
    context = {'objects': parts[0:100],
               'raphael': getResourceJs('raphael'),
               'maskRadius': mask_radius,
               'xdim':xdim
               }
    
    return render_to_response('wiz_particle_mask.html', context)

def wiz_volume_mask(protocol, request):
    
    vols = [vol for vol in protocol.input3DReferences.get()]
    
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
    
    vols = [vol for vol in protocol.input3DReferences.get()]
    
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
    
    parts = [p for p in protocol.inputParticles.get()] 
    
    for p in parts:
        index = p.getIndex()
        text = p.getFileName()
        p.basename = basename(text)
        if index:
            p.text = "%03d@%s" % (index, text)
            p.basename = "%03d@%s" % (index, basename(text))
             
    context = {'objects': parts[0:100],
               'lowFreq': protocol.lowFreq.get(),
               'highFreq': protocol.highFreq.get(),
               'decayFreq': protocol.freqDecay.get()
               }
    
    return render_to_response('wiz_bandpass.html', context)

def wiz_gaussian(protocol, request):
    
    parts = [p for p in protocol.inputParticles.get()] 
    
    for p in parts:
        index = p.getIndex()
        text = p.getFileName()
        p.basename = basename(text)
        if index:
            p.text = "%03d@%s" % (index, text)
            p.basename = "%03d@%s" % (index, basename(text))
             
    context = {'objects': parts[0:100],
               'freqSigma': protocol.freqSigma.get()
               }
    
    return render_to_response('wiz_gaussian.html', context)

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
