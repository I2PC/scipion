import os
import xmipp
from pyworkflow.web.app.views_util import * 
from pyworkflow.manager import Manager
from pyworkflow.project import Project
from django.shortcuts import render_to_response
from pyworkflow.gui import getImage, getPILImage
from django.http import HttpResponse

def wizard(request):
    _ , protocol = loadProtocolProject(request, 'GET')
    functionName = request.GET.get('function', None)
    function = globals().get(functionName, None)
    
    if function is None:
        pass  # redirect to error: wizard not found
    elif not callable(function):
        pass  # redirect to error: name is not a function
    else:
        return function(protocol)

def wiz_downsampling(protocol):
    mics = [mic for mic in protocol.inputMicrographs.get()]
    
    context = {'objects': mics}
    
    return render_to_response('wiz_downsampling.html', context)

def wiz_frequencies(protocol):
    mics = [mic for mic in protocol.inputMicrographs.get()]
    
    context = {'objects': mics}
    
    return render_to_response('wiz_frequencies.html', context)
    
def get_image_psd(request):
    
    imagePath = request.GET.get('image', None)
    downsample = request.GET.get('downsample', None)
    dim = request.GET.get('dim', None)
    
    # create a xmipp image empty
    imgXmipp = xmipp.Image()
    
    # compute the PSD image
    xmipp.fastEstimateEnhancedPSD(imgXmipp, str(imagePath), int(downsample), int(dim), 2)
        
    # from PIL import Image
    img = getPILImage(imgXmipp, dim)
        
    response = HttpResponse(mimetype="image/png")
    img.save(response, "PNG")
    return response
