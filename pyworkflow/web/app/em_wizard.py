import os
import xmipp
from pyworkflow.web.app.views_util import * 
from pyworkflow.manager import Manager
from pyworkflow.project import Project
from django.shortcuts import render_to_response

def wizard(request):
    _ , protocol = loadProtocolProject(request, 'GET')
    functionName = request.GET.get('function', None)
    function = globals().get(functionName, None)
    
    if function is None:
        pass # redirect to error: wizard not found
    elif not callable(function):
        pass # redirect to error: name is not a function
    else:
        return function(protocol)

def wiz_downsampling(protocol):
    mics = [mic for mic in protocol.inputMicrographs.get()]
    
    context = {'objects': mics}
    
    return render_to_response('wiz_downsampling.html', context)
    
    