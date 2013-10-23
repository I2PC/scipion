import json
import pyworkflow.gui.graph as gg
from pyworkflow.em import *
from views_util import *
from views_protocol import updateProtocolParams 
from pyworkflow.manager import Manager
from pyworkflow.gui.tree import TreeProvider, ProjectRunsTreeProvider
from django.http import HttpResponse

from pyworkflow.em.packages.xmipp3.convert import writeSetOfParticles
from pyworkflow.em.packages.xmipp3.plotter import XmippPlotter
from pyworkflow.viewer import WEB_DJANGO

############## 1ST STEP: LAUNCH VIEWER METHODS ##############
def launch_viewer(request):
    if request.is_ajax():
        projectName = request.session['projectName']
        project = loadProject(projectName)
        protId = request.GET.get('protocolId', None)
        protocol = project.mapper.selectById(int(protId))
        
        viewers = findViewers(protocol.getClassName(), WEB_DJANGO)
        
        viewer = viewers[0]()
        functionName = viewer.getView()
        function = globals().get(functionName, None)
        
        if function is None:
            pass  # redirect to error: viewer not found
        elif not callable(function):
            pass  # redirect to error: name is not a function
        else:
            ioDict = function(request, protocol, viewer)
       
        jsonStr = json.dumps(ioDict, ensure_ascii=False)
    return HttpResponse(jsonStr, mimetype='application/javascript')

def viewerXmipp(request, protocol, viewer):
    protId = request.GET.get('protocolId', None)
    
    if getattr(protocol, 'outputMicrographs', False):
        objId = protocol.outputMicrographs.getObjId()
    elif getattr(protocol, 'outputParticles', False):
        objId = protocol.outputParticles.getObjId()
    
    url_plotter = "/view_plot/?protocolId="+protId
    
    from views_showj import visualizeObject
    url_showj = "/visualize_object/?objectId="+str(objId)
    
    ioDict = {"url": url_showj , "plot" : url_plotter}
    return ioDict

def viewerForm(request, protocol, viewer):
    viewerClassName = viewer.__class__.__name__
    ioDict = {"url_form": "/form/?protocolClass="+ viewerClassName +"&action=visualize"}
    return ioDict  

############## 2ND STEP: VIEWER FUNCTION METHODS ##############
def viewer(request):
    project, protocol = loadProtocolProject(request)
    updateProtocolParams(request, protocol, project)
    
    functionName = protocol.getViewFunction()
    function = globals().get(functionName, None)    
    
    if function is None:
        pass  # redirect to error: wizard not found
    elif not callable(function):
        pass  # redirect to error: name is not a function
    else:
        return function(request, protocol, viewer)

def viewerML2D (request, protocol, viewer):
    

    return HttpResponse("<h1>Under Construction</h1>");

############## AUX METHODS ##############
def view_plot(request):
    projectName = request.session['projectName']
    project = loadProject(projectName)
    protId = request.GET.get('protocolId', None)
    protocol = project.mapper.selectById(int(protId))
    
    mdFn = getattr(protocol.outputParticles, '_xmippMd', None)
    if mdFn:
        fn = mdFn.get()
    else:
        fn = project.getTmpPath(protocol.outputParticles.getName() + '_images.xmd')
        writeSetOfParticles(protocol.outputParticles, fn, protocol.getTmpPath())
        
    md = xmipp.MetaData(fn)
    if md.containsLabel(xmipp.MDL_ZSCORE):
        
        from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
        from matplotlib.figure import Figure
        from matplotlib.dates import DateFormatter
        
#        print "MD contains ZSCORE"
        xplotter = XmippPlotter(windowTitle="Zscore particles sorting")
        xplotter.createSubPlot("Particle sorting", "Particle number", "Zscore")
        xplotter.plotMd(md, False, mdLabelY=xmipp.MDL_ZSCORE)
        figFn = fn.replace('.xmd', '.png')
        canvas = xplotter.getCanvas()
        
        response= HttpResponse(content_type='image/png')
        canvas.print_png(response)
        return response
