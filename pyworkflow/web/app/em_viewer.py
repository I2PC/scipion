import json
import pyworkflow.gui.graph as gg
from pyworkflow.em import *
from views_util import *
from views_protocol import updateProtocolParams
from views_showj import visualizeObject
from pyworkflow.manager import Manager
from pyworkflow.gui.tree import TreeProvider, ProjectRunsTreeProvider
from pyworkflow.em import emProtocolsDict
from django.http import HttpResponse

from pyworkflow.em.packages.xmipp3.convert import writeSetOfParticles
from pyworkflow.em.packages.xmipp3.plotter import XmippPlotter
from pyworkflow.viewer import WEB_DJANGO

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.dates import DateFormatter

############## 1ST STEP: LAUNCH VIEWER METHODS ##############
def launch_viewer(request):
    if request.is_ajax():
        projectName = request.session['projectName']
        project = loadProject(projectName)
        protId = request.GET.get('protocolId', None)
        protocol = project.mapper.selectById(int(protId))
        
        viewers = findViewers(protocol.getClassName(), WEB_DJANGO)
        print viewers
        
        viewer = viewers[0]()
        functionName = viewer.getView()
        function = globals().get(functionName, None)
        
        if function is None:
            pass  # redirect to error: viewer not found
        elif not callable(function):
            pass  # redirect to error: name is not a function
        else:
            ioDict = function(project, protocol, viewer)
       
        jsonStr = json.dumps(ioDict, ensure_ascii=False)
    return HttpResponse(jsonStr, mimetype='application/javascript')

def viewerXmipp(project, protocol, viewer):
    
    if getattr(protocol, 'outputMicrographs', False):
        objId = protocol.outputMicrographs.getObjId()
    elif getattr(protocol, 'outputParticles', False):
        objId = protocol.outputParticles.getObjId()
    
    from views_showj import visualizeObject
    url_showj = "/visualize_object/?objectId="+str(objId)

    protId = protocol.getObjId()
    url_plotter = "/view_plot_xmipp/?protocolId="+ str(protId)
    
    ioDict = {"url": url_showj , "plot" : url_plotter}
    return ioDict

def viewerForm(project, protocol, viewer):
    protId = protocol.getObjId()
    viewerClassName = viewer.getClassName();
    
    ioDict = {"url_form": "/form/?protocolClass="+ viewerClassName +
              "&protRunIdViewer="+ str(protId) +
              "&action=visualize"}
    
    return ioDict 

############## 2ND STEP: VIEWER FUNCTION METHODS ##############
def viewer(request):
    project, protocolViewer = loadProtocolProject(request)
    updateProtocolParams(request, protocolViewer, project)
    protId = request.POST.get('protRunIdViewer', None)
    protocol = project.mapper.selectById(int(protId))
    functionName = protocolViewer.getViewFunction()
    
    function = globals().get(functionName, None)
    
    if function is None:
        pass  # redirect to error: viewer not found
    elif not callable(function):
        pass  # redirect to error: name is not a function
    else:
        ioDict = function(request, protocol, protocolViewer)
    
    jsonStr = json.dumps(ioDict, ensure_ascii=False)
#    print jsonStr
    return HttpResponse(jsonStr, mimetype='application/javascript')

def viewerElement(request):
    project, protocolViewer = loadProtocolProject(request)
    updateProtocolParams(request, protocolViewer, project)
    protId = request.POST.get('protRunIdViewer', None)
    viewerParam = request.POST.get('viewerParam', None)
    protocol = project.mapper.selectById(int(protId))
    
    functionName = protocolViewer.getVisualizeDictWeb()[viewerParam]
    function = globals().get(functionName, None)
    
    if function is None:
        pass  # redirect to error: viewer not found
    elif not callable(function):
        pass  # redirect to error: name is not a function
    else:
        ioDict= {}
        typeUrl, url = function(request, protocol, protocolViewer)
        ioDict[typeUrl] = url
    
    jsonStr = json.dumps(ioDict, ensure_ascii=False)
    return HttpResponse(jsonStr, mimetype='application/javascript')

############## VIEWER XMIPP CL2D ##############
def viewerCL2D(request, protocol, protocolViewer):
    ioDict = {}
   
    text = viewLevelFiles(request, protocol, protocolViewer)
    ioDict["showj"] = text
    
    if protocolViewer.doShowClassHierarchy:
        text = viewClassHierarchy(request, protocol, protocolViewer)
        ioDict["html"] = text
                    
    return ioDict

def viewLevelFiles(request, protocol, protocolViewer):
    fnSubset = protocolViewer._getSubset()
    levelFiles = protocol._getLevelMdFiles(fnSubset)
    url=""
    if levelFiles:
        levelFiles.sort()
        lastLevelFile = levelFiles[-1]
        if protocolViewer.doShowLastLevel:
            url = "/showj/?path="+lastLevelFile
        else:
            if protocolViewer.showSeveralLevels.empty():
                print 'Please select the levels that you want to visualize.'
            else:
                listOfLevels = []
                try:
                    listOfLevels = protocolViewer._getListFromRangeString(protocolViewer.showSeveralLevels.get())
                except Exception:
                    print 'Invalid levels range.'
                    
                files = "";
                for level in listOfLevels:
                    fn = protocol._getExtraPath("level_%02d/level_classes%s.xmd"%(level,fnSubset))
                    if os.path.exists(fn):
                        files += "classes_sorted@"+fn+" "
                    else:
                        print 'Level %s does not exist.' % level
                if files != "":
                    print 'runShowJ(files)' 
        return url
        

def viewClassHierarchy(request, protocol, protocolViewer):
    fnSubset = protocolViewer._getSubset()
    fnHierarchy = protocol._getExtraPath("classes%s_hierarchy.txt" % fnSubset)
    if os.path.exists(fnHierarchy):
        html = textfileViewer(fnHierarchy, [fnHierarchy])
    return html

def textfileViewer(title, fileList):
    f = open(fileList[0], 'r')
        
    style = "background-color:black;color:white;font-family:Monospace;padding:1em;font-size:90%;"
    title = "<title>"+ title + "</title>"
    html = "<div style="+ style +">"+ title
    
    x = 0
    while 1:
        line = f.readline()
        if not line:
            break
        if len(line) > 1:
            x = x+1
            html = html +"<p><span style='color:cyan;'>" + str(x) + ":    </span>"+ line +" </p>"
            
    html = html + "</div>"
    return html

############## VIEWER XMIPP ML2D ##############
def viewerML2D(request, protocol, protocolViewer):
    ioDict = {}

    if protocolViewer.doShowClasses:
        typeUrl, url = doShowClasses(request, protocol, protocolViewer)
        ioDict[typeUrl]= url
    if protocolViewer.doShowPlots:
        typeUrl, url = doAllPlotsML2D(request, protocol, protocolViewer)
        ioDict[typeUrl]= url
    else:
        typeUrl, url = doSomePlotsML2D(protocol, protocolViewer)
        if url != None:
            ioDict[typeUrl]= url
        
    return ioDict

def doShowClasses(request, protocol, protocolViewer):
    objId = protocol.outputClasses.getObjId()
    return "showj","/visualize_object/?objectId="+str(objId)
            
def doAllPlotsML2D(request, protocol, protocolViewer):
    return "plots","/view_plots/?function=allPlotsML2D&protViewerClass="+ str(protocolViewer.getClassName())+ "&protId="+ str(protocol.getObjId())

def doSomePlotsML2D(protocol, protocolViewer):
    plots=""
    for p in protocolViewer._plotVars:
        if protocolViewer.getAttributeValue(p):
            plots = plots + p + "-"
    
    if plots != "":
        return "plots","/view_plots/?function=somePLotsML2D&plots="+str(plots)+"&protViewerClass="+ str(protocolViewer.getClassName())+ "&protId="+ str(protocol.getObjId())
    else:
        return "", None

def doShowLL(request, protocol, protocolViewer):
    return "plot","/view_plots/?function=somePLotsML2D&plots=doShowLL&protViewerClass="+ str(protocolViewer.getClassName())+ "&protId="+ str(protocol.getObjId())

def doShowPmax(request, protocol, protocolViewer):
    return "plot","/view_plots/?function=somePLotsML2D&plots=doShowPmax&protViewerClass="+ str(protocolViewer.getClassName())+ "&protId="+ str(protocol.getObjId())

def doShowSignalChange(request, protocol, protocolViewer):
    return "plot","/view_plots/?function=somePLotsML2D&plots=doShowSignalChange&protViewerClass="+ str(protocolViewer.getClassName())+ "&protId="+ str(protocol.getObjId())
    
def doShowMirror(request, protocol, protocolViewer):
    return "plot","/view_plots/?function=somePLotsML2D&plots=doShowMirror&protViewerClass="+ str(protocolViewer.getClassName())+ "&protId="+ str(protocol.getObjId())

def allPlotsML2D(request, protocolViewer, protocol):
    xplotter = protocolViewer.createPlots(protocol, protocolViewer._plotVars)
    return xplotter

def somePLotsML2D(request, protocolViewer, protocol):
    plots = request.GET.get('plots', None)
    plots = str(plots).split("-")
    if len(plots) > 1:
        plots.remove(plots[len(plots)-1])
    
    xplotter = protocolViewer.createPlots(protocol, plots)
    return xplotter

############## AUX METHODS ##############
def view_plots(request):
    projectName = request.session['projectName']
    project = loadProject(projectName)
    protId = request.GET.get('protId', None)
    protocol = project.mapper.selectById(int(protId))
    
    protViewerClass = request.GET.get('protViewerClass', None)
    protocolClass = emProtocolsDict.get(protViewerClass, None)
    protocolViewer = protocolClass()
    
    functionName = request.GET.get('function', None)
    function = globals().get(functionName, None)
    
    xplotter = function(request, protocolViewer, protocol)
    
    canvas = xplotter.getCanvas()
    response = HttpResponse(content_type='image/png')
    canvas.print_png(response)
    return response
    
def view_plot_xmipp(request):
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
               
#        print "MD contains ZSCORE"
        xplotter = XmippPlotter(windowTitle="Zscore particles sorting")
        xplotter.createSubPlot("Particle sorting", "Particle number", "Zscore")
        xplotter.plotMd(md, False, mdLabelY=xmipp.MDL_ZSCORE)
        figFn = fn.replace('.xmd', '.png')
        canvas = xplotter.getCanvas()
        
        response= HttpResponse(content_type='image/png')
        canvas.print_png(response)
        return response
