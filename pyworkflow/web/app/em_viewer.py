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

import json
from pyworkflow.em import *
from views_util import *
from views_protocol import updateProtocolParams
from pyworkflow.manager import Manager
from pyworkflow.em import emProtocolsDict
from django.http import HttpResponse

from pyworkflow.em.packages.xmipp3.convert import writeSetOfParticles
from pyworkflow.em.packages.xmipp3.plotter import XmippPlotter
from pyworkflow.viewer import WEB_DJANGO

# XMIPP
from viewers.xmipp_ml2d import *
from viewers.xmipp_cl2d import *
from viewers.xmipp_ml3d import *
from viewers.xmipp_nma import *
from viewers.xmipp_nma_align import *
# SPIDER
from viewers.spider_capca import *
from viewers.spider_ward import *

############## 1ST STEP: LAUNCH VIEWER METHODS ##############
def launch_viewer(request):
    if request.is_ajax():
        projectName = request.session['projectName']
        project = loadProject(projectName)
        protId = request.GET.get('protocolId', None)
        protocol = project.mapper.selectById(int(protId))
        
        viewers = findViewers(protocol.getClassName(), WEB_DJANGO)
        
        if len(viewers) == 0:
            msg = "There is not viewer for protocol: <strong>" + protocol.getClassName() +"</strong>"
            ioDict = {'error': msg}
        else:
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
    ioDict={}
    
    if getattr(protocol, 'outputMicrographs', False):
        objId = protocol.outputMicrographs.getObjId()
    elif getattr(protocol, 'outputClasses', False):
        objId = protocol.outputClasses.getObjId()
    elif getattr(protocol, 'outputParticles', False):
        objId = protocol.outputParticles.getObjId()
        protId = protocol.getObjId()
        ioDict["plot"] = "/view_plot_xmipp/?protocolId="+ str(protId)
        
    ioDict["url"] = "/visualize_object/?objectId="+str(objId)
    
    if isinstance(protocol, XmippProtKerdensom):
        ioDict['url'] += '&mode=gallery&colRowMode=On&cols=%d' % protocol.SomXdim.get()
    if isinstance(protocol, XmippProtRotSpectra):
        ioDict['url'] += '&classCount___renderable=True&classCount___renderFunc=getTestPlot'
    
    return ioDict

def viewerSpider(project, protocol, viewer):
    ioDict={}    
        
    if isinstance(protocol, PcaFile):
        print "PcaFile"
        html = textfileViewer("PCA files", [protocol.filename.get()])
        ioDict["html"] = html
        
    elif isinstance(protocol, SpiderProtFilter):
        particles = protocol.outputParticles
        url = "/visualize_object/?objectId="+str(particles.getObjId())
        ioDict["url"] = url
        
    elif isinstance(protocol, SpiderProtCustomMask):
        mask = protocol.outputMask
        url1 = "/visualize_object/?objectId="+str(mask.getObjId())
        url2 = "/visualize_object/?path="+str(mask.getFileName())    
            
        ioDict["urls"] = [url1, url2] 
        
    return ioDict

def viewerRelion(project, protocol, viewer):
    ioDict={}    
        
    if isinstance(protocol, Relion3DClassification):
        
        volumes = protocol.outputVolumes
        url1 = "/visualize_object/?objectId="+str(volumes.getObjId())
            
        ioDict["urls"] = [url1] 
        
    return ioDict


def viewerForm(project, protocol, viewer):
    protId = protocol.getObjId()
    viewerClassName = viewer.getClassName()
    
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
    protocolViewer.setProtocol(protocol)
    protocolViewer.showPlot = False # Get xplotter instead of show()
    functionName = protocolViewer.getViewFunction()
    
    function = globals().get(functionName, None)
    
    if function is None:
        pass  # redirect to error: viewer not found
    elif not callable(function):
        pass  # redirect to error: name is not a function
    else:
        ioDict = function(request, protocolViewer)
    
    jsonStr = json.dumps(ioDict, ensure_ascii=False)
#    print jsonStr
    return HttpResponse(jsonStr, mimetype='application/javascript')

def viewerElement(request):
    project, protocolViewer = loadProtocolProject(request)
    updateProtocolParams(request, protocolViewer, project)
    protId = request.POST.get('protRunIdViewer', None)
    viewerParam = request.POST.get('viewerParam', None)
    protocol = project.mapper.selectById(int(protId))
    protocolViewer.setProtocol(protocol)
    protocolViewer.showPlot = False # Get xplotter instead of show()
    functionName = protocolViewer.getVisualizeDictWeb()[viewerParam]
    function = globals().get(functionName, None)
    
    if function is None:
        pass  # redirect to error: viewer not found
    elif not callable(function):
        pass  # redirect to error: name is not a function
    else:
        ioDict= {}
        typeUrl, url = function(request, protocolViewer)
        ioDict[typeUrl] = url
    
    jsonStr = json.dumps(ioDict, ensure_ascii=False)
    return HttpResponse(jsonStr, mimetype='application/javascript')

############## AUX METHODS ##############
def view_plots(request):
    from pyworkflow.em.packages.xmipp3.plotter import XmippPlotter
    XmippPlotter.setInteractive(False)
#    print "plotter.interactive: ", XmippPlotter.interactive
    
    projectName = request.session['projectName']
    project = loadProject(projectName)
    protId = request.GET.get('protId', None)
    protocol = project.mapper.selectById(int(protId))
    
    protViewerClass = request.GET.get('protViewerClass', None)
    protocolClass = emProtocolsDict.get(protViewerClass, None)
    protocolViewer = protocolClass()
    protocolViewer.setProtocol(protocol)
    protocolViewer.showPlot = False # Get xplotter instead of show()
    
    updateProtocolParams(request, protocolViewer, project)
    
    functionName = request.GET.get('function', None)
    function = globals().get(functionName, None)
    
    xplotter = function(request, protocolViewer)
    
#    figure = xplotter.getFigure()
#    print "width:", figure.get_figwidth()*100
#    print "height:", figure.get_figheight()*100
        
    canvas = xplotter.getCanvas()
    
    # Adding this line the space between plots is fixed
#    try:
#        xplotter.show()
#    except:
#        pass
    
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


