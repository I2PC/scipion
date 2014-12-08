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
from django.http import HttpResponse
from pyworkflow.viewer import WEB_DJANGO, MessageView, ProtocolViewer, TextView

############## 1ST STEP: LAUNCH VIEWER METHODS ##############

def launch_viewer(request):
    from views_util import loadProject
    from pyworkflow.em import findViewers, EMProtocol, EMObject
    
    projectName = request.session['projectName']
    project = loadProject(projectName)
    objId = request.GET.get('objectId', None).split("::")[0]
    obj = project.getProtocol(int(objId))
    
    if obj.isPointer():
        obj = obj.get()
        
    viewers = findViewers(obj.getClassName(), WEB_DJANGO)
    
    if len(viewers) == 0:
        views = []
        
        if isinstance(obj, EMProtocol):
            for _, output in obj.iterOutputAttributes(EMObject):
                objViewers = findViewers(output.getClassName(), WEB_DJANGO)
                if objViewers:
                    views += objViewers[0](project=project)._visualize(output) or []

        if not views:
            views = [MessageView("No viewers found for object type: %s" % obj.getClassName())]                    
        
        urls = [viewToUrl(request, view) for view in views]
    else:
        viewer = viewers[0](project=project)
        # Lets assume if there is a viewer for the protocol
        # it will be a ProtocolViewer with a Form
        if isinstance(viewer, ProtocolViewer):
            urls = [viewerForm(project, obj, viewer)]
        else:
            views = viewer._visualize(obj)
            urls = [viewToUrl(request, v) for v in views]
            
    jsonStr = json.dumps(urls, ensure_ascii=False)
        
    return HttpResponse(jsonStr, mimetype='application/javascript')


def viewerForm(project, protocol, viewer):
    protId = protocol.getObjId()
    viewerClassName = viewer.getClassName()
    
    return "url::/form/?protocolClass=%s&protRunIdViewer=%s&action=visualize" % (viewerClassName, protId)
    

############## 2ND STEP: VIEWER FUNCTION METHODS ##############

def viewToUrl(request, view):
    from views_util import savePlot
    from pyworkflow.gui.plotter import Plotter
    from pyworkflow.em import TABLE_NAME, DataView, PATH

    # PLOT
    if isinstance(view, Plotter):
        url = 'url::/' + savePlot(request, view)
    
    # SHOWJ
    elif isinstance(view, DataView):
        url = "/showj/?%s=%s" % (PATH, view.getPath())
        showjParams = view.getShowJWebParams()
        
        if showjParams:
            for key in showjParams:
                url += '&%s=%s' % (key, showjParams[key])
        
        if view.getTableName():
            url += '&%s=%s' % (TABLE_NAME, view.getTableName())
        url = 'showj::'+ url
        
    # TEXT VIEWER
    elif isinstance(view, TextView):
        fn  = view.getFileList()[0]
        url = "file::/file_viewer/?%s=%s" % (PATH, fn)
        
    # MESSAGE
    elif isinstance(view, MessageView):
        url = 'error::' + view.getMessage()
    
    return url

def viewer_element(request):
    from views_util import loadProtocolProject
    from views_protocol import updateProtocolParams
    
    project, protocolViewer = loadProtocolProject(request)
    protId = request.POST.get('protRunIdViewer', None)
    viewerParam = request.POST.get('viewerParam', None)
    
    protocol = project.getProtocol(int(protId))
    protocolViewer.setProtocol(protocol)
    updateProtocolParams(request, protocolViewer, project)
    
    views = protocolViewer._getVisualizeDict()[viewerParam](viewerParam)
    
    urls = [viewToUrl(request, v) for v in views]
    
    jsonStr = json.dumps(urls, ensure_ascii=False)
    return HttpResponse(jsonStr, mimetype='application/javascript')
    
