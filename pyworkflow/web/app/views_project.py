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
import pyworkflow.gui.graph as gg
from pyworkflow.em import *
from views_util import * 
from views_protocol import updateParam 
from pyworkflow.utils.utils import prettyDate
from pyworkflow.manager import Manager
from pyworkflow.apps.pw_project_viewprotocols import STATUS_COLORS
from pyworkflow.gui.tree import TreeProvider, ProjectRunsTreeProvider
from django.http import HttpResponse, HttpRequest
from django.contrib.gis.shortcuts import render_to_text

from pyworkflow.em.packages.xmipp3.convert import writeSetOfParticles
from pyworkflow.em.packages.xmipp3.plotter import XmippPlotter
from pyworkflow.viewer import WEB_DJANGO

def projects(request):
    manager = Manager()
    
    projects = manager.listProjects()
    for p in projects:
        p.pTime = prettyDate(p.mTime)

    context = {'projects': projects,
               'css': getResourceCss('projects'),
               'messi_css': getResourceCss('messi'),
               'project_utils': getResourceJs('project_utils'),
               'contentConfig': 'full'}
    
    return render_to_response('projects.html', context)

def create_project(request):
    manager = Manager()
    
    if request.is_ajax():
        projectName = request.GET.get('projectName')
        manager.createProject(projectName)       
        
    return HttpResponse(mimetype='application/javascript')

def delete_project(request):
    
    manager = Manager()
    
    if request.is_ajax():
        projectName = request.GET.get('projectName')
        manager.deleteProject(projectName)       
        
    return HttpResponse(mimetype='application/javascript')

class WebNode(object):
    def __init__(self, text, x=0, y=0):
        self.text = text
        self.moveTo(x, y)
        self.width, self.height = 0, 0
        
    def getDimensions(self):
        return (self.width, self.height)
    
    def moveTo(self, x, y):
        self.x = x
        self.y = y

def createNode(canvas, node, y):
    try:
        item = WebNode(node.getName(), y=y)
        item.width = node.w
        item.height = node.h
    except Exception:
        print "Error with node: ", node.getName()
        raise
    return item
    
def createEdge(srcItem, dstItem):
    pass
    

def getNodeStateColor(node):
    color = '#ADD8E6';  # Lightblue
    status = ''
    if node.run:
        status = node.run.status.get(STATUS_FAILED)
        color = STATUS_COLORS[status]
        
    return status, color

def project_graph (request):
    if request.is_ajax():
        boxList = request.GET.get('list')
        # Project Id(or Name) should be stored in SESSION
        projectName = request.session['projectName']
        # projectName = request.GET.get('projectName')
        project = loadProject(projectName)
        g = project.getRunsGraph()
        root = g.getRoot()
        root.w = 100
        root.h = 40
        root.item = WebNode('project', x=0, y=0)
        
        
        for box in boxList.split(','):
            i, w, h = box.split('-')
            node = g.getNode(i)
#            print node.getName()
            node.w = float(w)
            node.h = float(h)
            
        lt = gg.LevelTree(g)
        lt.paint(createNode, createEdge)
        nodeList = []
        
#        nodeList = [{'id': node.getName(), 'x': node.item.x, 'y': node.item.y} 
#                    for node in g.getNodes()]
        for node in g.getNodes():
            try:
                hx = node.w / 2
                hy = node.h / 2
                childs = [c.getName() for c in node.getChilds()]
                status, color = getNodeStateColor(node)
                nodeList.append({'id': node.getName(), 'x': node.item.x - hx, 'y': node.item.y - hy,
                                 'color': color, 'status': status,
                                 'childs': childs})
            except Exception:
                print "Error with node: ", node.getName()
                raise
        
#        print nodeList
        jsonStr = json.dumps(nodeList, ensure_ascii=False)   
         
        return HttpResponse(jsonStr, mimetype='application/javascript')

class TreeItem():
    def __init__(self, name, tag, protClass=None):
        self.name = name
        self.tag = tag
        self.protClass = protClass
        self.childs = []
        
def populateTree(tree, obj):    
    for sub in obj:
        text = sub.text.get()
        value = sub.value.get(text)
        tag = sub.tag.get('')
        item = TreeItem(text, tag)
        tree.childs.append(item)
        # If have tag 'protocol_base', fill dynamically with protocol sub-classes
        protClassName = value.split('.')[-1]  # Take last part
        if sub.value.hasValue() and tag == 'protocol_base':
            prot = emProtocolsDict.get(protClassName, None)
            if prot is not None:
                for k, v in emProtocolsDict.iteritems():
                    if not v is prot and issubclass(v, prot):
                        protItem = TreeItem(k, 'protocol_class', protClassName)
                        item.childs.append(protItem)
        else:
            item.protClass = protClassName
            populateTree(item, sub)                


def loadProtTree(project):
    protCfg = project.getSettings().getCurrentProtocolMenu()
    root = TreeItem('root', 'root')
    populateTree(root, protCfg)
    return root    

def update_graph_view(request):
    status = request.GET.get('status', None)
    projectName = request.session['projectName']
    project = loadProject(projectName) 
#    status = project.getSettings().graphView.get()
    if status == "True":
        project.getSettings().graphView.set(True)
    else :
        project.getSettings().graphView.set(False)
    project.getSettings().write()
    return HttpResponse(mimetype='application/javascript')
    
def project_content(request):        
    projectName = request.GET.get('projectName', None)
    
    if projectName is None:
        projectName = request.POST.get('projectName', None)
        
    request.session['projectName'] = projectName
    manager = Manager()
    request.session['projectPath'] = manager.getProjectPath(projectName)
   
    project = loadProject(projectName)    
    provider = ProjectRunsTreeProvider(project)
    graphView = project.getSettings().graphView.get()
    
    root = loadProtTree(project)
    
    context = {'projectName': projectName,
               'editTool': getResourceIcon('edit_toolbar'),
               'copyTool': getResourceIcon('copy_toolbar'),
               'deleteTool': getResourceIcon('delete_toolbar'),
               'browseTool': getResourceIcon('browse_toolbar'),
               'stopTool': getResourceIcon('stop_toolbar'),
               'analyzeTool': getResourceIcon('analyze_toolbar'),
               'treeTool': getResourceIcon('tree_toolbar'),
               'listTool': getResourceIcon('list_toolbar'),
               'utils': getResourceJs('utils'),
               'graph_utils': getResourceJs('graph_utils'),
               'project_content_utils': getResourceJs('project_content_utils'),
               'jquery_cookie': getResourceJs('jquery_cookie'),
               'jquery_treeview': getResourceJs('jquery_treeview'),
               'tabs_config': getResourceJs('tabs_config'),
               'css':getResourceCss('project_content'),
               'jquery_ui':getResourceCss('jquery_ui'),
               'sections': root.childs,
               'provider':provider,
               'messi_css': getResourceCss('messi'),
               'view': 'protocols',
               'graphView': graphView,
               'contentConfig': 'divided'}
    
    return render_to_response('project_content.html', context)

def protocol_status(request):
    if request.is_ajax():
        projectName = request.session['projectName']
        project = loadProject(projectName)
        protId = request.GET.get('protocolId', None)
        protocol = project.mapper.selectById(int(protId))
        status = protocol.status.get()

#        print "======================= in protocol_status...."
#        print jsonStr        
    return HttpResponse(status, mimetype='application/javascript')

def protocol_io(request):
    if request.is_ajax():
        projectName = request.session['projectName']
        project = loadProject(projectName)
        protId = request.GET.get('protocolId', None)
        protocol = project.mapper.selectById(int(protId))
        ioDict = {'inputs': [{'name':n, 'id': attr.getObjId()} for n, attr in protocol.iterInputAttributes()],
                  'outputs': [{'name':n, 'id': attr.getObjId()} for n, attr in protocol.iterOutputAttributes(EMObject)]}
        jsonStr = json.dumps(ioDict, ensure_ascii=False)

#        print "======================= in protocol_io...."
#        print jsonStr        
    return HttpResponse(jsonStr, mimetype='application/javascript')

def protocol_summary(request):
    if request.is_ajax():
        projectName = request.session['projectName']
        project = loadProject(projectName)
        protId = request.GET.get('protocolId', None)
        protocol = project.mapper.selectById(int(protId))
        summary = protocol.summary()
        print summary
        jsonStr = json.dumps(summary, ensure_ascii=False)
        
#        print "======================= in protocol_summary...."
#        print jsonStr
    return HttpResponse(jsonStr, mimetype='application/javascript')    
    