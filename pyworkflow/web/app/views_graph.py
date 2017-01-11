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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import json
import pyworkflow.gui.graph as gg
from views_util import loadProject
from views_project import getNodeStateColor
from pyworkflow.gui.tree import ProjectRunsTreeProvider
from django.http import HttpResponse


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

#===============================================================================
# PROTOCOL GRAPH
#===============================================================================

def project_graph(request):
    if request.is_ajax():
        boxList = request.GET.get('list')
        # Project Id(or Name) should be stored in SESSION
        projectName = request.session['projectName']
        # projectName = request.GET.get('projectName')
        project = loadProject(request)  
        
        provider = ProjectRunsTreeProvider(project)
        
        g = project.getRunsGraph()
        root = g.getRoot()
        root.w = 100
        root.h = 40
        root.item = WebNode('project', x=0, y=0)
        
        # Assign the width and height
        for box in boxList.split(','):
            id, w, h = box.split('-')
            node = g.getNode(id)
            if node is None:
                print "Get NONE node: i=%s" % id
            else:
                node.id = id
                node.w = float(w)
                node.h = float(h)
            
        lt = gg.LevelTree(g)
        lt.paint(createNode, createEdge)
        nodeList = []
        
        for node in g.getNodes():
            try:
                hx = node.w / 2
                hy = node.h / 2
                childs = [c.getName() for c in node.getChilds()]
                status, color = getNodeStateColor(node)
                
                info = ""
                if str(node.id) != "PROJECT":
                    protocol = project.getProtocol(int(node.id))
                    info = provider.getObjectInfo(protocol)["values"][0]
                
                nodeList.append({'id': node.getName(),
                                 'x': node.item.x - hx, 
                                 'y': node.item.y - hy,
                                 'color': color, 
                                 'status': info,
                                 'childs': childs})
            except Exception:
                print "Error with node: ", node.getName()
                raise
        
#        print nodeList
        jsonStr = json.dumps(nodeList, ensure_ascii=False)   
        return HttpResponse(jsonStr, mimetype='application/javascript')


#===============================================================================
# OBJECT GRAPH
#===============================================================================
    
def elements_graph(request):
    if request.is_ajax():
        projectName = request.session['projectName']
        project = loadProject(request)
        g = project.getSourceGraph()
        
        elmList = []
        for node in g.getNodes():
            elmList.append({'id': node.getName(),'label': node.getLabel()})
    
        jsonStr = json.dumps(elmList, ensure_ascii=False)   
        return HttpResponse(jsonStr, mimetype='application/javascript')
    
    
def object_graph(request):
    if request.is_ajax():
        boxList = request.GET.get('list')
        projectName = request.session['projectName']
        project = loadProject(request)
        g = project.getSourceGraph()
        
        root = g.getRoot()
        root.w = 100
        root.h = 40
        root.item = WebNode('project', x=0, y=0)
        
        # Assign the width and height
        for box in boxList.split(','):
            id, w, h = box.split('-')
            node = g.getNode(id)
            if node is None:
                print "Get NONE node: i=%s" % id
            else:
                node.id = id
                node.w = float(w)
                node.h = float(h)
        
        lt = gg.LevelTree(g)
        lt.paint(createNode, createEdge)
        
        nodeList = []
        for node in g.getNodes():
            try:
                hx = node.w / 2
                hy = node.h / 2
                childs = [c.getName() for c in node.getChilds()]
            
                nodeList.append({'id': node.getName(),
                                 'x': node.item.x - hx, 
                                 'y': node.item.y - hy,
                                 'status': None,
                                 'color': '#ADD8E6', # Lightblue
                                 'childs': childs})
                
            except Exception:
                    print "Error with node: ", node.getName()
                    raise
        
        jsonStr = json.dumps(nodeList, ensure_ascii=False)   
        return HttpResponse(jsonStr, mimetype='application/javascript')
        
    