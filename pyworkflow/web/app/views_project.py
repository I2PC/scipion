import json
import pyworkflow.gui.graph as gg
from pyworkflow.em import *
from views_util import * 
from pyworkflow.utils.utils import prettyDate
from pyworkflow.manager import Manager
from pyworkflow.apps.pw_project_viewprotocols import STATUS_COLORS
from pyworkflow.gui.tree import TreeProvider, ProjectRunsTreeProvider
from django.http import HttpResponse, HttpRequest
from django.contrib.gis.shortcuts import render_to_text

from pyworkflow.em.packages.xmipp3.convert import writeSetOfParticles
from pyworkflow.em.packages.xmipp3.plotter import XmippPlotter

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

def createNode(node, y):
    try:
        item = gg.TNode(node.getName(), y=y)
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
        root.item = gg.TNode('project', x=0, y=0)
        
        
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
    
def project_content(request):        
    projectName = request.GET.get('projectName', None)
    
    if projectName is None:
        projectName = request.POST.get('projectName', None)
        
    request.session['projectName'] = projectName
    manager = Manager()
    request.session['projectPath'] = manager.getProjectPath(projectName)
        
    project = loadProject(projectName)    
    provider = ProjectRunsTreeProvider(project)
    
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
        jsonStr = json.dumps(summary, ensure_ascii=False)
        
#        print "======================= in protocol_summary...."
#        print jsonStr
    return HttpResponse(jsonStr, mimetype='application/javascript')

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

def viewer(request):    
    if request.is_ajax():
        projectPath= request.session['projectPath']
        projectName = request.session['projectName']
        project = loadProject(projectName)
        protId = request.GET.get('protocolId', None)
        protocol = project.mapper.selectById(int(protId))
        
        print "======================= in viewer...."
          
        request.GET = request.GET.copy()
        
        # TODO: This should be changed to a more general approach 
        # when the different protocol viewer be defined (in the same way of wizards.)

        if getattr(protocol, 'outputMicrographs', False):
            objId = protocol.outputMicrographs.getObjId()
            request.GET['objectId'] = protocol.outputMicrographs.getObjId()
        elif getattr(protocol, 'outputParticles', False):
            objId = protocol.outputParticles.getObjId()
            request.GET['objectId'] = protocol.outputParticles.getObjId()
        
        if isinstance(protocol, XmippProtML2D):
            ioDict = viewer_ML2D(protId)
        else:
            ioDict = viewer_default(protId, objId)
       
        jsonStr = json.dumps(ioDict, ensure_ascii=False)
                
    return HttpResponse(jsonStr, mimetype='application/javascript')

def viewer_default(protId, objId):
    #==XmippPloter Functionality============================================
    url_plotter = "/view_plot/?protocolId="+protId
    #==Showj visualizer=====================================================
    from views_showj import visualizeObject
#   response = visualizeObject(request)
    url_showj = "/visualize_object/?objectId="+str(objId)
                
    ioDict = {"url": url_showj , "plot" : url_plotter}
    
    return ioDict

def viewer_ML2D(protId):
    url_form = "/form/?protocolId=" + protId + "&action=visualize"
    ioDict = {"url_form": url_form}
    
    return ioDict
  
