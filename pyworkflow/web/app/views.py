# from scipion.models import *
import pyworkflow as pw
import os
import xmipp
from django.shortcuts import render_to_response
from django.core.context_processors import csrf
from django.template import RequestContext
import json
from pyworkflow.manager import Manager
from pyworkflow.project import Project
import pyworkflow.gui.graph as gg
from pyworkflow.gui.tree import TreeProvider, ProjectRunsTreeProvider
from pyworkflow.utils.path import findResource
from pyworkflow.utils.utils import prettyDate
from pyworkflow.web.pages import settings
from pyworkflow.apps.config import *
from pyworkflow.apps.pw_project_viewprotocols import STATUS_COLORS
from pyworkflow.em import *
from pyworkflow.hosts import HostMapper
from pyworkflow.tests import getInputPath 
from forms import HostForm, VolVisualizationForm
from django.core.urlresolvers import reverse
from django.http import HttpResponseRedirect, HttpResponse, HttpRequest
from commons import staticPath

iconDict = {
            'logoScipion': 'scipion_logo.png',
            'favicon': 'scipion_bn.png',
            }
def getResource(request):
#    return iconDict[request]
    
    if request == 'logoScipion':
        img = 'scipion_logo.png'
    elif request == 'favicon':
        img = 'scipion_bn.png'
    elif request == 'help':
        img = 'system_help24.png'
    elif request == 'browse':
        img = 'zoom.png'
    elif request == 'wizard':
        img = 'tools_wizard.png'
    elif request == 'edit_toolbar':
        img = 'edit.gif'
    elif request == 'copy_toolbar':
        img = 'copy.gif'
    elif request == 'delete_toolbar':
        img = 'delete.gif'
    elif request == 'browse_toolbar':
        img = 'run_steps.gif'
    elif request == 'tree_toolbar':
        img = 'tree2.gif'
    elif request == 'new_toolbar':
        img = 'new_object.gif'
        
    path = os.path.join(settings.MEDIA_URL, img)
    return path

# Resources #
new_tool_path = getResource('new_toolbar')
edit_tool_path = getResource('edit_toolbar')
copy_tool_path = getResource('copy_toolbar')
delete_tool_path = getResource('delete_toolbar')
browse_tool_path = getResource('browse_toolbar')
tree_tool_path = getResource('tree_toolbar')

######    Projects template    #####
def projects(request):
    # CSS #
    css_path = staticPath('css/projects_style.css')
    messi_css_path = staticPath('css/messi.css')
    
    # JS #
    projectForm_path = staticPath('js/projectForm.js')
    
    manager = Manager()
    
    projects = manager.listProjects()
    for p in projects:
        p.pTime = prettyDate(p.mTime)

    context = {'projects': projects,
               'css': css_path,
               'messi_css': messi_css_path,
               'projectForm':projectForm_path,
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

######    Project Content template    #####
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

    
def loadProject(projectName):
    manager = Manager()
    projPath = manager.getProjectPath(projectName)
    project = Project(projPath)
    project.load()
    return project
    
    
def project_content(request):        
    # CSS #
    css_path = staticPath('css/project_content_style.css')
    messi_css_path = staticPath('css/messi.css')

    # JS #
    jquery_cookie = staticPath('js/jquery.cookie.js')
    jquery_treeview = staticPath('js/jquery.treeview.js')
    launchTreeview = staticPath('js/launchTreeview.js')
    utils_path = staticPath('js/utils.js')
    tabs_config = staticPath('js/tabsConfig.js')
    
    projectName = request.GET.get('projectName', None)
    if projectName is None:
        projectName = request.POST.get('projectName', None)
        
    request.session['projectName'] = projectName
        
    project = loadProject(projectName)    
    provider = ProjectRunsTreeProvider(project)
    
    root = loadProtTree(project)
    
    context = {'projectName': projectName,
               'editTool': edit_tool_path,
               'copyTool': copy_tool_path,
               'deleteTool': delete_tool_path,
               'browseTool': browse_tool_path,
               'treeTool': tree_tool_path,
               'utils': utils_path,
               'jquery_cookie': jquery_cookie,
               'jquery_treeview': jquery_treeview,
               'launchTreeview': launchTreeview,
               'tabs_config': tabs_config,
               'css':css_path,
               'sections': root.childs,
               'provider':provider,
               'messi_css': messi_css_path,
               'view': 'protocols',
               'contentConfig': 'divided'}
    
    return render_to_response('project_content.html', context)

def delete_protocol(request):
    # Project Id(or Name) should be stored in SESSION
    if request.is_ajax():
        projectName = request.GET.get('projectName')
        project = loadProject(projectName)
        protId = request.GET.get('protocolId', None)
        protocol = project.mapper.selectById(int(protId))
     
        project.deleteProtocol(protocol)         
        
    return HttpResponse(mimetype='application/javascript')

def protocol_io(request):
    # Project Id(or Name) should be stored in SESSION
    if request.is_ajax():
        projectName = request.GET.get('projectName')
        project = loadProject(projectName)
        protId = request.GET.get('protocolId', None)
        protocol = project.mapper.selectById(int(protId))
#        print "======================= in protocol_io...."
        ioDict = {'inputs': [{'name':n, 'id': attr.getObjId()} for n, attr in protocol.iterInputAttributes()],
                  'outputs': [{'name':n, 'id': attr.getObjId()} for n, attr in protocol.iterOutputAttributes(EMObject)]}
        jsonStr = json.dumps(ioDict, ensure_ascii=False)
#        print jsonStr
        
    return HttpResponse(jsonStr, mimetype='application/javascript')

def protocol_summary(request):
    # Project Id(or Name) should be stored in SESSION
    if request.is_ajax():
        projectName = request.GET.get('projectName')
        project = loadProject(projectName)
        protId = request.GET.get('protocolId', None)
        protocol = project.mapper.selectById(int(protId))
        summary = protocol.summary()
#        print "======================= in protocol_summary...."
        jsonStr = json.dumps(summary, ensure_ascii=False)
#        print jsonStr
        
    return HttpResponse(jsonStr, mimetype='application/javascript')

######    Project Form template    #####
def form(request):
    
    # Resources #
    favicon_path = getResource('favicon')
    logo_help = getResource('help')
    logo_browse = getResource('browse')
    logo_wiz = getResource('wizard')
    logo_edit = getResource('edit_toolbar')
    
    # CSS #
    css_path = staticPath('css/form.css')
    messi_css_path = staticPath('css/messi.css')
    
    # JS #
    jquery_path = staticPath('js/jquery.js')
    jsForm_path = staticPath('js/form.js')
    utils_path = staticPath('js/utils.js')
    messi_path = staticPath('js/messi.js')
    
    project, protocol = loadProtocolProject(request, requestType='GET')
    
    action = request.GET.get('action', None)
    hosts = [host.getLabel() for host in project.getSettings().getHosts()]
    
    if action == 'copy':
        protocol = project.copyProtocol(protocol)
    
    # TODO: Add error page validation when protocol is None
    for section in protocol._definition.iterSections():
        for paramName, param in section.iterParams():
            protVar = getattr(protocol, paramName, None)
            if protVar is None:
                raise Exception("_fillSection: param '%s' not found in protocol" % paramName)
                # Create the label
            if protVar.isPointer():
                if protVar.hasValue():
                    param.htmlValue = protVar.get().getNameId()
                else:
                    param.htmlValue = ""
            else:
                param.htmlValue = protVar.get(param.default.get(""))
                if isinstance(protVar, Boolean):
                    if param.htmlValue:
                        param.htmlValue = 'true'
                    else:
                        param.htmlValue = 'false' 
            param.htmlCond = param.condition.get()
            param.htmlDepend = ','.join(param._dependants)
            param.htmlCondParams = ','.join(param._conditionParams)
#            param.htmlExpertLevel = param.expertLevel.get()   
    
    context = {'projectName':project.getName(),
               'protocol':protocol,
               'definition': protocol._definition,
               'favicon': favicon_path,
               'help': logo_help,
               'comment':logo_edit,
               'form': jsForm_path,
               'jquery': jquery_path,
               'browse': logo_browse,
               'wizard': logo_wiz,
               'utils': utils_path,
               'css':css_path,
               'messi': messi_path,
               'messi_css': messi_css_path,
               'hosts':hosts
               }
    # Update the context dictionary with the special params
    for paramName in SPECIAL_PARAMS:
        context[paramName] = protocol.getAttributeValue(paramName, '')
    
    # Cross Site Request Forgery protection is need it
    context.update(csrf(request))
    
    return render_to_response('form.html', context)

def wizard(request):
    action = request.GET.get('action', None)
    
    if(action=='downsampling'):
        response = 'wiz_downsampling.html'
    
    context={'action':action
             }
    
    return render_to_response(response, context)
    

def loadProtocolProject(request, requestType='POST'):
    """ Retrieve the project and protocol from this request.
    Return:
        (project, protocol) tuple
    """
    requestDict = getattr(request, requestType)
    projectName = request.session['projectName']
    protId = requestDict.get("protocolId")
    protClass = requestDict.get("protocolClass")
    
    # Load the project
    project = loadProject(projectName)
    
    # Create the protocol object
    if protId and protId != 'None':  # Case of new protocol
        protId = requestDict.get('protocolId', None)
        protocol = project.mapper.selectById(int(protId))
    else:
        protocolClass = emProtocolsDict.get(protClass, None)
        protocol = protocolClass()
        
    return (project, protocol)
    
def updateParam(request, project, protocol, paramName):
    """
    Params:
        request: current request handler
        project: current working project
        protocol: current protocol
        paramName: name of the attribute to be set in the protocol
            from the web form
    """
    attr = getattr(protocol, paramName)
    value = request.POST.get(paramName)
    if attr.isPointer():
        if len(value.strip()) > 0:
            objId = int(value.split('.')[-1])  # Get the id string for last part after .
            value = project.mapper.selectById(objId)  # Get the object from its id
            if attr.getObjId() == value.getObjId():
                raise Exception("Param: %s is autoreferencing with id: %d" % (paramName, objId))
        else:
            value = None
    attr.set(value)
    print "setting attr %s with value:" % paramName, value 
    
SPECIAL_PARAMS = ['runName', 'numberOfMpi', 'numberOfThreads', 'hostName','expertLevel','_useQueue']

def updateProtocolParams(request, protocol, project):
    """ Update the protocol values from the Web-form.
    This function will be used from save_protocol and execute_protocol.
    """
    for paramName, _ in protocol.iterDefinitionAttributes():
        updateParam(request, project, protocol, paramName)
        
    for paramName in SPECIAL_PARAMS:
        updateParam(request, project, protocol, paramName)
        
def save_protocol(request):
    project, protocol = loadProtocolProject(request)
    updateProtocolParams(request, protocol, project)            
    project.saveProtocol(protocol)
    
    return HttpResponse(mimetype='application/javascript')

# Method to launch a protocol #
def protocol(request):
    project, protocol = loadProtocolProject(request)
    updateProtocolParams(request, protocol, project)    
    errors = protocol.validate()

    if len(errors) == 0:
        # No errors 
        # Finally, launch the protocol
        project.launchProtocol(protocol)
    jsonStr = json.dumps({'errors' : errors}, ensure_ascii=False)
    
    return HttpResponse(jsonStr, mimetype='application/javascript')

def browse_objects(request):
    """ Browse objects from the database. """
    if request.is_ajax():
        objClass = request.GET.get('objClass')
        projectName = request.GET.get('projectName')
        project = loadProject(projectName)    
        
        objs = []
        for obj in project.mapper.selectByClass(objClass, iterate=True):
            objs.append(obj.getNameId())
        jsonStr = json.dumps({'objects' : objs},
                             ensure_ascii=False)
        return HttpResponse(jsonStr, mimetype='application/javascript')


######    Hosts template    #####
def getScipionHosts():
    from pyworkflow.apps.config import getSettingsPath
    defaultHosts = getSettingsPath()
    return HostMapper(defaultHosts).selectAll()

def viewHosts(request):  
    # Resources #
    css_path = staticPath('css/general_style.css')
    jquery_path = staticPath('js/jquery.js')
    utils_path = staticPath('js/utils.js')
    host_utils_path = staticPath('js/hostUtils.js')
    
    projectName = request.session['projectName']    
    project = loadProject(projectName)
    projectHosts = project.getSettings().getHosts()   
    scpnHostsChoices = []
    scpnHostsChoices.append(('', ''))

    message = request.GET.get("message")
    context = {'projectName' : projectName,
               'editTool': edit_tool_path,
               'newTool': new_tool_path,
               'deleteTool': delete_tool_path,
               'browseTool': browse_tool_path,
               'hosts': projectHosts,
               'jquery': jquery_path,
               'utils': utils_path,
               'hostUtils': host_utils_path,
               'css':css_path,
               'message': message,
               'view': 'hosts',
               'contentConfig': 'full'}    
      
    return render_to_response('hosts.html', context)

# def getHost(request):
#     from django.http import HttpResponse
#     import json
#     from django.utils import simplejson
#     
#     if request.is_ajax():
#         hostLabel = request.GET.get('hostLabel')
#         projectName = request.session['projectName']
#         project = loadProject(projectName)
#         hostsMapper = HostMapper(project.settingsPath)
#         HostConfig = hostsMapper.selectByLabel(hostLabel)
#         jsonStr = json.dumps({'host':HostConfig.getDictionary()})
#         return HttpResponse(jsonStr, mimetype='application/javascript')


def getHostFormContext(request, host=None, initialContext=None):
    css_path = staticPath('css/general_style.css')
    jquery_path = staticPath('js/jquery.js')
    utils_path = staticPath('js/utils.js')
    form = None
#     scpnHostsChoices = []
#     scpnHostsChoices.append(('', ''))
#     scipionHosts = getScipionHosts()
#     for hostConfig in scipionHosts:
#         scpnHostsChoices.append((hostConfig.getLabel(), hostConfig.getHostName()))
#     form.fields['scpnHosts'].choices = scpnHostsChoices        
    # We check if we are going to edit a host
    tittle = None
    if host is not None:
        form = initialContext['form']
        tittle = host.getLabel() + " host configuration"
    else:
        form = HostForm()
        tittle = "New host configuration"  
            
    context = {'tittle': tittle,
               'jquery': jquery_path,
               'utils': utils_path,
               'css':css_path,
               'form': form}
     
    if initialContext is not None:
        context.update(initialContext)
        
    return context

def hostForm(request):
    hostId = request.GET.get('hostId')
    context = None
    hostConfig = None
    if hostId != None and hostId != '':
        projectName = request.session['projectName']
        project = loadProject(projectName)
        hostsMapper = HostMapper(project.settingsPath)
        hostConfig = hostsMapper.selectById(hostId)
        form = HostForm()
        form.setFormHost(hostConfig)
        context = {'form': form}
    return render_to_response('host_form.html', RequestContext(request, getHostFormContext(request, hostConfig, context)))  # Form Django forms

def updateHostsConfig(request):    
    form = HostForm(request.POST, queueSystemConfCont=request.POST.get('queueSystemConfigCount'), queueConfCont=request.POST.get('queueConfigCount'))  # A form bound to the POST data
    context = {'form': form}
    if form.is_valid():  # All validation rules pass
        projectName = request.session['projectName']
        project = loadProject(projectName)
        hostId = request.POST.get('objId')
        hostsMapper = HostMapper(project.settingsPath)
        hostConfig = hostsMapper.selectById(hostId)
        form.host = hostConfig
        host = form.getFormHost()
        savedHost = project.getSettings().saveHost(host)
        form.setFormHost(savedHost)
        return render_to_response('host_form.html', RequestContext(request, getHostFormContext(request, host, context))) 
    else:   
        return render_to_response('host_form.html', RequestContext(request, getHostFormContext(request, None, context)))  # Form Django forms

def deleteHost(request):
    hostId = request.GET.get("hostId")    
    projectName = request.session['projectName']
    project = loadProject(projectName)
    project.deleteHost(hostId)
#     context = {'message': "Host succesfully deleted"}
    return HttpResponseRedirect('/view_hosts')#, RequestContext(request))

def visualizeObject(request):
    probandoCTFParam = True
    
    objectId = request.GET.get("objectId")    
    #projectName = request.session['projectName']
    projectName = request.GET.get("projectName")
    
#    project = loadProject(projectName)
    manager = Manager()
    projPath = manager.getProjectPath(projectName)
    request.session['projectPath'] = projPath
    project = Project(projPath)
    project.load()
    
    object = project.mapper.selectById(int(objectId))
    if object.isPointer():
        object = object.get()
        
    if probandoCTFParam:
        inputParameters = {'path': join(request.session['projectPath'], "Runs/XmippProtCTFMicrographs218/extra/BPV_1386/xmipp_ctf.ctfparam"),
               'allowRender': True,
               'mode': 'column',
               'zoom': '150px',
               'goto': 1,
               'colRowMode': 'Off'}    
    elif isinstance(object, SetOfMicrographs):
        fn = project.getTmpPath(object.getName() + '_micrographs.xmd')
        mics = XmippSetOfMicrographs.convert(object, fn)
        inputParameters = {'path': join(request.session['projectPath'], mics.getFileName()),
                       'allowRender': True,
                       'mode': 'gallery',
                       'zoom': '150px',
                       'goto': 1,
                       'colRowMode': 'Off'}
    elif isinstance(object, SetOfVolumes):
        print ("XXXXX", object.getObjId())
        inputParameters = {'setOfVolumes' : object, 'setOfVolumesId': object.getObjId()}  
    elif isinstance(object, SetOfImages):
        fn = project.getTmpPath(object.getName() + '_images.xmd')
        imgs = XmippSetOfImages.convert(object, fn)
        inputParameters = {'path': join(request.session['projectPath'], imgs.getFileName()),
               'allowRender': True,
               'mode': 'gallery',
               'zoom': '150px',
               'goto': 1,
               'colRowMode': 'Off'}

    elif isinstance(object, XmippClassification2D):
        mdPath = object.getClassesMdFileName()
        block, path = mdPath.split('@')
        inputParameters = {'path': join(request.session['projectPath'], path),
               'allowRender': True,
               'mode': 'gallery',
               'zoom': '150px',
               'goto': 1,
               'colRowMode': 'Off'}
#        runShowJ(obj.getClassesMdFileName())
    else:
        raise Exception('Showj Web visualizer: can not visualize class: %s' % object.getClassName())

    if isinstance(object, SetOfVolumes):
        return render_to_response('volume_visualization.html', inputParameters)
    else:
        from views_showj import showj 
        return showj(request, inputParameters)
    
            
    
#    url2 = reverse('app.views_showj.showj', kwargs={'path': path})
#    print "url"
#    print url2
#    return HttpResponseRedirect(url2)

#    from django.shortcuts import redirect

    # return redirect('/showj', args=inputParameters)

#    return HttpResponseRedirect('/showj', inputParameters)

def visualizeVolume(request):
    from django.http import HttpResponse
    import json
     
    if request.is_ajax():
        setOfVolumesId = int(request.GET.get('setOfVolumesId'))
        volumeId = int(request.GET.get('volumeId'))
        projectName = request.session['projectName']
        project = loadProject(projectName)
        setOfVolume = project.mapper.selectById(setOfVolumesId)
        volume = setOfVolume[volumeId]
        # Chimera 
        from subprocess import Popen, PIPE, STDOUT
#         inputVolume = join(os.getcwd(),volume.getFileName())
#         outputHtmlFile = join(os.getcwd(),project.getTmpPath("volume_" + str(volume.getObjId()) + '.html'))
        inputVolume = volume.getFileName()
        outputHtmlFile = project.getTmpPath("volume_" + str(volume.getObjId()) + '.html')
        if (request.GET.get('threshold') is None or request.GET.get('threshold') == ''):
            threshold = 1
        else:
            threshold = float(request.GET.get('threshold'))
        print ("THRESHOLD",threshold)
        # TODO: Get with Xmipp an approximate threshold
        p = Popen(['chimera', inputVolume], stdout=PIPE, stdin=PIPE, stderr=PIPE)
#         p = Popen(['chimera', '--start', 'ReadStdin', inputVolume], stdout=PIPE, stdin=PIPE, stderr=PIPE)        
        stdout_data = p.communicate(input='volume #0 level ' + str(threshold) + '; export format WebGL ' + outputHtmlFile + '; stop')[0]
        f = open(outputHtmlFile)
        volumeHtml = f.read().decode('string-escape').decode("utf-8").split("</html>")[1]
        jsonStr = json.dumps({'volumeHtml': volumeHtml})
        return HttpResponse(jsonStr, mimetype='application/javascript')
    
    
def showVolVisualization(request):
    form = None
    volLinkPath = None
    volLink = None
    chimeraHtml = None
    if (request.POST.get('operation') == 'visualize'):
        form = VolVisualizationForm(request.POST, request.FILES)
        if form.is_valid():
            volPath = form.cleaned_data['volPath']
            # Astex viewer            
            from random import randint
            linkName = 'test_link_' + str(randint(0, 10000)) + '.map'
            volLinkPath = os.path.join(pw.HOME, 'web', 'pages', 'resources', 'astex', 'tmp', linkName)
            from pyworkflow.utils.path import cleanPath, createLink
            cleanPath(volLinkPath)
            createLink(volPath, volLinkPath)
#             os.system("ln -s " + str(volPath) + " " + volLinkPath)
            volLink = os.path.join('/', 'static', 'astex', 'tmp', linkName)
            # Chimera 
            from subprocess import Popen, PIPE, STDOUT
#             p = Popen(['chimera', '--start', 'ReadStdin', volPath], stdout=PIPE, stdin=PIPE, stderr=PIPE)
            p = Popen(['chimera', volPath], stdout=PIPE, stdin=PIPE, stderr=PIPE)
            outputHtmlFile = '/home/antonio/test.html'
            threshold = form.cleaned_data['threshold']
            stdout_data = p.communicate(input='volume #0 level ' + str(threshold) + '; export format WebGL ' + outputHtmlFile + '; stop')[0]
            f = open(outputHtmlFile)
            chimeraHtml = f.read().decode('string-escape').decode("utf-8").split("</html>")[1]
    else:
        form = VolVisualizationForm()
    context = {'MEDIA_URL' : settings.MEDIA_URL, 'STATIC_URL' :settings.STATIC_URL, 'form': form, 'volLink': volLink, 'chimeraHtml': chimeraHtml}    
    return render_to_response('showVolVisualization.html',  RequestContext(request, context))   

        

if __name__ == '__main__':
    root = loadProtTree()    
    for s in root.childs:
        print s.name, '-', s.tag
        for p in s.childs:
            print p.name, '-', p.tag
