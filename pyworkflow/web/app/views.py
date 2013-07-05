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
from pyworkflow.gui.tree import TreeProvider, ProjectRunsTreeProvider
from pyworkflow.utils.path import findResource
from pyworkflow.utils.utils import prettyDate
from pyworkflow.web.pages import settings
from pyworkflow.apps.config import *
from pyworkflow.em import *
from pyworkflow.hosts import HostMapper
from pyworkflow.tests import getInputPath 
from forms import HostForm
from django.core.urlresolvers import reverse
from django.http import HttpResponseRedirect, HttpResponse, HttpRequest

def getResource(request):
    if request == 'logoScipion':
        img = 'scipion_logo.png'
    elif request == 'favicon':
        img = 'scipion_bn.png'
    elif request == 'help':
        img = 'contents24.png'
    elif request == 'browse':
        img = 'zoom.png'
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

def projects(request):
    # CSS #
    css_path = os.path.join(settings.STATIC_URL, 'css/projects_style.css')
    messi_css_path = os.path.join(settings.STATIC_URL, 'css/messi.css')
    
    # JS #
    projectForm_path = os.path.join(settings.STATIC_URL, 'js/projectForm.js')
    
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
        if sub.value.hasValue() and tag == 'protocol_base':
            protClassName = value.split('.')[-1]  # Take last part
            prot = emProtocolsDict.get(protClassName, None)
            if prot is not None:
                for k, v in emProtocolsDict.iteritems():
                    if not v is prot and issubclass(v, prot):
                        protItem = TreeItem(k, 'protocol_class', protClassName)
                        item.childs.append(protItem)
        else:
            populateTree(item, sub)                
       
def loadConfig(config, name):
    c = getattr(config, name) 
    fn = getConfigPath(c.get())
    if not os.path.exists(fn):
        raise Exception('loadMenuConfig: menu file "%s" not found' % fn)
    mapper = ConfigMapper(getConfigPath(fn), globals())
    menuConfig = mapper.getConfig()
    return menuConfig


def loadProtTree():
    configMapper = ConfigMapper(getConfigPath('configuration.xml'), globals())
    generalCfg = configMapper.getConfig()
    protCfg = loadConfig(generalCfg, 'protocols')    
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
    css_path = os.path.join(settings.STATIC_URL, 'css/project_content_style.css')
    messi_css_path = os.path.join(settings.STATIC_URL, 'css/messi.css')
    
    # JS #
    jquery_cookie = os.path.join(settings.STATIC_URL, 'js/jquery.cookie.js')
    jquery_treeview = os.path.join(settings.STATIC_URL, 'js/jquery.treeview.js')
    launchTreeview = os.path.join(settings.STATIC_URL, 'js/launchTreeview.js')
    utils_path = os.path.join(settings.STATIC_URL, 'js/utils.js')
    tabs_config = os.path.join(settings.STATIC_URL, 'js/tabsConfig.js')
    
    
    projectName = request.GET.get('projectName', None)
    if projectName is None:
        projectName = request.POST.get('projectName', None)
        
    request.session['projectName'] = projectName
        
    project = loadProject(projectName)    
    provider = ProjectRunsTreeProvider(project)
    
    root = loadProtTree()
    
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
        print "======================= in protocol_io...."
        ioDict = {'inputs': [{'name':n, 'id': attr.getObjId()} for n, attr in protocol.iterInputAttributes()],
                  'outputs': [{'name':n, 'id': attr.getObjId()} for n, attr in protocol.iterOutputAttributes(EMObject)]}
        jsonStr = json.dumps(ioDict, ensure_ascii=False)
        print jsonStr
        
    return HttpResponse(jsonStr, mimetype='application/javascript')

def protocol_summary(request):
    # Project Id(or Name) should be stored in SESSION
    if request.is_ajax():
        projectName = request.GET.get('projectName')
        project = loadProject(projectName)
        protId = request.GET.get('protocolId', None)
        protocol = project.mapper.selectById(int(protId))
        summary = protocol.summary()
        print "======================= in protocol_summary...."
        jsonStr = json.dumps(summary, ensure_ascii=False)
        print jsonStr
        
    return HttpResponse(jsonStr, mimetype='application/javascript')
    
def form(request):
    
    # Resources #
    favicon_path = getResource('favicon')
    logo_help = getResource('help')
    logo_browse = getResource('browse')
    
    # CSS #
    css_path = os.path.join(settings.STATIC_URL, 'css/form.css')
    messi_css_path = os.path.join(settings.STATIC_URL, 'css/messi.css')
    
    # JS #
    jquery_path = os.path.join(settings.STATIC_URL, 'js/jquery.js')
    jsForm_path = os.path.join(settings.STATIC_URL, 'js/form.js')
    utils_path = os.path.join(settings.STATIC_URL, 'js/utils.js')
    messi_path = os.path.join(settings.STATIC_URL, 'js/messi.js')
    
    # Project Id(or Name) should be stored in SESSION
    projectName = request.session['projectName']
    # projectName = request.GET.get('projectName')
    project = loadProject(projectName)        
    protocolName = request.GET.get('protocol', None)
    action = request.GET.get('action', None)
    
    if protocolName is None:
        protId = request.GET.get('protocolId', None)
        protocol = project.mapper.selectById(int(protId))
    else:
        protocolClass = emProtocolsDict.get(protocolName, None)
        protocol = protocolClass()
    
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
    
    context = {'projectName':projectName,
               'protocol':protocol,
               'definition': protocol._definition,
               'favicon': favicon_path,
               'help': logo_help,
               'form': jsForm_path,
               'jquery': jquery_path,
               'browse': logo_browse,
               'utils': utils_path,
               'css':css_path,
               'messi': messi_path,
               'messi_css': messi_css_path}
    
    # Cross Site Request Forgery protection is need it
    context.update(csrf(request))
    
    return render_to_response('form.html', context)

def save_protocol(request):
    projectName = request.POST.get('projectName')
    protId = request.POST.get("protocolId")
    protClass = request.POST.get("protocolClass")
    
    # Load the project
    project = loadProject(projectName)
    
    # Create the protocol object
    if protId != 'None':  # Case of new protocol
        protId = request.POST.get('protocolId', None)
        protocol = project.mapper.selectById(int(protId))
    else:
        protocolClass = emProtocolsDict.get(protClass, None)
        protocol = protocolClass() 
    
    # Update parameter set in the form
    for paramName, attr in protocol.iterDefinitionAttributes():
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
        
    project.saveProtocol(protocol)
    
    return HttpResponse(mimetype='application/javascript')


def protocol(request):
    projectName = request.POST.get('projectName')
    protId = request.POST.get("protocolId")
    protClass = request.POST.get("protocolClass")
    
    # Load the project
    project = loadProject(projectName)
    # Create the protocol object
    if protId != 'None':  # Case of new protocol
        protId = request.POST.get('protocolId', None)
        protocol = project.mapper.selectById(int(protId))
    else:
        protocolClass = emProtocolsDict.get(protClass, None)
        protocol = protocolClass() 
    # Update parameter set in the form
    for paramName, attr in protocol.iterDefinitionAttributes():
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
    
    errors = protocol.validate()

    if len(errors) == 0:
        # No errors 
        # Finally, launch the protocol
        project.launchProtocol(protocol)
    jsonStr = json.dumps({'errors' : errors},
                     ensure_ascii=False)
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

def getScipionHosts():
    from pyworkflow.apps.config import getSettingsPath
    defaultHosts = getSettingsPath()
    return HostMapper(defaultHosts).selectAll()

def viewHosts(request):  
    # Resources #
    css_path = os.path.join(settings.STATIC_URL, 'css/general_style.css')
    jquery_path = os.path.join(settings.STATIC_URL, 'js/jquery.js')
    utils_path = os.path.join(settings.STATIC_URL, 'js/utils.js')
    host_utils_path = os.path.join(settings.STATIC_URL, 'js/hostUtils.js')
    
    projectName = request.session['projectName']    
    project = loadProject(projectName)
    projectHosts = project.getHosts()   
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


def getHostFormContext(request, initialContext=None):
    css_path = os.path.join(settings.STATIC_URL, 'css/general_style.css')
    jquery_path = os.path.join(settings.STATIC_URL, 'js/jquery.js')
    utils_path = os.path.join(settings.STATIC_URL, 'js/utils.js')
    hostId = request.GET.get("hostId")
    if hostId is None or hostId == "":
        if initialContext is not None:
            hostId = initialContext['hostId']
    form = None
    if request.method == 'GET':        
        form = HostForm(auto_id=True)
    else:
        form = HostForm(request.POST, queueConfCont=request.POST.get('queueConfigCount'), queueSystemConfCont=request.POST.get('queueSystemConfigCount'))
    projectName = request.session['projectName']
    project = loadProject(projectName)
    hostsMapper = HostMapper(project.settingsPath)
#     scpnHostsChoices = []
#     scpnHostsChoices.append(('', ''))
#     scipionHosts = getScipionHosts()
#     for hostConfig in scipionHosts:
#         scpnHostsChoices.append((hostConfig.getLabel(), hostConfig.getHostName()))
#     form.fields['scpnHosts'].choices = scpnHostsChoices        
    # We check if we are going to edit a host
    tittle = None
    if hostId is not None and hostId != "":
        hostConfig = hostsMapper.selectById(hostId)
        form.setHost(hostConfig)
        tittle = hostConfig.getLabel() + " host configuration"
    else:
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
    return render_to_response('hostForm.html', RequestContext(request, getHostFormContext(request)))  # Form Django forms

def updateHostsConfig(request):
    form = HostForm(request.POST)  # A form bound to the POST data
    if form.is_valid():  # All validation rules pass
        projectName = request.session['projectName']
        project = loadProject(projectName)
        host = project.saveHost(form.getHost())      
        context = {'hostId' : host.getObjId(),
                   'message': "Project hosts config sucesfully updated"}
        return render_to_response('hostForm.html', RequestContext(request, getHostFormContext(request, context)))  # Form Django forms
    else:
        return render_to_response('hostForm.html', RequestContext(request, getHostFormContext(request)))  # Form Django forms

def deleteHost(request):
    hostId = request.GET.get("hostId")    
    projectName = request.session['projectName']
    project = loadProject(projectName)
    project.deleteHost(hostId)
#     context = {'message': "Host succesfully deleted"}
    return HttpResponseRedirect('/viewHosts')

def visualizeObject(request):
    objectId = request.GET.get("objectId")    
    projectName = request.session['projectName']
    
#    project = loadProject(projectName)
    manager = Manager()
    projPath = manager.getProjectPath(projectName)
    request.session['projectPath'] = projPath
    project = Project(projPath)
    project.load()
    
    object = project.mapper.selectById(int(objectId))
    if object.isPointer():
        object = object.get()
        
    if isinstance(object, SetOfMicrographs):
        fn = project.getTmpPath(object.getName() + '_micrographs.xmd')
        mics = XmippSetOfMicrographs.convert(object, fn)
        inputParameters = {'path': join(request.session['projectPath'], mics.getFileName()),
                       'allowRender': True,
                       'mode': 'gallery',
                       'zoom': 150,
                       'gotoContainer': 1}
  
    elif isinstance(object, SetOfImages):
        fn = project.getTmpPath(object.getName() + '_images.xmd')
        imgs = XmippSetOfImages.convert(object, fn)
        inputParameters = {'path': join(request.session['projectPath'], imgs.getFileName()),
               'allowRender': True,
               'mode': 'gallery',
               'zoom': 150,
               'gotoContainer': 1}

    elif isinstance(object, XmippClassification2D):
        mdPath = object.getClassesMdFileName()
        block, path = mdPath.split('@')
        inputParameters = {'path': join(request.session['projectPath'], path),
               'allowRender': True,
               'mode': 'gallery',
               'zoom': 150,
               'gotoContainer': 1}
#        runShowJ(obj.getClassesMdFileName())
    else:
        raise Exception('Showj Web visualizer: can not visualize class: %s' % object.getClassName())

    
    
    from views_showj import showj 
    return showj(request, inputParameters)
    
    
    
#    url2 = reverse('app.views_showj.showj', kwargs={'path': path})
#    print "url"
#    print url2
#    return HttpResponseRedirect(url2)

#    from django.shortcuts import redirect

    # return redirect('/showj', args=inputParameters)

#    return HttpResponseRedirect('/showj', inputParameters)
    
if __name__ == '__main__':
    root = loadProtTree()    
    for s in root.childs:
        print s.name, '-', s.tag
        for p in s.childs:
            print p.name, '-', p.tag
