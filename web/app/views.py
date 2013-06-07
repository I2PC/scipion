# from scipion.models import *
import os
from django.shortcuts import render_to_response
from django.core.context_processors import csrf
from pyworkflow.manager import Manager
from pyworkflow.project import Project
from pyworkflow.gui.tree import TreeProvider
from pyworkflow.utils.path import findResource
from pyworkflow.utils.utils import prettyDate
from pyworkflow.web.pages import settings
from pyworkflow.apps.config import *
from pyworkflow.em import *
from django.http.request import HttpRequest

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
        
    path = os.path.join(settings.MEDIA_URL, img)
    return path

def projects(request):
    manager = Manager()
#    logo_path = findResource('scipion_logo.png')

    # Resources #
    css_path = os.path.join(settings.STATIC_URL, 'css/projects_style.css')
    #############
    projectForm_path = os.path.join(settings.STATIC_URL, 'js/projectForm.js')
    jquery_path = os.path.join(settings.STATIC_URL, 'js/jquery.js')
    
    # Messi Plugin #
    messi_path = os.path.join(settings.STATIC_URL, 'js/messi.js')
    messi_css_path = os.path.join(settings.STATIC_URL, 'css/messi.css')
    #############
    
    projects = manager.listProjects()
    for p in projects:
        p.pTime = prettyDate(p.mTime)

    context = {'jquery':jquery_path,
               'projects': projects,
               'css': css_path,
               'messi': messi_path,
               'messi_css': messi_css_path,
               'projectForm':projectForm_path}
    
    return render_to_response('projects.html', context)

def create_project(request):
    from django.http import HttpResponse    
    
    manager = Manager()
    
    if request.is_ajax():
        projectName = request.GET.get('projectName')
        manager.createProject(projectName)       
        
    return HttpResponse(mimetype='application/javascript')

def delete_project(request):
    from django.http import HttpResponse    
    
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

# to do a module from pw_project
class RunsTreeProvider(TreeProvider):
    """Provide runs info to populate tree"""
    def __init__(self, mapper, actionFunc=None):
        self.actionFunc = actionFunc
        self.getObjects = lambda: mapper.selectAll()
        
    def getColumns(self):
        return [('Run', 250), ('State', 100), ('Modified', 100)]
    
    def getObjectInfo(self, obj):
        return {'key': obj.getObjId(),
                'text': '%s.%s' % (obj.getClassName(), obj.strId()),
                'values': (obj.status.get(), obj.endTime.get())}
      
    def getObjectActions(self, obj):
        prot = obj  # Object should be a protocol
        actionsList = [(ACTION_EDIT, 'Edit     '),
                       # (ACTION_COPY, 'Duplicate   '),
                       (ACTION_DELETE, 'Delete    '),
                       # (None, None),
                       # (ACTION_STOP, 'Stop'),
                       (ACTION_STEPS, 'Browse data')
                       ]
        status = prot.status.get()
        if status == STATUS_RUNNING:
            actionsList.insert(0, (ACTION_STOP, 'Stop execution'))
            actionsList.insert(1, None)
        elif status == STATUS_WAITING_APPROVAL:
            actionsList.insert(0, (ACTION_CONTINUE, 'Approve continue'))
            actionsList.insert(1, None)
        
        actions = []
        def appendAction(a):
            v = a
            if v is not None:
                action = a[0]
                text = a[1]
                v = (text, lambda: self.actionFunc(action), ActionIcons[action])
            actions.append(v)
            
        for a in actionsList:
            appendAction(a)
            
        return actions 
    
def loadProject(projectName):
    manager = Manager()
    projPath = manager.getProjectPath(projectName)
    project = Project(projPath)
    project.load()
    return project    
    
def project_content(request):    
    # Resources #
    edit_tool_path = getResource('edit_toolbar')
    copy_tool_path = getResource('copy_toolbar')
    delete_tool_path = getResource('delete_toolbar')
    browse_tool_path = getResource('browse_toolbar')
    
    css_path = os.path.join(settings.STATIC_URL, 'css/project_content_style.css')
    jquery_path = os.path.join(settings.STATIC_URL, 'js/jquery.js')
    jquery_cookie = os.path.join(settings.STATIC_URL, 'js/jquery.cookie.js')
    jquery_treeview = os.path.join(settings.STATIC_URL, 'js/jquery.treeview.js')
    launchTreeview = os.path.join(settings.STATIC_URL, 'js/launchTreeview.js')
    utils_path = os.path.join(settings.STATIC_URL, 'js/utils.js')
    #############
     # Messi Plugin
    messi_path = os.path.join(settings.STATIC_URL, 'js/messi.js')
    messi_css_path = os.path.join(settings.STATIC_URL, 'css/messi.css')
    #############
    
    projectName = request.GET.get('projectName', None)
    if projectName is None:
        projectName = request.POST.get('projectName', None)
        
    project = loadProject(projectName)    
    provider = RunsTreeProvider(project.mapper)
    
    root = loadProtTree()
    
    context = {'projectName':projectName,
               'editTool': edit_tool_path,
               'copyTool': copy_tool_path,
               'deleteTool': delete_tool_path,
               'browseTool': browse_tool_path,
               'jquery': jquery_path,
               'utils': utils_path,
               'jquery_cookie': jquery_cookie,
               'jquery_treeview': jquery_treeview,
               'launchTreeview': launchTreeview,
               'css':css_path,
               'sections': root.childs,
               'provider':provider,
               'messi': messi_path,
               'messi_css': messi_css_path}
    
    return render_to_response('project_content.html', context)

def delete_protocol(request):
    from django.http import HttpResponse   
    
    # Project Id(or Name) should be stored in SESSION
    if request.is_ajax():
        projectName = request.GET.get('projectName')
        project = loadProject(projectName)
        protId = request.GET.get('protocolId', None)
        protocol = project.mapper.selectById(int(protId))
     
        project.deleteProtocol(protocol)         
        
    return HttpResponse(mimetype='application/javascript')

def form(request):
    
    # Resources #
    favicon_path = getResource('favicon')
    logo_help = getResource('help')
    logo_browse = getResource('browse')
    jquery_path = os.path.join(settings.STATIC_URL, 'js/jquery.js')
    jsForm_path = os.path.join(settings.STATIC_URL, 'js/form.js')
    utils_path = os.path.join(settings.STATIC_URL, 'js/utils.js')
    css_path = os.path.join(settings.STATIC_URL, 'css/form.css')
    # Messi Plugin
    messi_path = os.path.join(settings.STATIC_URL, 'js/messi.js')
    messi_css_path = os.path.join(settings.STATIC_URL, 'css/messi.css')
    #############
    
    # # Project Id(or Name) should be stored in SESSION
    projectName = request.GET.get('projectName')
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
    # Finally, launch the protocol
    error = protocol.validate()
    if error == []:
        pass
    else:
        # Errors
        pass
    
    
    project.launchProtocol(protocol)
    
    return project_content(request)

def browse_objects(request):
    """ Browse objects from the database. """
    from django.http import HttpResponse
    import json
    
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

def hosts(request):
    # Resources #
    favicon_path = getResource('favicon')
    jquery_path = os.path.join(settings.STATIC_URL, 'js/jquery.js')
    # # Project Id(or Name) should be stored in SESSION
    projectName = request.GET.get('projectName')
    project = loadProject(projectName)
    context = {'projectName' : projectName,
               'project' : project,
               'hosts': project.getHosts(),
               'favicon': favicon_path,
               'jquery': jquery_path}
    
    return render_to_response('hosts.html', context)

def showj(request):
    # manager = Manager()
#    logo_path = findResource('scipion_logo.png')

    # Resources #
    css_path = os.path.join(settings.STATIC_URL, 'css/showj_style.css')
    favicon_path = getResource('favicon')
    jquery_path = os.path.join(settings.STATIC_URL, 'js/jquery.js')
    #############

    context = {'favicon': favicon_path,
               'jquery': jquery_path,
               'css':css_path}
    
    return render_to_response('showj.html', context)
    
    
if __name__ == '__main__':
    root = loadProtTree()    
    for s in root.childs:
        print s.name, '-', s.tag
        for p in s.childs:
            print p.name, '-', p.tag
