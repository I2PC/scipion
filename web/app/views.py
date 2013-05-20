# from scipion.models import *
import os
from django.shortcuts import render_to_response
from pyworkflow.manager import Manager
from pyworkflow.project import Project
from pyworkflow.gui.tree import TreeProvider
from pyworkflow.utils.path import findResource
from pyworkflow.utils.utils import prettyDate
from pyworkflow.web.pages import settings
from pyworkflow.apps.config import *
from pyworkflow.em import *

def getResource(request):
    if request == 'logoScipion':
        img = 'scipion_logo.png'
    elif request == 'favicon':
        img = 'scipion_bn.png'
    elif request == 'help':
        img = 'contents24.png'
    elif request == 'browse':
        img = 'zoom.png'
    path = os.path.join(settings.MEDIA_URL, img)
    return path

def projects(request):
    manager = Manager()
#    logo_path = findResource('scipion_logo.png')

    # Resources #
    css_path = os.path.join(settings.STATIC_URL, 'css/projects_style.css')
    #############
    
    projects = manager.listProjects()
    for p in projects:
        p.pTime = prettyDate(p.mTime)

    context = {'projects': projects,
               'css':css_path}
    
    return render_to_response('projects.html', context)

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
        return {'key': obj.getId(),
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
    
        
def project_content(request):
    
    # Resources #
    css_path = os.path.join(settings.STATIC_URL, 'css/project_content_style.css')
    jquery_path = os.path.join(settings.STATIC_URL, 'js/jquery.js')
    jquery_cookie = os.path.join(settings.STATIC_URL, 'js/jquery.cookie.js')
    jquery_treeview = os.path.join(settings.STATIC_URL, 'js/jquery.treeview.js')
    launchTreeview = os.path.join(settings.STATIC_URL, 'js/launchTreeview.js')
    popup_path = os.path.join(settings.STATIC_URL, 'js/popup.js')
    #############

    manager = Manager()
    project_name = request.GET.get('project_name')
    projPath = manager.getProjectPath(project_name)
    project = Project(projPath)
    project.load()
    provider = RunsTreeProvider(project.mapper)
    
    root = loadProtTree()
    
    context = {'project_name':project_name,
               'jquery': jquery_path,
               'popup': popup_path,
               'jquery_cookie': jquery_cookie,
               'jquery_treeview': jquery_treeview,
               'launchTreeview': launchTreeview,
               'css':css_path,
               'sections': root.childs,
               'provider':provider}
    
    return render_to_response('project_content.html', context)


def formTable(request):
    
    # Resources #
    favicon_path = getResource('favicon')
    logo_help = getResource('help')
    logo_browse = getResource('browse')
    jquery_path = os.path.join(settings.STATIC_URL, 'js/jquery.js')
    jsForm_path = os.path.join(settings.STATIC_URL, 'js/form.js')
    css_path = os.path.join(settings.STATIC_URL, 'css/formTable.css')
    #############
    
    ## Project Id(or Name) should be stored in SESSION
    manager = Manager()
    project_name = request.GET.get('project_name')
    projPath = manager.getProjectPath(project_name)
    project = Project(projPath)
    project.load()
    
    protocolName = request.GET.get('protocol', None)
    if protocolName is None:
        protId = request.GET.get('protocolId', None)
        protocol = project.mapper.selectById(int(protId))
    else:
        protocolClass = emProtocolsDict.get(protocolName, None)
        protocol = protocolClass()
    
    #TODO: Add error page validation when protocol is None
    for section in protocol._definition.iterSections():
        for paramName, param in section.iterParams():
            protVar = getattr(protocol, paramName, None)
            if protVar is None:
                raise Exception("_fillSection: param '%s' not found in protocol" % paramName)
                # Create the label
            param.htmlValue = protVar.get(param.default.get(""))
            param.htmlCond = param.condition.get()
            param.htmlDepend = ','.join(param._dependants)
            param.htmlCondParams =  ','.join(param._conditionParams)
#            if isinstance(param, EnumParam):
#                param.htmlChoicesSize = len(param.choices)
    
    
    context = {'protocol':protocol,
               'definition': protocol._definition,
               'favicon': favicon_path,
               'help': logo_help,
               'form': jsForm_path,
               'jquery': jquery_path,
               'browse': logo_browse,
               'css':css_path}
    
    return render_to_response('formTable.html', context)


if __name__ == '__main__':
    root = loadProtTree()    
    for s in root.childs:
        print s.name, '-', s.tag
        for p in s.childs:
            print p.name, '-', p.tag
