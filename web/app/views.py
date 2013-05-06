# from scipion.models import *
import os
from django.shortcuts import render_to_response
from pyworkflow.manager import Manager
from pyworkflow.utils.path import findResource
from pyworkflow.utils.utils import prettyDate
from pyworkflow.web.pages import settings

def getResource(request):
    if request == 'logoScipion':
        img = 'images/logo.png'
    elif request == 'favicon':
        img = 'images/scipion_bn.png'
    path = os.path.join(settings.MEDIA_URL, img)
    return path

def projects(request):
    manager = Manager()
#    logo_path = findResource('scipion_logo.png')

    # Resources #
    logo_path = getResource('logoScipion')
    favicon_path = getResource('favicon')
    css_path = os.path.join(settings.MEDIA_URL, 'css/projects_style.css')
    #############
    
    projects = manager.listProjects()
    for p in projects:
        p.pTime = prettyDate(p.mTime)

    context = {'projects': projects,
               'logo': logo_path,
               'favicon': favicon_path,
               'css':css_path}
    
    return render_to_response('projects.html', context)

class TreeItem():
    def __init__(self, name, type, protClass=None):
        self.name = name
        self.type = type
        self.protClass = protClass
        self.childs = []
        
def populateTree(self, tree, prefix, obj, level=0):
    
    for sub in obj:
        text = sub.text.get()
        if text:
            value = sub.value.get(text)
            tag = obj.tag.get('')
        item = TreeItem(text, tag)
        populateTree(self, item, sub, level+1)
        # If have tag 'protocol_base', fill dynamically with protocol sub-classes
        if obj.value.hasValue() and tag == 'protocol_base':
            protClassName = value.split('.')[-1] # Take last part
            prot = emProtocolsDict.get(protClassName, None)
            if prot is not None:
                for k, v in emProtocolsDict.iteritems():
                    if not v is prot and issubclass(v, prot):
                        prot = TreeItem(k, 'protocol_class', protClassName)
                        item.childs.append(prot)                        
        
        
def project_content(request):
    
    # Resources #
    logo_path = getResource('logoScipion')
    favicon_path = getResource('favicon')
    css_path = os.path.join(settings.MEDIA_URL, 'css/project_content_style.css')
    treeview_css_path = os.path.join(settings.MEDIA_URL, 'css/treeview.css')
    jquery_path = os.path.join(settings.MEDIA_URL, 'libs/jquery.js')
    jquery_cookie = os.path.join(settings.MEDIA_URL, 'libs/jquery.cookie.js')
    jquery_treeview = os.path.join(settings.MEDIA_URL, 'libs/jquery.treeview.js')
    launchTreeview = os.path.join(settings.MEDIA_URL, 'libs/launchTreeview.js')
    #############
    manager = Manager()
    projects = manager.listProjects()
    

    context = {'logo': logo_path,
               'favicon': favicon_path,
               'treeview': treeview_css_path,
               'jquery': jquery_path,
               'jquery_cokkie': jquery_cookie,
               'jquery_treeview': jquery_treeview,
               'launchTreeview': launchTreeview,
               'css':css_path,
               'projects': projects}
    
    return render_to_response('project_content.html', context)

