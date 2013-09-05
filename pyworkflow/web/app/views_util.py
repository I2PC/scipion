import os
import xmipp
import json
from pyworkflow.tests import getInputPath
from pyworkflow.web.pages import settings
from pyworkflow.manager import Manager
from pyworkflow.project import Project
from django.shortcuts import render_to_response
from django.http import HttpResponse

iconDict = {
            'logo_scipion': 'scipion_logo.png',
            'favicon': 'scipion_bn.png',
            'help': 'system_help24.png',
            'browse': 'zoom.png',
            'wizard': 'tools_wizard.png',
            'edit_toolbar': 'edit.gif',
            'copy_toolbar': 'copy.gif',
            'delete_toolbar': 'delete.gif',
            'browse_toolbar': 'run_steps.gif',
            'tree_toolbar': 'tree2.gif',
            'new_toolbar': 'new_object.gif'
            }

cssDict = {'project_content': 'project_content_style.css',
           'messi': 'messi.css',
           'projects': 'projects_style.css',
           'showj': 'showj_style.css',
           'general': 'general_style.css',
           'form': 'form.css',
           'ui_smoothness': 'jquery-ui_smoothness.css',
           'jquery_ui': 'jquery-ui.css',
           'showj_demo_table_jui': 'demo_table_jui.css'
           
           
           }

jsDict = {'jquery': 'jquery/jquery.js',
          'jquery_cookie': 'jquery/jquery.cookie.js',
          'jquery_treeview': 'jquery/jquery.treeview.js',
          'jquery_datatables': 'jquery/jquery.dataTables.js',
          'jquery_editable': 'jquery/jquery.jeditable.js',
          'jquery_ui': 'jquery/jquery-ui.js',
          
          'utils': 'templates_libs/utils.js',
          'host_utils': 'templates_libs/host_utils.js',
          'graph_utils': 'templates_libs/graph_utils.js',
          'project_content_utils':'templates_libs/project_content_utils.js',
          'project_utils': 'templates_libs/project_utils.js',
          'protocols_utils': 'templates_libs/protocols_utils.js',
          'protocol_form_utils': 'templates_libs/protocol_form_utils.js',
          'wizard_utils': 'templates_libs/wizard_utils.js',

          'tabs_config': 'tabs_config.js',
          'jquery_colreorder': 'ColReorder.js',
          'jquery_colreorder_resize': 'ColReorderWithResize.js',
          'jquery_waypoints': 'waypoints.min.js',
          'messi': 'messi.js'
          }

def getResourceIcon(icon):
    return os.path.join(settings.MEDIA_URL, iconDict[icon])

def getResourceCss(css):
    return os.path.join(settings.STATIC_URL, "css/", cssDict[css])

def getResourceJs(js):
    return os.path.join(settings.STATIC_URL, "js/", jsDict[js])

def loadProject(projectName):
    manager = Manager()
    projPath = manager.getProjectPath(projectName)
    project = Project(projPath)
    project.load()
    return project

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


def get_image(request):
#    from django.http import HttpResponse
    from pyworkflow.gui import getImage, getPILImage
#    print "request.session['projectPath']2", request.session['projectPath']
    
    imageNo = None
    imagePath = request.GET.get('image')
    imageDim = request.GET.get('dim', 150)
    
    # PAJM: Como vamos a gestionar lsa imagen    
    if imagePath.endswith('png') or imagePath.endswith('gif'):
        img = getImage(imagePath, tk=False)
    else:
        if '@' in imagePath:
            parts = imagePath.split('@')
            imageNo = parts[0]
            imagePath = parts[1]

        if 'projectPath' in request.session:
            imagePathTmp = os.path.join(request.session['projectPath'], imagePath)
            if not os.path.isfile(imagePathTmp):
                imagePath = getInputPath('showj', imagePath)      

#        imagePath = join(request.session['projectPath'],imagePath)
            
        if imageNo:
            imagePath = '%s@%s' % (imageNo, imagePath) 
            
        imgXmipp = xmipp.Image(imagePath)
        
        if ('mirrorY' in request.GET):
            imgXmipp.mirrorY()
        
        
        # from PIL import Image
        img = getPILImage(imgXmipp, imageDim)
         
    # response = HttpResponse(mimetype="image/png")    
    response = HttpResponse(mimetype="image/png")
    img.save(response, "PNG")
    return response

def get_image_dimensions(projectPath, imagePath):
    from django.http import HttpResponse
    from pyworkflow.gui import getImage
    imageNo = None
#    imagePath = request.GET.get('image')
    
    # PAJM: Como vamos a gestionar lsa imagen    
    if imagePath.endswith('png') or imagePath.endswith('gif'):
        img = getImage(imagePath, tk=False)
    else:
        if '@' in imagePath:
            parts = imagePath.split('@')
            imageNo = parts[0]
            imagePath = parts[1]
            
        if projectPath != '':
            imagePathTmp = os.path.join(projectPath, imagePath)
            if not os.path.isfile(imagePathTmp):
                imagePath = getInputPath('showj', imagePath)      
            

#        imagePath = join(request.session['projectPath'],imagePath)
        
        if imageNo:
            imagePath = '%s@%s' % (imageNo, imagePath) 
            
        imgXmipp = xmipp.Image(imagePath)
        
        return imgXmipp.getDimensions()
        

