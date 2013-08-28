from django.shortcuts import render_to_response
from pyworkflow.tests import getInputPath
import os
import xmipp
from pyworkflow.web.pages import settings
from pyworkflow.manager import Manager
from pyworkflow.project import Project

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
           'form': 'form.css'
           }

jsDict = {'jquery': 'jquery.js',
          'jquery_cookie': 'jquery.cookie.js',
          'jquery_treeview': 'jquery.treeview.js',
          'utils': 'utils.js',
          'host_util': 'host_utils.js',
          'tabs_config': 'tabs_config.js',
          'project_form': 'project_form.js',
          'jquery_datatables': 'jquery.dataTables.js',
          'jquery_colreorder': 'ColReorder.js',
          'jquery_editable': 'jquery.jeditable.js',
          'jquery_ui': 'jquery-ui.js',
          'jquery_waypoints': 'waypoints.min.js',
          'form': 'form.js',
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
            imagePathTmp = os.path.join(projectPath,imagePath)
            if not os.path.isfile(imagePathTmp):
                imagePath = getInputPath('showj', imagePath)      
            

#        imagePath = join(request.session['projectPath'],imagePath)
        
        if imageNo:
            imagePath = '%s@%s' % (imageNo, imagePath) 
            
        imgXmipp = xmipp.Image(imagePath)
        
        return imgXmipp.getDimensions()
        

def get_image(request):
    from django.http import HttpResponse
    from pyworkflow.gui import getImage, getPILImage
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
            imagePathTmp = os.path.join(request.session['projectPath'],imagePath)
            if not os.path.isfile(imagePathTmp):
                imagePath = getInputPath('showj', imagePath)      
            

#        imagePath = join(request.session['projectPath'],imagePath)
        
        if imageNo:
            imagePath = '%s@%s' % (imageNo, imagePath) 
            
        imgXmipp = xmipp.Image(imagePath)
        
        # from PIL import Image
        img = getPILImage(imgXmipp, imageDim)
        
        
        
    # response = HttpResponse(mimetype="image/png")    
    response = HttpResponse(mimetype="image/png")
    img.save(response, "PNG")
    return response