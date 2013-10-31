import os
import xmipp
import json
from pyworkflow.em import emProtocolsDict
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
            'stop_toolbar': 'stop.gif',
            'delete_toolbar': 'delete.gif',
            'browse_toolbar': 'run_steps.gif',
            'tree_toolbar': 'tree2.gif',
            'list_toolbar': 'md_view.gif',
            'analyze_toolbar': 'visualize.gif',
            'new_toolbar': 'new_object.gif',
            'no_image': 'no-image.png'
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
          'jquery_hover_intent': 'jquery/jquery.hoverIntent.minified.js',
          
          'utils': 'templates_libs/utils.js',
          'host_utils': 'templates_libs/host_utils.js',
          'graph_utils': 'templates_libs/graph_utils.js',
          'project_content_utils':'templates_libs/project_content_utils.js',
          'project_utils': 'templates_libs/project_utils.js',
          'protocol_form_utils': 'templates_libs/protocol_form_utils.js',
          'wizard_utils': 'templates_libs/wizard_utils.js',

          'tabs_config': 'tabs_config.js',
          'jquery_colreorder': 'ColReorder.js',
          'jquery_colreorder_resize': 'ColReorderWithResize.js',
          'jquery_waypoints': 'waypoints.min.js',
          'messi': 'messi/messi.js',
          'raphael': 'raphael/raphael.js'
          
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
    
    imageNo = None
    imagePath = request.GET.get('image')
    imageDim = request.GET.get('dim', 150)
    
    mirrorY = 'mirrorY' in request.GET
    
    applyTransformMatrix = 'applyTransformMatrix' in request.GET
    onlyShifts = 'onlyShifts' in request.GET
    wrap = 'wrap' in request.GET
    
     
    matrix = request.GET.get('matrix',None)
        
    try:
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
    
            if imageNo:
                imagePath = '%s@%s' % (imageNo, imagePath) 
                
            #imgXmipp = xmipp.Image(imagePath)
            imgXmipp = xmipp.Image()
    
            imgXmipp.readPreview(imagePath, int(imageDim))
            if applyTransformMatrix: 
                print "akitikiri"
                takarras=[tMatrix[0][0], tMatrix[0][1], tMatrix[0][2], x if x!=None else 0,
                tMatrix[1][0], tMatrix[1][1], tMatrix[1][2], y if y!=None else 0,
                tMatrix[2][0], tMatrix[2][1], tMatrix[2][2], z if z!=None else 0]
#                imgXmipp.applyTransforMatScipion(matrix, 
#                                                 onlyShifts,
#                                                 wrap)
                imgXmipp.applyTransforMatScipion(takarras, 
                                                 onlyShifts,
                                                 wrap)
            
            if mirrorY: 
                imgXmipp.mirrorY()
            
            # from PIL import Image
            img = getPILImage(imgXmipp, None)
    except Exception:
        from pyworkflow import findResource
        img = getImage(findResource(getResourceIcon("no_image")), tk=False)


    response = HttpResponse(mimetype="image/png")
    img.save(response, "PNG")
    return response

def get_slice(request):
#    from django.http import HttpResponse
    from pyworkflow.gui import getImage, getPILImage
    
    sliceNo = None
    imagePath = request.GET.get('image')
    imageDim = request.GET.get('dim', 150)
    mirrorY = 'mirrorY' in request.GET
#    applyTransformMatrix = 'applyTransformMatrix' in request.GET
#    onlyApplyShifts = request.GET.get('onlyApplyShifts',False)
#    wrap = request.GET.get('wrap',False)
#    transformMatrix = request.GET.get('transformMatrix',None)
#    
      
    try:
            # PAJM: Como vamos a gestionar lsa imagen    
        if not '@' in imagePath:
            raise Exception('Slice number required.')
        
        parts = imagePath.split('@')
        sliceNo = parts[0]
        imagePath = parts[1]
    
        if 'projectPath' in request.session:
            imagePathTmp = os.path.join(request.session['projectPath'], imagePath)
            if not os.path.isfile(imagePathTmp):
                imagePath = getInputPath('showj', imagePath)

        imgXmipp = xmipp.Image()
        imgXmipp.readPreview(imagePath, int(imageDim), int(sliceNo))
                
    #        if applyTransformMatrix and transformMatrix != None: 
    #            imgXmipp.applyTransforMatScipion(transformMatrix, onlyApplyShifts, wrap)
    #        
        if mirrorY: 
            imgXmipp.mirrorY()
        
        
        # from PIL import Image
        img = getPILImage(imgXmipp, None, False)
    except Exception:
        from pyworkflow import findResource
        img = getImage(findResource(getResourceIcon("no_image")), tk=False)         
    
    response = HttpResponse(mimetype="image/png")
    img.save(response, "PNG")
    return response

def getImageXdim(request, imagePath):
    return getImageDim(request, imagePath)[0]

def getImageDim(request, imagePath):
    img = xmipp.Image()
    imgFn = os.path.join(request.session['projectPath'], imagePath)
    img.read(str(imgFn), xmipp.HEADER)
    return img.getDimensions()


def readImageVolume(request, path, getImageDim, convert, dataType, reslice, axis, getStats):
    _newPath=path
    _imageDim=None
    _stats=None
    
    img = xmipp.Image()
    imgFn = os.path.join(request.session['projectPath'], path)
    
    if not convert and not reslice and not getStats and not getImageDim:
        img.read(str(imgFn), xmipp.HEADER)
    else:
        img.read(str(imgFn))
        
    if getImageDim:
         _imageDim = img.getDimensions()   
        
    if convert:
         img.convert2DataType(dataType, xmipp.CW_ADJUST)
         
    if reslice:
        if axis !=xmipp.VIEW_Z_NEG:
            img.reslice(axis)    
    
    if getStats:
        _stats=img.computeStats()
    
    if convert or reslice:
        fileName, fileExtension = os.path.splitext(path)
        _newPath = '%s_tmp%s' % (fileName, '.mrc')
        img.write(str(_newPath))
    
    return _newPath, _imageDim, _stats
    
    
