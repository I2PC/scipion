# **************************************************************************
# *
# * Authors:    Jose Gutierrez (jose.gutierrez@cnb.csic.es)
# *             Adrian Quintana (aquintana@cnb.csic.es)   
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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************

import os
import xmipp
import json
import mimetypes
from django.shortcuts import render_to_response
from django.http import HttpResponse, HttpResponseForbidden, HttpResponseNotFound
from django.core.servers.basehttp import FileWrapper

from pyworkflow.em import emProtocolsDict
from pyworkflow.web.pages import settings as django_settings
from pyworkflow.manager import Manager
from pyworkflow.project import Project
from pyworkflow.utils import *
from pyworkflow.gui import getImage, getPILImage
from pyworkflow.dataset import COL_RENDER_IMAGE, COL_RENDER_VOLUME

iconDict = {
            'logo_scipion': 'scipion_logo_small_web.png',
            'logo_scipion_small': 'scipion_logo.png',
            'logo_scipion_normal': 'scipion_logo_normal.png',
            'logo_scipion_transparent': 'scipion_logo_transparent.png',
            'favicon': 'favicon.png',
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
            'no_image': 'no-image.png',
            'loading' : 'loading.gif',
            'error_page' : 'error_page.jpg',
            
            #Extensions file
            'folder': 'fa-folder-open.png',
            'file_normal': 'fa-file-o.png',
            'file_text':'file_text.gif',
            'file_image':'file_image.gif',
            'file_python': 'file_python.gif',
            'file_java':'file_java.gif',
            'file_md':'file_md.gif',
            'file_sqlite':'file_sqlite.gif',
            'file_vol':'file_vol.gif',
            'file_stack':'file_stack.gif',
            
            #MyFirstMap
            '1st': 'myfirstmap/1st.png',
            
            }

cssDict = {'project_content': 'project_content_style.css',
           'messi': 'messi.css',
           'projects': 'projects_style.css',
           'showj': 'showj_style.css',
           'general': 'general_style.css',
           'general_flex': 'general_style_base_flex.css',
           'general_grid': 'general_style_base_grid.css',
           'form': 'form.css',
           'ui_smoothness': 'jquery-ui_smoothness.css',
           'jquery_ui': 'jquery-ui.css',
           'showj_demo_table_jui': 'demo_table_jui.css',
           'wizard': 'wizard_style.css',
           'font_awesome': 'font-awesome/css/font-awesome.min.css'
           }

jsDict = {'jquery': 'jquery/jquery.js',
          'jquery_cookie': 'jquery/jquery.cookie.js',
          'jquery_treeview': 'jquery/jquery.treeview.js',
          'jquery_datatables': 'jquery/jquery.dataTables.js',
          'jquery_editable': 'jquery/jquery.jeditable.js',
          'jquery_sizes': 'jquery/jquery.sizes.js',
          'jquery_layout': 'jquery/jquery.jlayout.js',
          'jlayout_border': 'jquery/jlayout.border.js',
          'jquery_ui': 'jquery/jquery-ui.js',
          'jquery_ui_touch': 'jquery/jquery.ui.touch-punch.min.js',
          'jquery_hover_intent': 'jquery/jquery.hoverIntent.minified.js',
          'jquery_browser':'jquery/jquery.serverBrowser.js',
          
          'config': 'templates_libs/config.js',
          'utils': 'templates_libs/utils.js',
          'host_utils': 'templates_libs/host_utils.js',
          'graph_utils': 'templates_libs/graph_utils.js',
          'project_content_utils':'templates_libs/project_content_utils.js',
          'data_content_utils':'templates_libs/data_content_utils.js',
          'project_utils': 'templates_libs/project_utils.js',
          'protocol_form_utils': 'templates_libs/protocol_form_utils.js',
          'wizard_utils': 'templates_libs/wizard_utils.js',
          'upload_utils': 'templates_libs/upload_utils.js',

#          'tabs_config': 'tabs_config.js',
          'jquery_colreorder': 'showj_libs/colReorder.js',
          'jquery_colreorder_resize': 'showj_libs/colReorderWithResize.js',
          'jquery_waypoints': 'showj_libs/waypoints.min.js',
          'transpose': 'showj_libs/transpose.js',
          
          'showj_utils': 'showj_libs/showj_utils.js',
          'showj_menu_utils': 'showj_libs/showj_menu_utils.js',
          'showj_gallery_utils': 'showj_libs/showj_gallery_utils.js',
          'showj_table_utils': 'showj_libs/showj_table_utils.js',
          'showj_column_utils': 'showj_libs/showj_column_utils.js',
          
          #Utils libs
          'jsplumb': 'jsPlumb.js',
          'messi': 'messi.js',
          'raphael': 'raphael.js',
          'philogl': 'PhiloGL-1.3.0.js',
          
          #JSmol
          'jsmol': 'jsmol/JSmol.min.js',
          'jsmolFolder': 'jsmol/j2s',
          
          }

def getResourceIcon(icon):
    return os.path.join(django_settings.MEDIA_URL, iconDict[icon])

def getResourceLogo(logo):
    return os.path.join(django_settings.MEDIA_URL, logo)

def getResourceCss(css):
    return os.path.join(django_settings.STATIC_URL, "css/", cssDict[css])

def getResourceJs(js):
    return os.path.join(django_settings.STATIC_URL, "js/", jsDict[js])

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
        protocol = project.getProtocol(int(protId))
        
        # Remove this and create loadProtocol method in project class
        protocol.setProject(project) 
    else:
        protocolClass = emProtocolsDict.get(protClass, None)
        protocol = project.newProtocol(protocolClass)
        
    return (project, protocol)
#===============================================================================
# Browse to relations objects
#===============================================================================

def browse_relations(request):
    """ Browse relation objects from the database. """
    if request.is_ajax():
        projectName = request.session['projectName']
        project = loadProject(projectName)

        relationName = request.GET.get('relationName')
        attributeName = request.GET.get('attributeName')
        protId = request.GET.get('protId')
        direction = request.GET.get('direction')
        
        protocol = project.getProtocol(int(protId))
        item = protocol.getAttributeValue(attributeName)

        objs = {}
        for obj in project.getRelatedObjects(relationName, item, direction):
            objs[obj.getObjId()] = {"nameId":obj.getNameId(), "info": str(obj)}


        jsonStr = json.dumps(objs, ensure_ascii=False)
        return HttpResponse(jsonStr, mimetype='application/javascript')
    
#===============================================================================
# Browse objects
#===============================================================================

def browse_objects(request):
    """ Browse objects from the database. """
    if request.is_ajax():
        objClassList = request.GET.get('objClass')
        projectName = request.session['projectName']
        
        objFilterParam = request.GET.get('objFilter', None)
        filterObject = FilterObject(objFilterParam)
        
        project = loadProject(projectName)

        objs = {}
        for obj in project.iterSubclasses(objClassList, filterObject.objFilter):
                objs[obj.getObjId()] = {"nameId":obj.getNameId(), "info": str(obj)} 
        
        jsonStr = json.dumps(objs, ensure_ascii=False)
        return HttpResponse(jsonStr, mimetype='application/javascript')
   

class FilterObject():
    def __init__(self, condition):
        self.condition = None if condition == 'None' else condition
        
    def objFilter(self, obj):
        result = True
        if self.condition:
            result = obj.evalCondition(self.condition)
        return result
    
#===============================================================================
# Browse protocols like objects
#===============================================================================
    
def browse_protocol_class(request):
    if request.is_ajax():
        protClassName = request.GET.get('protClassName')
        from pyworkflow.em import findSubClasses, emProtocolsDict
        objs = findSubClasses(emProtocolsDict, protClassName).keys()
        
        jsonStr = json.dumps({'objects' : objs},ensure_ascii=False)
        return HttpResponse(jsonStr, mimetype='application/javascript')

def get_attributes(request):
    if request.is_ajax():
        projectName = request.session['projectName']
        project = loadProject(projectName)
        objId = request.GET.get('objId', None)
        
        obj = project.getProtocol(int(objId))
        if obj is None:
            obj = project.getProtocol(int(objId)).get()
            
        res = obj.getObjLabel() + "_-_" + obj.getObjComment()
        return HttpResponse(res, mimetype='application/javascript')
    
def set_attributes(request):
    if request.is_ajax():
        id = request.GET.get('id', None)
        
        # New values
        label = request.GET.get('label', None)
        comment = request.GET.get('comment', None)

        projectName = request.session['projectName']
        project = loadProject(projectName)

        obj = project.getProtocol(int(id))
        if obj is None:
            obj = project.getProtocol(int(id)).get()
                
        obj.setObjLabel(label)
        obj.setObjComment(comment)
        
        # Save the protocol 
        project._storeProtocol(obj)
        
    return HttpResponse(mimetype='application/javascript')

def file_viewer(request):
    file = request.GET.get("path")
    html = textfileViewer('title', file)
    return HttpResponse(html, mimetype='application/javascript')

def textfileViewer(title, file):
    f = open(file, 'r')
        
    style = "background-color:black;color:white;font-family:Monospace;padding:1em;font-size:90%;"
    title = "<title>"+ title + "</title>"
    html = "<div style="+ style +">"+ title
    
    x = 0
    while 1:
        line = f.readline()
        
        if not line:
            break
        if len(line) > 1:
            x = x+1
            html = html +"<p><span style='color:cyan;'>" + str(x) + ":    </span>"+ line +" </p>"
            
    html = html + "</div>"
    
    return html

def get_log(request):
    "Return a response with the content of the file mentioned in ?path=fname"
    # Got the idea from here:
    # https://stackoverflow.com/questions/8600843
    path = request.GET.get("path")

    # First some simple security: only allow to serve log files
    if not any(path.endswith(x) for x in ['.log', '.stdout', '.stderr']):
        return HttpResponseForbidden('Forbidden: Sorry, invalid path requested.')

    if not os.path.exists(path):
        return HttpResponseNotFound('Path not found: %s' % path)

    response = HttpResponse(FileWrapper(open(path)),
                            content_type=mimetypes.guess_type(path)[0])
    response['Content-Length'] = os.path.getsize(path)
    response['Content-Disposition'] = 'attachment; filename=%s' % path
    return response

def get_file(request):
    "Return a response with the content of the file mentioned in ?path=fname"
    path = request.GET.get("path")
    filename = request.GET.get("filename", path)
    

    if not os.path.exists(path):
        return HttpResponseNotFound('Path not found: %s' % path)

    response = HttpResponse(FileWrapper(open(path)),
                            content_type=mimetypes.guess_type(path)[0])
    response['Content-Length'] = os.path.getsize(path)
    response['Content-Disposition'] = 'attachment; filename=%s' % filename
    return response


def download_output(request):
    if request.is_ajax():
        projectName = request.session['projectName']
        project = loadProject(projectName)
        objId = request.GET.get('objId', None)
        
        obj = project.getProtocol(int(objId))
        if obj is None:
            obj = project.getProtocol(int(objId)).get()
            
        files = obj.getFiles()
        
        import zipfile
        z = zipfile.ZipFile("output.zip", "w")
        
        for f in files:
            z.write(f, arcname=os.path.basename(f))
        z.close()
        
        pathFile = os.path.join(request.session['projectPath'], "output.zip")
        
        return HttpResponse(pathFile, mimetype='application/javascript')
    
    
def render_column(request):
    
    renderFunction = request.GET.get("renderFunc")
    
    #PAJM: No se puede llamar a una funcion con reflex sino pertenece auna clase
    if renderFunction == "get_image":
        return get_image(request)
    elif renderFunction == "get_slice":
        return get_slice(request)
    elif renderFunction == "get_image_psd":
        from pyworkflow.web.app.em_wizard import get_image_psd
        return get_image_psd(request)
    elif renderFunction == "getTestPlot":
        return getTestPlot(request)
#    return getattr(self, renderFunction)

      
def get_image_plot(request):
    from PIL import Image
    imagePath = os.path.join(request.GET.get('image'))
    img = Image.open(imagePath)
    
    response = HttpResponse(mimetype="image/png")
    
    # Create and save the image
    img.save(response, "PNG")
    # after the image is removed from the file system
    os.remove(imagePath)
    
    return response   
    
def get_image(request):
    imageNo = None
    
    # TO DO: Change the way to obtain the separate string of the imagePath
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
            img = getImage(imagePath, tkImage=False)
        else:
            if '@' in imagePath:
                parts = imagePath.split('@')
                imageNo = parts[0]
                imagePath = parts[1]
    
            if 'projectPath' in request.session:
                imagePathTmp = os.path.join(request.session['projectPath'], imagePath)
                if not os.path.isfile(imagePathTmp):
                    raise Exception('should not use getInputPath')
                    #imagePath = getInputPath('showj', imagePath)      
    
            if imageNo:
                imagePath = '%s@%s' % (imageNo, imagePath) 
                
            imgXmipp = xmipp.Image()
            imgXmipp.readPreview(imagePath, int(imageDim))
            
            #===================================================================
            # Transform Matrix
            if applyTransformMatrix: 
                takarras=[tMatrix[0][0], tMatrix[0][1], tMatrix[0][2], x if x!=None else 0,
                tMatrix[1][0], tMatrix[1][1], tMatrix[1][2], y if y!=None else 0,
                tMatrix[2][0], tMatrix[2][1], tMatrix[2][2], z if z!=None else 0]
#               imgXmipp.applyTransforMatScipion(matrix, onlyShifts, wrap)
                imgXmipp.applyTransforMatScipion(takarras, onlyShifts, wrap)
            #===================================================================
            
            #===================================================================
            # Invert Y axis
            if mirrorY: 
                imgXmipp.mirrorY()
            #===================================================================
            
            #TO DO: PSD FIX
            if imagePath.endswith('.psd'):
                imgXmipp.convertPSD()
            
            # from PIL import Image
            img = getPILImage(imgXmipp, None)
    except Exception:
        img = getImage(findResource(getResourceIcon("no_image")), tkImage=False)


    response = HttpResponse(mimetype="image/png")
    img.save(response, "PNG")
    return response

def get_slice(request):
    imageNo = None
    sliceNo = None
    imagePath = request.GET.get('image')
    imageDim = request.GET.get('dim', 150)
    mirrorY = 'mirrorY' in request.GET
    
#    applyTransformMatrix = 'applyTransformMatrix' in request.GET
#    onlyApplyShifts = request.GET.get('onlyApplyShifts',False)
#    wrap = request.GET.get('wrap',False)
#    transformMatrix = request.GET.get('transformMatrix',None)

    # This validations not works for volumes into a stack
    if '__slice__' in imagePath:
        parts = imagePath.split('__slice__', 1)
        sliceNo = int(parts[0])
        imagePath = parts[1]
        
#     if '@' in imagePath:
#             parts = imagePath.split('@')
#             imageNo = parts[0]
#             imagePath = parts[1]
            
    if 'projectPath' in request.session:
        imagePathTmp = os.path.join(request.session['projectPath'], imagePath)
#         if not os.path.isfile(imagePathTmp):
#             raise Exception('should not use getInputPath')
            #imagePath = getInputPath('showj', imagePath)
            
#     if imageNo:
#         imagePath = '%s@%s' % (imageNo, imagePath)

    imgXmipp = xmipp.Image()
    
    if sliceNo is None:
        imgXmipp.readPreview(imagePath, int(imageDim))
    else:
        imgXmipp.readPreview(imagePath, int(imageDim), sliceNo)
        
#        if applyTransformMatrix and transformMatrix != None: 
#            imgXmipp.applyTransforMatScipion(transformMatrix, onlyApplyShifts, wrap)
#        
    if mirrorY: 
        imgXmipp.mirrorY()
    
    # from PIL import Image
#   img = getPILImage(imgXmipp, None, False)
    img = getPILImage(imgXmipp)

    
    response = HttpResponse(mimetype="image/png")
    
    if img.mode != 'RGB':
        img = img.convert('RGB')
    img.save(response, "PNG")
    return response

def getImageXdim(request, imagePath):
    return getImageDim(request, imagePath)[0]

def getImageDim(request, imagePath):
    projectPath = request.session['projectPath']
    img = xmipp.Image()
    imgFn = os.path.join(projectPath, imagePath)
    img.read(str(imgFn), xmipp.HEADER)
    return img.getDimensions()

def readDimensions(request, path, typeOfColumn):
    if (typeOfColumn == COL_RENDER_IMAGE or
        typeOfColumn == COL_RENDER_VOLUME):
        return getImageDim(request, path)
    return (300,300,1,1) 

def readImageVolume(request, path, convert, dataType, reslice, axis, getStats):
    _newPath = path
    _stats = None
    
    img = xmipp.Image()
    imgFn = os.path.join(request.session['projectPath'], path)
    
    if not convert and not reslice and not getStats:
        img.read(str(imgFn), xmipp.HEADER)
    else:
        img.read(str(imgFn))
    
    if getStats:
        _stats = img.computeStats()
        
    if convert:
        img.convert2DataType(dataType, xmipp.CW_ADJUST)
         
    if reslice:
        if axis != xmipp.VIEW_Z_NEG:
            img.reslice(axis)    
    
    if convert or reslice:
        fileName, _ = os.path.splitext(path)
        _newPath = '%s_tmp%s' % (fileName, '.mrc')
        img.write(str(_newPath))
    
    return _newPath, _stats

def getTestPlot(request):
    """ Just a test of a custom render function. """
    from pyworkflow.gui.plotter import Plotter
    xplotter = Plotter()
    xplotter.createSubPlot("Particle sorting", "Particle number", "Zscore")
    x = range(100)
    xplotter.plot(x)
    
    canvas = xplotter.getCanvas()
    response = HttpResponse(content_type='image/png')
    canvas.print_png(response)
    return response   

def replacePattern(m, mode):
    g1 = m.group(mode)
    if mode == HYPER_BOLD:
        text = " <b>%s</b> " % g1
    elif mode == HYPER_ITALIC:
        text = " <i>%s</i> " % g1
    elif mode == HYPER_LINK1:
        text = " <a href='%s' target='_blank'>%s</a> " % (g1, g1)
    elif mode == HYPER_LINK2:
        text = " <a href='%s' target='_blank'>%s</a> " % (g1, m.group('link2_label'))
    else:
        raise Exception("Unrecognized pattern mode: " + mode)
    
    return text

def parseText(text, func=replacePattern):
    """ Parse the text adding some basic tags for html.
    Params:
        text: can be string or list, if it is a list, a <br> tag will be generated.
    """
    parsedText = ""
    if isinstance(text, list):
        for itemText in text:
            splitLines=itemText.splitlines(True)
            if len(splitLines) == 0:
                parsedText += '<br />'
            else:    
                for lineText in splitLines:
                    parsedText += parseHyperText(lineText, func)+'<br />'
    else:
        splitLines=text.splitlines(True)
        for lineText in splitLines:
            parsedText += parseHyperText(lineText, func)+'<br />'
#        parsedText = parseHyperText(text, func)
    return parsedText[:-6]    


def savePlot(request, plot):
    projectPath = request.session['projectPath']
    
    name_img = 'image%s.png' % id(plot)
    fn = os.path.join(projectPath,'Tmp', name_img)
    plot.savefig(fn)
    url_plot = "get_image_plot/?image=" + fn
        
    return url_plot

#===============================================================================
# ERROR PAGE
#===============================================================================

def error(request):
    from views_base import base_grid
    context = {"logoScipionNormal": getResourceIcon("logo_scipion_normal"),
               "logoErrorPage":getResourceIcon("error_page"),
               }
    context = base_grid(request, context)
    return render_to_response('error.html', context)

