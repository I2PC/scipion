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
from pyworkflow.em import emProtocolsDict
from pyworkflow.tests import getInputPath
from pyworkflow.web.pages import settings
from pyworkflow.manager import Manager
from pyworkflow.project import Project
from django.shortcuts import render_to_response
from django.http import HttpResponse
from pyworkflow.utils import *


iconDict = {
            'logo_scipion': 'scipion_logo_small.gif',
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
            'no_image': 'no-image.png'
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
          'jquery_hover_intent': 'jquery/jquery.hoverIntent.minified.js',
          'jsplumb': 'jsPlumb/jquery.jsPlumb.js',
          
          'utils': 'templates_libs/utils.js',
          'host_utils': 'templates_libs/host_utils.js',
          'graph_utils': 'templates_libs/graph_utils.js',
          'project_content_utils':'templates_libs/project_content_utils.js',
          'project_utils': 'templates_libs/project_utils.js',
          'protocol_form_utils': 'templates_libs/protocol_form_utils.js',
          'wizard_utils': 'templates_libs/wizard_utils.js',
          'showj_utils': 'showj_libs/showj_utils.js',

#          'tabs_config': 'tabs_config.js',
          'jquery_colreorder': 'showj_libs/colReorder.js',
          'jquery_colreorder_resize': 'showj_libs/colReorderWithResize.js',
          'jquery_waypoints': 'showj_libs/waypoints.min.js',
          'transpose': 'showj_libs/transpose.js',
          'messi': 'messi/messi.js',
          'raphael': 'raphael/raphael.js'
          
          }

def getResourceIcon(icon):
    return os.path.join(settings.MEDIA_URL, iconDict[icon])

def getResourceLogo(logo):
    return os.path.join(settings.MEDIA_URL, logo)

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
        objClassList = request.GET.get('objClass')
        projectName = request.GET.get('projectName')
        
        objFilterParam = request.GET.get('objFilter',None)
        filterObject = FilterObject(objFilterParam)
        
        project = loadProject(projectName)    

        objs = []
        for objClass in objClassList.split(","):
            for obj in project.mapper.selectByClass(objClass, objectFilter=filterObject.objFilter, iterate=True):
                objs.append(obj.getNameId())
        jsonStr = json.dumps({'objects' : objs},
                             ensure_ascii=False)
        return HttpResponse(jsonStr, mimetype='application/javascript')

class FilterObject():
    def __init__(self, condition):
        self.condition = None if condition == 'None' else condition
        
    def objFilter(self, obj):
        result = True
        if self.condition:
            result = obj.evalCondition(self.condition)
        return result     

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
        obj = project.mapper.selectById(int(objId)).get()
        res = obj.getObjLabel() + "_-_" + obj.getObjComment()
        return HttpResponse(res, mimetype='application/javascript')
    
def set_attributes(request):
    if request.is_ajax():
        id = request.GET.get('id', None)
        label = request.GET.get('label', None)
        comment = request.GET.get('comment', None)
        
        typeObj = request.GET.get('typeObj', None)

        projectName = request.session['projectName']
        project = loadProject(projectName)
        
        if typeObj=='object':
            obj = project.mapper.selectById(int(id)).get()
        elif typeObj=='protocol':
            obj = project.mapper.selectById(int(id))
        
        obj.setObjLabel(label)
        obj.setObjComment(comment)
        
        if typeObj=='object':
            project._storeProtocol(obj)
        elif typeObj=='protocol':
            project.saveProtocol(obj)
#            project.mapper.store(protocol)
        
    return HttpResponse(mimetype='application/javascript')

def getSizePlotter(plots):
    figsize = (800, 600)
    
    if plots == -1:
        figsize = (800, 400)
    if plots == 1:
        figsize = (600, 500)
    elif plots == 2:
        figsize = (600, 400)
    elif plots == 3 or plots == 4:
        figsize = (800, 600)
    
    return figsize

def textfileViewer(title, fileList):
    f = open(fileList[0], 'r')
        
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

def convertTktoHtml(text):
    text = text.replace('\n', '<br/>')
    return text

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

def get_image(request):
#    from django.http import HttpResponse
    from pyworkflow.gui import getImage, getPILImage
    
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
    
    imageNo = None
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
        
        parts = imagePath.split('@',1)
        sliceNo = parts[0]
        imagePath = parts[1]

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


def readDimensions(request, path, typeOfColumn):
    if typeOfColumn=="image":
        img = xmipp.Image()
        imgFn = os.path.join(request.session['projectPath'], path)
        img.read(str(imgFn), xmipp.HEADER)
        return img.getDimensions()
    return (300,300,1,1) 

def readImageVolume(request, path, convert, dataType, reslice, axis, getStats):
    _newPath=path
    _stats=None
    
    img = xmipp.Image()
    imgFn = os.path.join(request.session['projectPath'], path)
    
    if not convert and not reslice and not getStats:
        img.read(str(imgFn), xmipp.HEADER)
    else:
        img.read(str(imgFn))
        
    if convert:
        img.convert2DataType(dataType, xmipp.CW_ADJUST)
         
    if reslice:
        if axis !=xmipp.VIEW_Z_NEG:
            img.reslice(axis)    
    
    if getStats:
        _stats=img.computeStats()
    
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
