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
import zipfile

import xmipp
import json
import mimetypes
from django.shortcuts import render_to_response
from django.http import HttpResponse, HttpResponseForbidden, HttpResponseNotFound

# Depending on DJANGO version (first is for DJANGO 1.9) second for 1.5.5
try:
    from wsgiref.util import FileWrapper
except ImportError:
    from django.core.servers.basehttp import FileWrapper

import pyworkflow.em as em
from pyworkflow.web.pages import settings as django_settings
from pyworkflow.manager import Manager
from pyworkflow.project import Project
from pyworkflow.utils import *
from pyworkflow.gui import getImage, getPILImage
from pyworkflow.dataset import COL_RENDER_IMAGE, COL_RENDER_VOLUME
from pyworkflow.em.convert import ImageHandler

# CONSTANTS
# Requests methods names
GET = 'GET'
POST = 'POST'

# CONTEXT KEYS
CTX_PROJECT_NAME = 'projectName'
CTX_SERVICE_NAME = 'serviceName'
CTX_PROJECT_PATH = 'projectPath'

# REQUEST PARAMS
PROJECT_NAME = 'p'
SERVICE_NAME = 's'

iconDict = {
    'logo_scipion': 'scipion_logo_small_web.png',
    'logo_scipion_small': 'scipion_logo.png',
    'logo_scipion_normal': 'scipion_logo_normal.png',
    'logo_scipion_transparent': 'scipion_logo_transparent.png',
    'favicon': 'favicon.ico',
    'scipion_mail': 'scipion_mail.png',
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
    'loading': 'loading.gif',
    'error_page': 'error_page.jpg',

    # Extensions file
    'folder': 'fa-folder-open.png',
    'file_normal': 'fa-file-o.png',
    'file_text': 'file_text.gif',
    'file_image': 'file_image.gif',
    'file_python': 'file_python.gif',
    'file_java': 'file_java.gif',
    'file_md': 'file_md.gif',
    'file_sqlite': 'file_sqlite.gif',
    'file_vol': 'file_vol.gif',
    'file_stack': 'file_stack.gif',

    # MyFirstMap
    'importAverages': 'myfirstmap/importAverages.png',
    'useProtocols': 'myfirstmap/useProtocols.png',
    'protForm': 'myfirstmap/protForm.png',
    'summary': 'myfirstmap/summary.png',
    'showj': 'myfirstmap/showj.png',
    'alignVol': 'myfirstmap/alignVol.png',
    'download': 'myfirstmap/download.png',

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
           'font_awesome': 'font-awesome/css/font-awesome.min.css',
           }

jsDict = {'jquery': 'jquery/jquery.js',
          'jquery_cookie': 'jquery/jquery.cookie.js',
          'jquery_treeview': 'jquery/jquery.treeview.js',
          'jquery_datatables': 'DataTables-1.10.6/media/js/jquery.dataTables.js',
          'jquery_editable': 'jquery/jquery.jeditable.js',
          'jquery_sizes': 'jquery/jquery.sizes.js',
          'jquery_layout': 'jquery/jquery.jlayout.js',
          'jlayout_border': 'jquery/jlayout.border.js',
          'jquery_ui': 'jquery/jquery-ui.js',
          'jquery_ui_touch': 'jquery/jquery.ui.touch-punch.min.js',
          'jquery_hover_intent': 'jquery/jquery.hoverIntent.minified.js',
          'jquery_browser': 'jquery/jquery.serverBrowser.js',

          'config': 'config.js',
          'utils': 'utils.js',
          'host_utils': 'host_utils.js',
          'graph_utils': 'graph_utils.js',
          'project_content_utils': 'project_content_utils.js',
          'data_content_utils': 'data_content_utils.js',
          'project_utils': 'project_utils.js',
          'protocol_form_utils': 'protocol_form_utils.js',
          'wizard_utils': 'wizard_utils.js',
          'upload_utils': 'upload_utils.js',
          'download_utils': 'download_utils.js',

          #          'tabs_config': 'tabs_config.js',
          'jquery_colreorder': 'DataTables-1.10.6/extensions/ColReorder/js/dataTables.colReorder.js',
          # 'jquery_colreorder_resize': 'showj_libs/colReorderWithResize.js',
          'jquery_waypoints': 'showj_libs/waypoints.min.js',
          'transpose': 'showj_libs/transpose.js',

          'showj_utils': 'showj_libs/showj_utils.js',
          'showj_menu_utils': 'showj_libs/showj_menu_utils.js',
          'showj_gallery_utils': 'showj_libs/showj_gallery_utils.js',
          'showj_table_utils': 'showj_libs/showj_table_utils.js',
          'showj_column_utils': 'showj_libs/showj_column_utils.js',

          # Utils libs
          'jsplumb': 'ext_libs/jsPlumb.js',
          'messi': 'ext_libs/messi.js',
          'raphael': 'ext_libs/raphael.js',
          'philogl': 'ext_libs/PhiloGL-1.3.0.js',

          # JSmol
          'jsmol': 'jsmol/JSmol.min.js',
          'jsmolFolder': 'jsmol/j2s',

          }

# To cache service managers
SERVICE_MANAGERS = {}


def getResourceIcon(icon):
    return getMedia(iconDict[icon])


def getResourceLogo(logo):
    return getMedia(logo)


def getMedia(resource):
    return os.path.join('/' + django_settings.MEDIA_URL, resource)


def getResourceCss(css=None):

    cssFile = ''

    if css is not None:
        cssFile = cssDict[css]

    return getResource(os.path.join("css/", cssFile))


def getResourceJs(js=None):

    jsFile = ''

    if js is not None:
        jsFile = jsDict[js]

    return getResource(os.path.join("js/", jsFile))


def getResource(resource=None):

    if resource is None:
        resource = ''

    return "/" + os.path.join(django_settings.STATIC_URL, resource)

def getVarFromRequest(request, varName):
    value = None

    if request.method == GET and varName in request.GET:
        value = request.GET.get(varName)

    if request.method == POST and varName in request.POST:
        value = request.POST.get(varName)

    # Convert "None" to None
    if value == str(None):
        value = None

    return value


def getProjectContextToGET(request):
    getString = "&" + SERVICE_NAME + "=" + str(getVarFromRequest(request, SERVICE_NAME))
    getString = getString + "&" + PROJECT_NAME + "=" + getVarFromRequest(request, PROJECT_NAME)

    return getString


def getProjectPathFromRequest(request):
    """ Retrieve the project path from the request,
    try first to get is from the GET, then POST and if not from the session.
    Return:
        String: project absolute path

    Parameters
    ----------
    request : request
    """
    projectName = getVarFromRequest(request, PROJECT_NAME)

    if projectName is None:
        # TODO: Remove this case after the refactoring of passing always
        # projectName and serviceName in the request (either in GET or POST)
        projectPath = request.session[CTX_PROJECT_PATH]
        print "WARNING: Project not in request!!! Reading path from SESSION!!! path: " + request.path \
              + ", Referrer: " + request.environ.get('HTTP_REFERER', "No referrer")
    else:

        serviceName = getVarFromRequest(request, SERVICE_NAME)
        # In many cases we already have a "p" as the projectName

        manager = getServiceManager(serviceName)
        projectPath = manager.getProjectPath(projectName)

    return projectPath


def loadProject(request):
    # Get the project path using a function
    projectPath = getProjectPathFromRequest(request)
    return loadProjectFromPath(projectPath)


def loadProjectFromPath(projectPath):
    project = Project(projectPath)
    project.load(chdir=False)
    return project


# def check_project_id(request):
#     result = 0
#     projectName = request.GET.get('code', None)
#     serviceName = request.GET.get('serviceName', None)
#
#     try:
#         if serviceName is None:
#             project = Manager().loadProject(projectName)
#         else:
#             manager = getServiceManager(serviceName)
#             project = manager.loadProject(projectName,
#                               protocolsConf=manager.protocols,
#                               hostsConf=manager.hosts, chdir=False)
#         if project is not None:
#             result = 1
#     except Exception:
#         pass
#
#     return HttpResponse(result, content_type='application/javascript')


def loadProtocolProject(request, requestType=POST):
    """ Retrieve the project and protocol from this request.
    Return:
        (project, protocol) tuple

    Parameters
    ----------
    requestType
    request
    """
    requestDict = getattr(request, requestType)

    protId = requestDict.get("protocolId")
    protClass = requestDict.get("protocolClass")

    # Load the project
    project = loadProject(request)

    # Create the protocol object
    if protId and protId != 'None':  # Case of existing protocol
        protId = requestDict.get('protocolId', None)
        protocol = project.getProtocol(int(protId))

        # Remove this and create loadProtocol method in project class
        protocol.setProject(project)

        # Create a new protocol (added from the menu)
    else:
        protocolClass = em.getProtocols().get(protClass, None)
        protocol = project.newProtocol(protocolClass)
        loadProtocolConf(protocol)

    return project, protocol


def getServiceManager(serviceName):
    global SERVICE_MANAGERS

    manager = SERVICE_MANAGERS.get(serviceName, None)

    if manager is None:
        print "...Creating new cached Manager"
        scipionUserData = os.environ['SCIPION_USER_DATA']
        servicePath = '' if serviceName is None else serviceName
        scipionUserData = os.path.join(scipionUserData, servicePath)

        manager = Manager(SCIPION_USER_DATA=scipionUserData)

        serviceConf = os.path.join(os.environ['HOME'], '.config', 'scipion', servicePath)
        manager.config = os.path.join(serviceConf, 'scipion.conf')
        manager.protocols = os.path.join(serviceConf, 'protocols.conf')
        manager.hosts = os.path.join(serviceConf, 'hosts.conf')
        SERVICE_MANAGERS[serviceName] = manager

    return manager


# ===============================================================================
# Browse to relations objects
# ===============================================================================

def browse_relations(request):
    """ Browse relation objects from the database.

    Parameters
    ----------
    request
    """
    if request.is_ajax():
        project = loadProject(request)

        relationName = request.GET.get('relationName')
        attributeName = request.GET.get('attributeName')
        protId = request.GET.get('protId')
        direction = request.GET.get('direction')

        protocol = project.getProtocol(int(protId))
        item = protocol.getAttributeValue(attributeName)

        objs = {}
        for obj in project.getRelatedObjects(relationName, item, direction):
            objs[obj.getObjId()] = {"nameId": obj.getNameId(), "info": str(obj)}

        jsonStr = json.dumps(objs, ensure_ascii=False)
        return HttpResponse(jsonStr, content_type='application/javascript')


# ===============================================================================
# Browse objects
# ===============================================================================

def browse_objects(request):
    """ Browse objects from the database.

    Parameters
    ----------
    request
    """
    if request.is_ajax():
        objClassList = request.GET.get('objClass')

        objFilterParam = request.GET.get('objFilter', None)
        filterObject = FilterObject(objFilterParam, objClassList)

        project = loadProject(request)
        objs = {}

        # Object Filter
        for obj in project.iterSubclasses(objClassList, filterObject.objFilter):
            objParent = project.mapper.selectById(obj.getObjParentId())

            objs[obj.getObjId()] = {"type": "obj",
                                    "nameId": obj.getNameId(),
                                    "objParentName": objParent.getRunName(),
                                    "objId": obj.getObjId(),
                                    "info": str(obj)
                                    }
            # Class Filter
        for obj in project.iterSubclasses("Set", filterObject.classFilter):
            objParent = project.mapper.selectById(obj.getObjParentId())

            context = {"type": "set",
                       "nameId": obj.getNameId(),
                       "objParentName": objParent.getRunName(),
                       "objId": obj.getObjId(),
                       "info": str(obj),
                       "objects": []}
            # Let's set manually now the projectPath
            # Quick fix to have absolute paths for the objects
            obj.projectPath = project.getPath()
            for child in obj.iterItems():
                obj_context = {"nameId": child.getNameId(),
                               "objId": child.getObjId(),
                               "info": str(child)}
                context["objects"].append(obj_context)
            objs[obj.getObjId()] = context

        jsonStr = json.dumps(objs, ensure_ascii=False)
        return HttpResponse(jsonStr, content_type='application/javascript')


class FilterObject:
    def __init__(self, condition, className):
        self.condition = None if condition == 'None' else condition
        self.className = className

    def objFilter(self, obj):
        result = True
        if self.condition:
            result = obj.evalCondition(self.condition)
        return result

    def classFilter(self, obj):
        """ Filter Set with the specified class name.

        Parameters
        ----------
        obj
        """
        itemType = getattr(obj, 'ITEM_TYPE', None)
        return itemType.__name__ == self.className


# ===============================================================================
# Browse protocols like objects
# ===============================================================================

def browse_protocol_class(request):
    if request.is_ajax():
        protClassName = request.GET.get('protClassName')
        from pyworkflow.em import findSubClasses, getProtocols
        objs = findSubClasses(getProtocols(), protClassName).keys()

        jsonStr = json.dumps({'objects': objs}, ensure_ascii=False)
        return HttpResponse(jsonStr, content_type='application/javascript')


def get_attributes(request):
    if request.is_ajax():
        project = loadProject(request)
        objId = request.GET.get('objId', None)

        obj = project.getObject(int(objId))
        if obj is None:
            obj = project.getObject(int(objId)).get()

        res = obj.getObjLabel() + "_-_" + obj.getObjComment()
        return HttpResponse(res, content_type='application/javascript')


def set_attributes(request):
    if request.is_ajax():
        objId = request.GET.get('id', None)

        # New values
        label = request.GET.get('label', None)
        comment = request.GET.get('comment', None)

        project = loadProject(request)

        obj = project.getObject(int(objId))
        if obj is None:
            obj = project.getObject(int(objId)).get()

        obj.setObjLabel(label)
        obj.setObjComment(comment)

        # Save the protocol 
        project._storeProtocol(obj)

    return HttpResponse(content_type='application/javascript')


def file_viewer(request):
    fileToOpen = request.GET.get("path")
    html = textFileViewer('title', fileToOpen)
    return HttpResponse(html, content_type='application/javascript')


def textFileViewer(title, fileToOpen):
    f = open(fileToOpen, 'r')

    style = "background-color:black;color:white;font-family:Monospace;padding:1em;font-size:90%;"
    title = "<title>" + title + "</title>"
    html = "<div style=" + style + ">" + title

    x = 0
    while 1:
        line = f.readline()
        # ? --> code
        if not line:
            break
        if len(line) > 1:
            x += 1
            html = html + "<p><span style='color:cyan;'>" + str(x) + ":    </span>" + line + " </p>"

    html += "</div>"

    return html


def get_log(request):
    """Return a response with the content of the file mentioned in ?path=fname

    Parameters
    ----------
    request: http request
    """
    # Got the idea from here:
    # https://stackoverflow.com/questions/8600843
    pathToLog = request.GET.get("path")

    # First some simple security: only allow to serve log files
    if not any(pathToLog.endswith(x) for x in ['.log', '.stdout', '.stderr']):
        return HttpResponseForbidden('Forbidden: Sorry, invalid path requested.')

    if not os.path.exists(pathToLog):
        return HttpResponseNotFound('Path not found: %s' % pathToLog)

    response = HttpResponse(FileWrapper(open(pathToLog)),
                            content_type=mimetypes.guess_type(pathToLog)[0])
    response['Content-Length'] = os.path.getsize(pathToLog)
    response['Content-Disposition'] = 'attachment; filename=%s' % pathToLog
    return response


def get_file(request):
    """Return a response with the content of the file mentioned in ?path=fname

    Parameters
    ----------
    request: http request
    """
    pathToFile = request.GET.get("path")
    filename = request.GET.get("filename", pathToFile)

    # If pathToFile is not absolute
    if not os.path.isabs(pathToFile):
        projectPath = getProjectPathFromRequest(request)
        pathToFile = os.path.join(projectPath, pathToFile)

    if not os.path.exists(pathToFile):
        return HttpResponseNotFound('Path not found: %s' % pathToFile)

    response = HttpResponse(FileWrapper(open(pathToFile)),
                            content_type=mimetypes.guess_type(pathToFile)[0])
    response['Content-Length'] = os.path.getsize(pathToFile)
    response['Content-Disposition'] = 'attachment; filename=%s' % os.path.basename(filename)
    return response


def delete_file(request):
    """Return a response with the results of the deletion of a file

    Parameters
    ----------
    request: http request
    """
    partialPath = request.GET.get("partialPath")
    projectPath = getProjectPathFromRequest(request)
    fileToDelete = os.path.join(projectPath, partialPath)

    if not os.path.exists(fileToDelete):
        return HttpResponseNotFound('File to delete not found: %s' % partialPath)

    os.remove(fileToDelete)

    return HttpResponse('File deleted %s' % partialPath)


def download_output(request):
    if request.is_ajax():
        project = loadProject(request)
        # This objId is a protocol
        objId = request.GET.get('objId', None)

        protocol = project.getProtocol(int(objId))
        if protocol is None:
            protocol = project.getProtocol(int(objId)).get()

        import zipfile

        # Use absolute path
        # z = zipfile.ZipFile("output.zip", "w")
        outputFile = project.getAbsPath(project.getTmpPath("output.zip"))
        z = zipfile.ZipFile(outputFile, "w")

        # This are relative
        files = protocol.getOutputFiles()

        if (files is not None) and (len(files) > 0):
            for f in files:
                fileName = project.getAbsPath(f)
                z.write(fileName, arcname=os.path.basename(f))
        else:
            print "No output defined for the protocol, returning the whole folder of the protocol instead: ", protocol.getRunName()

            f = project.getAbsPath(protocol.getWorkingDir())

            if os.path.exists(f):

                zipdir(f, outputFile)
            else:

                return HttpResponseNotFound('Output not found or ready.')

        z.close()

        return HttpResponse(outputFile, content_type='application/javascript')


def render_column(request):
    renderFunction = request.GET.get("renderFunc")

    # PAJM: No se puede llamar a una funcion con reflex sino pertenece auna clase
    if renderFunction == "get_image":
        return get_image(request)
    elif renderFunction == "get_slice":
        return get_slice(request)


# Seems to be there is no cases where renderFunction is get_image_psd(used instead get_image_psd url) or getTestPlot
#     elif renderFunction == "get_image_psd":
#         from pyworkflow.web.app.em_wizard import get_image_psd
#         return get_image_psd(request)
#     elif renderFunction == "getTestPlot":
#         return getTestPlot(request)
#    return getattr(self, renderFunction)


def get_image_plot(request):
    from PIL import Image
    imagePath = getImageFullPathFromRequest(request, request.GET.get('image'))
    img = Image.open(imagePath)
    response = HttpResponse(content_type="image/png")
    # Create and save the image
    img.save(response, "PNG")
    # after the image is removed from the file system
    os.remove(imagePath)
    return response


def get_image_path(request):
    from PIL import Image
    imagePath = getImageFullPathFromRequest(request, request.GET.get('image'))
    img = Image.open(imagePath)
    response = HttpResponse(content_type="image/png")
    img.save(response, "PNG")
    return response


def get_image(request):
    imageNo = None

    # TO DO: Change the way to obtain the separate string of the imagePath
    imagePath = request.GET.get('image')
    imageDim = request.GET.get('dim', 150)

    # This prefix can be passed to avoid that image is not refresh when cached by browser (name does not change)
    prefix = request.GET.get('prefix', "")

    mirrorY = 'mirrorY' in request.GET

    applyTransformMatrix = 'applyTransformMatrix' in request.GET
    onlyShifts = 'onlyShifts' in request.GET
    wrap = 'wrap' in request.GET

    matrix = request.GET.get('matrix', None)

    try:
        # PAJM: Como vamos a gestionar las imagen
        if imagePath.endswith('png') or imagePath.endswith('gif'):
            imagePathTmp = getImageFullPathFromRequest(request, prefix + imagePath)
            img = getImage(imagePathTmp, tkImage=False)
        else:
            if '@' in imagePath:
                parts = imagePath.split('@')
                imageNo = parts[0]
                imagePath = parts[1]

            # Get the image full path
            imagePath = getImageFullPathFromRequest(request, imagePath)

            if imageNo:
                imagePath = '%s@%s' % (imageNo, imagePath)

            imgXmipp = xmipp.Image()
            imgXmipp.readPreview(imagePath, int(imageDim))

            # ===================================================================
            # Transform Matrix
            if applyTransformMatrix:
                takarras = [tMatrix[0][0], tMatrix[0][1], tMatrix[0][2], x if x != None else 0,
                            tMatrix[1][0], tMatrix[1][1], tMatrix[1][2], y if y != None else 0,
                            tMatrix[2][0], tMatrix[2][1], tMatrix[2][2], z if z != None else 0]
                #               imgXmipp.applyTransforMatScipion(matrix, onlyShifts, wrap)
                imgXmipp.applyTransforMatScipion(takarras, onlyShifts, wrap)
            # ===================================================================

            # ===================================================================
            # Invert Y axis
            if mirrorY:
                imgXmipp.mirrorY()
            # ===================================================================

            # TO DO: PSD FIX
            if imagePath.endswith('.psd'):
                imgXmipp.convertPSD()

            # from PIL import Image
            img = getPILImage(imgXmipp, None)
    except Exception:
        img = getImage(findResource(getResourceIcon("no_image")), tkImage=False)

    response = HttpResponse(content_type="image/png")
    img.save(response, "PNG")
    return response


def get_slice(request):
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

    # if '@' in imagePath:
    #             parts = imagePath.split('@')
    #             imageNo = parts[0]
    #             imagePath = parts[1]

    imagePath = convertVolume(request, imagePath)
    imgXmipp = xmipp.Image()
    if sliceNo is None:
        imgXmipp.readPreview(imagePath, int(imageDim))
    else:
        imgXmipp.readPreview(imagePath, int(imageDim), sliceNo)

    # if applyTransformMatrix and transformMatrix != None:
    #            imgXmipp.applyTransforMatScipion(transformMatrix, onlyApplyShifts, wrap)
    #
    if mirrorY:
        imgXmipp.mirrorY()

        # from PIL import Image
    #   img = getPILImage(imgXmipp, None, False)
    img = getPILImage(imgXmipp, normalize=False)
    response = HttpResponse(content_type="image/png")

    if img.mode != 'RGB':
        img = img.convert('RGB')
    img.save(response, "PNG")

    return response


def get_image_dim(request):
    pathToImage = request.GET.get('path', None)
    x, y, _, _ = getImageDim(request, pathToImage)

    dimMax = 1024
    dimMin = 256

    xdim = min(max(dimMin, x), dimMax)
    rate = float(x) / xdim
    ydim = int(y / rate)

    objs = [xdim, ydim]

    #     print "x, y = %s:%s" % (x,y)
    #     print "rate: ", rate
    #     print "xdim, ydim = %s:%s" % (xdim,ydim)

    jsonStr = json.dumps(objs, ensure_ascii=False)
    return HttpResponse(jsonStr, content_type='application/javascript')


def getImageDim(request, imagePath):
    imagePath = getImageFullPathFromRequest(request, imagePath)
    from pyworkflow.em.packages.xmipp3.convert import xmippToLocation
    location = xmippToLocation(imagePath)
    x, y, z, n = ImageHandler().getDimensions(location)
    return x, y, z, n


def getImageFullPathFromRequest(request, imagePath):
    if not os.path.isabs(imagePath):

        projectPath = getProjectPathFromRequest(request)

        return getImageFullPath(projectPath, imagePath)
    else:
        return imagePath


def getImageFullPathFromProject(project, imagePath):
    return getImageFullPath(project.getPath(), imagePath)


def getImageFullPath(projectPath, imagePath):
    imageNo = None
    if '@' in imagePath:
        parts = imagePath.split('@')
        imageNo = parts[0]
        imagePath = parts[1]
    if not os.path.isabs(imagePath):
        imagePath = os.path.join(projectPath, imagePath)
    if imageNo:
        imagePath = '%s@%s' % (imageNo, imagePath)
    return imagePath


def getImageXdim(request, imagePath):
    xdim = getImageDim(request, imagePath)[0]
    #     print "x dimension: ", xdim
    return xdim


def readDimensions(request, path, typeOfColumn):
    if (typeOfColumn == COL_RENDER_IMAGE or
                typeOfColumn == COL_RENDER_VOLUME):
        return getImageDim(request, path)
    return 300, 300, 1, 1


def getTmpVolumePath(fileName):
    """ Return the temporarly filename of converted volumes
    to be rendered over the web.

    Parameters
    ----------
    fileName
    """
    baseName, _ = os.path.splitext(fileName)
    return '%s_tmp%s' % (baseName, '.mrc')

    # This convert is only used in table mode


def convertVolume(request, pathToVolume):
    imgFn = getImageFullPathFromRequest(request, pathToVolume)
    imgConvertedFn = getTmpVolumePath(imgFn.replace(':mrc', ''))
    # For rendering volume slices over the web, PIL need that
    # the volume is stored with a specific datatype
    # So, we write a temporarly volume the first time
    if not os.path.exists(imgConvertedFn):
        img = xmipp.Image()
        img.read(str(imgFn))
        img.convert2DataType(xmipp.DT_UCHAR, xmipp.CW_ADJUST)
        img.write(imgConvertedFn)
    return imgConvertedFn


def readImageVolume(request, pathToVolume, convert, dataType, reslice, axis, getStats):
    _stats = None
    # _newPath = path

    img = xmipp.Image()
    imgFn = os.path.join(getProjectPathFromRequest(request), pathToVolume)

    if not convert and not reslice and not getStats:
        #         img.read(str(imgFn), xmipp.HEADER)
        pass
    else:
        img.read(str(imgFn))

    if getStats:
        _stats = img.computeStats()

    if convert:
        img.convert2DataType(dataType, xmipp.CW_ADJUST)

    if reslice:
        if axis != xmipp.VIEW_Z_NEG:
            img.reslice(axis)
    _newPath = getTmpVolumePath(imgFn.replace(':mrc', ''))
    if (convert and not os.path.exists(_newPath)) or reslice:
        img.write(_newPath)
    return _newPath, _stats


def getTestPlot():
    """ Just a test of a custom render function.

    """
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
        text = " <a href='%s' target='_blank' style='color:firebrick;'>%s</a> " % (g1, g1)
    elif mode == HYPER_LINK2:
        if g1.startswith("sci-open:"):
            url = 'javascript:launchViewer(%s)' % g1[len("sci-open:"):]
        else:
            url = g1
        text = " <a href='%s' target='_blank' style='color:firebrick;'>%s</a> " % (url, m.group('link2_label'))
    else:
        raise Exception("Unrecognized pattern mode: " + mode)

    return text


def parseText(text, func=replacePattern):
    """ Parse the text adding some basic tags for html.

    Parameters
    ----------
    func
    text: can be string or list, if it is a list, a <br> tag will be generated.
    """
    parsedText = ""
    if isinstance(text, list):
        for itemText in text:
            splitLines = itemText.splitlines(True)
            if len(splitLines) == 0:
                parsedText += '<br />'
            else:
                for lineText in splitLines:
                    parsedText += parseHyperText(lineText, func) + '<br />'
    else:
        splitLines = text.splitlines(True)
        for lineText in splitLines:
            parsedText += parseHyperText(lineText, func) + '<br />'
            #        parsedText = parseHyperText(text, func)
    return parsedText[:-6]


def getImageUrl(filename):
    abs_url = django_settings.ABSOLUTE_URL
    url_plot = "/"
    if len(abs_url) != 0:
        url_plot = getAbsoluteURL()
    url_plot = url_plot + "get_image_path/?image=" + filename
    return url_plot


def savePlot(request, plot, close=False):
    projectPath = getProjectPathFromRequest(request)

    name_img = 'image%s.png' % id(plot)
    fn = os.path.join(projectPath, 'Tmp', name_img)
    plot.savefig(fn)

    if close: plot.close()

    return getImageUrl(fn)


def getAbsoluteURL(additionalPath=None):

    if additionalPath is None:
        additionalPath = ''

    return '/' + django_settings.ABSOLUTE_URL + additionalPath

# ===============================================================================
# ERROR PAGE
# ===============================================================================

def handle404error(request):
    from views_base import base_grid
    context = {"logoScipionNormal": getResourceIcon("logo_scipion_normal"),
               "logoErrorPage": getResourceIcon("error_page"),
               "abs_url": getAbsoluteURL()
               }
    context = base_grid(request, context)
    return render_to_response('error.html', context)

def handle500error(request):

    # So far use the same error page.
    return handle404error(request);


def loadProtocolConf(protocol):
    """ Load some properties of the protocol object if defined in WEB_PROTOCOLS.
    Properties to be read:
        - NumberOfMpi
        - NumberOfThreads
        - submitToQueue
        - queueName

    Parameters
    ----------
    protocol
    """
    from pyworkflow.web.pages.settings import WEB_CONF
    protDict = WEB_CONF['PROTOCOLS'].get(protocol.getClassName(), None)

    if protDict:
        if 'numberOfMpi' in protDict:
            protocol.numberOfMpi.set(protDict.get('numberOfMpi'))

        if 'numberOfThreads' in protDict:
            protocol.numberOfThreads.set(protDict.get('numberOfThreads'))

        if 'hostName' in protDict:
            protocol.hostName.set(protDict.get('hostName'))

        if 'useQueue' in protDict:
            protocol._useQueue.set(protDict.get('useQueue'))

        if 'queueParams' in protDict:
            protocol.setQueueParams(protDict.get('queueParams'))


def zipdir(dirPath=None, zipFilePath=None, includeDirInZip=True):
    """Create a zip archive from a directory.
    Note that this function is designed to put files in the zip archive with
    either no parent directory or just one parent directory, so it will trim any
    leading directories in the filesystem paths and not include them inside the
    zip archive paths. This is generally the case when you want to just take a
    directory and make it into a zip file that can be extracted in different
    locations.

    Keyword arguments:

    dirPath -- string path to the directory to archive. This is the only
    required argument. It can be absolute or relative, but only one or zero
    leading directories will be included in the zip archive.

    zipFilePath -- string path to the output zip file. This can be an absolute
    or relative path. If the zip file already exists, it will be updated. If
    not, it will be created. If you want to replace it from scratch, delete it
    prior to calling this function. (default is computed as dirPath + ".zip")

    includeDirInZip -- boolean indicating whether the top level directory should
    be included in the archive or omitted. (default True)


    taken from http://peterlyons.com/problog/2009/04/zip-dir-python
"""
    if not zipFilePath:
        zipFilePath = dirPath + ".zip"
    if not os.path.isdir(dirPath):
        raise OSError("dirPath argument must point to a directory. "
                      "'%s' does not." % dirPath)
    parentDir, dirToZip = os.path.split(dirPath)

    # Little nested function to prepare the proper archive path
    def trimPath(path):
        archivePath = path.replace(parentDir, "", 1)
        if parentDir:
            archivePath = archivePath.replace(os.path.sep, "", 1)
        if not includeDirInZip:
            archivePath = archivePath.replace(dirToZip + os.path.sep, "", 1)
        return os.path.normcase(archivePath)

    outFile = zipfile.ZipFile(zipFilePath, "w",
                              compression=zipfile.ZIP_DEFLATED)
    for (archiveDirPath, dirNames, fileNames) in os.walk(dirPath):
        for fileName in fileNames:
            filePath = os.path.join(archiveDirPath, fileName)
            outFile.write(filePath, trimPath(filePath))
        # Make sure we get empty directories as well
        if not fileNames and not dirNames:
            zipInfo = zipfile.ZipInfo(trimPath(archiveDirPath) + "/")
            # some web sites suggest doing
            # zipInfo.external_attr = 16
            # or
            # zipInfo.external_attr = 48
            # Here to allow for inserting an empty directory.  Still TBD/TODO.
            outFile.writestr(zipInfo, "")
    outFile.close()
