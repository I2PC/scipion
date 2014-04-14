# **************************************************************************
# *
# * Authors:    Adrian Quintana (aquintana@cnb.csic.es)
# *             Jose Gutierrez (jose.gutierrez@cnb.csic.es)
# *                
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
from django.http import HttpResponse
from pyworkflow.web.pages import settings
from django.shortcuts import render_to_response
from django.template import RequestContext
from pyworkflow.web.app.views_util import *
from forms import VolVisualizationForm, ShowjForm
from pyworkflow.em import *
from pyworkflow.em.packages.xmipp3.convert import *
from layout_configuration import *

from views_base import * 

# Define some constants
PATH = 'path'
DATASET = 'dataset'

TABLE_NAME = 'blockComboBox'
COLS_CONFIG = 'tableLayoutConfiguration'
COLS_CONFIG_DEFAULT = 'defaultColumnsLayoutProperties'
LABEL_SELECTED = 'labelsToRenderComboBox'

MODE = 'mode'
MODE_GALLERY = 'gallery'
MODE_TABLE = 'table'
MODE_VOL_ASTEX = 'volume_astex'
MODE_VOL_CHIMERA = 'volume_chimera'

ALLOW_RENDER = 'allowRender'
VOL_VALUES = 'volumesToRenderComboBox'
IMG_DIMS = 'imageDimensions'
IMG_ZOOM = 'zoom'
IMG_ZOOM_DEFAULT = 'defaultZoom'

PROJECT_NAME = 'projectName'
PROJECT_PATH = 'projectPath'


def loadDataSet(request, filename, firstTime):
    """ Load the DataSet from the session or from file. Also load some table.
    Params:
        request: web request variable, where dataset can be in session.
        filename: the path from where to load the dataset.
        firstTime: if True, for to load from file
    """
    if firstTime or not DATASET in request.session:
        dataset = loadDatasetXmipp(filename)
    else:
        dataset = request.session[DATASET]
        
    return dataset


def loadColumnsConfig(request, dataset, table, inputParams, extraParams, firstTime):
    """ Load table layout configuration. How to display columns and attributes (visible, render, editable) """ 
    
    if firstTime:
        columns_properties = getExtraParameters(extraParams, table)
        request.session[COLS_CONFIG_DEFAULT] = columns_properties

        print "loadColumnsConfig: creating ColumnsConfig"
        layoutConfig = ColumnsConfig(dataset, table, inputParams[ALLOW_RENDER], columns_properties)
        request.session[COLS_CONFIG] = layoutConfig
        
    
    else:
        layoutConfig = request.session[COLS_CONFIG] 
            
            
    return layoutConfig


def setLabelToRender(request, table, inputParams, extraParams, firstTime):
    """ If no label is set to render, set the first one if exists """
    if (not inputParams.get(LABEL_SELECTED, False) or 
        request.session.get(TABLE_NAME, None) != inputParams.get(TABLE_NAME, None)):
        print "-"*50, "HERE"
        inputParams[COLS_CONFIG].printColumns()
        labelsToRender = inputParams[COLS_CONFIG].getRenderableColumns()
        
        if labelsToRender:
            inputParams[LABEL_SELECTED] = labelsToRender[0]
        else:
            # If there is no image to display and it is initial load, switch to table mode 
            if firstTime and inputParams[MODE] != MODE_TABLE:
                inputParams[MODE] = MODE_TABLE
                showj(request, inputParams, extraParams)
            inputParams[LABEL_SELECTED] = None
    
    print "label: ", inputParams[LABEL_SELECTED]
    table.setLabelToRender(inputParams[LABEL_SELECTED]) 
    

def setRenderingOptions(request, dataset, tableDataset, inputParams):
    
    #Setting the _imageVolName
    _imageVolName = inputParams[VOL_VALUES] if (VOL_VALUES in inputParams and inputParams[VOL_VALUES] != '') else tableDataset.getElementById(0,inputParams[LABEL_SELECTED])
 
    #Setting the _typeOfColumnToRender
    label = inputParams[LABEL_SELECTED]
    _typeOfColumnToRender = inputParams[COLS_CONFIG].getColumnProperty(label, 'columnType')
    
    #Setting the _imageDimensions
    _imageDimensions = readDimensions(request, _imageVolName, _typeOfColumnToRender)
    
    dataset.setNumberSlices(_imageDimensions[2])
    
    if _typeOfColumnToRender == "image":
        isVol = dataset.getNumberSlices() > 1
        #Setting the _convert 
        _convert = isVol and (inputParams[MODE]==MODE_GALLERY or inputParams[MODE]==MODE_VOL_ASTEX or inputParams[MODE]==MODE_VOL_CHIMERA)
        #Setting the _reslice 
        _reslice = isVol and inputParams[MODE]==MODE_GALLERY
        #Setting the _getStats 
        _getStats = isVol and (inputParams[MODE]==MODE_VOL_ASTEX or inputParams[MODE]==MODE_VOL_CHIMERA)
        #Setting the _dataType 
        _dataType = xmipp.DT_FLOAT if isVol and inputParams[MODE]==MODE_VOL_ASTEX else xmipp.DT_UCHAR
        #Setting the _imageVolName and _stats     
        _imageVolName, _stats = readImageVolume(request, _imageVolName, _convert, _dataType, _reslice, int(inputParams['resliceComboBox']), _getStats)
        
        if isVol:
            inputParams[COLS_CONFIG].configColumn(label, renderFunc="get_slice")
            dataset.setVolumeName(_imageVolName)
            
    volPath = os.path.join(request.session[PROJECT_PATH], _imageVolName)
    
    return dataset, volPath, _stats, _imageDimensions
    

def showj(request, inputParams=None, extraParams=None, firstTime=False):
    
    #=TIME CONTROL==============================================================
#    from datetime import datetime
#    start = datetime.now()
#    print "INIT SHOWJ: ", datetime.now()-start
    #===========================================================================
    
    #############
    # WEB INPUT PARAMETERS
    _imageDimensions = '' 
    _imageVolName = ''
    _stats = None
    dataset = None
    table = None
    volPath = None
    
    if not firstTime:
        updateParams = request.POST.copy()
        if inputParams is None:
            inputParams = updateParams
        else:
            inputParams.update(updateParams)

    if inputParams["typeVolume"] != 'pdb':
        # Load the initial dataset from file or session
        dataset =  loadDataSet(request, inputParams[PATH], firstTime)
        # Load the requested table (or the first if no specified)
        table = dataset.getTable(inputParams[TABLE_NAME])
        # Update inputParams to make sure have a valid table name (if using first table)
        inputParams[TABLE_NAME] = dataset.currentTable()
        
        # Load columns configuration. How to display columns and attributes (visible, render, editable)  
        inputParams[COLS_CONFIG] = loadColumnsConfig(request, dataset, table, inputParams, extraParams, firstTime)
    
#        for col in inputParams[COLS_CONFIG].columnsLayout.values():
#            print "val!! ", col.columnLayoutProperties.getValues() 
    
        #If no label is set to render, set the first one if exists
        
        setLabelToRender(request, table, inputParams, extraParams, firstTime)
        
        if not inputParams[LABEL_SELECTED]:
            inputParams[IMG_ZOOM] = 0
            _imageDimensions = None
            dataset.setNumberSlices(0)
        else:
            dataset, volPath, _stats, _imageDimensions = setRenderingOptions(request, dataset, table, inputParams)
    
        #Store variables into session 
        request = storeToSession(request, inputParams, dataset, _imageDimensions)
        
    else:
        inputParams[COLS_CONFIG] = None
        volPath = inputParams[PATH]

    context, return_page = createContextShowj(request, inputParams, dataset, table, _stats, volPath)

    render_var = render_to_response(return_page, RequestContext(request, context))
    
    #=TIME CONTROL==============================================================
#    print "FINISH SHOWJ: ", datetime.now()-start    
    #===========================================================================
    
    return render_var


def storeToSession(request, inputParams, dataset, _imageDimensions):
    #Store dataset and labelsToRender in session 
    request.session[DATASET] = dataset
    request.session[LABEL_SELECTED] = inputParams[LABEL_SELECTED]
    request.session[TABLE_NAME] = inputParams[TABLE_NAME]
    request.session[IMG_DIMS] = _imageDimensions
    
    return request


def createContextShowj(request, inputParams, dataset, table, paramStats, volPath=None):
    showjForm = ShowjForm(dataset,
                          inputParams[COLS_CONFIG],
                          inputParams) # A form bound for the POST data and unbound for the GET
        
    if showjForm.is_valid() is False:
        print showjForm.errors

    context = createContext(dataset, table, inputParams[COLS_CONFIG], request, showjForm)

    if inputParams[MODE]==MODE_VOL_ASTEX or inputParams[MODE]==MODE_VOL_CHIMERA:
        context.update(create_context_volume(request, inputParams, volPath, paramStats))
               
    elif inputParams[MODE]==MODE_GALLERY or inputParams[MODE]==MODE_TABLE or inputParams[MODE]=='column':
        context.update({"showj_alt_js": getResourceJs('showj_' + inputParams[MODE] + '_utils')})
        
    return_page = 'showj/%s%s%s' % ('showj_', showjForm.data[MODE], '.html')
#    return_page = 'showj/showj_base.html'
        
    return context, return_page
    

def createContext(dataset, table, columnsConfig, request, showjForm):
    #Create context to be send
    
    context = {
            IMG_DIMS: request.session[IMG_DIMS] if IMG_DIMS in request.session else 0,
            IMG_ZOOM_DEFAULT: request.session[IMG_ZOOM_DEFAULT] if IMG_ZOOM_DEFAULT in request.session else 0,
            PROJECT_NAME: request.session[PROJECT_NAME],
            'form': showjForm
            }
    
    
    if dataset is not None:
        context.update({DATASET: dataset})
    if columnsConfig is not None:
        context.update({COLS_CONFIG: json.dumps({'columnsLayout': columnsConfig._columnsDict,
                                                               #'colsOrder': columnsConfig.colsOrder
                                                               },
                                                               ensure_ascii=False,
                                                               cls=ColumnPropertiesEncoder)})
    if table is not None:
        context.update({'tableDataset': table})
    
    # showj_base context
    context = base_showj(request, context)
    context.update(context)

    return context


def getExtraParameters(extraParams, table):
#    print "extraParams",extraParams 
    defaultColumnsLayoutProperties = None
    if extraParams != None and extraParams != {}:
        defaultColumnsLayoutProperties = {k.getName(): {} for k in table.iterColumns()}
        for key, value in extraParams.iteritems():

            try:
                label, attribute = key.rsplit('___')
            except Exception, ex:
                raise Exception('Showj Web visualizer: Incorrect extra parameter')

            if table.hasColumn(label):
                defaultColumnsLayoutProperties[label].update({attribute:value})
                
    return defaultColumnsLayoutProperties

###################################################################
########################### BEGIN SAVE & LOAD #####################    
###################################################################
#### Load an Xmipp Dataset ###
def loadDatasetXmipp(path):
    """ Create a table from a metadata. """
    from pyworkflow.em.packages.xmipp3 import XmippDataSet
    if path.endswith('.star'):
        from pyworkflow.em.packages.relion import addRelionLabels
        addRelionLabels(extended=True)
    return XmippDataSet(str(path))

    
def save_showj_table(request):
    if request.is_ajax():
        changes = request.POST.get('changes')
        jsonChanges = json.loads(changes)
        
#        print "jsonChanges",jsonChanges 
        
        dataset=request.session[DATASET]
        blockComboBox=request.session[TABLE_NAME]
#        columnsConfig=request.session[COLS_CONFIG]
        
        table=dataset.getTable(blockComboBox)
        
        for key in jsonChanges:
            element_split = key.rsplit('___')
            if len(element_split)!=2: 
                print "this fails and sth has to be done"
            
            #NAPA de LUXE ahora mismo se realiza una conversion a int pero habria que ver de que tipo de datos se trata 
            #columnsConfig.columnsLayout[element_split[0]].typeOfColumn

            dictelement = {element_split[0]:jsonChanges[key]}
            table.updateRow(int(element_split[1]),**dictelement)
        
        dataset.writeTable(blockComboBox, table)
        
        return HttpResponse(json.dumps({'message':'Ok'}), mimetype='application/javascript')
###################################################################
########################### END SAVE & LOAD #######################    
###################################################################    



###################################################################
######################## DISPATCHER ###############################    
###################################################################  
def visualizeObject(request):
    #Initialize default values
    inputParams = {'objectId': '',
                       PATH: '',
                       ALLOW_RENDER: True,                 # Image can be displayed, depending on column layout 
                       MODE: MODE_GALLERY,                   # Mode Options: gallery, table, column, volume_astex, volume_chimera
                       IMG_ZOOM: '128px',                     # Zoom set by default
                       TABLE_NAME: None,                    # Table name to display. If None the first one will be displayed
                       LABEL_SELECTED: None,        # Column to be displayed in gallery mode. If None the first one will be displayed
                       VOL_VALUES: '',       # If 3D, Volume to be displayed in gallery, volume_astex and volume_chimera mode. If None the first one will be displayed
#                       'dims': '2d',                        # Object Dimensions
                       'goto': 1,                           # Element selected (metadata record) by default. It can be a row in table mode or an image in gallery mode
                       'colRowMode': 'Off',                 # In gallery mode 'On' means columns can be adjust manually by the user. When 'Off' columns are adjusted automatically to screen width.
                       'cols':  '',                         # In gallery mode (and colRowMode set to 'On') cols define number of columns to be displayed
                       'rows': '',                          # In gallery mode (and colRowMode set to 'On') rows define number of columns to be displayed
                       'mirrorY': False,                    # When 'True' image are mirrored in Y Axis 
                       'applyTransformMatrix': False,       # When 'True' if there is transform matrix, it will be applied
                       'onlyShifts': False,                 # When 'True' if there is transform matrix, only shifts will be applied
                       'wrap': False,                       # When 'True' if there is transform matrix, only shifts will be applied
                       'resliceComboBox': xmipp.VIEW_Z_NEG, # If 3D, axis to slice volume 
                       'imageMaxWidth': 512,                # Maximum image width (in pixels)
                       'imageMinWidth': 64,                 # Minimum image width (in pixels)
                       'imageMaxHeight': 512,               # Maximum image height (in pixels)
                       'imageMinHeight': 64,                # Minimum image height (in pixels)
                       'typeVolume': 'map'}                 # If map, it will be displayed normally, else if pdb only astexViewer and chimera display will be available


    # Extra parameters can be used to configure table layout and set render function for a column
    # Default layout configuration is set in ColumnLayoutProperties method in layout_configuration.py
    # 
    # Parameters are formed by: [label]___[property]: [value]. E.g.: id___visible:True or micrograph___renderFunc:"get_image_psd"
    # Properties to be configured are:
    #    visible: Defines if this column is displayed
    #    allowSetVisible: Defines if user can change visible property (show/hide this column).
    #    editable: Defines if this column is editable, ie user can change field value.
    #    allowSetEditable: Defines if user can change editable property (allow editing this column).
    #    renderable: Defines if this column is renderizable, ie it renders data column using renderFunc
    #    allowSetRenderable: Defines if user can change renderable property.
    #    renderFunc: Function to be used when this field is rendered. (it has to be inserted in render_column method)
    #    extraRenderFunc: Any extra parameters needed for rendering. Parameters are passed like in a url ie downsample=2&lowPass=3.5
    #
    # Example:
    # extraParams["id___visible"]=True
    # extraParams["micrograph___renderFunc"]="get_image_psd"
    # extraParams["micrograph___extraRenderFunc"]="downsample=2"
            
    extraParams = {}
    
# NAPA DE LUXE: PARA probarlo sin sesion iniciada    
#    request.session[PROJECT_PATH] = "/home/adrian/Scipion/projects/TestXmippWorkflow/"
#    request.session[PROJECT_NAME] = "Juanitpo"
    
    if PROJECT_PATH in request.session:
        projectPath = request.session[PROJECT_PATH]
    else:         
        raise Exception('Showj Web visualizer: No project loaded')
    
    for key, value in request.GET.iteritems():
        if key in inputParams:
            inputParams[key]=value
        else:
            extraParams[key]=value
    
    
    if "objectId" in request.GET: 
        project = Project(projectPath)
        project.load()
        
        obj = project.mapper.selectById(int(inputParams["objectId"]))
        if obj.isPointer():
            obj = obj.get()

        if isinstance(obj, SetOfMicrographs):
            fn = project.getTmpPath(obj.getName() + '_micrographs.xmd')
            fn = createXmippInputMicrographs(None, obj, micsFn=os.path.join(projectPath, fn))
            inputParams[PATH] = os.path.join(projectPath, fn)
#            writeSetOfMicrographs(obj, fn)
#            inputParams[PATH]= os.path.join(projectPath, fn)
        elif isinstance(obj, SetOfMovies):
            fn = project.getTmpPath(obj.getName() + '_movies.xmd')
            writeSetOfMovies(obj, fn)
            inputParams[PATH]= os.path.join(projectPath, fn)
        elif isinstance(obj, SetOfVolumes):
            fn = project.getTmpPath(obj.getName()+ '_volumes.xmd')
            inputParams[PATH] = os.path.join(projectPath, createXmippInputVolumes(None, obj, volsFn=os.path.join(projectPath, fn))) 
            inputParams[MODE]= MODE_TABLE
            
#            writeSetOfVolumes(obj, fn)
#            inputParams[PATH]= os.path.join(projectPath, fn)

        elif isinstance(obj, PdbFile):
            inputParams[PATH] = obj.getFileName()
            inputParams[MODE] = MODE_VOL_ASTEX
            inputParams['typeVolume'] = 'pdb'

        elif isinstance(obj, SetOfImages):
            fn = project.getTmpPath(obj.getName() + '_images.xmd')
            inputParams[PATH] = os.path.join(projectPath, createXmippInputImages(None, obj, imagesFn=os.path.join(projectPath, fn)))
#            writeSetOfParticles(obj, fn)
#            inputParams[PATH]= os.path.join(projectPath, fn)
        elif isinstance(obj, Image):  
            fn = project.getTmpPath(obj.getName() + '_image.xmd')
            
            md = xmipp.MetaData()
            imgFn = getImageLocation(obj)
            if isinstance(obj, Volume) and imgFn.endswith('.mrc'):
                imgFn += ':mrc'
            md.setValue(xmipp.MDL_IMAGE, imgFn, md.addObject())
            md.write(fn)
            
            inputParams[PATH] = os.path.join(projectPath, fn)
            
        elif isinstance(obj, SetOfClasses2D):
            fn = project.getTmpPath(obj.getName() + '_classes.xmd')
#            writeSetOfClasses2D(obj, fn)
#            inputParams[PATH]= os.path.join(projectPath, fn)
            inputParams[PATH] = os.path.join(projectPath, createXmippInputClasses2D(None, obj, classFn=os.path.join(projectPath, fn)))
        
        elif isinstance(obj, SetOfCTF):
            fn = project.getTmpPath(obj.getName() + '_ctfs.xmd')
            inputParams[PATH] = os.path.join(projectPath, createXmippInputCTF(None, obj, ctfFn=os.path.join(projectPath, fn)))
            
            extraParams["itemId___visible"]= False
#            extraParams["psd___visible"]= False
            extraParams["psdEnhanced___renderable"]= True
            extraParams["micrograph___renderable"]= True
            extraParams["image1___renderable"]= True
            extraParams["image2___renderable"]= True
            
#            writeSetOfCTFs(obj, fn)
#            inputParams[PATH]= os.path.join(projectPath, fn)
        elif isinstance(obj, NormalModes):
            inputParams[PATH] = os.path.join(projectPath, obj.getFileName())
        else:
            raise Exception('Showj Web visualizer: can not visualize class: %s' % obj.getClassName())
    
    elif PATH in request.GET:
        inputParams.update({PATH:os.path.join(projectPath, request.GET.get(PATH))})
        
    else:
        raise Exception('Showj Web visualizer: No object identifier or path found')         

    request.session[IMG_ZOOM_DEFAULT] = inputParams[IMG_ZOOM]
    
    # Clean the dataset if exist 
    if DATASET in request.session:
        request.session.__delitem__(DATASET)
    
    return showj(request, inputParams, extraParams, firstTime=True)  


def testingSSH(request):
    context = {}
#    return render_to_response("testing_ssh.html", RequestContext(request, context))
    return render_to_response("scipion.html", RequestContext(request, context))


def create_context_volume(request, inputParams, volPath, param_stats):
#        volPath = os.path.join(request.session[PROJECT_PATH], _imageVolName)

    threshold = calculateThreshold(param_stats)
    
    context = {"threshold":threshold,
                    'minStats':param_stats[2] if param_stats != None else 1,
                    'maxStats':param_stats[3] if param_stats != None else 1 }
    
    if inputParams[MODE] == MODE_VOL_ASTEX:
        context.update(create_context_astex(request, inputParams['typeVolume'], volPath))

#   'volType': 2, #0->byte, 1 ->Integer, 2-> Float
        
    elif inputParams[MODE] == MODE_VOL_CHIMERA:   
        context.update(create_context_chimera(volPath))
        
    return context
        

def create_context_astex(request, typeVolume, volPath):
#   volPath='/home/adrian/Scipion/tests/input/showj/emd_1042.map'        
#   fileName, fileExtension = os.path.splitext(_imageVolName)
            
#    linkName = 'test_link_' + request.session._session_key + '.' + inputParams["typeVolume"]

    linkName = 'test_link_' + request.session._session_key + '.' + typeVolume
    volLinkPath = os.path.join(pw.WEB_RESOURCES, 'astex', 'tmp', linkName)
    
    cleanPath(volLinkPath)
    createLink(volPath, volLinkPath)
    volLink = os.path.join('/', settings.STATIC_ROOT, 'astex', 'tmp', linkName)
    
    return {"volLink":volLink, 
            "jquery_ui_css": getResourceCss("jquery_ui")}
    
    
def create_context_chimera(volPath):
    from subprocess import Popen, PIPE, STDOUT
    
    p = Popen([os.environ.get('CHIMERA_HEADLESS'), volPath], stdout=PIPE, stdin=PIPE, stderr=PIPE)
    outputHtmlFile = os.path.join(pw.WEB_RESOURCES, 'astex', 'tmp', 'test.html')
    #chimeraCommand= 'volume #0 level ' + str(threshold) + '; export format WebGL ' + outputHtmlFile + '; stop'
    chimeraCommand= 'export format WebGL ' + outputHtmlFile + '; stop'
    stdout_data, stderr_data = p.communicate(input=chimeraCommand)
    f = open(outputHtmlFile)
    chimeraHtml = f.read().decode('string-escape').decode("utf-8").split("</html>")[1]
    
    return {"chimeraHtml":chimeraHtml}
    
    
def calculateThreshold(params):
    threshold = (params[2] + params[3])/2 if params != None else 1
        
    return threshold


def updateSessionTable(request):
    label = request.GET.get('label', None)
    type = request.GET.get('type', None)
    option = request.GET.get('option', None)
    
    if type == "renderable":
        request.session["columnsConfig"].columnsLayout[label].renderable = option
    elif type == "editable":    
        request.session["columnsConfig"].columnsLayout[label].editable = option
        
#    print request.session["columnsConfig"].columnsLayout[label].renderable
     
    return HttpResponse(mimetype='application/javascript')

