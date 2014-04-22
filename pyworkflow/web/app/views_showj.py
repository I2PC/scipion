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
from pyworkflow.em.viewer import *
from layout_configuration import *
from views_base import * 


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
        dataset = request.session[filename][DATASET]
        
    return dataset


def hasTableChanged(request, inputParams):
    return request.session.get(inputParams[PATH], {}).get(TABLE_NAME, None) != inputParams.get(TABLE_NAME, None)    
    
    
def loadColumnsConfig(request, dataset, table, inputParams, extraParams, firstTime):
    """ Load table layout configuration. How to display columns and attributes (visible, render, editable) """ 
    tableChanged = hasTableChanged(request, inputParams)
    
    # Clear table name and selected volume after a table change
    if not firstTime and tableChanged: 
        print "setting VOL_SELECTED = None"
        inputParams[VOL_SELECTED] = None
        
    if firstTime or tableChanged:
        columns_properties = getExtraParameters(extraParams, table)
        
        request.session[COLS_CONFIG_DEFAULT] = columns_properties

        columnsConfig = ColumnsConfig(dataset, table, inputParams[ALLOW_RENDER], columns_properties)
        
        request.session[COLS_CONFIG] = columnsConfig
    else:
        columnsConfig = request.session[COLS_CONFIG]
            
    inputParams[COLS_CONFIG] = columnsConfig
    
    setLabelToRender(request, table, inputParams, extraParams, firstTime)
     
     
def addProjectPrefix(request, fn):
    """ Split path in block and filename and add the project path to filename. """
    if PROJECT_PATH in request.session:
        projectPath = request.session[PROJECT_PATH]
    else:         
        raise Exception('Showj Web visualizer: No project loaded')
    
    return join(projectPath, fn)


def setLabelToRender(request, table, inputParams, extraParams, firstTime):
    """ If no label is set to render, set the first one if exists """
    if (not inputParams.get(LABEL_SELECTED, False) or 
        hasTableChanged(request, inputParams)):
        labelsToRender = inputParams[COLS_CONFIG].getRenderableColumns()
        
        if labelsToRender:
            inputParams[LABEL_SELECTED] = labelsToRender[0]
        else:
            # If there is no image to display and it is initial load, switch to table mode 
            if firstTime:
                inputParams[MODE] = MODE_TABLE
                # FIXME: we really need this other call to showj???
                #showj(request, inputParams, extraParams)
            inputParams[LABEL_SELECTED] = None
    
    table.setLabelToRender(inputParams[LABEL_SELECTED]) 
    

def setRenderingOptions(request, dataset, table, inputParams):
    """ Read the first renderizable item and setup some variables.
    For example, if we are in volume mode or not.
    """
    #Setting the _typeOfColumnToRender
    label = inputParams[LABEL_SELECTED]
    
    if not label:
        volPath = None
        _imageDimensions = None
        _stats = None
        inputParams[IMG_ZOOM] = 0
        inputParams[MODE] = MODE_TABLE
        dataset.setNumberSlices(0)
    else:
        #Setting the _imageVolName
        _imageVolName = inputParams.get(VOL_SELECTED, None) or table.getElementById(0, label)
     
        _typeOfColumnToRender = inputParams[COLS_CONFIG].getColumnProperty(label, 'columnType')
        
        #Setting the _imageDimensions
        _imageDimensions = readDimensions(request, _imageVolName, _typeOfColumnToRender)
        
        dataset.setNumberSlices(_imageDimensions[2])
        
        if _typeOfColumnToRender == "image":
            isVol = dataset.getNumberSlices() > 1
            is3D = inputParams[MODE]==MODE_VOL_ASTEX or inputParams[MODE]==MODE_VOL_CHIMERA
            #Setting the _convert 
            _convert = isVol and (inputParams[MODE]==MODE_GALLERY or is3D)
            #Setting the _reslice 
            _reslice = isVol and inputParams[MODE]==MODE_GALLERY
            #Setting the _getStats 
            _getStats = isVol and is3D
            #Setting the _dataType 
            _dataType = xmipp.DT_FLOAT if isVol and inputParams[MODE]==MODE_VOL_ASTEX else xmipp.DT_UCHAR
            #Setting the _imageVolName and _stats     
            _imageVolName, _stats = readImageVolume(request, _imageVolName, _convert, _dataType, _reslice, int(inputParams[VOL_VIEW]), _getStats)
            
            if isVol:
                inputParams[COLS_CONFIG].configColumn(label, renderFunc="get_slice")
                dataset.setVolumeName(_imageVolName)
            else:
                if inputParams[MODE] != MODE_TABLE:
                    inputParams[MODE] = MODE_GALLERY
                inputParams[COLS_CONFIG].configColumn(label, renderFunc="get_image")
                #dataset.setVolumeName(None)  
                #dataset.setNumberSlices(0)         
                
        volPath = os.path.join(request.session[PROJECT_PATH], _imageVolName)
    
    return volPath, _stats, _imageDimensions
    
#Initialize default values
DEFAULT_PARAMS = {
               PATH: None,
               ALLOW_RENDER: True,                 # Image can be displayed, depending on column layout 
               MODE: MODE_GALLERY,                   # Mode Options: gallery, table, column, volume_astex, volume_chimera
               TABLE_NAME: None,                    # Table name to display. If None the first one will be displayed
               LABEL_SELECTED: None,        # Column to be displayed in gallery mode. If None the first one will be displayed
               GOTO: 1,                           # Element selected (metadata record) by default. It can be a row in table mode or an image in gallery mode
               MANUAL_ADJUST: 'Off',                 # In gallery mode 'On' means columns can be adjust manually by the user. When 'Off' columns are adjusted automatically to screen width.
               COLS: None,                         # In gallery mode (and colRowMode set to 'On') cols define number of columns to be displayed
               ROWS: None,                          # In gallery mode (and colRowMode set to 'On') rows define number of columns to be displayed
               
               IMG_ZOOM: '128px',                     # Zoom set by default
               IMG_MIRRORY: False,                    # When 'True' image are mirrored in Y Axis 
               IMG_APPLY_TRANSFORM: False,       # When 'True' if there is transform matrix, it will be applied
               IMG_ONLY_SHIFTS: False,                 # When 'True' if there is transform matrix, only shifts will be applied
               IMG_WRAP: False,                       # When 'True' if there is transform matrix, only shifts will be applied
               IMG_MAX_WIDTH: 512,                # Maximum image width (in pixels)
               IMG_MIN_WIDTH: 64,                 # Minimum image width (in pixels)
               IMG_MAX_HEIGHT: 512,               # Maximum image height (in pixels)
               IMG_MIN_HEIGHT: 64,                # Minimum image height (in pixels)
               
               VOL_SELECTED: None,       # If 3D, Volume to be displayed in gallery, volume_astex and volume_chimera mode. If None the first one will be displayed
               VOL_VIEW: xmipp.VIEW_Z_NEG,      # If 3D, axis to slice volume 
               VOL_TYPE: 'map'}                 # If map, it will be displayed normally, else if pdb only astexViewer and chimera display will be available


def showj(request):
    """ This method will visualize data objects in table or gallery mode.
    Columns can be customized to be rendered, visible or more...
    This method can be called in two modes:
    GET: the first time the parameters will be parsed from GET request.
    POST: next call to the method will be done throgh POST, some paramters will
         be store also in SESSION
    """
    firstTime = request.method == 'GET'
    
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
    
    inputParams = DEFAULT_PARAMS.copy()
    extraParams = {}
    
    if firstTime:
        for key, value in request.GET.iteritems():
            if key in inputParams:
                inputParams[key] = value
            else:
                extraParams[key] = value   
        inputParams[PATH] = addProjectPrefix(request, inputParams[PATH])
        
        cleanSession(request, inputParams[PATH])
    else:
        for key, value in request.POST.iteritems():
            if key in inputParams:
                inputParams[key] = value
        # extraParams will be readed from SESSION
        
        
    request.session[IMG_ZOOM_DEFAULT] = inputParams[IMG_ZOOM]    
    
    #=DEBUG=====================================================================
#    from pprint import pprint
#    pprint(inputParams)
#    pprint(extraParams)
    #===========================================================================

    if inputParams[VOL_TYPE] != 'pdb':
        # Load the initial dataset from file or session
        dataset = loadDataSet(request, inputParams[PATH], firstTime)
        # Load the requested table (or the first if no specified)
        table = dataset.getTable(inputParams[TABLE_NAME])
        # Update inputParams to make sure have a valid table name (if using first table)
        inputParams[TABLE_NAME] = dataset.currentTable()
        
        # Load columns configuration. How to display columns and attributes (visible, render, editable)  
        loadColumnsConfig(request, dataset, table, inputParams, extraParams, firstTime)
    
        volPath, _stats, _imageDimensions = setRenderingOptions(request, dataset, table, inputParams)
    
        #Store variables into session 
        storeToSession(request, inputParams, dataset, _imageDimensions)
        
    else:
        inputParams[COLS_CONFIG] = None
        volPath = inputParams[PATH]

    context, return_page = createContextShowj(request, inputParams, dataset, table, _stats, volPath)

    render_var = render_to_response(return_page, RequestContext(request, context))
    
    #=TIME CONTROL==============================================================
#    print "FINISH SHOWJ: ", datetime.now()-start    
    #===========================================================================
    
    return render_var


def cleanSession(request, filename):
    """ Clean data stored in session for a new visualization. """
    #Store dataset and labelsToRender in session 
    if filename in request.session:
        del request.session[filename]
        
#     for key in [DATASET, LABEL_SELECTED, TABLE_NAME, IMG_DIMS]:
#         if key in request.session:
#             del request.session[key]
       
def storeToSession(request, inputParams, dataset, _imageDimensions):
    #Store dataset and labelsToRender in session 
    datasetDict = {}
    datasetDict[DATASET] = dataset
    datasetDict[LABEL_SELECTED] = inputParams[LABEL_SELECTED]
    datasetDict[TABLE_NAME] = inputParams[TABLE_NAME]
    datasetDict[IMG_DIMS] = _imageDimensions
    
    request.session[inputParams[PATH]] = datasetDict


def createContextShowj(request, inputParams, dataset, table, paramStats, volPath=None):
    showjForm = ShowjForm(dataset,
                          inputParams[COLS_CONFIG],
                          inputParams) # A form bound for the POST data and unbound for the GET
        
    if showjForm.is_valid() is False:
        print showjForm.errors

    context = createContext(dataset, table, inputParams[COLS_CONFIG], request, showjForm, inputParams)

    if inputParams[MODE]==MODE_VOL_ASTEX or inputParams[MODE]==MODE_VOL_CHIMERA:
        context.update(create_context_volume(request, inputParams, volPath, paramStats))
               
    elif inputParams[MODE]==MODE_GALLERY or inputParams[MODE]==MODE_TABLE or inputParams[MODE]=='column':
        context.update({"showj_alt_js": getResourceJs('showj_' + inputParams[MODE] + '_utils')})
        
    return_page = 'showj/%s%s%s' % ('showj_', showjForm.data[MODE], '.html')
#    return_page = 'showj/showj_base.html'
        
    return context, return_page
    

def createContext(dataset, table, columnsConfig, request, showjForm, inputParams):
    #Create context to be send
    
    context = {
            IMG_DIMS: request.session[inputParams[PATH]].get(IMG_DIMS, 0),
            IMG_ZOOM_DEFAULT: request.session.get(IMG_ZOOM_DEFAULT, 0),
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
        
        dataset = request.session[DATASET]
        blockComboBox = request.session[TABLE_NAME]
        
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
        context.update(create_context_astex(request, inputParams[VOL_TYPE], volPath))

#   'volType': 2, #0->byte, 1 ->Integer, 2-> Float
        
    elif inputParams[MODE] == MODE_VOL_CHIMERA:   
        context.update(create_context_chimera(volPath))
        
    return context
        

def create_context_astex(request, typeVolume, volPath):
#   volPath='/home/adrian/Scipion/tests/input/showj/emd_1042.map'        
#   fileName, fileExtension = os.path.splitext(_imageVolName)
            
#    linkName = 'test_link_' + request.session._session_key + '.' + inputParams[VOL_TYPE]

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
        
    return HttpResponse(mimetype='application/javascript')

