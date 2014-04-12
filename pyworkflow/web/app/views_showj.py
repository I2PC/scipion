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


def getTableLayoutConfig(request, dataset, tableDataset, inputParams, extraParams, firstTime):
    """ Load table layout configuration. How to display columns and attributes (visible, render, editable) """ 
    
    if firstTime:
        columns_properties = getExtraParameters(extraParams, tableDataset)
        request.session['defaultColumnsLayoutProperties'] = columns_properties

        layoutConfig = TableLayoutConfiguration(dataset, tableDataset, inputParams['allowRender'], columns_properties)
        request.session['tableLayoutConfiguration'] = layoutConfig
        
    
    else:
        layoutConfig = request.session['tableLayoutConfiguration'] 
            
            
    return layoutConfig


def setLabelToRender(request, dataset, tableDataset, inputParams, extraParams, firstTime):
    """ If no label is set to render, set the first one if exists """
    
    if 'labelsToRenderComboBox' not in inputParams or inputParams['labelsToRenderComboBox'] == '' or request.session[TABLE_NAME] != inputParams[TABLE_NAME]:
        labelsToRenderComboBoxValues = inputParams['tableLayoutConfiguration'].getLabelsToRenderComboBoxValues()
        if len(labelsToRenderComboBoxValues) > 0:
            inputParams['labelsToRenderComboBox'] = labelsToRenderComboBoxValues[0][0]
        else:
            # If there is no image to display and it is initial load, switch to table mode 
            if firstTime and inputParams['mode']!='table':
                inputParams['mode']='table'
                showj(request, inputParams, extraParams) 
            inputParams['labelsToRenderComboBox'] = ''
    
    tableDataset.setLabelToRender(inputParams['labelsToRenderComboBox']) 
    
    return dataset, tableDataset


def setRenderingOptions(request, dataset, tableDataset, inputParams):
    
    #Setting the _imageVolName
    _imageVolName = inputParams['volumesToRenderComboBox'] if ("volumesToRenderComboBox" in inputParams and inputParams['volumesToRenderComboBox'] != '') else tableDataset.getElementById(0,inputParams['labelsToRenderComboBox'])
 
    #Setting the _typeOfColumnToRender
    _typeOfColumnToRender = inputParams['tableLayoutConfiguration'].columnsLayout[inputParams['labelsToRenderComboBox']].typeOfColumn
    
    #Setting the _imageDimensions
    _imageDimensions = readDimensions(request, _imageVolName, _typeOfColumnToRender)
    
    dataset.setNumberSlices(_imageDimensions[2])
    
    if _typeOfColumnToRender == "image":
        
        #Setting the _convert 
        _convert = dataset.getNumberSlices()>1 and (inputParams['mode']=='gallery' or inputParams['mode']=='volume_astex' or inputParams['mode']=='volume_chimera')
        #Setting the _reslice 
        _reslice = dataset.getNumberSlices()>1 and inputParams['mode']=='gallery'
        #Setting the _getStats 
        _getStats = dataset.getNumberSlices()>1 and (inputParams['mode']=='volume_astex' or inputParams['mode']=='volume_chimera')
        #Setting the _dataType 
        _dataType = xmipp.DT_FLOAT if dataset.getNumberSlices()>1 and inputParams['mode']=='volume_astex' else xmipp.DT_UCHAR
        #Setting the _imageVolName and _stats     
        _imageVolName, _stats = readImageVolume(request, _imageVolName, _convert, _dataType, _reslice, int(inputParams['resliceComboBox']), _getStats)
        
        if dataset.getNumberSlices()>1:
            inputParams['tableLayoutConfiguration'].columnsLayout[inputParams['labelsToRenderComboBox']].columnLayoutProperties.renderFunc = "get_slice"
            dataset.setVolumeName(_imageVolName)
            
    volPath = os.path.join(request.session['projectPath'], _imageVolName)
    
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
        
        #Load table layout configuration. How to display columns and attributes (visible, render, editable)  
        inputParams['tableLayoutConfiguration'] = getTableLayoutConfig(request, dataset, table, inputParams, extraParams, firstTime)
    
#        for col in inputParams['tableLayoutConfiguration'].columnsLayout.values():
#            print "val!! ", col.columnLayoutProperties.getValues() 
    
        #If no label is set to render, set the first one if exists
        
        dataset, table = setLabelToRender(request, dataset, table, inputParams, extraParams, firstTime)
        
        if inputParams['labelsToRenderComboBox'] == '':
            inputParams['zoom'] = 0
            _imageDimensions = None
            dataset.setNumberSlices(0)
        else:
            dataset, volPath, _stats, _imageDimensions = setRenderingOptions(request, dataset, table, inputParams)
    
        #Store variables into session 
        request = storeToSession(request, inputParams, dataset, _imageDimensions)
        
    else:
        inputParams['tableLayoutConfiguration'] = None
        volPath = inputParams['path']

    context, return_page = createContextShowj(request, inputParams, dataset, table, _stats, volPath)

    render_var = render_to_response(return_page, RequestContext(request, context))
    
    #=TIME CONTROL==============================================================
#    print "FINISH SHOWJ: ", datetime.now()-start    
    #===========================================================================
    
    return render_var


def storeToSession(request, inputParams, dataset, _imageDimensions):
    #Store dataset and labelsToRender in session 
    request.session['dataset'] = dataset
    request.session['labelsToRenderComboBox'] = inputParams['labelsToRenderComboBox']
    request.session[TABLE_NAME] = inputParams[TABLE_NAME]
    request.session['imageDimensions'] = _imageDimensions
    
    return request


def createContextShowj(request, inputParams, dataset, table, paramStats, volPath=None):
    showjForm = ShowjForm(dataset,
                          inputParams['tableLayoutConfiguration'],
                          inputParams) # A form bound for the POST data and unbound for the GET
        
    if showjForm.is_valid() is False:
        print showjForm.errors

    context = createContext(dataset, table, inputParams['tableLayoutConfiguration'], request, showjForm)

    if inputParams['mode']=='volume_astex' or inputParams['mode']=='volume_chimera':
        context.update(create_context_volume(request, inputParams, volPath, paramStats))
               
    elif inputParams['mode']=='gallery' or inputParams['mode']=='table' or inputParams['mode']=='column':
        context.update({"showj_alt_js": getResourceJs('showj_' + inputParams['mode'] + '_utils')})
        
    return_page = 'showj/%s%s%s' % ('showj_', showjForm.data['mode'], '.html')
#    return_page = 'showj/showj_base.html'
        
    return context, return_page
    

def createContext(dataset, table, tableLayoutConfiguration, request, showjForm):
    #Create context to be send
    
    context = {
            'imageDimensions': request.session['imageDimensions'] if 'imageDimensions' in request.session else 0,
            'defaultZoom': request.session['defaultZoom'] if 'defaultZoom' in request.session else 0,
            'projectName': request.session['projectName'],
            'form': showjForm
            }
    
    
    if dataset is not None:
        context.update({'dataset': dataset})
    if tableLayoutConfiguration is not None:
        context.update({'tableLayoutConfiguration': json.dumps({'columnsLayout': tableLayoutConfiguration.columnsLayout,
                                                               #'colsOrder': tableLayoutConfiguration.colsOrder
                                                               },
                                                               ensure_ascii=False,
                                                               cls=ColumnLayoutConfigurationEncoder)})
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
    return XmippDataSet(str(path))

    
def save_showj_table(request):
    if request.is_ajax():
        changes = request.POST.get('changes')
        jsonChanges = json.loads(changes)
        
#        print "jsonChanges",jsonChanges 
        
        dataset=request.session['dataset']
        blockComboBox=request.session[TABLE_NAME]
#        tableLayoutConfiguration=request.session['tableLayoutConfiguration']
        
        table=dataset.getTable(blockComboBox)
        
        for key in jsonChanges:
            element_split = key.rsplit('___')
            if len(element_split)!=2: 
                print "this fails and sth has to be done"
            
            #NAPA de LUXE ahora mismo se realiza una conversion a int pero habria que ver de que tipo de datos se trata 
            #tableLayoutConfiguration.columnsLayout[element_split[0]].typeOfColumn

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
                       'path': '',
                       'allowRender': True,                 # Image can be displayed, depending on column layout 
                       'mode': 'gallery',                   # Mode Options: gallery, table, column, volume_astex, volume_chimera
                       'zoom': '128px',                     # Zoom set by default
                       TABLE_NAME: None,                    # Table name to display. If None the first one will be displayed
                       'labelsToRenderComboBox': '',        # Column to be displayed in gallery mode. If None the first one will be displayed
                       'volumesToRenderComboBox': '',       # If 3D, Volume to be displayed in gallery, volume_astex and volume_chimera mode. If None the first one will be displayed
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
#    request.session['projectPath'] = "/home/adrian/Scipion/projects/TestXmippWorkflow/"
#    request.session['projectName'] = "Juanitpo"
    
    if "projectPath" in request.session:
        projectPath = request.session['projectPath']
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
            inputParams['path'] = os.path.join(projectPath, fn)
#            writeSetOfMicrographs(obj, fn)
#            inputParams['path']= os.path.join(projectPath, fn)
        elif isinstance(obj, SetOfMovies):
            fn = project.getTmpPath(obj.getName() + '_movies.xmd')
            writeSetOfMovies(obj, fn)
            inputParams['path']= os.path.join(projectPath, fn)
        elif isinstance(obj, SetOfVolumes):
            fn = project.getTmpPath(obj.getName()+ '_volumes.xmd')
            inputParams['path'] = os.path.join(projectPath, createXmippInputVolumes(None, obj, volsFn=os.path.join(projectPath, fn))) 
            inputParams['mode']= 'table'
            
#            writeSetOfVolumes(obj, fn)
#            inputParams['path']= os.path.join(projectPath, fn)

        elif isinstance(obj, PdbFile):
            inputParams['path'] = obj.getFileName()
            inputParams['mode'] = 'volume_astex'
            inputParams['typeVolume'] = 'pdb'

        elif isinstance(obj, SetOfImages):
            fn = project.getTmpPath(obj.getName() + '_images.xmd')
            inputParams['path'] = os.path.join(projectPath, createXmippInputImages(None, obj, imagesFn=os.path.join(projectPath, fn)))
#            writeSetOfParticles(obj, fn)
#            inputParams['path']= os.path.join(projectPath, fn)
        elif isinstance(obj, Image):  
            fn = project.getTmpPath(obj.getName() + '_image.xmd')
            
            md = xmipp.MetaData()
            imgFn = getImageLocation(obj)
            if isinstance(obj, Volume) and imgFn.endswith('.mrc'):
                imgFn += ':mrc'
            md.setValue(xmipp.MDL_IMAGE, imgFn, md.addObject())
            md.write(fn)
            
            inputParams['path'] = os.path.join(projectPath, fn)
            
        elif isinstance(obj, SetOfClasses2D):
            fn = project.getTmpPath(obj.getName() + '_classes.xmd')
#            writeSetOfClasses2D(obj, fn)
#            inputParams['path']= os.path.join(projectPath, fn)
            inputParams['path'] = os.path.join(projectPath, createXmippInputClasses2D(None, obj, classFn=os.path.join(projectPath, fn)))
        
        elif isinstance(obj, SetOfCTF):
            fn = project.getTmpPath(obj.getName() + '_ctfs.xmd')
            inputParams['path'] = os.path.join(projectPath, createXmippInputCTF(None, obj, ctfFn=os.path.join(projectPath, fn)))
            
            extraParams["itemId___visible"]= False
#            extraParams["psd___visible"]= False
            extraParams["psdEnhanced___renderable"]= True
            extraParams["micrograph___renderable"]= True
            extraParams["image1___renderable"]= True
            extraParams["image2___renderable"]= True
            
#            writeSetOfCTFs(obj, fn)
#            inputParams['path']= os.path.join(projectPath, fn)
        elif isinstance(obj, NormalModes):
            inputParams['path'] = os.path.join(projectPath, obj.getFileName())
        else:
            raise Exception('Showj Web visualizer: can not visualize class: %s' % obj.getClassName())
    
    elif "path" in request.GET:
        inputParams.update({'path':os.path.join(projectPath, request.GET.get("path"))})
        
    else:
        raise Exception('Showj Web visualizer: No object identifier or path found')         

    request.session['defaultZoom'] = inputParams['zoom']
    
    # Clean the dataset if exist 
    if 'dataset' in request.session:
        request.session.__delitem__('dataset')
    
    return showj(request, inputParams, extraParams, firstTime=True)  


def testingSSH(request):
    context = {}
#    return render_to_response("testing_ssh.html", RequestContext(request, context))
    return render_to_response("scipion.html", RequestContext(request, context))



def create_context_volume(request, inputParams, volPath, param_stats):
#        volPath = os.path.join(request.session['projectPath'], _imageVolName)

    threshold = calculateThreshold(param_stats)
    
    context = {"threshold":threshold,
                    'minStats':param_stats[2] if param_stats != None else 1,
                    'maxStats':param_stats[3] if param_stats != None else 1 }
    
    if inputParams['mode'] == 'volume_astex':
        context.update(create_context_astex(request, inputParams['typeVolume'], volPath))

#   'volType': 2, #0->byte, 1 ->Integer, 2-> Float
        
    elif inputParams['mode'] == 'volume_chimera':   
        context.update(create_context_chimera(volPath))
        
    return context
        

def create_context_astex(request, typeVolume, volPath):
#   volPath='/home/adrian/Scipion/tests/input/showj/emd_1042.map'        
#   fileName, fileExtension = os.path.splitext(_imageVolName)
            
#    linkName = 'test_link_' + request.session._session_key + '.' + inputParams["typeVolume"]

    linkName = 'test_link_' + request.session._session_key + '.' + typeVolume
    volLinkPath = os.path.join(pw.WEB_RESOURCES, 'astex', 'tmp', linkName)
    
    from pyworkflow.utils.path import cleanPath, createLink
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
#            print "stdout_data",stdout_data
#            print "stderr_data",stderr_data
    f = open(outputHtmlFile)
    chimeraHtml = f.read().decode('string-escape').decode("utf-8").split("</html>")[1]
    
    return {"chimeraHtml":chimeraHtml}
    
    
def calculateThreshold(params):
    threshold = (params[2] + params[3])/2 if params != None else 1
        
#    print "Threshold:", threshold
#    print "stats:",params
#    print "minStats:", params[2]
#    print "maxStats:", params[3]
    
    return threshold


def updateSessionTable(request):
    label = request.GET.get('label', None)
    type = request.GET.get('type', None)
    option = request.GET.get('option', None)
    
    if type == "renderable":
        request.session["tableLayoutConfiguration"].columnsLayout[label].renderable = option
    elif type == "editable":    
        request.session["tableLayoutConfiguration"].columnsLayout[label].editable = option
        
#    print request.session["tableLayoutConfiguration"].columnsLayout[label].renderable
     
    return HttpResponse(mimetype='application/javascript')

