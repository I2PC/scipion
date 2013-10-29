import os
import xmipp
from django.http import HttpResponse
from pyworkflow.web.pages import settings
from django.shortcuts import render_to_response
from pyworkflow.tests import getInputPath
from pyworkflow.web.app.forms import ShowjForm, getLabelsToRenderComboBoxValues
from django.template import RequestContext
from os.path import join    
from collections import OrderedDict
import json
from pyworkflow.web.app.views_util import *
from forms import VolVisualizationForm
from pyworkflow.em import *
from pyworkflow.em.packages.xmipp3.convert import *

    
def showj(request, inputParameters=None):
    #############
    # WEB INPUT PARAMETERS
    _imageDimensions = ''
    _imageVolName=''
    _stats=None
    
    
    if request.method == 'POST': # If the form has been submitted... Post method
        inputParameters=request.POST.copy()
    else:
        request.session['defaultZoom'] = inputParameters['zoom'] if 'zoom' in inputParameters else '150px'
    
    #Init Dataset
    #NAPA DE LUXE: Check type of Dataset 
    dataset = loadDatasetXmipp(inputParameters['path']) 
    if inputParameters['blockComboBox'] == '':
        inputParameters['blockComboBox'] = dataset.listTables()[0]
    
    #Get table from block name
    dataset.setTableName(inputParameters['blockComboBox'])
    tableDataset=dataset.getTable()
    
    #Load table layout configuration. How to display columns and attributes (visible, render, editable)  
    inputParameters['tableLayoutConfiguration']= TableLayoutConfiguration(dataset, tableDataset, inputParameters['allowRender'])

    #If no label is set to render, set the first one if exists
    if inputParameters['labelsToRenderComboBox'] == '':
        labelsToRenderComboBoxValues = getLabelsToRenderComboBoxValues(inputParameters['tableLayoutConfiguration'].columnsLayout)
        inputParameters['labelsToRenderComboBox']=labelsToRenderComboBoxValues[0][0] if len(labelsToRenderComboBoxValues) > 0 else ''

    
    if inputParameters['labelsToRenderComboBox'] == '':
        inputParameters['zoom']=0
        _imageDimensions = None
    else:
        _imageVolName = inputParameters['volumesToRenderComboBox'] if ("volumesToRenderComboBox" in inputParameters and inputParameters['volumesToRenderComboBox'] != '') else tableDataset.getElementById(0,inputParameters['labelsToRenderComboBox'])
        
        _getStats=False
        _convert=True
        _dataType=xmipp.DT_UCHAR
        _reslice=False 
        
        if inputParameters['dims']=='3d':
            if inputParameters['mode']=='column' or inputParameters['mode']=='table':
                _convert=False
            elif inputParameters['mode']=='gallery':
                _reslice=True 
            elif inputParameters['mode']=='volume_astex' or inputParameters['mode']=='volume_chimera':
                _getStats=True
                if inputParameters['mode']=='volume_astex':
                    _dataType = xmipp.DT_FLOAT
            
        _imageVolName, _imageDimensions, _stats = readVolume(request, _imageVolName, True, _convert, _dataType, _reslice, inputParameters['resliceComboBox'], _getStats)

    dataset.setLabelToRender(inputParameters['labelsToRenderComboBox'])    
    if _imageDimensions != None and _imageDimensions[2]>1:
        dataset.setNumberSlices(_imageDimensions[2])
        dataset.setVolumeName(_imageVolName)

    showjForm = ShowjForm(dataset,
                          inputParameters['tableLayoutConfiguration'],
                          request.POST if request.method == 'POST' else inputParameters) # A form bound for the POST data and unbound for the GET
        
    if showjForm.is_valid() is False:
        print showjForm.errors

    threshold = (_stats[2] + _stats[3])/2 if _stats != None else 1
    volLink = ""
    chimeraHtml = ""
    volPath = os.path.join(request.session['projectPath'], _imageVolName)
    
    if showjForm.data['mode']=='volume_astex':
        
        fileName, fileExtension = os.path.splitext(_imageVolName)
        
#        volPath='/home/adrian/Scipion/tests/input/showj/emd_1042.map'
        
        # Astex viewer            
        #Hay qye ver como gestionamos el tema de la extension (con .map no me lei un map.gz)
        linkName = 'test_link_' + request.session._session_key + '.map'
        volLinkPath = os.path.join(pw.WEB_RESOURCES, 'astex', 'tmp', linkName)
        from pyworkflow.utils.path import cleanPath, createLink
        cleanPath(volLinkPath)
        print "volPath",volPath
        print "volLinkPath",volLinkPath
        createLink(volPath, volLinkPath)
        volLink = os.path.join('/', settings.STATIC_ROOT, 'astex', 'tmp', linkName)
     
    elif showjForm.data['mode']=='volume_chimera':   
        from subprocess import Popen, PIPE, STDOUT
#        p = Popen(['/home/adrian/.local/UCSF-Chimera64-2013-10-16/bin/chimera', volPath], stdout=PIPE, stdin=PIPE, stderr=PIPE)
        p = Popen([os.environ.get('CHIMERA_HEADLESS'), volPath], stdout=PIPE, stdin=PIPE, stderr=PIPE)
        outputHtmlFile = os.path.join(pw.WEB_RESOURCES, 'astex', 'tmp', 'test.html')
        #chimeraCommand= 'volume #0 level ' + str(threshold) + '; export format WebGL ' + outputHtmlFile + '; stop'
        chimeraCommand= 'export format WebGL ' + outputHtmlFile + '; stop'
        stdout_data, stderr_data = p.communicate(input=chimeraCommand)
        print "stdout_data",stdout_data
        print "stderr_data",stderr_data
        f = open(outputHtmlFile)
        chimeraHtml = f.read().decode('string-escape').decode("utf-8").split("</html>")[1]
        
        
    #Store dataset and labelsToRender in session 
    request.session['dataset'] = dataset
    request.session['labelsToRenderComboBox'] = inputParameters['labelsToRenderComboBox']
    request.session['blockComboBox'] = inputParameters['blockComboBox']
#    request.session['tableLayoutConfiguration'] = tableLayoutConfiguration
    if (_imageDimensions != ''):
        request.session['imageDimensions'] = _imageDimensions

    #Create context to be send
    context = {'css': getResourceCss('showj'),
               'smoothness': getResourceCss('ui_smoothness'),
               'demo_table_jui': getResourceCss('showj_demo_table_jui'),
               
               'favicon': getResourceIcon('favicon'),
               'logo': getResourceIcon('logo_scipion'),
               
               'jquery': getResourceJs('jquery'), #Configuration variables
               'jquery_datatable': getResourceJs('jquery_datatables'),
               'jquerydataTables_colreorder': getResourceJs('jquery_colreorder'),
               'jquerydataTables_colreorder_resize': getResourceJs('jquery_colreorder_resize'),
               'jeditable': getResourceJs('jquery_editable'),
               'jquery_waypoints':getResourceJs('jquery_waypoints'),
               'jquery_hover_intent':getResourceJs('jquery_hover_intent'),
               
               'dataset': dataset,
               'tableLayoutConfiguration': json.dumps({'columnsLayout': inputParameters['tableLayoutConfiguration'].columnsLayout, 'colsOrder': inputParameters['tableLayoutConfiguration'].colsOrder}, ensure_ascii=False, cls=ColumnLayoutConfigurationEncoder), #Data variables
               'tableDataset': tableDataset,
               
               'imageDimensions': request.session['imageDimensions'],
               'defaultZoom': request.session['defaultZoom'],
               'projectName': request.session['projectName'],
               
               'form': showjForm,
               
               
               'volLink': volLink, #Astex Viewer
#               'volType': 2, #0->byte, 1 ->Integer, 2-> Float
               'minStats':_stats[2] if _stats != None else 1,
               'maxStats':_stats[3] if _stats != None else 1, 
               'threshold': threshold,
               'chimeraHtml':chimeraHtml} #Form
    
    return_page = '%s%s%s' % ('showj_', showjForm.data['mode'], '.html')
    return render_to_response(return_page, RequestContext(request, context))

###################################################################
################################ BEGIN LAYOUT #####################  
###################################################################
class ColumnLayoutConfigurationEncoder(json.JSONEncoder):
    def default(self, columnLayoutConfiguration):
        columnLayoutConfigurationCoded={}
        columnLayoutConfigurationCoded={"typeOfColumn":columnLayoutConfiguration.typeOfColumn,
                                        "columnLayoutProperties":{"visible":columnLayoutConfiguration.columnLayoutProperties.visible,
                                                                  "allowSetVisible":columnLayoutConfiguration.columnLayoutProperties.allowSetVisible,
                                                                  "editable":columnLayoutConfiguration.columnLayoutProperties.editable,
                                                                  "allowSetEditable":columnLayoutConfiguration.columnLayoutProperties.allowSetEditable,
                                                                  "renderable":columnLayoutConfiguration.columnLayoutProperties.renderable,
                                                                  "allowSetRenderable":columnLayoutConfiguration.columnLayoutProperties.allowSetRenderable,
                                                                  "renderFunc":columnLayoutConfiguration.columnLayoutProperties.renderFunc}
                                        }
        return columnLayoutConfigurationCoded                     
                 
            
class TableLayoutConfiguration():
    def __init__(self, ds, tableDataset, allowRender=True):
        
        self.columnsLayout = OrderedDict() 
         
        for col in tableDataset.iterColumns():
            self.columnsLayout[col.getName()]=ColumnLayoutConfiguration(col, ds.getTypeOfColumn(col.getName()), allowRender)
            
        self.colsOrder = defineColsLayout(self.columnsLayout.keys())
        
            
class ColumnLayoutConfiguration():
    def __init__(self, col, typeOfColumn, allowRender):
        self.columns = col
        
        self.label = col.getName()
        self.typeOfColumn = typeOfColumn
        
        self.columnLayoutProperties = ColumnLayoutProperties(self.typeOfColumn, allowRender)

class ColumnLayoutProperties():
    def __init__(self, typeOfColumn, allowRender=True):
#        self.layoutPropertiesDict = {}
        self.visible = not(typeOfColumn == 'id')
        self.allowSetVisible = True 
        
        self.editable = (typeOfColumn == 'text')
        self.allowSetEditable = self.editable
        
        self.renderable = False
        self.allowSetRenderable = True if (typeOfColumn == 'image' and allowRender) else False

        self.renderFunc = "taka"
        
def defineColsLayout(labels):
    colsOrder = range(len(labels))
    if 'enabled' in labels:
        colsOrder.insert(0, colsOrder.pop(labels.index('enabled')))
    return colsOrder

###################################################################
################################ END LAYOUT #######################    
###################################################################

###################################################################
########################### BEGIN SAVE & LOAD #####################    
###################################################################
#### Load an Xmipp Dataset ###
def loadDatasetXmipp(path):
    """ Create a table from a metadata. """
    from pyworkflow.em.packages.xmipp3 import XmippDataSet
    mdPath = getInputPath('showj', path)
    return XmippDataSet(str(mdPath))

    
def save_showj_table(request):
    if request.is_ajax():
        changes = request.POST.get('changes')
        jsonChanges = json.loads(changes)
        
        print "jsonChanges",jsonChanges 
        
        dataset=request.session['dataset']
        blockComboBox=request.session['blockComboBox']
#        tableLayoutConfiguration=request.session['tableLayoutConfiguration']
        
        tableDataset=dataset.getTable(blockComboBox)
        
        for key in jsonChanges:
            element_split = key.rsplit('___')
            if len(element_split)!=2: 
                print "this fails and sth has to be done"
            
            #NAPA de LUXE ahora mismo se realiza una conversion a int pero habria que ver de que tipo de datos se trata 
            #tableLayoutConfiguration.columnsLayout[element_split[0]].typeOfColumn

            dictelement = {element_split[0]:jsonChanges[key]}
            tableDataset.updateRow(int(element_split[1]),**dictelement)
        
        dataset.writeTable(blockComboBox, tableDataset)
        
        return HttpResponse(json.dumps({'message':'Ok'}), mimetype='application/javascript')
###################################################################
########################### END SAVE & LOAD #######################    
###################################################################    


###################################################################
######################## BEGIN VOLUME DISPLAY #####################    
###################################################################    

def visualizeVolume(request):
    from django.http import HttpResponse
    import json
     
    if request.is_ajax():
        setOfVolumesId = int(request.GET.get('setOfVolumesId'))
        volumeId = int(request.GET.get('volumeId'))
        projectName = request.session['projectName']
        project = loadProject(projectName)
        setOfVolume = project.mapper.selectById(setOfVolumesId)
        volume = setOfVolume[volumeId]
        # Chimera 
        from subprocess import Popen, PIPE, STDOUT
#         inputVolume = join(os.getcwd(),volume.getFileName())
#         outputHtmlFile = join(os.getcwd(),project.getTmpPath("volume_" + str(volume.getObjId()) + '.html'))
        inputVolume = volume.getFileName()
        outputHtmlFile = project.getTmpPath("volume_" + str(volume.getObjId()) + '.html')
        if (request.GET.get('threshold') is None or request.GET.get('threshold') == ''):
            threshold = 1
        else:
            threshold = float(request.GET.get('threshold'))
        print ("THRESHOLD",threshold)
        # TODO: Get with Xmipp an approximate threshold
        p = Popen(['chimera', inputVolume], stdout=PIPE, stdin=PIPE, stderr=PIPE)
#         p = Popen(['chimera', '--start', 'ReadStdin', inputVolume], stdout=PIPE, stdin=PIPE, stderr=PIPE)        
        stdout_data = p.communicate(input='volume #0 level ' + str(threshold) + '; export format WebGL ' + outputHtmlFile + '; stop')[0]
        f = open(outputHtmlFile)
        volumeHtml = f.read().decode('string-escape').decode("utf-8").split("</html>")[1]
        jsonStr = json.dumps({'volumeHtml': volumeHtml})
        return HttpResponse(jsonStr, mimetype='application/javascript')
    

def showVolVisualization(request):
    form = None
    volLinkPath = None
    volLink = None
    chimeraHtml = None
    if (request.POST.get('operation') == 'visualize'):
        form = VolVisualizationForm(request.POST, request.FILES)
        if form.is_valid():
            volPath = form.cleaned_data['volPath']
            
            fileName, fileExtension = os.path.splitext(volPath)
            
            print "fileExtension",fileExtension
            if (fileExtension != '.map' and fileExtension != '.map.gz') :
                
                img = xmipp.Image()
    #            imgFn = os.path.join(request.session['projectPath'], _imageVolName)
                img.read(str(volPath))
                #img.convert2DataType(xmipp.DT_UCHAR, xmipp.CW_CAST)
                img.convert2DataType(xmipp.DT_FLOAT, xmipp.CW_CAST)
                
                _imageVolName2 = '%s_tmp%s' % (fileName, '.map')
                print "_imageVolName2",_imageVolName2 
                print "_imageVolName22", os.path.basename(_imageVolName2)
                volPathTmp=os.path.join(pw.HOME, 'web', 'pages', 'resources', 'astex', 'tmp', os.path.basename(_imageVolName2))
                print "volPathTmp",volPathTmp
                img.write(str(volPathTmp))
                volPath=volPathTmp
                
            
            # Astex viewer            
            from random import randint
            #Hay qye ver como gestionamos el tema de la extension (con .map no me lei un map.gz)
            linkName = 'test_link_' + str(randint(0, 10000)) +'.map'
            volLinkPath = os.path.join(pw.HOME, 'web', 'pages', 'resources', 'astex', 'tmp', linkName)
            from pyworkflow.utils.path import cleanPath, createLink
            cleanPath(volLinkPath)
            print "volPath",volPath
            print "volLinkPath",volLinkPath
            createLink(volPath, volLinkPath)
#             os.system("ln -s " + str(volPath) + " " + volLinkPath)
            volLink = os.path.join('/', 'static', 'astex', 'tmp', linkName)
           
           
            # Chimera 
            from subprocess import Popen, PIPE, STDOUT
##             p = Popen(['chimera', '--start', 'ReadStdin', volPath], stdout=PIPE, stdin=PIPE, stderr=PIPE)
#            p = Popen(['/home/adrian/.local/UCSF-Chimera64-2013-10-16/bin/chimera', volPath], stdout=PIPE, stdin=PIPE, stderr=PIPE)
            p = Popen(['$CHIMERA_HEADLESS', volPath], stdout=PIPE, stdin=PIPE, stderr=PIPE)

            outputHtmlFile = '/home/adrian/test.html'
            threshold = form.cleaned_data['threshold']
            stdout_data, stderr_data = p.communicate(input='volume #0 level ' + str(threshold) + '; export format WebGL ' + outputHtmlFile + '; stop')
            print "stdout_data",stdout_data
            print "stderr_data",stderr_data 
            f = open(outputHtmlFile)
            chimeraHtml = f.read().decode('string-escape').decode("utf-8").split("</html>")[1]
#            chimeraHtml = ""
    else:
        form = VolVisualizationForm()
        
    print "brokenpipe"    
    context = {'MEDIA_URL' : settings.MEDIA_URL, 'STATIC_URL' :settings.STATIC_URL, 'form': form, 'volLink': volLink, 'chimeraHtml': chimeraHtml}
#    print "context",context    
    return render_to_response('showVolVisualization.html',  RequestContext(request, context))   

###################################################################
######################## END VOLUME DISPLAY #######################    
###################################################################  


###################################################################
######################## DISPATCHER ###############################    
###################################################################  
def visualizeObject(request):
    #Initialize default values
    inputParameters = {'allowRender': True,                 # Image can be displayed, depending on column layout 
                       'mode': 'gallery',                   # Mode Options: gallery, table, column, volume_astex, volume_chimera
                       'zoom': '150px',                     # Zoom set by default
                       'blockComboBox': '',                 # Metadata Block to display. If None the first one will be displayed
                       'labelsToRenderComboBox': '',        # Column to be displayed in gallery mode. If None the first one will be displayed
                       'volumesToRenderComboBox': '',       # If 3D, Volume to be displayed in gallery, volume_astex and volume_chimera mode. If None the first one will be displayed
                       'dims': '2d',                        # Object Dimensions
                       'goto': 1,                           # Element selected (metadata record) by default. It can be a row in table mode or an image in gallery mode
                       'colRowMode': 'Off',                 # In gallery mode 'On' means columns can be adjust manually by the user. When 'Off' columns are adjusted automatically to screen width. 
                       'mirrorY': False,                    # When 'True' image are mirrored in Y Axis 
                       'applyTransformMatrix': False,       # When 'True' if there is transform matrix, it will be applied
                       'onlyShifts': False,                 # When 'True' if there is transform matrix, only shifts will be applied
                       'wrap': False,                       # When 'True' if there is transform matrix, only shifts will be applied
                       'resliceComboBox': xmipp.VIEW_Z_NEG, # If 3D, axis to slice volume 
                       'imageMaxWidth': 512,                # Maximum image width (in pixels)
                       'imageMinWidth': 30,                 # Minimum image width (in pixels)
                       'imageMaxHeight': 512,               # Maximum image height (in pixels)
                       'imageMinHeight': 30}                # Minimum image height (in pixels)

# NAPA DE LUXE: PARA probarlo sin sesion iniciada    
#    request.session['projectPath'] = "/home/adrian/Scipion/projects/TestXmippWorkflow/"
#    request.session['projectName'] = "Juanitpo"
    
    if "projectPath" in request.session:
        projectPath = request.session['projectPath']
    else:         
        raise Exception('Showj Web visualizer: No project loaded')
    
    if "objectId" in request.GET: 
        objectId = request.GET.get("objectId")

        project = Project(projectPath)
        project.load()
        
        obj = project.mapper.selectById(int(objectId))
        if obj.isPointer():
            obj = obj.get()

        if isinstance(obj, SetOfMicrographs):
            fn = project.getTmpPath(obj.getName() + '_micrographs.xmd')
            writeSetOfMicrographs(obj, fn)
            inputParameters['path']= join(projectPath, fn)
        elif isinstance(obj, SetOfVolumes):
            fn = project.getTmpPath(obj.getName()+ '_volumes.xmd')
            writeSetOfVolumes(obj, fn)
            inputParameters['path']= join(projectPath, fn)
#            inputParameters['setOfVolumes']= obj
#            inputParameters['setOfVolumesId']= obj.getObjId()
            inputParameters['dims']= '3d'
            inputParameters['mode']= 'table'
        elif isinstance(obj, SetOfImages):
    #        PAJM aqui falla para el cl2d align y se esta perdiendo la matrix de transformacion en la conversion
            fn = project.getTmpPath(obj.getName() + '_images.xmd')
            writeSetOfParticles(obj, fn)
            inputParameters['path']= join(projectPath, fn)
        elif isinstance(obj, Image):
            fn = project.getTmpPath(obj.getName() + '_image.xmd')
            fnSet = fn.replace('.xmd', '.sqlite')
            cleanPath(fn, fnSet)
            
            imgSet = SetOfImages()
            imgSet.setFileName(fnSet)
            img = Image()
            #img.copyInfo(obj)
            img.copyLocation(obj)
            imgSet.append(img)
            writeSetOfParticles(imgSet, fn)
            inputParameters['path']= join(projectPath, fn)
            
        elif isinstance(obj, SetOfClasses2D):
            fn = project.getTmpPath(obj.getName() + '_classes.xmd')
            writeSetOfClasses2D(obj, fn)
            inputParameters['path']= join(projectPath, fn)
        else:
            raise Exception('Showj Web visualizer: can not visualize class: %s' % obj.getClassName())
    
    elif "path" in request.GET:
        inputParameters.update(request.GET.items())
        inputParameters.update({'path':join(projectPath, request.GET.get("path"))})
    else:
        raise Exception('Showj Web visualizer: No object identifier or path found')         

    return showj(request, inputParameters)  


def testingSSH(request):
    context = {}
    return render_to_response("testing_ssh.html", RequestContext(request, context))
