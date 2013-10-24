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
    
    if request.method == 'POST': # If the form has been submitted... Post method
        _path = request.POST.get('path')
        _blockComboBox = request.POST.get('blockComboBox')
        _render = request.POST.get('allowRender')
        _labelsToRenderComboBox = request.POST.get('labelsToRenderComboBox')
         
    else: #If the form is called by get 
        _path = inputParameters['path']
        _blockComboBox = inputParameters['blockComboBox'] if 'blockComboBox' in inputParameters else '' 
        _render = inputParameters['allowRender']
        _labelsToRenderComboBox = inputParameters['labelsToRenderComboBox'] if 'labelsToRenderComboBox' in inputParameters else ''
        request.session['defaultZoom'] = inputParameters['zoom'] if 'zoom' in inputParameters else '150px' 
    
    print("path",_path)
    
    
    #Init Dataset
    #NAPA DE LUXE: Check type of Dataset 
    dataset = loadDatasetXmipp(_path) 
    
    if _blockComboBox == '':
        _blockComboBox = dataset.listTables()[0]
    
    #Get table from block name
    dataset.setTableName(_blockComboBox)
    tableDataset=dataset.getTable()
    
    #Load table layout configuration. How to display columns and attributes (visible, render, editable)  
    #tableLayoutConfiguration = TableLayoutConfiguration(tableDataset, _render) if 'blockComboBox' not in request.session or _blockComboBox != request.session['blockComboBox'] or request.method == 'GET' else None
    tableLayoutConfiguration = TableLayoutConfiguration(dataset, tableDataset, _render)

    #Initialize Showj Form (button toolbar)
    if request.method == 'POST': # If the form has been submitted... Post method
        if _labelsToRenderComboBox != '':
#            if _labelsToRenderComboBox != request.session['labelsToRenderComboBox']:
#            print request.POST.get('volumesToRenderComboBox') if ('volumesToRenderComboBox' in request.POST and request.POST.get('volumesToRenderComboBox') != '') else "taka"
            _imageVolName = request.POST.get('volumesToRenderComboBox') if ('volumesToRenderComboBox' in request.POST and request.POST.get('volumesToRenderComboBox') != '') else tableDataset.getElementById(0,_labelsToRenderComboBox)
            if request.POST.get('dims')=='3d':
                if request.POST.get('mode')=='volume_astex':
                    _imageVolName = readVolumeAndReslice(request.session['projectPath'], _imageVolName, int(request.POST.get('resliceComboBox')), xmipp.DT_FLOAT)
                elif request.POST.get('mode')=='column' or request.POST.get('mode')=='gallery':        
                    _imageVolName = readVolumeAndReslice(request.session['projectPath'], _imageVolName, int(request.POST.get('resliceComboBox')), xmipp.DT_UCHAR)
                
                
                
            _imageDimensions = getImageDim(request, _imageVolName)
        else:
            _imageDimensions = None

    else:
        if _labelsToRenderComboBox == '':
            labelsToRenderComboBoxValues = getLabelsToRenderComboBoxValues(tableLayoutConfiguration.columnsLayout)
            _labelsToRenderComboBox=labelsToRenderComboBoxValues[0][0] if len(labelsToRenderComboBoxValues) > 0 else ''
            inputParameters['labelsToRenderComboBox']=_labelsToRenderComboBox
            
        if inputParameters['labelsToRenderComboBox'] == '':
            inputParameters['zoom']=0
            _imageDimensions = None
        else:    
            _imageVolName = tableDataset.getElementById(0,inputParameters['labelsToRenderComboBox'])
            if inputParameters['dims']=='3d':
                if inputParameters['mode']=='volume_astex':
                    _imageVolName = readVolumeAndReslice(request.session['projectPath'], _imageVolName, inputParameters['resliceComboBox'], xmipp.DT_FLOAT)
                elif inputParameters['mode']=='column' or inputParameters['mode']=='gallery':    
                    _imageVolName = readVolumeAndReslice(request.session['projectPath'], _imageVolName, inputParameters['resliceComboBox'], xmipp.DT_UCHAR)
                
            _imageDimensions = getImageDim(request, _imageVolName)
            
        inputParameters['blockComboBox']=_blockComboBox
        inputParameters['tableLayoutConfiguration']=tableLayoutConfiguration

    if _imageDimensions != None and _imageDimensions[2]>1:
        dataset.setNumberSlices(_imageDimensions[2])
        dataset.setVolumeName(_imageVolName)

    
                    
    showjForm = ShowjForm(dataset,
                          tableLayoutConfiguration,
                          request.POST if request.method == 'POST' else inputParameters) # A form bound for the POST data and unbound for the GET
        
    if showjForm.is_valid() is False:
        print showjForm.errors

    dataset.setLabelToRender(_labelsToRenderComboBox)
    
    volLink=''
#    threshold =0.285 #Para emd_1042.map
    #threshold =11580 #Para hand.vol
    threshold = 1
    chimeraHtml=''
    
    
    volPath = os.path.join(request.session['projectPath'], _imageVolName)
    print "path",volPath
    if showjForm.data['mode']=='volume_astex':
        
        fileName, fileExtension = os.path.splitext(_imageVolName)
        
#        volPath='/home/adrian/Scipion/tests/input/showj/emd_1042.map'
        
        # Astex viewer            
        from random import randint
        #Hay qye ver como gestionamos el tema de la extension (con .map no me lei un map.gz)
        linkName = 'test_link_' + str(randint(0, 10000)) + '.map'
        volLinkPath = os.path.join(pw.HOME, 'web', 'pages', 'resources', 'astex', 'tmp', linkName)
        from pyworkflow.utils.path import cleanPath, createLink
        cleanPath(volLinkPath)
        createLink(volPath, volLinkPath)
        volLink = os.path.join('/', 'static', 'astex', 'tmp', linkName)
     
    elif showjForm.data['mode']=='volume_chimera':   
        from subprocess import Popen, PIPE, STDOUT
#        p = Popen(['/home/adrian/.local/UCSF-Chimera64-2013-10-16/bin/chimera', volPath], stdout=PIPE, stdin=PIPE, stderr=PIPE)
        p = Popen([os.environ.get('CHIMERA_HEADLESS'), volPath], stdout=PIPE, stdin=PIPE, stderr=PIPE)
        outputHtmlFile = os.path.join(pw.HOME, 'web', 'pages', 'resources', 'astex', 'tmp', 'test.html')
        #chimeraCommand= 'volume #0 level ' + str(threshold) + '; export format WebGL ' + outputHtmlFile + '; stop'
        chimeraCommand= 'export format WebGL ' + outputHtmlFile + '; stop'
        stdout_data, stderr_data = p.communicate(input=chimeraCommand)
        print "stdout_data",stdout_data
        print "stderr_data",stderr_data
        f = open(outputHtmlFile)
        chimeraHtml = f.read().decode('string-escape').decode("utf-8").split("</html>")[1]
        
        
    #Store dataset and labelsToRender in session 
    request.session['dataset'] = dataset
    request.session['labelsToRenderComboBox'] = _labelsToRenderComboBox
    request.session['blockComboBox'] = _blockComboBox
#    request.session['tableLayoutConfiguration'] = tableLayoutConfiguration
    #Esto falla y creo que es por el tamano
    #request.session['table'] = tableDataset
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
               'tableLayoutConfiguration': json.dumps({'columnsLayout': tableLayoutConfiguration.columnsLayout, 'colsOrder': tableLayoutConfiguration.colsOrder}, ensure_ascii=False, cls=ColumnLayoutConfigurationEncoder), #Data variables
               'tableDataset': tableDataset,
               'imageDimensions': request.session['imageDimensions'],
               'defaultZoom': request.session['defaultZoom'],
               'projectName': request.session['projectName'],
               'form': showjForm,
               
#               NAPA DE LUXE: ESTO es un poco cutre
               'STATIC_URL' :settings.STATIC_URL,
               'volLink': volLink,
               'volType': 2, #0->byte, 1 ->Integer, 2-> Float
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
    #Napa de Luxe
    probandoCTFParam = False
    probandoVolume = False
    
    objectId = request.GET.get("objectId")    
    projectPath = request.session['projectPath']
    project = Project(projectPath)
    project.load()
    
    obj = project.mapper.selectById(int(objectId))
    if obj.isPointer():
        obj = obj.get()

    #Initialize default values
    inputParameters = {'allowRender': True,
                       'mode': 'gallery', #Gallery,table,column,volume_astex, volume_chimera
                       'zoom': '150px',
                       'dims': '2d',
                       'goto': 1,
                       'colRowMode': 'Off',
                       'mirrorY': False,
                       'applyTransformMatrix': False,
                       'onlyShifts': False,
                       'wrap': False,
                       'resliceComboBox': xmipp.VIEW_Z_NEG,
                       'imageMaxWidth': 512,
                       'imageMinWidth': 30,
                       'imageMaxHeight': 512,
                       'imageMinHeight': 30}      
    
    if probandoCTFParam:
        inputParameters['path']= join(projectPath, "Runs/XmippProtCTFMicrographs175/extra/BPV_1386/xmipp_ctf.ctfparam")
        inputParameters['mode']= 'column'
    elif isinstance(obj, SetOfMicrographs):
        fn = project.getTmpPath(obj.getName() + '_micrographs.xmd')
        writeSetOfMicrographs(obj, fn)
        inputParameters['path']= join(projectPath, fn)
    elif isinstance(obj, SetOfVolumes):
        fn = project.getTmpPath(obj.getName()+ '_volumes.xmd')
        writeSetOfVolumes(obj, fn)
        inputParameters['path']= join(projectPath, fn)
        inputParameters['setOfVolumes']= obj
        inputParameters['setOfVolumesId']= obj.getObjId()
        inputParameters['dims']= '3d'
        inputParameters['mode']= 'table'
    elif isinstance(obj, SetOfImages):
#        PAJM aqui falla para el cl2d align y se esta perdiendo la matrix de transformacion en la conversion
        fn = project.getTmpPath(obj.getName() + '_images.xmd')
        writeSetOfParticles(obj, fn)
        inputParameters['path']= join(projectPath, fn)
    elif isinstance(obj, Image):
        fn = project.getTmpPath(obj.getName() + '_image.xmd')
        writeSetOfParticles(obj, fn)
        inputParameters['path']= join(projectPath, fn)
        
        
    elif isinstance(obj, SetOfClasses2D):
        fn = project.getTmpPath(obj.getName() + '_classes.xmd')
        writeSetOfClasses2D(obj, fn)
        inputParameters['path']= join(projectPath, fn)
#        runShowJ(obj.getClassesMdFileName())
    else:
        raise Exception('Showj Web visualizer: can not visualize class: %s' % obj.getClassName())

    if isinstance(obj, SetOfVolumes) and probandoVolume:
        return render_to_response('volume_visualization.html', inputParameters)
    else:
        return showj(request, inputParameters)  
