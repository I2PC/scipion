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
    
    #Init Dataset
    #Todo: Check type of Dataset 
    dataset = loadDatasetXmipp(_path) 
    
    if _blockComboBox == '':
        _blockComboBox = dataset.listTables()[0]
    
    #Get table from block name
    tableDataset=dataset.getTable(_blockComboBox)
    
    #Load table layout configuration. How to display columns and attributes (visible, render, editable)    
    tableLayoutConfiguration = TableLayoutConfiguration(tableDataset, _render)

    #Initialize Showj Form (button toolbar)
    if request.method == 'POST': # If the form has been submitted... Post method  
        if _labelsToRenderComboBox != '':
            if _labelsToRenderComboBox != request.session['labelsToRenderComboBox']:
                _imageDimensions = get_image_dimensions(request.session['projectPath'], tableDataset.getElementById(0,_labelsToRenderComboBox))
        else:
            _imageDimensions = None
        
        showjForm = ShowjForm(dataset,
                              tableLayoutConfiguration,
                              request.POST) # A form bound to the POST data
    else:
        if _labelsToRenderComboBox == '':
            labelsToRenderComboBoxValues = getLabelsToRenderComboBoxValues(tableLayoutConfiguration.columnsLayout)
            _labelsToRenderComboBox=labelsToRenderComboBoxValues[0][0] if len(labelsToRenderComboBoxValues) > 0 else ''
            inputParameters['labelsToRenderComboBox']=_labelsToRenderComboBox
            
        if inputParameters['labelsToRenderComboBox'] == '':
            inputParameters['zoom']=0
            _imageDimensions = None
        else:    
            _imageDimensions = get_image_dimensions(request.session['projectPath'], tableDataset.getElementById(0,inputParameters['labelsToRenderComboBox']))  
            
        inputParameters['blockComboBox']=_blockComboBox
        if 'mirrorY' not in inputParameters:
            inputParameters['mirrorY']=False
        
        showjForm = ShowjForm(dataset,
                              tableLayoutConfiguration,
                              inputParameters) # An unbound form

    if showjForm.is_valid() is False:
        print showjForm.errors
        
    #Store dataset and labelsToRender in session 
    request.session['dataset'] = dataset
    request.session['labelsToRenderComboBox'] = _labelsToRenderComboBox
    request.session['blockComboBox'] = _blockComboBox
    request.session['tableLayoutConfiguration'] = tableLayoutConfiguration
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
               
               'tableLayoutConfiguration' : tableLayoutConfiguration if (showjForm.data['mode']=='gallery') else json.dumps({'columnsLayout': tableLayoutConfiguration.columnsLayout, 'colsOrder': tableLayoutConfiguration.colsOrder}, ensure_ascii=False, cls=ColumnLayoutConfigurationEncoder), #Data variables
               'tableDataset': tableDataset,
               'imageDimensions': request.session['imageDimensions'],
               'defaultZoom': request.session['defaultZoom'],
               'projectName': request.session['projectName'],
               'form': showjForm} #Form
    
    return_page = '%s%s%s' % ('showj_', showjForm.data['mode'], '.html')
    return render_to_response(return_page, RequestContext(request, context))

                 
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
    def __init__(self, tableDataset, allowRender=True):
        
        self.columnsLayout = OrderedDict() 
         
        for col in tableDataset.iterColumns():
            self.columnsLayout[col.getName()]=ColumnLayoutConfiguration(col, allowRender)
            
        self.colsOrder = defineColsLayout(self.columnsLayout.keys())
        
            
class ColumnLayoutConfiguration():
    def __init__(self, col, allowRender):
        self.columns = col
        
        self.label = col.getName()
        self.typeOfColumn = getTypeOfColumn(col.getName())
        
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
        

def getTypeOfColumn(label):
    if (label == "id"):
        return "id"
    elif (label!='image_transformationMatrix' and xmipp.labelIsImage(str(label))):
        return "image"
    elif (label == "enabled"):
        return "checkbox"
    else:
        return "text"    
    
def defineColsLayout(labels):
    colsOrder = range(len(labels))
    if 'enabled' in labels:
        colsOrder.insert(0, colsOrder.pop(labels.index('enabled')))
    return colsOrder

################################### END LAYOUT ##########################    

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
        tableLayoutConfiguration=request.session['tableLayoutConfiguration']
        
        tableDataset=dataset.getTable(blockComboBox)
        
        for key in jsonChanges:
            element_split = key.rsplit('___')
            if len(element_split)!=2: 
                print "this fails and sth has to be done"
            
            #NAPA de LUXE ahora mismo se realiza una conversion a int pero habria que ver de que tipo de datos se trata 
            #tableLayoutConfiguration.columnsLayout[element_split[0]].typeOfColumn

            dictelement = {element_split[0]:int(jsonChanges[key])}
            tableDataset.updateRow(int(element_split[1]),**dictelement)
        
        dataset.writeTable(blockComboBox, tableDataset)
        
        return HttpResponse(json.dumps({'message':'Ok'}), mimetype='application/javascript')


def visualizeObject(request):
    probandoCTFParam = False
    
    objectId = request.GET.get("objectId")    
    #projectName = request.session['projectName']
    projectName = request.GET.get("projectName")
    
#    project = loadProject(projectName)
    project = Project(request.session['projectPath'])
    project.load()
    
    obj = project.mapper.selectById(int(objectId))
    if obj.isPointer():
        obj = obj.get()
        
    if probandoCTFParam:
        inputParameters = {'path': join(request.session['projectPath'], "Runs/XmippProtCTFMicrographs218/extra/BPV_1386/xmipp_ctf.ctfparam"),
               'allowRender': True,
               'mode': 'column',
               'zoom': '150px',
               'goto': 1,
               'colRowMode': 'Off'}    
    elif isinstance(obj, SetOfMicrographs):
        fn = project.getTmpPath(obj.getName() + '_micrographs.xmd')
        writeSetOfMicrographs(obj, fn)
        inputParameters = {'path': join(request.session['projectPath'], fn),
                       'allowRender': True,
                       'mode': 'gallery',
                       'zoom': '150px',
                       'goto': 1,
                       'colRowMode': 'Off'}
    elif isinstance(obj, SetOfVolumes):
        fn = project.getTmpPath(obj.getName()+ '_volumes.xmd')
        vols = XmippSetOfVolumes.convert(obj,fn)
        inputParameters = {'path': join(request.session['projectPath'], vols.getFileName()),
                           'setOfVolumes' : obj,
                           'setOfVolumesId': obj.getObjId(),
                           'dims': '3d',
                           'mode': 'table'}  
    elif isinstance(obj, SetOfImages):
        fn = project.getTmpPath(obj.getName() + '_images.xmd')
        imgs = XmippSetOfImages.convert(obj, fn)
        inputParameters = {'path': join(request.session['projectPath'], imgs.getFileName()),
               'allowRender': True,
               'mode': 'gallery',
               'zoom': '150px',
               'goto': 1,
               'colRowMode': 'Off'}

    elif isinstance(obj, XmippClassification2D):
        mdPath = obj.getClassesMdFileName()
        block, path = mdPath.split('@')
        inputParameters = {'path': join(request.session['projectPath'], path),
               'allowRender': True,
               'mode': 'gallery',
               'zoom': '150px',
               'goto': 1,
               'colRowMode': 'Off'}
#        runShowJ(obj.getClassesMdFileName())
    else:
        raise Exception('Showj Web visualizer: can not visualize class: %s' % obj.getClassName())

    if isinstance(obj, SetOfVolumes):
        return render_to_response('volume_visualization.html', inputParameters)
    else:
        return showj(request, inputParameters)
    
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
            # Astex viewer            
            from random import randint
            linkName = 'test_link_' + str(randint(0, 10000)) + '.map'
            volLinkPath = os.path.join(pw.HOME, 'web', 'pages', 'resources', 'astex', 'tmp', linkName)
            from pyworkflow.utils.path import cleanPath, createLink
            cleanPath(volLinkPath)
            createLink(volPath, volLinkPath)
#             os.system("ln -s " + str(volPath) + " " + volLinkPath)
            volLink = os.path.join('/', 'static', 'astex', 'tmp', linkName)
            # Chimera 
            from subprocess import Popen, PIPE, STDOUT
#             p = Popen(['chimera', '--start', 'ReadStdin', volPath], stdout=PIPE, stdin=PIPE, stderr=PIPE)
            p = Popen(['chimera', volPath], stdout=PIPE, stdin=PIPE, stderr=PIPE)
            outputHtmlFile = '/home/antonio/test.html'
            threshold = form.cleaned_data['threshold']
            stdout_data = p.communicate(input='volume #0 level ' + str(threshold) + '; export format WebGL ' + outputHtmlFile + '; stop')[0]
            f = open(outputHtmlFile)
            chimeraHtml = f.read().decode('string-escape').decode("utf-8").split("</html>")[1]
    else:
        form = VolVisualizationForm()
    context = {'MEDIA_URL' : settings.MEDIA_URL, 'STATIC_URL' :settings.STATIC_URL, 'form': form, 'volLink': volLink, 'chimeraHtml': chimeraHtml}    
    return render_to_response('showVolVisualization.html',  RequestContext(request, context))   

