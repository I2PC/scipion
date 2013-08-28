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
    
def showj(request, inputParameters=None):
    
    # Resources #
    # Style Sheets
    css_path = join(settings.STATIC_URL, 'css/showj_style.css')

    #Favicon
#    favicon_path = getResource('favicon')
    
    #General jquery libs
    jquery_path = join(settings.STATIC_URL, 'js/jquery.js')
    jquery_cookie = join(settings.STATIC_URL, 'js/jquery.cookie.js')
    jquery_treeview = join(settings.STATIC_URL, 'js/jquery.treeview.js')
    launchTreeview = join(settings.STATIC_URL, 'js/launchTreeview.js')
    utils_path = join(settings.STATIC_URL, 'js/utils.js')
    
    #Table View jquery libs
    jquerydataTables_path = join(settings.STATIC_URL, 'js/jquery.dataTables.js')
    jquerydataTables_colreorder_path = join(settings.STATIC_URL, 'js/ColReorder.js')
    jeditable_path = join(settings.STATIC_URL, 'js/jquery.jeditable.js')
    jquery_ui_path = join(settings.STATIC_URL, 'js/jquery-ui.js')
    jquery_waypoints_path = join(settings.STATIC_URL, 'js/waypoints.min.js')    
#    jquery_tabletools_path = join(settings.STATIC_URL, 'js/TableTools.min.js')
    
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
            inputParameters['labelsToRenderComboBox']=labelsToRenderComboBoxValues[0][0] if len(labelsToRenderComboBoxValues) > 0 else ''
       
        _imageDimensions = get_image_dimensions(request.session['projectPath'], tableDataset.getElementById(0,inputParameters['labelsToRenderComboBox']))  if inputParameters['labelsToRenderComboBox']!='' else None
            
        inputParameters['blockComboBox']=_blockComboBox
        
        showjForm = ShowjForm(dataset,
                              tableLayoutConfiguration,
                              inputParameters) # An unbound form

    if showjForm.is_valid() is False:
        print showjForm.errors
        
    #Store dataset and labelsToRender in session 
    request.session['dataset'] = dataset
    request.session['labelsToRenderComboBox'] = _labelsToRenderComboBox
    if (_imageDimensions != ''):
        request.session['imageDimensions'] = _imageDimensions

    #Create context to be send
    context = {'jquery': jquery_path, #Configuration variables
               'utils': utils_path,
               'jquery_cookie': jquery_cookie,
               'jquery_treeview': jquery_treeview,
               'launchTreeview': launchTreeview,
               'jquery_datatable': jquerydataTables_path,
               'jquerydataTables_colreorder': jquerydataTables_colreorder_path,
               'jeditable': jeditable_path,
               'jquery_ui': jquery_ui_path,
               'jquery_waypoints':jquery_waypoints_path,
               'css': css_path,
               'tableLayoutConfiguration' : tableLayoutConfiguration if (showjForm.data['mode']=='gallery') else json.dumps({'columnsLayout': tableLayoutConfiguration.columnsLayout, 'colsOrder': tableLayoutConfiguration.colsOrder}, ensure_ascii=False, cls=ColumnLayoutConfigurationEncoder), #Data variables
#               'tableLayoutConfiguration' : json.dumps({'columnsLayout': tableLayoutConfiguration.columnsLayout, 'colsOrder': tableLayoutConfiguration.colsOrder}, ensure_ascii=False, cls=ColumnLayoutConfigurationEncoder),
               'tableDataset': tableDataset,
               'imageDimensions': request.session['imageDimensions'],
               'defaultZoom': request.session['defaultZoom'], 
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

        
#        self.labels = [xmipp.label2Str(l) for l in labels]
#        self.typeOfColumns = getTypeOfColumns(labels, allowRender)
#        self.colsOrder = defineColsLayout(self.labels)
#        #Esto es un napeidus que habria que arreglar
#        self.labels_typeOfColumns= zip(self.labels,self.typeOfColumns)

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
        

#PAJM         
def getTypeOfColumn(label):
    if (label == "id"):
        return "id"
    elif (xmipp.labelIsImage(str(label))):
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
    
    from django.http import HttpResponse
    import json
    from django.utils import simplejson
    print "akiiiiiiii"
    
    if request.is_ajax():
        changes = request.POST.get('changes')
        print "changes", changes
        print "type", type(changes)
        request.session['dataset']
#        NAPA DE LUXE
        
#        element_id = request.GET.get('element_id')
#        try:
#            label, idRow = element_id.split("___")
#        except ValueError:
#            return HttpResponse(json.dumps({'message':'Error'}), mimetype='application/javascript')
#        
##        if len(element_id_split)!=2: 
##            print "esto peto y hay que hacer alguna movidita"
#        element_value= request.GET.get('element_value')
#        #conversion for checkbox element
#        if (element_value == 'true'): element_value = 1
#        else: element_value = 0
#        
#        md = request.session['md']
#
#        for index, mdObject in enumerate(md.objects[int(idRow)].values):
#            if label in mdObject.label:
#                md.objects[int(idRow)].values[index].strValue=element_value
#        
#        request.session['md']=md
##        mdXmipp.setValue(element_id_split[0], False, mdXmipp[element_id_split[1]])
#        return HttpResponse(json.dumps({'message':'Ok'}), mimetype='application/javascript')

#    request.get.get('value')

def get_image_dimensions(projectPath, imagePath):
    from django.http import HttpResponse
    from pyworkflow.gui import getImage
    imageNo = None
#    imagePath = request.GET.get('image')

    
    # PAJM: Como vamos a gestionar lsa imagen    
    if imagePath.endswith('png') or imagePath.endswith('gif'):
        img = getImage(imagePath, tk=False)
    else:
        if '@' in imagePath:
            parts = imagePath.split('@')
            imageNo = parts[0]
            imagePath = parts[1]
            
        if projectPath != '':
            imagePathTmp = join(projectPath,imagePath)
            if not os.path.isfile(imagePathTmp):
                imagePath = getInputPath('showj', imagePath)      
            

#        imagePath = join(request.session['projectPath'],imagePath)
        
        if imageNo:
            imagePath = '%s@%s' % (imageNo, imagePath) 
            
        imgXmipp = xmipp.Image(imagePath)
        
        return imgXmipp.getDimensions()
        

def get_image(request):
    from django.http import HttpResponse
    from pyworkflow.gui import getImage, getPILImage
    imageNo = None
    imagePath = request.GET.get('image')
    imageDim = request.GET.get('dim', 150)
    
    # PAJM: Como vamos a gestionar lsa imagen    
    if imagePath.endswith('png') or imagePath.endswith('gif'):
        img = getImage(imagePath, tk=False)
    else:
        if '@' in imagePath:
            parts = imagePath.split('@')
            imageNo = parts[0]
            imagePath = parts[1]
            
        if 'projectPath' in request.session:
            imagePathTmp = join(request.session['projectPath'],imagePath)
            if not os.path.isfile(imagePathTmp):
                imagePath = getInputPath('showj', imagePath)      
            

#        imagePath = join(request.session['projectPath'],imagePath)
        
        if imageNo:
            imagePath = '%s@%s' % (imageNo, imagePath) 
            
        imgXmipp = xmipp.Image(imagePath)
        
        # from PIL import Image
        img = getPILImage(imgXmipp, imageDim)
        
        
        
    # response = HttpResponse(mimetype="image/png")    
    response = HttpResponse(mimetype="image/png")
    img.save(response, "PNG")
    return response
