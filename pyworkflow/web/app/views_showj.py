import os
import xmipp
from django.http import HttpResponse
from pyworkflow.web.pages import settings
from django.shortcuts import render_to_response
from pyworkflow.tests import getInputPath
from pyworkflow.web.app.forms import ShowjForm, getBlockComboBoxValues, getMetadataComboBoxValues
from django.template import RequestContext
from os.path import join    
    
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
    jquery_tabletools_path = join(settings.STATIC_URL, 'js/TableTools.min.js')
        
    
    #############
    # WEB INPUT PARAMETERS
    if request.method == 'POST': # If the form has been submitted...
        print "POST METHOD"
        mdXmipp = loadMetaDataXmipp(request.POST.get('path'), request.POST.get('blockComboBox'))
        showjForm = ShowjForm(mdXmipp, request.POST) # A form bound to the POST data
    else:
        print "GET METHOD"
        
        if inputParameters == None:
            mdXmipp = loadMetaDataXmipp(request.GET.get('path', 'tux_vol.xmd'), request.GET.get('block',''))
            inputParameters = {'path': request.GET.get('path', 'tux_vol.xmd'),
                         'blockComboBox': request.GET.get('block', getBlockComboBoxValues(request.GET.get('path', 'tux_vol.xmd'))[0][0]),
                         'metadataComboBox': getMetadataComboBoxValues(mdXmipp, 'render' in request.GET)[0][0],
                         'allowRender': 'render' in request.GET,
                         'zoom' : request.GET.get('dim', 150),
                         'mode': request.GET.get('mode', 'gallery'),
                         'goto': 1,
                         'colRowMode': request.GET.get('colRowMode', 'Off')}
        else:
            mdXmipp = loadMetaDataXmipp(inputParameters['path'], inputParameters['blockComboBox'] if ('blockComboBox' in inputParameters) else '')
            if 'blockComboBox' not in inputParameters:
                inputParameters['blockComboBox']=getBlockComboBoxValues(inputParameters['path'])[0][0]
            if 'metadataComboBox' not in inputParameters:
                inputParameters['metadataComboBox']=getMetadataComboBoxValues(mdXmipp, inputParameters['allowRender'])[0][0]

        showjForm = ShowjForm(mdXmipp, inputParameters) # An unbound form
        
    if showjForm.is_valid() is False:
        print showjForm.errors
    
    md = MdData(mdXmipp, showjForm.data['allowRender'], showjForm.data['zoom'])
    request.session['md'] = md

#    menuLayoutConfig = MenuLayoutConfig(showjForm.data['mode'], showjForm.data['path'], showjForm.data['blockComboBox'], showjForm.data['allowRender'], showjForm.data['zoom'])
    
    context = {'jquery': jquery_path,
               'utils': utils_path,
               'jquery_cookie': jquery_cookie,
               'jquery_treeview': jquery_treeview,
               'launchTreeview': launchTreeview,
               'jquery_datatable': jquerydataTables_path,
               'jquerydataTables_colreorder': jquerydataTables_colreorder_path,
               'jeditable': jeditable_path,
               'jquery_ui': jquery_ui_path,
               'jquery_waypoints':jquery_waypoints_path,
#               'jquery_tabletools':jquery_tabletools_path,
               'css': css_path,
               
               'metadata': md,
               
#               'inputParameters': inputParameters,
#               'menuLayoutConfig': menuLayoutConfig,
               'form': showjForm}
    
    return_page = '%s%s%s' % ('showj_', showjForm.data['mode'], '.html')

    return render_to_response(return_page, RequestContext(request, context))

#class InitialValuesShowj():
#    def __init__(self, md, path, allowRender

AT = '__at__'

class MdObj():
    pass
    
class MdValue():
    def __init__(self, md, label, objId, typeOfColumn):
             

        self.label = xmipp.label2Str(label)
        
#        self.allowRender = allowRender

        # check if enabled label
#        self.displayCheckbox = (label == xmipp.MDL_ENABLED)

        
        self.strValue = str(md.getValue(label, objId))   
        
        # Prepare path for image
        self.imgValue = self.strValue
        
        self.typeOfColumn = typeOfColumn
        
        if typeOfColumn=="image" and '@' in self.strValue:
            self.imgValue = self.imgValue.replace('@', AT)
#            if imageDim:
#                self.imgValue += '&dim=%s' % imageDim

class MdData():
    def __init__(self, md, allowRender=True, imageDim=None):        
        labels = md.getActiveLabels()
        
        self.tableLayoutConfiguration = TableLayoutConfiguration(labels, allowRender)
        
        self.objects = []
        for objId in md:
            obj = MdObj()
            #PAJM que es este objId
            obj.id = objId
            obj.values = [MdValue(md, l, objId, typeOfColumn) for l, typeOfColumn in zip(labels, self.tableLayoutConfiguration.typeOfColumns)]
            self.objects.append(obj)
            
class TableLayoutConfiguration():
    def __init__(self, labels, allowRender=True): 
        self.labels = [xmipp.label2Str(l) for l in labels]
        self.typeOfColumns = getTypeOfColumns(labels, allowRender)
        self.colsOrder = defineColsLayout(self.labels)
        #Esto es un napeidus que habria que arreglar
        self.labels_typeOfColumns= zip(self.labels,self.typeOfColumns)

  
         
def getTypeOfColumns(label, allowRender):
    typeOfColumns = []
    for l in label:
        if (xmipp.labelIsImage(l) and allowRender):
            typeOfColumns.append("image")
        elif (l == xmipp.MDL_ENABLED):
            typeOfColumns.append("checkbox")
        else:
            typeOfColumns.append("text")    
    
    return typeOfColumns
        
def defineColsLayout(labels):
    colsOrder = range(len(labels))
    if 'enabled' in labels:
        colsOrder.insert(0, colsOrder.pop(labels.index('enabled')))
    return colsOrder    

def loadMetaDataXmipp(path, block):
    path= getInputPath('showj', path)
    if len(block):
        path = '%s@%s' % (block, path)
    return xmipp.MetaData(path)
        
    # path2 = 'Volumes@' + path1
#    return MdData(path, allowRender, imageDim)   

def save_showj_metadata(request):    
    from django.http import HttpResponse
    import json
    from django.utils import simplejson
    
    if request.is_ajax():
        path = request.GET.get('path')
        
#        md = request.session['md']
#        mdXmipp = xmipp.MetaData()
#        mdXmipp.write(path)
        
        return HttpResponse(json.dumps({'message':'Ok'}), mimetype='application/javascript')

def save_showj_table(request):
    
    from django.http import HttpResponse
    import json
    from django.utils import simplejson
    
    
    if request.is_ajax():
        element_id = request.GET.get('element_id')
        try:
            label, idRow = element_id.split("___")
        except ValueError:
            return HttpResponse(json.dumps({'message':'Error'}), mimetype='application/javascript')
        
#        if len(element_id_split)!=2: 
#            print "esto peto y hay que hacer alguna movidita"
        element_value= request.GET.get('element_value')
        #conversion for checkbox element
        if (element_value == 'true'): element_value = 1
        else: element_value = 0
        
        md = request.session['md']

        for index, mdObject in enumerate(md.objects[int(idRow)].values):
            if label in mdObject.label:
                md.objects[int(idRow)].values[index].strValue=element_value
        
        request.session['md']=md
#        mdXmipp.setValue(element_id_split[0], False, mdXmipp[element_id_split[1]])
        return HttpResponse(json.dumps({'message':'Ok'}), mimetype='application/javascript')

#    request.get.get('value')

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
        if AT in imagePath:
            parts = imagePath.split(AT)
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
