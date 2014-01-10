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
from django.http import HttpResponse
from pyworkflow.web.pages import settings
from django.shortcuts import render_to_response
from pyworkflow.tests import getInputPath
from django.template import RequestContext
from pyworkflow.web.app.views_util import *
from forms import VolVisualizationForm, ShowjForm
from pyworkflow.em import *
from pyworkflow.em.packages.xmipp3.convert import *
from layout_configuration import *



def showj(request, inputParameters=None, extraParameters=None):
    #############
    # WEB INPUT PARAMETERS
    _imageDimensions = ''
    _imageVolName=''
    _stats=None
    
    if request.method == 'POST': # If the form has been submitted... Post method
        inputParameters=request.POST.copy()
        
    #Init Dataset
    #NAPA DE LUXE: Check type of Dataset 
    dataset = loadDatasetXmipp(inputParameters['path']) 
    if inputParameters['blockComboBox'] == '':
        inputParameters['blockComboBox'] = dataset.listTables()[0]
    
    #Get table from block name
    dataset.setTableName(inputParameters['blockComboBox'])
    tableDataset=dataset.getTable()

    if request.method == 'GET':
        request.session['defaultColumnsLayoutProperties'] = getExtraParameters(extraParameters, tableDataset)
    
    #Load table layout configuration. How to display columns and attributes (visible, render, editable)  
    inputParameters['tableLayoutConfiguration']= TableLayoutConfiguration(dataset, tableDataset, inputParameters['allowRender'], request.session['defaultColumnsLayoutProperties'])

    #If no label is set to render, set the first one if exists
    if 'labelsToRenderComboBox' not in inputParameters or inputParameters['labelsToRenderComboBox'] == '' or request.session['blockComboBox'] != inputParameters['blockComboBox']:
        labelsToRenderComboBoxValues = inputParameters['tableLayoutConfiguration'].getLabelsToRenderComboBoxValues()
        if len(labelsToRenderComboBoxValues) > 0:
            inputParameters['labelsToRenderComboBox'] = labelsToRenderComboBoxValues[0][0]
        else:
            # If there is no image to display and it is initial load, switch to table mode 
            if request.method == 'GET' and inputParameters['mode']!='table':
                inputParameters['mode']='table'
                showj(request, inputParameters, extraParameters) 
            inputParameters['labelsToRenderComboBox'] = ''

    dataset.setLabelToRender(inputParameters['labelsToRenderComboBox']) 

    if inputParameters['labelsToRenderComboBox'] == '':
        inputParameters['zoom']=0
        _imageDimensions = None
        dataset.setNumberSlices(0)
    else:
        
        _imageVolName = inputParameters['volumesToRenderComboBox'] if ("volumesToRenderComboBox" in inputParameters and inputParameters['volumesToRenderComboBox'] != '') else tableDataset.getElementById(0,inputParameters['labelsToRenderComboBox'])
        _typeOfColumnToRender = inputParameters['tableLayoutConfiguration'].columnsLayout[inputParameters['labelsToRenderComboBox']].typeOfColumn
        _imageDimensions = readDimensions(request,
                                          _imageVolName,
                                          _typeOfColumnToRender)
        dataset.setNumberSlices(_imageDimensions[2])
        
        if _typeOfColumnToRender == "image":
            _convert = dataset.getNumberSlices()>1 and (inputParameters['mode']=='gallery' or inputParameters['mode']=='volume_astex' or inputParameters['mode']=='volume_chimera')
            _reslice = dataset.getNumberSlices()>1 and inputParameters['mode']=='gallery'
            _getStats = dataset.getNumberSlices()>1 and (inputParameters['mode']=='volume_astex' or inputParameters['mode']=='volume_chimera')
            _dataType = xmipp.DT_FLOAT if dataset.getNumberSlices()>1 and inputParameters['mode']=='volume_astex' else xmipp.DT_UCHAR
                
            _imageVolName, _stats = readImageVolume(request, _imageVolName, _convert, _dataType, _reslice, int(inputParameters['resliceComboBox']), _getStats)
            
            if dataset.getNumberSlices()>1:
                inputParameters['tableLayoutConfiguration'].columnsLayout[inputParameters['labelsToRenderComboBox']].columnLayoutProperties.renderFunc = "get_slice"
                dataset.setVolumeName(_imageVolName)
                

    #Store dataset and labelsToRender in session 
    request.session['dataset'] = dataset
    request.session['labelsToRenderComboBox'] = inputParameters['labelsToRenderComboBox']
    request.session['blockComboBox'] = inputParameters['blockComboBox']
    request.session['imageDimensions'] = _imageDimensions

    showjForm = ShowjForm(dataset,
                          inputParameters['tableLayoutConfiguration'],
                          #request.POST if request.method == 'POST' else inputParameters) # A form bound for the POST data and unbound for the GET
                          inputParameters) # A form bound for the POST data and unbound for the GET
        
    if showjForm.is_valid() is False:
        print showjForm.errors

    context = createContext(dataset, tableDataset, inputParameters['tableLayoutConfiguration'], request, showjForm)

    if inputParameters['mode']=='volume_astex' or inputParameters['mode']=='volume_chimera':

        threshold = (_stats[2] + _stats[3])/2 if _stats != None else 1
        volPath = os.path.join(request.session['projectPath'], _imageVolName)
        
        context.update({"threshold":threshold,
                        'minStats':_stats[2] if _stats != None else 1,
                        'maxStats':_stats[3] if _stats != None else 1, })
        
        if inputParameters['mode']=='volume_astex':
    
            #        volPath='/home/adrian/Scipion/tests/input/showj/emd_1042.map'        
#            fileName, fileExtension = os.path.splitext(_imageVolName)
            
            linkName = 'test_link_' + request.session._session_key + '.map'
            volLinkPath = os.path.join(pw.WEB_RESOURCES, 'astex', 'tmp', linkName)
            from pyworkflow.utils.path import cleanPath, createLink
            cleanPath(volLinkPath)
            createLink(volPath, volLinkPath)
            volLink = os.path.join('/', settings.STATIC_ROOT, 'astex', 'tmp', linkName)
            
            context.update({"volLink":volLink})
         
        elif inputParameters['mode']=='volume_chimera':   
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
            
            context.update({"chimeraHtml":chimeraHtml})
               
#               'volType': 2, #0->byte, 1 ->Integer, 2-> Float
    
    return_page = '%s%s%s' % ('showj_', showjForm.data['mode'], '.html')
    return render_to_response(return_page, RequestContext(request, context))

def createContext(dataset, tableDataset, tableLayoutConfiguration, request, showjForm):
    #Create context to be send
    return {'showj_css': getResourceCss('showj'),
            'general_css': getResourceCss('general'),
            'smoothness': getResourceCss('ui_smoothness'),
            'demo_table_jui': getResourceCss('showj_demo_table_jui'),
           
            'favicon': getResourceIcon('favicon'),
            'logo': getResourceIcon('logo_scipion'),
            'logo_transparent': getResourceIcon('logo_scipion_transparent'),
           
            'jquery': getResourceJs('jquery'), #Configuration variables
            'jquery_datatable': getResourceJs('jquery_datatables'),
            'jquerydataTables_colreorder': getResourceJs('jquery_colreorder'),
            'jquerydataTables_colreorder_resize': getResourceJs('jquery_colreorder_resize'),
            'jeditable': getResourceJs('jquery_editable'),
            'jquery_waypoints':getResourceJs('jquery_waypoints'),
            'jquery_hover_intent':getResourceJs('jquery_hover_intent'),
            
            'transpose_lib':getResourceJs('transpose'),
           
            'dataset': dataset,
            'tableLayoutConfiguration': json.dumps({'columnsLayout': tableLayoutConfiguration.columnsLayout,
                                                   'colsOrder': tableLayoutConfiguration.colsOrder},
                                                   ensure_ascii=False,
                                                   cls=ColumnLayoutConfigurationEncoder), #Data variables
            'tableDataset': tableDataset,
           
            'imageDimensions': request.session['imageDimensions'],
            'defaultZoom': request.session['defaultZoom'],
            'projectName': request.session['projectName'],
           
            'form': showjForm
            }


def getExtraParameters(extraParameters, tableDataset):
    print "extraParameters",extraParameters 
    defaultColumnsLayoutProperties = None
    if extraParameters != None and extraParameters != {}:
        defaultColumnsLayoutProperties = {k.getName(): {} for k in tableDataset.iterColumns()}
        for key, value in extraParameters.iteritems():

            try:
                label, attribute = key.rsplit('___')
            except Exception, ex:
                raise Exception('Showj Web visualizer: Incorrect extra parameter')

            if tableDataset.hasColumn(label):
                defaultColumnsLayoutProperties[label].update({attribute:value})
                
    return defaultColumnsLayoutProperties

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
#         inputVolume = os.path.join(os.getcwd(),volume.getFileName())
#         outputHtmlFile = os.path.join(os.getcwd(),project.getTmpPath("volume_" + str(volume.getObjId()) + '.html'))
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
    inputParameters = {'objectId': '',
                       'path': '',
                       'allowRender': True,                 # Image can be displayed, depending on column layout 
                       'mode': 'gallery',                   # Mode Options: gallery, table, column, volume_astex, volume_chimera
                       'zoom': '150px',                     # Zoom set by default
                       'blockComboBox': '',                 # Metadata Block to display. If None the first one will be displayed
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
                       'imageMinWidth': 30,                 # Minimum image width (in pixels)
                       'imageMaxHeight': 512,               # Maximum image height (in pixels)
                       'imageMinHeight': 30}                # Minimum image height (in pixels)


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
    # extraParameters["id___visible"]=True
    # extraParameters["micrograph___renderFunc"]="get_image_psd"
    # extraParameters["micrograph___extraRenderFunc"]="downsample=2"
            
    extraParameters = {}
    
# NAPA DE LUXE: PARA probarlo sin sesion iniciada    
#    request.session['projectPath'] = "/home/adrian/Scipion/projects/TestXmippWorkflow/"
#    request.session['projectName'] = "Juanitpo"
    
    if "projectPath" in request.session:
        projectPath = request.session['projectPath']
    else:         
        raise Exception('Showj Web visualizer: No project loaded')
    
    for key, value in request.GET.iteritems():
        if key in inputParameters:
            inputParameters[key]=value
        else:
            extraParameters[key]=value
    
    
    if "objectId" in request.GET: 
        project = Project(projectPath)
        project.load()
        
        obj = project.mapper.selectById(int(inputParameters["objectId"]))
        if obj.isPointer():
            obj = obj.get()

        if isinstance(obj, SetOfMicrographs):
            fn = project.getTmpPath(obj.getName() + '_micrographs.xmd')
            inputParameters['path'] = os.path.join(projectPath, createXmippInputMicrographs(None, obj, micsFn=os.path.join(projectPath, fn)))
#            writeSetOfMicrographs(obj, fn)
#            inputParameters['path']= os.path.join(projectPath, fn)
        elif isinstance(obj, SetOfVolumes):
            fn = project.getTmpPath(obj.getName()+ '_volumes.xmd')
            inputParameters['path'] = os.path.join(projectPath, createXmippInputVolumes(None, obj, volsFn=os.path.join(projectPath, fn))) 
            inputParameters['mode']= 'table'
            
#            writeSetOfVolumes(obj, fn)
#            inputParameters['path']= os.path.join(projectPath, fn)
        elif isinstance(obj, SetOfImages):
            fn = project.getTmpPath(obj.getName() + '_images.xmd')
            inputParameters['path'] = os.path.join(projectPath, createXmippInputImages(None, obj, imagesFn=os.path.join(projectPath, fn)))
            
#            writeSetOfParticles(obj, fn)
#            inputParameters['path']= os.path.join(projectPath, fn)
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
            inputParameters['path'] = os.path.join(projectPath, createXmippInputImages(None, imgSet, imagesFn=os.path.join(projectPath, fn)))
#            writeSetOfParticles(imgSet, fn)
#            inputParameters['path']= os.path.join(projectPath, fn)
            
        elif isinstance(obj, SetOfClasses2D):
            fn = project.getTmpPath(obj.getName() + '_classes.xmd')
#            writeSetOfClasses2D(obj, fn)
#            inputParameters['path']= os.path.join(projectPath, fn)
            inputParameters['path'] = os.path.join(projectPath, createXmippInputClasses2D(None, obj, classFn=os.path.join(projectPath, fn)))
        elif isinstance(obj, SetOfCTF):
            fn = project.getTmpPath(obj.getName() + '_ctfs.xmd')
            inputParameters['path'] = os.path.join(projectPath, createXmippInputCTF(None, obj, ctfFn=os.path.join(projectPath, fn)))
#            writeSetOfCTFs(obj, fn)
#            inputParameters['path']= os.path.join(projectPath, fn)
        else:
            raise Exception('Showj Web visualizer: can not visualize class: %s' % obj.getClassName())
    
    elif "path" in request.GET:
        inputParameters.update({'path':os.path.join(projectPath, request.GET.get("path"))})
    else:
        raise Exception('Showj Web visualizer: No object identifier or path found')         

    request.session['defaultZoom'] = inputParameters['zoom']
    return showj(request, inputParameters, extraParameters)  


def testingSSH(request):
    context = {}
#    return render_to_response("testing_ssh.html", RequestContext(request, context))
    return render_to_response("scipion.html", RequestContext(request, context))
