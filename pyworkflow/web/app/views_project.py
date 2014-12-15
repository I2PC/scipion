# **************************************************************************
# *
# * Authors:    Jose Gutierrez (jose.gutierrez@cnb.csic.es)
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

from os.path import exists, join, basename
import json
from views_base import base_grid, base_flex
from views_util import loadProject, getResourceCss, getResourceIcon, getResourceJs
from views_tree import loadProtTree

from pyworkflow.manager import Manager
from pyworkflow.utils.path import copyFile

from django.http import HttpResponse, HttpRequest
from django.shortcuts import render_to_response
from pyworkflow.tests.tests import DataSet

def projects(request):
    from pyworkflow.utils.utils import prettyDate
    
    manager = Manager()
    projects = manager.listProjects()
    for p in projects:
        p.pTime = prettyDate(p.mTime)

    if 'projectName' in request.session: request.session['projectName'] = ""
    if 'projectPath' in request.session: request.session['projectPath'] = ""

    context = {'projects': projects,
               'projects_css': getResourceCss('projects'),
               'project_utils_js': getResourceJs('project_utils'),
               }
    
    context = base_grid(request, context)
    
    return render_to_response('projects.html', context)

def create_project(request):
    manager = Manager()
    
    if request.is_ajax():
        projectName = request.GET.get('projectName')
        manager.createProject(projectName)       
        
    return HttpResponse(mimetype='application/javascript')

def delete_project(request):
    manager = Manager()
    
    if request.is_ajax():
        projectName = request.GET.get('projectName')
        manager.deleteProject(projectName)       
        
    return HttpResponse(mimetype='application/javascript')


def getNodeStateColor(node):
    from pyworkflow.gui.project.viewprotocols import STATUS_COLORS
    from pyworkflow.protocol.constants import STATUS_FAILED
    
    color = '#ADD8E6'  # Lightblue
    status = ''
    if node.run:
        status = node.run.status.get(STATUS_FAILED)
        color = STATUS_COLORS[status]
        
    return status, color


def update_prot_tree(request):
    projectName = request.session['projectName']
    project = loadProject(projectName)
    project_settings = project.getSettings()
    index = request.GET.get('index', None)

    # set the new protocol tree chosen
    project_settings.setCurrentProtocolMenu(index)
    project_settings.write()
        
    return HttpResponse(mimetype='application/javascript')


def update_graph_view(request):
    status = request.GET.get('status', None)
    projectName = request.session['projectName']
    project = loadProject(projectName)
    project_settings = project.getSettings()

    if status == "True":
        project_settings.graphView.set(True)
    else :
        project_settings.graphView.set(False)
    project_settings.write()
    return HttpResponse(mimetype='application/javascript')


def save_selection(request):
    if request.is_ajax():
        mark = request.GET.get('mark', None)
        
        projectName = request.session['projectName']
        project = loadProject(projectName)
        project_settings = project.getSettings()
        
        # Set the selected runs stored in BD    
        project_settings.runSelection.set(mark)
        
        project_settings.write()
        
    return HttpResponse(mimetype='application/javascript')


def tree_prot_view(request):
    projectName = request.session['projectName'] 
    project = loadProject(projectName)
     
    # load the protocol tree current active
    htmlTree = loadProtTree(project)
    
    return render_to_response('project_content/tree_prot_view.html', {'protTreeHtml': htmlTree})
    
def run_table_graph(request):
    from pyworkflow.gui.tree import ProjectRunsTreeProvider
    
    try:
        projectName = request.session['projectName']
        project = loadProject(projectName)
        project_settings = project.getSettings()
        provider = ProjectRunsTreeProvider(project)
        
        runs = request.session['runs']
        runsNew = formatProvider(provider, "runs")
        
        refresh = False
        listNewElm = []
        
        if len(runs) != len(runsNew):
            print 'Change detected, different size'
            refresh = True
        else:
            for kx, vx in runs:
                for ky, vy in runsNew:
                    if kx == ky and vx != vy:
                        print 'Change detected', vx, vy
    #                   refresh = True
                        listNewElm.append(vy)
        if refresh:
            request.session['runs'] = runsNew
            
            # Get the selected runs stored in BD    
            selectedRuns = project_settings.runSelection
    
            # Get the mode view (list or graph) stored in BD
            graphView = project_settings.graphView.get()
            
            context = {'runs': runsNew,
                       'columns': provider.getColumns(),
                       'graphView': graphView, 
                       'selectedRuns' : selectedRuns}
            
            return render_to_response('project_content/run_table_graph.html', context)
        
        elif listNewElm:
            request.session['runs'] = runsNew
            jsonStr = json.dumps(listNewElm, ensure_ascii=False)
            return HttpResponse(jsonStr, mimetype='application/json')
        
        else:
            print "No changes detected"
            return HttpResponse("ok", mimetype='text/plain')
        
    except Exception:
        print "Stopped script"
        return HttpResponse("stop", mimetype='text/plain')


def formatProvider(provider, mode):
    objs = []
    for obj in provider.getObjects():
        objInfo = provider.getObjectInfo(obj)
        
        id = objInfo["key"]
        name = objInfo["text"]
        info = objInfo["values"]
        
        if mode == "runs":
            status = info[0]
            time = info[1]
            objs.append((id, [id, name, status, time]))
        
        elif mode == "objects":
            objs.append((id, [id, name, info]))
        
    return objs

def project_content(request):        
    from pyworkflow.gui.tree import ProjectRunsTreeProvider
    
    projectName = request.GET.get('projectName', None)
    mode = request.GET.get('mode', None)
    
    if projectName is None:
        projectName = request.POST.get('projectName', None)
        
    request.session['projectName'] = projectName
    manager = Manager()
    request.session['projectPath'] = manager.getProjectPath(projectName)
   
    project = loadProject(projectName)
    project_settings = project.getSettings()
    
    provider = ProjectRunsTreeProvider(project)
    runs = formatProvider(provider, "runs")
    request.session['runs'] = runs

    # Get the selected runs stored in BD    
    selectedRuns = project_settings.runSelection

    # Get the mode view (list or graph) stored in BD
    graphView = project_settings.graphView.get()
    
    # load the protocol tree current active
    htmlTree = loadProtTree(project)
    
    # get the choices to load protocol trees
    choices = [pm.text.get() for pm in project_settings.protMenuList]

    # get the choice current 
    choiceSelected =  project_settings.protMenuList.getIndex()
    
    # show the project name in the header.html
    projectNameHeader = 'Project '+ str(projectName)
    
    context = {'projectName': projectName,
               'mode': mode,
               'view':'protocols',
               'editTool': getResourceIcon('edit_toolbar'),
               'copyTool': getResourceIcon('copy_toolbar'),
               'deleteTool': getResourceIcon('delete_toolbar'),
               'browseTool': getResourceIcon('browse_toolbar'),
               'stopTool': getResourceIcon('stop_toolbar'),
               'analyzeTool': getResourceIcon('analyze_toolbar'),
               'treeTool': getResourceIcon('tree_toolbar'),
               'listTool': getResourceIcon('list_toolbar'),
               'graph_utils': getResourceJs('graph_utils'),
               'project_content_utils': getResourceJs('project_content_utils'),
               'jquery_cookie': getResourceJs('jquery_cookie'),
               'jquery_treeview': getResourceJs('jquery_treeview'),
               'project_content_css':getResourceCss('project_content'),
               'protTreeHtml': htmlTree,
               'choices':choices,
               'choiceSelected': choiceSelected,
               'runs': runs,
               'columns': provider.getColumns(),
               'projectNameHeader': projectNameHeader, 
               'graphView': graphView,
               'selectedRuns': selectedRuns
               }
    
    if mode == 'service':
        context = service_context(context)
    
    context = base_flex(request, context)
    
    return render_to_response('project_content/project_content.html', context)

def service_context(context):
    context.update({'importAverages': getResourceIcon('importAverages'),
                    'useProtocols': getResourceIcon('useProtocols'),
                    'protForm': getResourceIcon('protForm'),
                    'summary': getResourceIcon('summary'),
                    'showj': getResourceIcon('showj'),
                    'download': getResourceIcon('download'),
                    })
    return context

def protocol_info(request):
    from pyworkflow.web.app.views_util import parseText
    from pyworkflow.em.data import EMObject  
    
#    print "ENTER IN PROTOCOL INFO METHOD"

    if request.is_ajax():
        jsonStr = ''
        projectName = request.session['projectName']
        protId = request.GET.get('protocolId', None)

        project = loadProject(projectName)
        
        if len(protId) > 0: 
            protocol = project.getProtocol(int(protId))
            
            # PROTOCOL IO
            input_obj = [{'name':name, 
                          'nameId': attr.getNameId(), 
                          'id': attr.getObjId(), 
                          'info': str(attr.get()) if attr.isPointer() else str(attr)} 
                         for name, attr in protocol.iterInputAttributes()]
            
            output_obj = [{'name':name, 
                           'nameId': attr.getNameId(), 
                           'id': attr.getObjId(), 
                           'info': str(attr)} 
                          for name, attr in protocol.iterOutputAttributes(EMObject)]
    
            # PROTOCOL SUMMARY
            summary = parseText(protocol.summary())
            
            # PROTOCOL METHODS
            methods = parseText(protocol.methods())
    
            # STATUS
            status = protocol.status.get()
            
            # LOGS (ERROR & OUTPUT)
            fOutString, fErrString, fScpnString = protocol.getLogsAsStrings()
    
            ioDict = {'inputs': input_obj,
                      'outputs': output_obj,
                      'summary': summary,
                      'methods': methods, 
                      'status': status,
                      'logs_out': parseText(fOutString),
                      'logs_error': parseText(fErrString),
                      'logs_scipion': parseText(fScpnString)
                      }
            
    #        print "ioDict: ", ioDict
            jsonStr = json.dumps(ioDict, ensure_ascii=False)
        
    return HttpResponse(jsonStr, mimetype='application/javascript')

#===============================================================================
# SERVICE PROJECTS
#===============================================================================

def service_projects(request):
    #Example Projects to be showed
    exProjects = {'TestXmippWorkflow'}
    
    if 'projectName' in request.session: request.session['projectName'] = ""
    if 'projectPath' in request.session: request.session['projectPath'] = ""

    context = {'projects_css': getResourceCss('projects'),
               'project_utils_js': getResourceJs('project_utils'),
               'hiddenTreeProt': True,
               }
    
    context = base_grid(request, context)
    return render_to_response('service_projects.html', context)


def writeCustomMenu(customMenu):
    if not exists(customMenu):
        f = open(customMenu, 'w+')
        f.write('''
[PROTOCOLS]

Initial_Volume = [
    {"tag": "section", "text": "1. Upload data", "children": [
        {"tag": "protocol", "value": "ProtImportAverages",     "text": "Import averages", "icon": "bookmark.png"}]},
    {"tag": "section", "text": "2. Create a 3D volume", "children": [
        {"tag": "protocol", "value": "XmippProtRansac", "text": "xmipp3 - ransac"},
        {"tag": "protocol", "value": "EmanProtInitModel", "text": "eman2 - Initial volume"},
        {"tag": "protocol", "value": "XmippProtReconstructSignificant", "text": "xmipp3 - significant"},
        {"tag": "protocol", "value": "ProtPrime", "text": "simple - prime"}]},
    {"tag": "section", "text": "3. Align volumes.", "children": [
        {"tag": "protocol", "value": "XmippProtAlignVolume", "text": "xmipp3 - align volumes"}]}]
            ''')
        f.close()
        
def create_service_project(request):
    if request.is_ajax():
        import os
        from pyworkflow.object import Pointer
        from pyworkflow.em.protocol import ProtUnionSet, ProtImportAverages
        from pyworkflow.em.packages.xmipp3 import XmippProtRansac
        from pyworkflow.em.packages.eman2 import EmanProtInitModel
        
        # Create a new project
        manager = Manager()
        projectName = request.GET.get('projectName')
        
        # Filename to use as test data 
        testDataKey = request.GET.get('testData')
        
        #customMenu = os.path.join(os.path.dirname(os.environ['SCIPION_MENU']), 'menu_initvolume.conf')
        customMenu = os.path.join(os.environ['HOME'], '.config/scipion/menu_initvolume.conf')
        writeCustomMenu(customMenu)
        confs = {'protocols': customMenu}
        project = manager.createProject(projectName, confs, graphView=True)   
        
        # 1. Import averages
        protImport = project.newProtocol(ProtImportAverages,
                                         objLabel='import averages')
        
        # If using test data execute the import averages run
        # options are set in 'project_utils.js'
        dsMDA = DataSet.getDataSet('initial_volume')
        
        if testDataKey :
            fn = dsMDA.getFile(testDataKey)
            newFn = join(project.uploadPath, basename(fn))
            copyFile(fn, newFn)
            
            protImport.filesPath.set(newFn)
            protImport.samplingRate.set(1.)
            
            project.launchProtocol(protImport, wait=True)
        else:
            project.saveProtocol(protImport)
            
        
        # 2a. Ransac 
        protRansac = project.newProtocol(XmippProtRansac)
        protRansac.setObjLabel('xmipp - ransac')
        protRansac.inputSet.set(protImport)
        protRansac.inputSet.setExtendedAttribute('outputAverages')
        project.saveProtocol(protRansac)
        
        # 2b. Eman 
        protEmanInitVol = project.newProtocol(EmanProtInitModel)
        protEmanInitVol.setObjLabel('eman - initial vol')
        protEmanInitVol.inputSet.set(protImport)
        protEmanInitVol.inputSet.setExtendedAttribute('outputAverages')
        project.saveProtocol(protEmanInitVol)
        
        # 3. Join result volumes
        p1 = Pointer()
        p1.set(protRansac)
        p1.setExtendedAttribute('outputVolumes')
        
        p2 = Pointer()
        p2.set(protEmanInitVol)
        p2.setExtendedAttribute('outputVolumes')
        
        protJoin = project.newProtocol(ProtUnionSet)
        protJoin.setObjLabel('merge all volumes')
        protJoin.inputSets.append(p1)
        protJoin.inputSets.append(p2)
        project.saveProtocol(protJoin)
        
        
    return HttpResponse(mimetype='application/javascript')

def get_testdata(request):
    # Filename to use as test data 
    testDataKey = request.GET.get('testData')
    dsMDA = DataSet.getDataSet('initial_volume')
    fn = dsMDA.getFile(testDataKey)
    return HttpResponse(fn, mimetype='application/javascript')

def check_project_id(request):
    result = 0
    projectName = request.GET.get('code', None)
    
    try:
        manager = Manager()
        project = loadProject(projectName)
        result = 1
    except Exception:
        pass
    
    return HttpResponse(result, mimetype='application/javascript')
 

