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
from views_util import loadProject, loadProjectFromPath, getResourceCss, getResourceIcon, getResourceJs, \
    getServiceManager, PROJECT_NAME, SERVICE_NAME, \
    getVarFromRequest, CTX_PROJECT_PATH, CTX_PROJECT_NAME
from views_tree import loadProtTree
from pyworkflow.manager import Manager
from pyworkflow.utils.path import copyFile
from django.http import HttpResponse, HttpRequest
from django.shortcuts import render_to_response
import pyworkflow
from pyworkflow.em import getProtocols


def projects(request):
    from pyworkflow.utils.utils import prettyDate

    manager = Manager()
    projectsList = manager.listProjects()
    for p in projectsList:
        p.pTime = prettyDate(p.mTime)

    if CTX_PROJECT_NAME in request.session:
        request.session[CTX_PROJECT_NAME] = ""
    if CTX_PROJECT_PATH in request.session:
        request.session[CTX_PROJECT_PATH] = ""

    context = {'projects': projectsList,
               'projects_css': getResourceCss('projects'),
               'project_utils_js': getResourceJs('project_utils'),
               }

    context = base_grid(request, context)

    return render_to_response('projects.html', context)


def create_project(request):
    manager = Manager()

    if request.is_ajax():
        projectName = getVarFromRequest(request, PROJECT_NAME)
        manager.createProject(projectName, chdir=False)

    return HttpResponse(content_type='application/javascript')


def delete_project(request):
    manager = Manager()

    if request.is_ajax():
        projectName = getVarFromRequest(request, PROJECT_NAME)
        manager.deleteProject(projectName)

    return HttpResponse(content_type='application/javascript')


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
    project = loadProject(request)
    project_settings = project.getSettings()
    index = request.GET.get('index', None)

    # set the new protocol tree chosen
    views = project.getProtocolViews()
    project_settings.setProtocolView(views[int(index)])
    project_settings.write()

    return HttpResponse(content_type='application/javascript')


def update_graph_view(request):
    status = request.GET.get('status', None)
    project = loadProject(request)
    project_settings = project.getSettings()

    project_settings.runsView.set(int(status))

    project_settings.write()
    return HttpResponse(content_type='application/javascript')


def save_selection(request):
    if request.is_ajax():
        mark = request.GET.get('mark', None)

        project = loadProject(request)
        project_settings = project.getSettings()

        # Set the selected runs stored in BD    
        project_settings.runSelection.set(mark)

        project_settings.write()

    return HttpResponse(content_type='application/javascript')


def tree_prot_view(request):
    project = loadProject(request)

    # load the protocol tree current active
    htmlTree = loadProtTree(project)

    return render_to_response('project_content/tree_prot_view.html', {'protTreeHtml': htmlTree})


def run_table_graph(request):
    from pyworkflow.gui.tree import ProjectRunsTreeProvider

    try:
        project = loadProject(request)
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

            # Get the run mode view (0:list / 1:graph / 2:small graph) stored in BD
            runsView = project_settings.runsView.get()

            context = {'runs': runsNew,
                       'columns': provider.getColumns(),
                       'runsView': runsView,
                       'selectedRuns': selectedRuns}

            return render_to_response('project_content/run_table_graph.html', context)

        elif listNewElm:
            request.session['runs'] = runsNew
            jsonStr = json.dumps(listNewElm, ensure_ascii=False)
            return HttpResponse(jsonStr, content_type='application/json')

        else:
            return HttpResponse("ok", content_type='text/plain')

    except Exception:
        print "Stopped script"
        return HttpResponse("stop", content_type='text/plain')


def formatProvider(provider, mode):
    objs = []
    for obj in provider.getObjects():
        objInfo = provider.getObjectInfo(obj)

        objId = objInfo["key"]
        name = objInfo["text"]
        info = objInfo["values"]

        if mode == "runs":
            status = info[0]
            time = info[1]
            objs.append((objId, [objId, name, status, time]))

        elif mode == "objects":
            objs.append((objId, [objId, name, info]))

    return objs


def project_content(request):
    projectName = getVarFromRequest(request, PROJECT_NAME)
    project = Manager().loadProject(projectName, chdir=False)
    context = contentContext(request, project)
    context.update({'mode': None,
                    'formUrl': 'form'})
    return render_to_response('project_content/project_content.html', context)


def search_protocol(request):
    context = base_flex(request, {})
    return render_to_response('project_content/search_protocol.html', context)


def get_protocols(request):
    search = request.GET.get('search', None)
    result = []
    emProtocolsDict = getProtocols()
    for key, prot in emProtocolsDict.iteritems():
        label = prot.getClassLabel()
        className = prot.__name__
        if search is None or search in label:
            result.append((label, className))

    jsonStr = json.dumps(result, ensure_ascii=False)
    return HttpResponse(jsonStr, content_type='application/javascript')


def contentContext(request, project, serviceName=None):
    from pyworkflow.gui.tree import ProjectRunsTreeProvider

    projectName = project.getShortName()
    if serviceName is None:
        serviceName = getVarFromRequest(request, SERVICE_NAME)

    project_settings = project.getSettings()
    request.session[CTX_PROJECT_PATH] = project.getPath()
    request.session[PROJECT_NAME] = projectName

    provider = ProjectRunsTreeProvider(project)
    runs = formatProvider(provider, "runs")
    request.session['runs'] = runs

    # Get the selected runs stored in BD    
    selectedRuns = project_settings.runSelection

    # Get the run mode view (0:list / 1:graph / 2:small graph) stored in BD
    runsView = project_settings.runsView.get()

    # load the protocol tree current active
    htmlTree = loadProtTree(project, serviceName)

    # get the choices to load protocol trees
    choices = project.getProtocolViews()

    # get the choice current 
    choiceSelected = choices.index(project.getCurrentProtocolView().text.get())

    context = {CTX_PROJECT_NAME: projectName,
               'view': 'protocols',
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
               'project_content_css': getResourceCss('project_content'),
               'protTreeHtml': htmlTree,
               'choices': choices,
               'choiceSelected': choiceSelected,
               'runs': runs,
               'columns': provider.getColumns(),
               'runsView': runsView,
               'selectedRuns': selectedRuns,
               SERVICE_NAME: serviceName
               }

    context = base_flex(request, context)

    return context


def protocol_info(request):
    from pyworkflow.web.app.views_util import parseText
    from pyworkflow.em.data import EMObject

    jsonStr = ''

    if request.is_ajax():

        protId = request.GET.get('protocolId', None)

        project = loadProject(request)

        if len(protId) > 0:
            protocol = project.getProtocol(int(protId))

            # PROTOCOL IO
            input_obj = [{'name': name,
                          'nameId': attr.getNameId(),
                          'id': attr.getObjId(),
                          'info': str(attr.get()) if attr.isPointer() else str(attr)}
                         for name, attr in protocol.iterInputAttributes()]

            output_obj = [{'name': name,
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
            fOutString, fErrString, fScpnString = protocol.getLogsAsStrings(project.path)

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

    return HttpResponse(jsonStr, content_type='application/javascript')


def webservice_projects(request):
    if CTX_PROJECT_NAME in request.session: request.session[CTX_PROJECT_NAME] = ""
    if CTX_PROJECT_PATH in request.session: request.session[CTX_PROJECT_PATH] = ""
    context = {'projects_css': getResourceCss('projects'),
               }

    context = base_grid(request, context)
    return render_to_response('webservice_projects.html', context)
