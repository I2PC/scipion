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
import urlparse
from os.path import exists, join, basename
from pyworkflow.web.app.views_util import getResourceCss, getResourceJs, getResourceIcon, getServiceManager, \
    loadProtocolConf, SERVICE_NAME, getVarFromRequest, PROJECT_NAME, CTX_PROJECT_PATH, CTX_PROJECT_NAME
from pyworkflow.web.app.views_base import base_grid, base_flex, base_form
from pyworkflow.web.app.views_project import contentContext
from pyworkflow.web.app.views_protocol import contextForm
from django.shortcuts import render_to_response
from pyworkflow.web.pages import settings as django_settings
from pyworkflow.manager import Manager
from django.http import HttpResponse
from pyworkflow.tests.tests import DataSet
from pyworkflow.utils import copyFile
import pyworkflow.utils as pwutils
from pyworkflow.utils.utils import prettyDelta
from django.contrib.sites.models import Site

MOVIES_SERVICE = 'movies'


def service_movies(request):
    if CTX_PROJECT_NAME in request.session: request.session[CTX_PROJECT_NAME] = ""
    if CTX_PROJECT_PATH in request.session: request.session[CTX_PROJECT_PATH] = ""

    movies_utils = join(django_settings.STATIC_URL, "js/", "movies_utils.js")

    context = {'projects_css': getResourceCss('projects'),
               'project_utils_js': getResourceJs('project_utils'),
               'scipion_mail': getResourceIcon('scipion_mail'),
               'movies_utils': movies_utils,
               }

    context = base_grid(request, context)
    return render_to_response('movies_projects.html', context)


def writeCustomMenu(customMenu):
    if not exists(customMenu):
        f = open(customMenu, 'w+')
        f.write('''
[PROTOCOLS]

Movies_Alignment = [
    {"tag": "section", "text": "1. Upload data", "children": [
        {"tag": "url", "value": "/upload_movies/", "text":"upload data", "icon": "fa-upload.png"}]},
    {"tag": "section", "text": "2. Import your data", "children": [
        {"tag": "protocol", "value": "ProtImportMovies", "text": "import Movies", "icon": "bookmark.png"}]},
    {"tag": "section", "text": "3. Align your Movies", "children": [
        {"tag": "protocol", "value": "ProtMovieAlignment", "text": "xmipp3 - movie alignment"}]}]
        ''')
        f.close()


def create_movies_project(request):
    if request.is_ajax():

        import os
        from pyworkflow.em.protocol import ProtImportMovies
        from pyworkflow.em.packages.xmipp3 import ProtMovieAlignment

        # Create a new project
        projectName = getVarFromRequest(request, PROJECT_NAME)

        # Filename to use as test data 
        testDataKey = request.GET.get('testData')

        manager = getServiceManager(MOVIES_SERVICE)
        writeCustomMenu(manager.protocols)
        project = manager.createProject(projectName, runsView=1,
                                        hostsConf=manager.hosts,
                                        protocolsConf=manager.protocols,
                                        chdir=False)

        project.getSettings().setLifeTime(336)  # 14 days * 24 hours
        project.saveSettings()

        #         copyFile(customMenu, project.getPath('.config', 'protocols.conf'))

        # Create symbolic link for uploads
        projectPath = manager.getProjectPath(projectName)
        dest = os.path.join(projectPath, 'Uploads')
        os.rmdir(dest)  # in movies uploads is created as a link
        # @todo: this path to uploads dir should be configurable outside the code...
        source = "/services/scipion/data/uploads/" + projectName
        pwutils.path.makePath(source)
        pwutils.createLink(source, dest)

        # 1. Import movies
        if testDataKey:
            attr = getAttrTestFile(testDataKey)
            path_test = attr['path']

            for f in os.listdir(path_test):
                # Create a symbolic link for each file
                file_path = os.path.join(path_test, f)
                source_file = os.path.join(source, f)
                pwutils.createAbsLink(file_path, source_file)

            label_import = "import movies (" + testDataKey + ")"
            protImport = project.newProtocol(ProtImportMovies, objLabel=label_import)

            protImport.filesPath.set(attr["filesPath"])
            protImport.voltage.set(attr['voltage'])
            protImport.sphericalAberration.set(attr['sphericalAberration'])
            protImport.amplitudeContrast.set(attr['amplitudeContrast'])
            protImport.magnification.set(attr['magnification'])
            protImport.samplingRate.set(attr['samplingRate'])

            project.launchProtocol(protImport, wait=True, chdir=False)
        else:
            protImport = project.newProtocol(ProtImportMovies, objLabel='import movies')
            project.saveProtocol(protImport)

        # 2. Movie Alignment 
        protMovAlign = project.newProtocol(ProtMovieAlignment)
        protMovAlign.setObjLabel('xmipp - movie alignment')
        protMovAlign.inputMovies.set(protImport)
        protMovAlign.inputMovies.setExtended('outputMovies')
        loadProtocolConf(protMovAlign)
        project.saveProtocol(protMovAlign)

    return HttpResponse(mimetype='application/javascript')


def getAttrTestFile(key):
    if key == "ribosome":
        riboDataset = DataSet.getDataSet('riboMovies')
        riboFiles = riboDataset.getFile("allMovies")
        attr = {"path": riboDataset.getPath(),
                "filesPath": riboFiles,
                "voltage": 300.0,
                "sphericalAberration": 2.7,
                "amplitudeContrast": 0.1,
                "magnification": 59000,
                "samplingRate": 1.77}
    if key == "falcon":
        jmbFalconDataset = DataSet.getDataSet('jmbFalconMovies')
        jmbFalconFiles = jmbFalconDataset.getFile("allMovies")
        attr = {"path": jmbFalconDataset.getPath(),
                "filesPath": jmbFalconFiles,
                "voltage": 300.0,
                "sphericalAberration": 2.7,
                "amplitudeContrast": 0.1,
                "magnification": 59000,
                "samplingRate": 1.34}

    return attr


def movies_content(request):
    projectName = request.GET.get('p', None)
    path_files = django_settings.ABSOLUTE_URL + '/resources_movies/img/'
    command = getSyncCommand(request)

    manager = getServiceManager(MOVIES_SERVICE)
    project = manager.loadProject(projectName,
                                  protocolsConf=manager.protocols,
                                  hostsConf=manager.hosts,
                                  chdir=False)
    daysLeft = prettyDelta(project.getLeftTime())

    context = contentContext(request, project, serviceName=MOVIES_SERVICE)
    context.update({
        # MODE
        'formUrl': 'mov_form',
        'mode': 'service',
        # IMAGES
        'importMovies': path_files + 'importMovies.png',
        'movieAlignment': path_files + 'movieAlignment.png',
        'protMovieAlign': path_files + 'protMovieAlign.png',
        'summary': path_files + 'summary.png',
        'showj': path_files + 'showj.png',
        'download': path_files + 'download.png',
        'command': command,
        'daysLeft': daysLeft
    })

    return render_to_response('movies_content.html', context)


def movies_form(request):
    from django.shortcuts import render_to_response
    context = contextForm(request)
    context.update({'path_mode': 'select',
                    'formUrl': 'mov_form',
                    'showHost': False,
                    'showParallel': True,
                    'hostSelected': 'localhost'})
    return render_to_response('form/form.html', context)


def upload_movies(request):
    command = getSyncCommand(request)
    context = {'command': command,
               'logo_scipion_small': getResourceIcon('logo_scipion_small'),
               }

    context = base_form(request, context)

    return render_to_response('upload_movies.html', context)


def getSyncCommand(request):
    domain = django_settings.SITE_URL
    domain = domain.split(":")[0]
    projectName = getVarFromRequest(request, PROJECT_NAME)
    command = "rsync -av --port 3333 USER_FOLDER/ %s::mws/%s" % (domain, projectName)
    return command
