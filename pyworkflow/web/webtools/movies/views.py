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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

from os.path import exists, join, basename
from pyworkflow.web.app.views_util import getResourceCss, getResourceJs, getResourceIcon, getServiceManager
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

def service_movies(request):

    if 'projectName' in request.session: request.session['projectName'] = ""
    if 'projectPath' in request.session: request.session['projectPath'] = ""

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
        {"tag": "url", "value": "/upload_movies/", "text":"Upload Data", "icon": "fa-upload.png"}]},
    {"tag": "section", "text": "2. Import your data", "children": [
        {"tag": "protocol", "value": "ProtImportMovies", "text": "Import Movies", "icon": "bookmark.png"}]},
    {"tag": "section", "text": "3. Align your Movies", "children": [
        {"tag": "protocol", "value": "ProtMovieAlignment", "text": "xmipp3 - movie alignment"}]}]
        ''')
        f.close()


def create_movies_project(request):
    
    if request.is_ajax():
        
        import os
        from pyworkflow.object import Pointer
        from pyworkflow.em.protocol import ProtImportMovies
        from pyworkflow.em.packages.xmipp3 import ProtMovieAlignment
        
        # Create a new project
        projectName = request.GET.get('projectName')
        
        # Filename to use as test data 
        testDataKey = request.GET.get('testData')
        
        manager = getServiceManager('movies')
        writeCustomMenu(manager.protocols)
        project = manager.createProject(projectName, runsView=1, 
                                        hostsConf=manager.hosts,
                                        protocolsConf=manager.protocols
                                        )   
        
        project.getSettings().setLifeTime(336) # 14 days * 24 hours
        project.saveSettings()
        
        
#         copyFile(customMenu, project.getPath('.config', 'protocols.conf'))
        
        # Create symbolic link for uploads
        projectPath = manager.getProjectPath(projectName)
        dest = os.path.join(projectPath,'Uploads')
        os.rmdir(dest)#in movies uploads is created as a link
        # @todo: this path to uploads dir should be configurable outside the code...
        source = "/services/scipion/data/uploads/"+ projectName
        pwutils.path.makePath(source)
        pwutils.createLink(source, dest)
        
        # 1. Import movies
        if testDataKey :
            attr = getAttrTestFile(testDataKey)
            path_test = attr['path']

            
            for f in os.listdir(path_test):           
                # Create a symbolic link for each file
                file_path = os.path.join(path_test, f)
                source_file = os.path.join(source, f)
                pwutils.createAbsLink(file_path, source_file)
            
            label_import = "import movies ("+ testDataKey +")" 
            protImport = project.newProtocol(ProtImportMovies, objLabel=label_import)

            protImport.filesPath.set(attr["filesPath"])
            protImport.voltage.set(attr['voltage'])
            protImport.sphericalAberration.set(attr['sphericalAberration'])
            protImport.amplitudeContrast.set(attr['amplitudeContrast'])
            protImport.magnification.set(attr['magnification'])
            protImport.samplingRate.set(attr['samplingRate'])
            
            project.launchProtocol(protImport, wait=True)
        else:
            protImport = project.newProtocol(ProtImportMovies, objLabel='import movies')
            project.saveProtocol(protImport)
        
        # 2. Movie Alignment 
        protMovAlign = project.newProtocol(ProtMovieAlignment)
        protMovAlign.setObjLabel('xmipp - movie alignment')
        protMovAlign.inputMovies.set(protImport)
        protMovAlign.inputMovies.setExtended('outputMovies')
        project.saveProtocol(protMovAlign)
        
    return HttpResponse(mimetype='application/javascript')

def getAttrTestFile(key):
    if(key == "ribosome"):
        riboDataset = DataSet.getDataSet('riboMovies')
        riboFiles = riboDataset.getFile("allMovies")
        attr = {#"filesPath" : "/services/scipion/data/scipionweb/movies_testdata/80S_ribosome/",
                "path": riboDataset.getPath(),
                "filesPath" : riboFiles,
                "voltage" : 300.0,
                "sphericalAberration" : 2.7,
                "amplitudeContrast" : 0.1,
                "magnification" : 59000,
                "samplingRate": 1.77}
    if(key == "falcon"):
        jmbFalconDataset = DataSet.getDataSet('jmbFalconMovies')
        jmbFalconFiles = jmbFalconDataset.getFile("allMovies")
        attr = {#"path" : "/services/scipion/data/scipionweb/movies_testdata/JMB_2015/",
                "path": jmbFalconDataset.getPath(),
                "filesPath" : jmbFalconFiles,
                "voltage" : 300.0,
                "sphericalAberration" : 2.7,
                "amplitudeContrast" : 0.1,
                "magnification" : 59000,
                "samplingRate": 1.34}
        
    return attr

 
def movies_content(request):
    
    domain = django_settings.SITE_URL
    projectName = request.GET.get('p', None)
    path_files = django_settings.ABSOLUTE_URL + '/resources_movies/img/'
    command = "rsync -av --port 3333 USER_FOLDER/ %s::mws/%s"%(domain, projectName)
    
    manager = getServiceManager('movies')
    project = manager.loadProject(projectName, 
                                  protocolsConf=manager.protocols,
                                  hostsConf=manager.hosts)
    daysLeft = prettyDelta(project.getLeftTime())
    
    context = contentContext(request, project)
    context.update({
                    # MODE
                    'formUrl': 'mov_form',
                    'mode':'service',
                    # IMAGES
                    'importMovies': path_files + 'importMovies.png',
                    'movieAlignment': path_files + 'movieAlignment.png',
                    'protMovieAlign': path_files + 'protMovieAlign.png',
                    'summary': path_files + 'summary.png',
                    'showj': path_files + 'showj.png',
                    'download': path_files + 'download.png',
                    'command' : command,
                    'daysLeft': daysLeft,
                    })
    
    return render_to_response('movies_content.html', context)


def movies_form(request):
    from django.shortcuts import render_to_response
    context = contextForm(request)
    context.update({'path_mode':'select',
                    'formUrl': 'mov_form',
                    'showHost': False,
                    'showParallel': False,
                    'hostSelected': 'localhost'})
    return render_to_response('form/form.html', context)


def upload_movies(request):
    domain = django_settings.SITE_URL
    projectName = request.session['projectName']
    
    command = "rsync -av --port 3333 USER_FOLDER/ %s::mws/%s"%(domain, projectName)

    context = {'command': command,
               'logo_scipion_small': getResourceIcon('logo_scipion_small'),
               }

    context = base_form(request, context)
    
    return render_to_response('upload_movies.html', context)

