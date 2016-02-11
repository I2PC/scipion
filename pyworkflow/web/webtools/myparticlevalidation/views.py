# **************************************************************************
# *
# * Authors:    Pablo Conesa (pconesa@cnb.csic.es)
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
from os.path import exists, join, basename, dirname

from django.http import HttpResponse
from django.shortcuts import render_to_response

import pyworkflow.utils as pwutils
from pyworkflow.em.packages.xmipp3.protocol_validate_overfitting import XmippProtValidateOverfitting
from pyworkflow.tests.tests import DataSet
from pyworkflow.utils.utils import prettyDelta
from pyworkflow.web.app.views_base import base_grid
from pyworkflow.web.app.views_project import contentContext
from pyworkflow.web.app.views_protocol import contextForm
from pyworkflow.web.app.views_util import ( getResourceCss, getResourceJs, getResourceIcon, getServiceManager,
                                            loadProtocolConf, CTX_PROJECT_PATH, CTX_PROJECT_NAME, PROJECT_NAME,
                                            getResource, getAbsoluteURL)
from pyworkflow.web.pages import settings as django_settings

MYPVAL_SERVICE = 'mypval'
MYPVAL_FORM_URL = 'p_form'


def particlevalidation_projects(request):
    if CTX_PROJECT_NAME in request.session: request.session[CTX_PROJECT_NAME] = ""
    if CTX_PROJECT_PATH in request.session: request.session[CTX_PROJECT_PATH] = ""

    mypval_utils = getResource("js/mypval_utils.js")

    context = {'projects_css': getResourceCss('projects'),
               'project_utils_js': getResourceJs('project_utils'),
               'scipion_mail': getResourceIcon('scipion_mail'),
               'mypval_utils': mypval_utils,
               'hiddenTreeProt': True,
               }

    context = base_grid(request, context)
    return render_to_response('pval_projects.html', context)


def writeCustomMenu(customMenu):
    if not exists(customMenu):
        # Make the parent path if it doesn't exists
        pwutils.makePath(dirname(customMenu))

        f = open(customMenu, 'w+')
        f.write('''
[PROTOCOLS]

Local_Resolution = [
    {"tag": "section", "text": "2. Import your data", "children": [
        {"tag": "protocol", "value": "ProtImportVolumes", "text": "import volumes", "icon": "bookmark.png"},
        {"tag": "protocol", "value": "ProtImportParticles", "text": "import particles", "icon": "bookmark.png"}]
    },
    {"tag": "section", "text": "3. Validation", "children": [
        {"tag": "protocol", "value": "XmippProtValidateOverfitting", "text": "xmipp3 - validate overfitting"}]
    }]
    ''')
        f.close()


def create_particlevalidation_project(request):
    if request.is_ajax():
        from pyworkflow.em.protocol import ProtImportVolumes
        from pyworkflow.em.protocol import ProtImportParticles
        from pyworkflow.em.packages.resmap.protocol_resmap import ProtResMap

        # Create a new project
        projectName = request.GET.get(PROJECT_NAME)

        # Filename to use as test data 
        testDataKey = request.GET.get('testData')

        manager = getServiceManager(MYPVAL_SERVICE)
        writeCustomMenu(manager.protocols)
        project = manager.createProject(projectName, runsView=1,
                                        hostsConf=manager.hosts,
                                        protocolsConf=manager.protocols,
                                        chdir=False
                                        )

        project.getSettings().setLifeTime(336)  # 14 days * 24 hours
        project.saveSettings()

        projectPath = manager.getProjectPath(projectName)

        # If we need to import test data...
        if testDataKey:

            # Get test data attributes
            attr = getAttrTestFile(testDataKey)

            # 1. Import volumes
            source = attr['volume']
            dest = os.path.join(projectPath, 'Uploads', basename(source))
            pwutils.createLink(source, dest)

            label_import = "import volumes (" + testDataKey + ")"
            protImportVol = project.newProtocol(ProtImportVolumes, objLabel=label_import)

            protImportVol.filesPath.set(dest)
            protImportVol.samplingRate.set(attr['samplingRate'])
            project.launchProtocol(protImportVol, wait=True, chdir=False)

            # 2. Import particles
            source = attr['particles']
            dest = os.path.join(projectPath, 'Uploads', basename(source))
            pwutils.createLink(source, dest)

            label_import = "import particles (" + testDataKey + ")"
            protImportParticles = project.newProtocol(ProtImportParticles, objLabel=label_import)

            protImportParticles.filesPath.set(dest)

            # Set import particle attributes
            protImportParticles.voltage.set(attr["microscopeVoltage"])
            protImportParticles.sphericalAberration.set(attr["sphericalAberration"])
            protImportParticles.amplitudeContrast.set(attr["amplitudeContrast"])
            protImportParticles.magnification.set(attr["magnificationRate"])
            protImportParticles.samplingRate.set(attr["particlesSamplingRate"])

            project.launchProtocol(protImportParticles, wait=True, chdir=False)

        else:

            # Empty import volumes protocol
            protImportVol = project.newProtocol(ProtImportVolumes, objLabel='import volumes')
            project.saveProtocol(protImportVol)

            # Empty import particles protocol
            protImportParticles = project.newProtocol(ProtImportParticles, objLabel='import particles')
            project.saveProtocol(protImportParticles)

        # 3. Validation
        protValidation = project.newProtocol(XmippProtValidateOverfitting)
        protValidation.setObjLabel('xmipp3 - validate overfitting')

        # Input volumes
        protValidation.input3DReference.set(protImportVol)
        protValidation.input3DReference.setExtended('outputVolume')

        # Input particles
        protValidation.inputParticles.set(protImportParticles)
        protValidation.inputParticles.setExtended('outputParticles')

        # Load additional configuration
        loadProtocolConf(protValidation)
        project.saveProtocol(protValidation)

    return HttpResponse(mimetype='application/javascript')


def getAttrTestFile(key):
    pval = DataSet.getDataSet('pval')

    if key == "pval":
        attr = {"volume": pval.getFile("pval_vol"),
                "samplingRate": 3.54,
                "particles": pval.getFile("pval_part"),
                "microscopeVoltage": 300,
                "sphericalAberration": 2,
                "amplitudeContrast": 0.1,
                "magnificationRate": 50000,
                "particlesSamplingRate": 3.54
                }

    return attr


def particlevalidation_form(request):
    from django.shortcuts import render_to_response
    context = contextForm(request)
    context.update({'path_mode': 'upload',
                    'formUrl': MYPVAL_FORM_URL,
                    'showHost': False,
                    'showParallel': True})
    return render_to_response('form/form.html', context)


def particlevalidation_content(request):
    projectName = request.GET.get('p', None)
    path_files = getAbsoluteURL('resources_mypval/img/')

    # Get info about when the project was created
    manager = getServiceManager(MYPVAL_SERVICE)
    project = manager.loadProject(projectName,
                                  protocolsConf=manager.protocols,
                                  hostsConf=manager.hosts,
                                  chdir=False)

    daysLeft = prettyDelta(project.getLeftTime())

    context = contentContext(request, project, serviceName=MYPVAL_SERVICE)

    # Resources for the help - guide, to be done.
    context.update({'formUrl': MYPVAL_FORM_URL,
                    'mode': 'service',
                    'daysLeft': daysLeft
                    })

    return render_to_response('pval_content.html', context)
