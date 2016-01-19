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
from pyworkflow.web.app.views_util import (getResourceCss, getResourceJs, getResourceIcon,
                                           getServiceManager, getImageFullPath, loadProtocolConf, SERVICE_NAME,
                                           CTX_PROJECT_PATH, CTX_PROJECT_NAME, PROJECT_NAME, getVarFromRequest,
                                           getResource, getAbsoluteURL)
from pyworkflow.web.app.views_base import base_grid
from pyworkflow.web.app.views_project import contentContext
from pyworkflow.web.app.views_protocol import contextForm
from django.shortcuts import render_to_response
from pyworkflow.web.pages import settings as django_settings
from django.http import HttpResponse
from pyworkflow.tests.tests import DataSet
from pyworkflow.utils import copyFile
from pyworkflow.utils.utils import prettyDelta
from pyworkflow.em.packages.xmipp3 import XmippProtRansac
from pyworkflow.em.packages.eman2.protocol_initialmodel import EmanProtInitModel
from pyworkflow.em.packages.xmipp3.protocol_reconstruct_significant import XmippProtReconstructSignificant
from pyworkflow.object import Pointer
from pyworkflow.em.protocol import ProtImportAverages
from pyworkflow.em.packages.xmipp3 import XmippProtAlignVolumeForWeb

MYFIRSTMAP_SERVICE = 'myfirstmap'


def service_projects(request):
   
    if CTX_PROJECT_NAME in request.session: request.session[CTX_PROJECT_NAME] = ""
    if CTX_PROJECT_PATH in request.session: request.session[CTX_PROJECT_PATH] = ""

    myfirstmap_utils = getResource("js/myfirstmap_utils.js")

    context = {'projects_css': getResourceCss('projects'),
               'project_utils_js': getResourceJs('project_utils'),
               'scipion_mail': getResourceIcon('scipion_mail'),
               'myfirstmap_utils': myfirstmap_utils,
               'hiddenTreeProt': True,
               SERVICE_NAME: MYFIRSTMAP_SERVICE
               }
    
    context = base_grid(request, context)
    return render_to_response('service_projects.html', context)


def writeCustomMenu(customMenu):
    print
    if not exists(customMenu):
        f = open(customMenu, 'w+')
        f.write('''
[PROTOCOLS]

Initial_Volume = [
    {"tag": "section", "text": "1. Upload data", "children": [
        {"tag": "protocol", "value": "ProtImportAverages",     "text": "import averages", "icon": "bookmark.png"}]},
    {"tag": "section", "text": "2. Create a 3D volume", "children": [
        {"tag": "protocol", "value": "XmippProtRansac", "text": "xmipp3 - ransac"},
        {"tag": "protocol", "value": "EmanProtInitModel", "text": "eman2 - initial volume"},
        {"tag": "protocol", "value": "XmippProtReconstructSignificant", "text": "xmipp3 - significant"}]},
    {"tag": "section", "text": "3. Align volumes.", "children": [
        {"tag": "protocol", "value": "XmippProtAlignVolumeForWeb", "text": "xmipp3 - align volumes"}
        ]}]
        ''')
        f.close()
        

def create_service_project(request):
    if request.is_ajax():
        
        # Create a new project
        projectName = getVarFromRequest(request, PROJECT_NAME)
        
        # Filename to use as test data 
        testDataKey = request.GET.get('testData')
        
        manager = getServiceManager(MYFIRSTMAP_SERVICE)
        writeCustomMenu(manager.protocols)
        project = manager.createProject(projectName, runsView=1, 
                                        hostsConf=manager.hosts,
                                        protocolsConf=manager.protocols,
                                        chdir=False)
        
        project.getSettings().setLifeTime(336) # 14 days * 24 hours
        project.saveSettings()
        #copyFile(customMenu, project.getPath('.config', 'protocols.conf'))
        
        # 1. Import averages
        
        # If using test data execute the import averages run
        # options are set in 'project_utils.js'
        dsMDA = DataSet.getDataSet('initial_volume')
        
        if testDataKey :
            fn = dsMDA.getFile(testDataKey)
            newFn = getImageFullPath(project.path, join(project.uploadPath, basename(fn)))
            copyFile(fn, newFn)
            
            label_import = 'import averages ('+ testDataKey +')'
            protImport = project.newProtocol(ProtImportAverages, objLabel=label_import)
            protImport.filesPath.set(newFn)
            protImport.samplingRate.set(1.)
            project.launchProtocol(protImport, wait=True, chdir=False)
        else:
            protImport = project.newProtocol(ProtImportAverages, objLabel='import averages')
            project.saveProtocol(protImport)
        
        # 2a. Ransac 
        protRansac = project.newProtocol(XmippProtRansac)
        protRansac.setObjLabel('xmipp - ransac')
        protRansac.inputSet.set(protImport)
        protRansac.inputSet.setExtended('outputAverages')
        setProtocolParams(protRansac, testDataKey)
        project.saveProtocol(protRansac)
        
        # 2b. Eman 
        protEmanInitVol = project.newProtocol(EmanProtInitModel)
        protEmanInitVol.setObjLabel('eman - initial vol')
        protEmanInitVol.inputSet.set(protImport)
        protEmanInitVol.inputSet.setExtended('outputAverages')
        setProtocolParams(protEmanInitVol, testDataKey)
        project.saveProtocol(protEmanInitVol)
        
        # 2c. Significant 
        protSignificant = project.newProtocol(XmippProtReconstructSignificant)
        protSignificant.setObjLabel('xmipp - significant')
        protSignificant.inputSet.set(protImport)
        protSignificant.inputSet.setExtended('outputAverages')
        setProtocolParams(protSignificant, testDataKey)
        project.saveProtocol(protSignificant)
        
        # 3. Join result volumes
        p1 = Pointer()
        p1.set(protRansac)
        p1.setExtended('outputVolumes')
        
        p2 = Pointer()
        p2.set(protEmanInitVol)
        p2.setExtended('outputVolumes')
        
        p3 = Pointer()
        p3.set(protSignificant)
        p3.setExtended('outputVolume')
        
        protJoin = project.newProtocol(XmippProtAlignVolumeForWeb)
        protJoin.setObjLabel('align volumes')
        protJoin.inputVolumes.append(p1)
        protJoin.inputVolumes.append(p2)
        protJoin.inputVolumes.append(p3)
#         protJoin.inputVolumes.append(p4)
        project.saveProtocol(protJoin)
        
#         protValidate = project.newProtocol(XmippProtValidateNonTilt)
#         protValidate.setObjLabel('validate nontilt')
#         protValidate.inputVolumes.set(protJoin)
#         protValidate.inputVolumes.setExtended('outputVolumes')
#         protValidate.inputParticles.set(protImport)
#         protValidate.inputParticles.setExtended('outputAverages')
#         protValidate.numberOfThreads.set(8)
#         if testDataKey :
#             setProtocolParams(protValidate, testDataKey)
# #         protJoin.inputVolumes.append(p4)
#         project.saveProtocol(protValidate)
        
        
    return HttpResponse(mimetype='application/javascript')

def get_testdata(request):
    # Filename to use as test data 
    testDataKey = request.GET.get('testData')
    dsMDA = DataSet.getDataSet('initial_volume')
    fn = dsMDA.getFile(testDataKey)
    return HttpResponse(fn, mimetype='application/javascript')

 
def myfirstmap_form(request):
    from django.shortcuts import render_to_response
    context = contextForm(request)
    context.update({'path_mode':'upload',
                    'formUrl': 'my_form',
                    'showHost': False,
                    'showParallel': True})
    return render_to_response('form/form.html', context)

 
def service_content(request):

    projectName = request.GET.get('p', None)
    path_files = getAbsoluteURL('resources_myfirstmap/img/')
    
    # Get info about when the project was created
    manager = getServiceManager(MYFIRSTMAP_SERVICE)
    project = manager.loadProject(projectName, 
                                  protocolsConf=manager.protocols,
                                  hostsConf=manager.hosts,
                                  chdir=False)
    
    daysLeft = prettyDelta(project.getLeftTime())
    
    context = contentContext(request, project, serviceName=MYFIRSTMAP_SERVICE)
    context.update({'importAverages': path_files + 'importAverages.png',
                    'useProtocols': path_files + 'useProtocols.png',
                    'protForm': path_files + 'protForm.png',
                    'summary': path_files + 'summary.png',
                    'showj': path_files + 'showj.png',
                    'alignVol': path_files + 'alignVol.png',
                    'download': path_files + 'download.png',
                    'formUrl': 'my_form',
                    'mode':'service',
                    'daysLeft': daysLeft
                    })
    
    return render_to_response('service_content.html', context)

def setProtocolParams(protocol, key):
    #Here we set protocol parameters for each test data
    if key:
        cls = type(protocol)
        
        if issubclass(cls, XmippProtRansac):
            if(key == "bpv"):
                attrs = {"symmetryGroup" : "i1",
                        "dimRed": True,
                        "numGrids": 2}
            if(key == "groel"):
                attrs = {"symmetryGroup" : "d7",
                        "dimRed": True, 
                        "numGrids": 3}
            if(key == "ribosome"):
                attrs = {"symmetryGroup" : "c1",
                        "dimRed": True,
                        "numGrids": 3}
        
        elif issubclass(cls, EmanProtInitModel):
            if(key == "bpv"):
                attrs = {"symmetry" : "icos"}
            if(key == "groel"):
                attrs = {"symmetry" : "d7"}
            if(key == "ribosome"):
                attrs = {"symmetry" : "c1"}
        
        elif issubclass(cls, XmippProtReconstructSignificant):
            if(key == "bpv"):
                attrs = {"symmetryGroup" : "i1",
                         "alpha0": 99.5, 
                         "alphaF": 99.5}
            if(key == "groel"):
                attrs = {"symmetryGroup" : "d7",
                         "alpha0": 98, 
                         "alphaF": 99.9}
            if(key == "ribosome"):
                attrs = {"symmetryGroup" : "c1", 
                         "alpha0": 80, 
                         "alphaF": 99.5}
        # if issubclass(cls, XmippProtValidateNonTilt):
        #     if(key == "bpv"):
        #         attrs = {"symmetryGroup" : "i1"}
        #     if(key == "groel"):
        #         attrs = {"symmetryGroup" : "d7"}
        #     if(key == "ribosome"):
        #         attrs = {"symmetryGroup" : "c1"}
        for key,value in attrs.iteritems():
            getattr(protocol, key).set(value)
    
    loadProtocolConf(protocol)
    