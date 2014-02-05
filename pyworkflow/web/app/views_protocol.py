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

import os
import xmipp
import json
from pyworkflow.viewer import WEB_DJANGO
from pyworkflow.em import *
from views_base import * 
from views_util import * 
from pyworkflow.manager import Manager
from pyworkflow.project import Project
from django.core.context_processors import csrf
from django.shortcuts import render_to_response
from django.http import HttpResponse

SPECIAL_PARAMS = ['numberOfMpi', 'numberOfThreads', 'hostName', 'expertLevel', '_useQueue']
OBJ_PARAMS =['runName', 'comment']

def getPointerHtml(protVar):
    if protVar.hasValue():
        return protVar.get().getNameId(), protVar.get().getObjId()
    return '',''
            
def form(request):
    project, protocol = loadProtocolProject(request, requestType='GET')
    action = request.GET.get('action', None)
    paramProt = request.GET.get('paramProt', None)
    protRunIdViewer = request.GET.get('protRunIdViewer', None)
    hosts = [host.getLabel() for host in project.getSettings().getHosts()]
    
    visualize = 0
    
    if action == 'visualize':
        visualize = 1
        viewerDict = protocol.getVisualizeDictWeb()
    elif action == 'protSimple':
        visualize = 2
    elif action == 'copy':
        protocol = project.copyProtocol(protocol)
        
    package = protocol._package
    logoPath = ''
    path = getattr(package, '_logo', '')
    if path != '': 
        logoPath = getResourceLogo(path)
            
    wizards = findWizards(protocol, WEB_DJANGO)
    
    protocol.htmlCitations = parseText(protocol.citations())
    protocol.htmlDoc = parseText(protocol.getDoc())
    
    for section in protocol._definition.iterSections():
        for paramName, param in section.iterParams():
            protVar = getattr(protocol, paramName, None)
            if protVar is None:
                raise Exception("_fillSection: param '%s' not found in protocol" % paramName)
                # Create the label
            
            if isinstance(param, MultiPointerParam):
                htmlValueList = []
                htmlIdValueList = []
                for pointer in protVar:
                    htmlValue, htmlIdValue = getPointerHtml(pointer)
                    htmlValueList.append(htmlValue)
                    htmlIdValueList.append(htmlIdValue)
                    
                param.htmlValueIdZip = zip(htmlValueList,htmlIdValueList)    
            elif isinstance(param, PointerParam):
                param.htmlValue, param.htmlIdValue = getPointerHtml(protVar)
            else:
                param.htmlValue = protVar.get(param.default.get(""))
                if isinstance(protVar, Boolean):
                    if param.htmlValue:
                        param.htmlValue = 'true'
                    else:
                        param.htmlValue = 'false'
            
            if paramName in wizards:
                param.hasWizard = True
                param.wizardName = wizards[paramName].getView()
#                print "param: ", paramName, " has wizard", " view: "
            
            if visualize == 1:
                if paramName in viewerDict:
                    param.hasViewer = True
#                    print "param: ", paramName, " has viewer"
                
            param.htmlCond = param.condition.get()
            param.htmlDepend = ','.join(param._dependants)
            param.htmlCondParams = ','.join(param._conditionParams)
            
            if not param.help.empty():
                param.htmlHelp = parseText(param.help.get())
#            param.htmlExpertLevel = param.expertLevel.get()   

            """ Workflow Addon """
            valueURL = request.GET.get(paramName, None)            
            if valueURL is not None:
                if valueURL is not param.htmlValue:
                    param.htmlValue = valueURL
    
    context = {'projectName':project.getName(),
               'package_logo': logoPath,
               'protocol':protocol,
               'protRunIdViewer':protRunIdViewer,
               'definition':protocol._definition,
               'visualize':visualize,
               'paramProt':paramProt,
               'favicon': getResourceIcon('favicon'),
               'comment':getResourceIcon('edit_toolbar'),
               'jquery_ui_css': getResourceCss('jquery_ui'),
               'wizard': getResourceIcon('wizard'),
               'protocol_form_utils': getResourceJs('protocol_form_utils'),
               'hosts':hosts,
               'comment': parseText(protocol.getObjComment())
               }
    # Update the context dictionary with the special params
    for paramName in SPECIAL_PARAMS:
        context[paramName] = protocol.getAttributeValue(paramName, '')
    
    # Cross Site Request Forgery protection is need it
    context.update(csrf(request))
    
    context = base_form(request, context)
    
    return render_to_response('form.html', context)
    
# Method to launch a protocol #
def protocol(request):
    project, protocol = loadProtocolProject(request)
    updateProtocolParams(request, protocol, project)    
    errors = protocol.validate()

    if len(errors) == 0:
        # No errors, launch the protocol
        try:
            project.launchProtocol(protocol)
            
        except Exception, ex:
            errors = [convertTktoHtml(str(ex))]
            
    jsonStr = json.dumps({'errors' : parseText(errors)}, ensure_ascii=False)
    
    return HttpResponse(jsonStr, mimetype='application/javascript')   

def updateProtocolParams(request, protocol, project):
    """ Update the protocol values from the Web-form.
    This function will be used from save_protocol and execute_protocol.
    """
    # PARAMS
    for paramName, _ in protocol.iterDefinitionAttributes():
        updateParam(request, project, protocol, paramName)
    
    # SPECIAL_PARAMS
    for paramName in SPECIAL_PARAMS:
        updateParam(request, project, protocol, paramName)
        
    # OBJ_PARAMS
    protocol.setObjLabel(request.POST.get('runName'))
    protocol.setObjComment(request.POST.get('comment'))

def getPointerValue(project, attr, value, paramName):
    if len(value.strip()) > 0:
      
        obj = project.mapper.selectById(value)  # Get the object from its id
        id1 = attr.getObjId()
        id2 = obj.getObjId()
        
        if id1 == id2 :
            raise Exception("Param: %s is autoreferencing with id: %d" % (paramName, value))
    else:
        obj = None
        
    return obj
 
def updateParam(request, project, protocol, paramName):
    """
    Params:
        request: current request handler
        project: current working project
        protocol: current protocol
        paramName: name of the attribute to be set in the protocol
            from the web form
    """
    attr = getattr(protocol, paramName)
        
    if isinstance(attr, PointerList):
        attr.clear()
        valueList = request.POST.getlist(paramName)
        for v in valueList:
            attr.append(Pointer(value=getPointerValue(project, attr, v, paramName)))
        value = None
            
    else:
        value = request.POST.get(paramName)
        if isinstance(attr, Pointer):
            value = getPointerValue(project, attr, value, paramName)
            
#    print "setting attr %s with value:" % paramName, value 
    attr.set(value)
        
def save_protocol(request):
    project, protocol = loadProtocolProject(request)
    updateProtocolParams(request, protocol, project)
    try:
        project.saveProtocol(protocol)
        protId = protocol.getObjId()
        res = {'success' : protId}
    except Exception, ex:
        errors = [convertTktoHtml(str(ex))]
        res = {'errors' : errors}
    
    jsonStr = json.dumps(res, ensure_ascii=False)
    return HttpResponse(jsonStr, mimetype='application/javascript')

# Method to delete a protocol #
def delete_protocol(request):
    # Project Id(or Name) should be stored in SESSION
    if request.is_ajax():
        projectName = request.session['projectName']
        project = loadProject(projectName)
        protId = request.GET.get('protocolId', None)
        protocol = project.mapper.selectById(int(protId))
        try:
            project.deleteProtocol(protocol)
            res = {'success' : 'Protocol deleted successful'}
        except Exception, ex:
            errors = [convertTktoHtml(str(ex))]
            res = {'errors' : errors}
            
    jsonStr = json.dumps(res, ensure_ascii=False)
    return HttpResponse(jsonStr, mimetype='application/javascript')   

# Method to stop a protocol #
def stop_protocol(request):
    # Project Id(or Name) should be stored in SESSION
    if request.is_ajax():
        projectName = request.session['projectName']
        project = loadProject(projectName)
        protId = request.GET.get('protocolId', None)
        protocol = project.mapper.selectById(int(protId))
     
        project.stopProtocol(protocol)
        
    return HttpResponse(mimetype='application/javascript') 
