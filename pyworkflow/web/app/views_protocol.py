import os
import xmipp
import json
from pyworkflow.viewer import WEB_DJANGO
from pyworkflow.em import *
from views_util import * 
from pyworkflow.manager import Manager
from pyworkflow.project import Project
from django.core.context_processors import csrf
from django.shortcuts import render_to_response
from django.http import HttpResponse

def form(request):
    project, protocol = loadProtocolProject(request, requestType='GET')
    action = request.GET.get('action', None)
    hosts = [host.getLabel() for host in project.getSettings().getHosts()]
    
    if action == 'copy':
        protocol = project.copyProtocol(protocol)
        
    wizards = findWizards(protocol.getDefinition(), WEB_DJANGO)
    
    # TODO: Add error page validation when protocol is None
    for section in protocol._definition.iterSections():
        for paramName, param in section.iterParams():
            protVar = getattr(protocol, paramName, None)
            if protVar is None:
                raise Exception("_fillSection: param '%s' not found in protocol" % paramName)
                # Create the label
            if protVar.isPointer():
                if protVar.hasValue():
                    param.htmlValue = protVar.get().getNameId()
                else:
                    param.htmlValue = ""
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
                print "param: ", paramName, " has wizard", " view: "
                
                
            param.htmlCond = param.condition.get()
            param.htmlDepend = ','.join(param._dependants)
            param.htmlCondParams = ','.join(param._conditionParams)
#            param.htmlExpertLevel = param.expertLevel.get()   
    
    context = {'projectName':project.getName(),
               'protocol':protocol,
               'definition': protocol._definition,
               'favicon': getResourceIcon('favicon'),
               'help': getResourceIcon('help'),
               'comment':getResourceIcon('edit_toolbar'),
               'jquery': getResourceJs('jquery'),
               'jquery_ui': getResourceJs('jquery_ui'),
               'jquery_ui_css': getResourceCss('jquery_ui'),
               'browse': getResourceIcon('browse'),
               'wizard': getResourceIcon('wizard'),
               'utils': getResourceJs('utils'),
               'protocol_form_utils': getResourceJs('protocol_form_utils'),
               'css':getResourceCss('form'),
               'messi': getResourceJs('messi'),
               'messi_css': getResourceCss('messi'),
               'hosts':hosts
               }
    # Update the context dictionary with the special params
    for paramName in SPECIAL_PARAMS:
        context[paramName] = protocol.getAttributeValue(paramName, '')
    
    # Cross Site Request Forgery protection is need it
    context.update(csrf(request))
    
    return render_to_response('form.html', context)
    
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
    value = request.POST.get(paramName)
    if attr.isPointer():
        if len(value.strip()) > 0:
            objId = int(value.split('.')[-1])  # Get the id string for last part after .
            value = project.mapper.selectById(objId)  # Get the object from its id
            if attr.getObjId() == value.getObjId():
                raise Exception("Param: %s is autoreferencing with id: %d" % (paramName, objId))
        else:
            value = None
    attr.set(value)
    print "setting attr %s with value:" % paramName, value 
    
SPECIAL_PARAMS = ['runName', 'numberOfMpi', 'numberOfThreads', 'hostName', 'expertLevel', '_useQueue']

def updateProtocolParams(request, protocol, project):
    """ Update the protocol values from the Web-form.
    This function will be used from save_protocol and execute_protocol.
    """
    for paramName, _ in protocol.iterDefinitionAttributes():
        updateParam(request, project, protocol, paramName)
        
    for paramName in SPECIAL_PARAMS:
        updateParam(request, project, protocol, paramName)
        
def save_protocol(request):
    project, protocol = loadProtocolProject(request)
    updateProtocolParams(request, protocol, project)
    project.saveProtocol(protocol)
    
    return HttpResponse(mimetype='application/javascript')

# Method to launch a protocol #
def protocol(request):
    project, protocol = loadProtocolProject(request)
    updateProtocolParams(request, protocol, project)    
    errors = protocol.validate()

    if len(errors) == 0:
        # No errors 
        # Finally, launch the protocol
        project.launchProtocol(protocol)
    jsonStr = json.dumps({'errors' : errors}, ensure_ascii=False)
    
    return HttpResponse(jsonStr, mimetype='application/javascript')   

def delete_protocol(request):
    # Project Id(or Name) should be stored in SESSION
    if request.is_ajax():
        projectName = request.GET.get('projectName')
        project = loadProject(projectName)
        protId = request.GET.get('protocolId', None)
        protocol = project.mapper.selectById(int(protId))
     
        project.deleteProtocol(protocol)         
        
    return HttpResponse(mimetype='application/javascript') 
