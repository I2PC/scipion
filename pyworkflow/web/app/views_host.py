import xmipp
from django.template import RequestContext
from forms import HostForm
from pyworkflow.hosts import HostMapper
from pyworkflow.web.app.views_util import *
from django.http import HttpResponseRedirect

######    Hosts template    #####
def getScipionHosts():
    from pyworkflow.apps.config import getSettingsPath
    defaultHosts = getSettingsPath()
    return HostMapper(defaultHosts).selectAll()

def viewHosts(request):      
    projectName = request.session['projectName']    
    project = loadProject(projectName)
    projectHosts = project.getSettings().getHosts()   
    scpnHostsChoices = []
    scpnHostsChoices.append(('', ''))

    message = request.GET.get("message")
    context = {'projectName' : projectName,
               'editTool': getResourceIcon('edit_toolbar'),
               'newTool': getResourceIcon('new_toolbar'),
               'deleteTool': getResourceIcon('delete_toolbar'),
               'browseTool': getResourceIcon('browse_toolbar'),
               'hosts': projectHosts,
               'jquery': getResourceJs('jquery'),
               'utils': getResourceJs('utils'),
               'host_utils': getResourceJs('host_utils'),
               'css':getResourceCss('general'),
               'message': message,
               'view': 'hosts',
               'contentConfig': 'full'}    
      
    return render_to_response('hosts.html', context)

# def getHost(request):
#     from django.http import HttpResponse
#     import json
#     from django.utils import simplejson
#     
#     if request.is_ajax():
#         hostLabel = request.GET.get('hostLabel')
#         projectName = request.session['projectName']
#         project = loadProject(projectName)
#         hostsMapper = HostMapper(project.settingsPath)
#         HostConfig = hostsMapper.selectByLabel(hostLabel)
#         jsonStr = json.dumps({'host':HostConfig.getObjDict()})
#         return HttpResponse(jsonStr, mimetype='application/javascript')


def getHostFormContext(request, host=None, initialContext=None):
#     scpnHostsChoices = []
#     scpnHostsChoices.append(('', ''))
#     scipionHosts = getScipionHosts()
#     for hostConfig in scipionHosts:
#         scpnHostsChoices.append((hostConfig.getLabel(), hostConfig.getHostName()))
#     form.fields['scpnHosts'].choices = scpnHostsChoices        
    # We check if we are going to edit a host
    tittle = None
    if host is not None:
        form = initialContext['form']
        tittle = host.getLabel() + " host configuration"
    else:
        form = HostForm()
        tittle = "New host configuration"  
            
    context = {'tittle': tittle,
               'jquery': getResourceJs('jquery'),
               'utils': getResourceJs('utils'),
               'css':getResourceCss('general'),
               'form': form}
     
    if initialContext is not None:
        context.update(initialContext)
        
    return context

def hostForm(request):
    hostId = request.GET.get('hostId')
    context = None
    hostConfig = None
    if hostId != None and hostId != '':
        projectName = request.session['projectName']
        project = loadProject(projectName)
        hostsMapper = HostMapper(project.settingsPath)
        hostConfig = hostsMapper.selectById(hostId)
        form = HostForm()
        form.setFormHost(hostConfig)
        context = {'form': form}
    return render_to_response('host_form.html', RequestContext(request, getHostFormContext(request, hostConfig, context)))  # Form Django forms

def updateHostsConfig(request):    
    form = HostForm(request.POST, queueSystemConfCont=request.POST.get('queueSystemConfigCount'), queueConfCont=request.POST.get('queueConfigCount'))  # A form bound to the POST data
    context = {'form': form}
    if form.is_valid():  # All validation rules pass
        projectName = request.session['projectName']
        project = loadProject(projectName)
        hostId = request.POST.get('objId')
        hostsMapper = HostMapper(project.settingsPath)
        hostConfig = hostsMapper.selectById(hostId)
        form.host = hostConfig
        host = form.getFormHost()
        savedHost = project.getSettings().saveHost(host)
        form.setFormHost(savedHost)
        return render_to_response('host_form.html', RequestContext(request, getHostFormContext(request, host, context))) 
    else:   
        return render_to_response('host_form.html', RequestContext(request, getHostFormContext(request, None, context)))  # Form Django forms

def deleteHost(request):
    hostId = request.GET.get("hostId")    
    projectName = request.session['projectName']
    project = loadProject(projectName)
    project.deleteHost(hostId)
#     context = {'message': "Host succesfully deleted"}
    return HttpResponseRedirect('/view_hosts')#, RequestContext(request))
