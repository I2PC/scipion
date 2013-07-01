'''
Created on Jun 7, 2013

@author: antonio
'''
from django import forms
from pyworkflow.hosts import HostConfig

class HostForm(forms.Form):
#     scpnHosts = forms.ChoiceField(label='Scipion hosts', widget = forms.Select(), required = False,)
#     scpnHosts.widget.attrs.update({'onchange' : 'changeScpnHostSelection()'})
    
    objId = forms.CharField(widget=forms.HiddenInput(), required = False)
    label = forms.CharField(label='Label', 
                            required=True,
                            error_messages={'required': 'Please, enter label for the configuration'},
                            widget=forms.TextInput(attrs={'class' : 'generalInput', 'size' : 20}))
    hostName = forms.CharField(label='Host name',
                               required=True, 
                               error_messages={'required': 'Please, enter host name'},
                               widget=forms.TextInput(attrs={'class' : 'generalInput', 'size' : 30}))
    userName = forms.CharField(label='User name', required=False,
                                widget=forms.TextInput(attrs={'class' : 'generalInput', 'size' : 20}))
    password = forms.CharField(label='Password', required=False,
                                widget=forms.PasswordInput())
    password.widget.attrs.update({'class' : 'generalInput', 'size' : 20})
    hostPath = forms.CharField(label='Host path', required=True, error_messages={'required': 'Please, enter your host path'},
                                widget=forms.TextInput(attrs={'class' : 'generalInput', 'size' : 30}))
    mpiCommand = forms.CharField(label='MPI command', required=False,
                                 widget=forms.Textarea(attrs={'cols': 35, 'rows': 5}))
    
    #label.widget.attrs.update({'class' : 'generalInput'})
    
    # Queue system
    queueSystemName = forms.CharField(label='Name', 
                                      required=False,
                                      error_messages={'required': 'Please, enter label for the configuration'},
                                      widget=forms.TextInput(attrs={'class' : 'generalInput', 'size' : 20}))
    queusSystemMandatory = forms.BooleanField(label='Mandatory',
                                              required=False,
                                              error_messages={'required': 'Please, select if queue is mandatory'})
    queusSystemSubmitCommand = forms.CharField(label='Submit command', 
                                      required=False,
                                      error_messages={'required': 'Please, insert submit command'},
                                      widget=forms.TextInput(attrs={'class' : 'generalInput', 'size' : 20}))
    queueSystemSubmitTemplate = forms.CharField(label='Submit template', 
                                      required=False,
                                      error_messages={'required': 'Please, insert submit template'},
                                      widget=forms.TextInput(attrs={'class' : 'generalInput', 'size' : 20}))
    queueSystemSubmitTemplate = forms.CharField(label='Submit template', 
                                      required=False,
                                      error_messages={'required': 'Please, insert submit template'},
                                      widget=forms.TextInput(attrs={'class' : 'generalInput', 'size' : 20}))
    queueSystemCheckCommand = forms.CharField(label='Check command', 
                                      required=False,
                                      error_messages={'required': 'Please, insert check command'},
                                      widget=forms.TextInput(attrs={'class' : 'generalInput', 'size' : 20}))
    queueSystemCheckCommand = forms.CharField(label='Cancel command', 
                                      required=False,
                                      error_messages={'required': 'Please, insert cancel command'},
                                      widget=forms.TextInput(attrs={'class' : 'generalInput', 'size' : 20}))
    
    
    def getHost(self):
        host = HostConfig()
        if self.cleaned_data['objId'] == '':
            host.setObjId(None)
        else:
            host.setObjId(self.cleaned_data['objId'])
        host.setLabel(self.cleaned_data['label'])
        host.setHostName(self.cleaned_data['hostName'])
        host.setUserName(self.cleaned_data['userName'])
        host.setHostPath(self.cleaned_data['hostPath'])
        host.setPassword(self.cleaned_data['password'])
        host.setPassword(self.cleaned_data['mpiCommand'])
        
        return host
    
    def setHost(self, host):
        self.fields['objId'].initial = host.getObjId()
        self.fields['label'].initial = host.getLabel()
        self.fields['hostName'].initial = host.getHostName()
        self.fields['userName'].initial = host.getUserName()
        self.fields['hostPath'].initial = host.getHostPath()
        self.fields['password'].initial = host.getPassword()
        self.fields['mpiCommand'].initial = host.getPassword()   
    
        