'''
Created on Jun 7, 2013

@author: antonio
'''
from django import forms
from pyworkflow.hosts import ExecutionHostConfig

class HostForm(forms.Form):
    scpnHosts = forms.ChoiceField(label='Scipion hosts', widget = forms.Select(), required = False,)
    scpnHosts.widget.attrs.update({'onchange' : 'changeScpnHostSelection()'})
    
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
    hostPath = forms.CharField(label='Host path', required=True, error_messages={'required': 'Please, enter your host path'},
                                widget=forms.TextInput(attrs={'class' : 'generalInput', 'size' : 30}))
    
    #label.widget.attrs.update({'class' : 'generalInput'})
    
    def getHost(self):
        host = ExecutionHostConfig()
        host.setLabel(self.cleaned_data['label'])
        host.setHostName(self.cleaned_data['hostName'])
        host.setUserName(self.cleaned_data['userName'])
        host.setHostPath(self.cleaned_data['hostPath'])
        return host
    
    def setHost(self, host):
        self.fields['label'].initial = host.getLabel()
        self.fields['hostName'].initial = host.getHostName()
        self.fields['userName'].initial = host.getUserName()
        self.fields['hostPath'].initial = host.getHostPath()