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
    
    #label.widget.attrs.update({'class' : 'generalInput'})
    
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
        return host
    
    def setHost(self, host):
        self.fields['objId'].initial = host.getObjId()
        self.fields['label'].initial = host.getLabel()
        self.fields['hostName'].initial = host.getHostName()
        self.fields['userName'].initial = host.getUserName()
        self.fields['hostPath'].initial = host.getHostPath()
        self.fields['password'].initial = host.getPassword()
        
        
class ShowjForm(forms.Form):
    
    zoom = forms.IntegerField(required=True,
                              max_value=512,
                              min_value=10,
                              localize=False,
                              widget=forms.TextInput(attrs={'class' : 'menuInputNumber'}))
    gotoContainer = forms.IntegerField(required=True,
                              max_value=100,
                              min_value=1,
                              localize=False,
                              widget=forms.TextInput(attrs={'class' : 'menuInputNumber'}))
    cols = forms.IntegerField(required=False,
                              max_value=100,
                              min_value=1,
                              localize=False,
                              widget=forms.TextInput(attrs={'class' : 'menuInputNumber'}))

    rows = forms.IntegerField(required=False,
                              max_value=100,
                              min_value=1,
                              localize=False,
                              widget=forms.TextInput(attrs={'class' : 'menuInputNumber'}))
    blockComboBox = forms.ChoiceField(required=False)
    
    metadataComBox = forms.ChoiceField(required=False)
    
    path = forms.CharField(widget=forms.HiddenInput())
    block = forms.CharField(required=False, widget=forms.HiddenInput())
    allowRender = forms.BooleanField(widget=forms.HiddenInput())
    #imageDim = forms.IntegerField(widget=forms.HiddenInput())#Se puede quitar
    mode = forms.CharField(widget=forms.HiddenInput())

    
    def setShowj(self, path, block, render, dim, mode):
#        self.fields['path'].initial = path
#        self.fields['block'].initial = block
#        self.fields['allowRender'].initial = render
#        self.fields['zoom'].initial = dim
#        self.fields['mode'].initial = mode
        
        self.fields['gotoContainer'].initial = 1
#        
#        print "bound"
#        print self.is_bound
#        print "valid"
#        print self.is_valid()
#        print "error"
#        print self.errors
#        print "cleaned"
#        print self.cleaned_data
    
    