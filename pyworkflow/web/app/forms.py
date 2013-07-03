'''
Created on Jun 7, 2013

@author: antonio
'''
from django import forms
from pyworkflow.hosts import HostConfig, QueueSystemConfig, QueueConfig

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
    queueSystemConfigCount = forms.CharField(widget=forms.HiddenInput())
    queueConfigCount = forms.CharField(widget=forms.HiddenInput())
    
    #label.widget.attrs.update({'class' : 'generalInput'})
    
    # Queue system
    '''
    queueSystemName = forms.CharField(label='Name', 
                                      required=True,
                                      error_messages={'required': 'Please, enter label for the queue system'},
                                      widget=forms.TextInput(attrs={'class' : 'generalInput', 'size' : 20}))
    queueSystemMandatory = forms.BooleanField(label='Mandatory',
                                              required=False,
                                              error_messages={'required': 'Please, select if queue is mandatory'})
    queueSystemSubmitTemplate = forms.CharField(label='Submit template', 
                                      required=False,
                                      error_messages={'required': 'Please, insert submit template'},
                                      widget=forms.TextInput(attrs={'class' : 'generalInput', 'size' : 20}))
    queueSystemSubmitCommand = forms.CharField(label='Submit command', 
                                      required=False,
                                      error_messages={'required': 'Please, insert submit command'},
                                      widget=forms.TextInput(attrs={'class' : 'generalInput', 'size' : 20}))
    queueSystemCheckCommand = forms.CharField(label='Check command', 
                                      required=False,
                                      error_messages={'required': 'Please, insert check command'},
                                      widget=forms.TextInput(attrs={'class' : 'generalInput', 'size' : 20}))
    queueSystemCancelCommand = forms.CharField(label='Cancel command', 
                                      required=False,
                                      error_messages={'required': 'Please, insert cancel command'},
                                      widget=forms.TextInput(attrs={'class' : 'generalInput', 'size' : 20}))
    '''
    def __init__(self, *args, **kwargs):
        extra_queueSystemConfigCount = int(kwargs.pop('queueSystemConfCont', 0))
        extra_queueConfigCount = int(kwargs.pop('queueConfCont', 0))
        super(HostForm, self).__init__(*args, **kwargs)
        self.fields['queueSystemConfigCount'].initial = extra_queueSystemConfigCount
        self.fields['queueConfigCount'].initial = extra_queueConfigCount
        
        if extra_queueSystemConfigCount > 0:
            
            self.fields['name'] = forms.CharField(label='Name', 
                                      required=True,
                                      error_messages={'required': 'Please, enter name for the queue system'},
                                      widget=forms.TextInput(attrs={'class' : 'generalInput', 'size' : 20}))
            self.fields['mandatory'] = forms.BooleanField(label='Mandatory',
                                                      required=False,
                                                      error_messages={'required': 'Please, select if queue is mandatory'})
            self.fields['submitTemplate'] = forms.CharField(label='Submit template', 
                                              required=False,
                                              error_messages={'required': 'Please, insert submit template'},
                                              widget=forms.TextInput(attrs={'class' : 'generalInput', 'size' : 20}))
            self.fields['submitCommand'] = forms.CharField(label='Submit command', 
                                              required=False,
                                              error_messages={'required': 'Please, insert submit command'},
                                              widget=forms.TextInput(attrs={'class' : 'generalInput', 'size' : 20}))
            self.fields['checkCommand'] = forms.CharField(label='Check command', 
                                              required=False,
                                              error_messages={'required': 'Please, insert check command'},
                                              widget=forms.TextInput(attrs={'class' : 'generalInput', 'size' : 20}))
            self.fields['cancelCommand'] = forms.CharField(label='Cancel command', 
                                              required=False,
                                              error_messages={'required': 'Please, insert cancel command'},
                                              widget=forms.TextInput(attrs={'class' : 'generalInput', 'size' : 20}))
        
            for index in range(extra_queueConfigCount):
                self.fields['name_{index}'.format(index=index)] = forms.CharField(required=True,
                                                                                  error_messages={'required': 'Please, enter name for the queue configuration'})
                self.fields['maxCores_{index}'.format(index=index)] = forms.IntegerField(required=False)
                self.fields['allowMPI_{index}'.format(index=index)] = forms.BooleanField(required=False)
                self.fields['allowThreads_{index}'.format(index=index)] = forms.BooleanField(required=False)
                self.fields['maxHours_{index}'.format(index=index)] = forms.IntegerField(required=False)
            
    
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
        host.setMpiCommand(self.cleaned_data['mpiCommand'])    
        if self.cleaned_data['queueSystemConfigCount'] >0:
            queueSystemConfig = QueueSystemConfig()
            queueSystemConfig.setName(self.cleaned_data['name'])
            queueSystemConfig.setMandatory(self.cleaned_data['mandatory'])
            queueSystemConfig.setSubmitTemplate(self.cleaned_data['submitTemplate'])
            queueSystemConfig.setSubmitCommand(self.cleaned_data['submitCommand'])
            queueSystemConfig.setCheckCommand(self.cleaned_data['checkCommand'])
            queueSystemConfig.setCancelCommand(self.cleaned_data['cancelCommand'])
            if self.cleaned_data['queueConfigCount'] >0:
                queuesList = []
                for index in range(self.cleaned_data['queueConfigCount']):
                    queueConfig = QueueConfig()
                    queueSystemConfig.setName(self.cleaned_data['name_{index}'.format(index=index)])
                    queueSystemConfig.setMaxCores(self.cleaned_data['maxCores_{index}'.format(index=index)])
                    queueSystemConfig.setAllowMPI(self.cleaned_data['allowMPI_{index}'.format(index=index)])
                    queueSystemConfig.setAllowThreads(self.cleaned_data['allowThreads_{index}'.format(index=index)])
                    queueSystemConfig.setMaxHours(self.cleaned_data['maxHours_{index}'.format(index=index)])
                    queuesList.append(queueConfig)
                queueSystemConfig.setQueues(queuesList)
            host.setQueueSystem(queueSystemConfig)                    
        return host
    
    def setHost(self, host):
        self.fields['objId'].initial = host.getObjId()
        self.fields['label'].initial = host.getLabel()
        self.fields['hostName'].initial = host.getHostName()
        self.fields['userName'].initial = host.getUserName()
        self.fields['hostPath'].initial = host.getHostPath()
        self.fields['password'].initial = host.getPassword()
        self.fields['mpiCommand'].initial = host.getMpiCommand() 
        
    
        
        
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
    
    
    cols = forms.IntegerField(label='Cols',
                              required=False,
                              max_value=100,
                              min_value=1,
                              localize=False,
                              widget=forms.TextInput(attrs={'class' : 'menuInputNumber'}))

    rows = forms.IntegerField(label='Rows',
                              required=False,
                              max_value=100,
                              min_value=1,
                              localize=False,
                              widget=forms.TextInput(attrs={'class' : 'menuInputNumber'}))
    
    
    path = forms.CharField(widget=forms.HiddenInput())
    allowRender = forms.BooleanField(widget=forms.HiddenInput())
    mode = forms.CharField(widget=forms.HiddenInput())
    
#    blockComboBox = forms.ChoiceField(required=False)
#    blockComboBox = forms.ChoiceField(required=False, choices=[('1','1'),('2','2'),('3','3')], initial = ('3','3'))

#    metadataComboBox = forms.ChoiceField(required=False)
    
    def __init__(self, mdXmipp, *args, **kwargs):
        super(ShowjForm, self).__init__(*args, **kwargs)
        
        blockComboBoxValues = getBlockComboBoxValues(self.data["path"])
        self.fields['blockComboBox'] = forms.ChoiceField(label='Select Block',
                                                         required=False,
                                                         choices = blockComboBoxValues)

        
        metadataComboBoxValues = getMetadataComboBoxValues(mdXmipp, self.data["allowRender"])
        if len(metadataComboBoxValues) > 0:
            self.fields['metadataComboBox'] = forms.ChoiceField(label='Select Metadata',
                                                            required=False,
                                                            choices = metadataComboBoxValues)
            if self.data['mode'] != 'gallery':
                self.fields['metadataComboBox'].widget=forms.HiddenInput()
                self.fields['cols'].widget.attrs['readonly'] = True
                self.fields['rows'].widget.attrs['readonly'] = True
                    

def getBlockComboBoxValues(path):    
    import xmipp
    from pyworkflow.tests import getInputPath
    blocks = xmipp.getBlocksInMetaDataFile(str(getInputPath('showj', path)))
    return tuple(zip(blocks, blocks))

def getMetadataComboBoxValues(mdXmipp, allowRender):
    import xmipp
    from pyworkflow.web.app.views_showj import getTypeOfColumns
    labels = mdXmipp.getActiveLabels()
    labelsToRender = [xmipp.label2Str(l) for l in labels if (xmipp.labelIsImage(l) and allowRender)]
    #self.fields['metadataComboBox'].choices = zip(labelsToRender,labelsToRender)
    return tuple(zip(labelsToRender,labelsToRender))
