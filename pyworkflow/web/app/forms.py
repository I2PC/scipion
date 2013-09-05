'''
Created on Jun 7, 2013

@author: antonio
'''
from django import forms
from pyworkflow.hosts import HostConfig, QueueSystemConfig, QueueConfig
from django.forms.forms import BoundField
from pyworkflow.object import List
#from pyworkflow.web.app.views_showj import get_image_dimensions 

class HostForm(forms.Form):
#     scpnHosts = forms.ChoiceField(label='Scipion hosts', widget = forms.Select(), required = False,)
#     scpnHosts.widget.attrs.update({'onchange' : 'changeScpnHostSelection()'})
    host = None
    
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
                                widget=forms.PasswordInput(render_value = True))
    password.widget.attrs.update({'class' : 'generalInput', 'size' : 20})
    hostPath = forms.CharField(label='Host path', required=True, error_messages={'required': 'Please, enter your host path'},
                                widget=forms.TextInput(attrs={'class' : 'generalInput', 'size' : 30}))
    mpiCommand = forms.CharField(label='MPI command', required=False,
                                 widget=forms.Textarea(attrs={'cols': 35, 'rows': 5}))
    queueSystemConfigCount = forms.CharField(widget=forms.HiddenInput(), initial=0)
    queueConfigCount = forms.CharField(widget=forms.HiddenInput(), initial=0)
    
    #label.widget.attrs.update({'class' : 'generalInput'})

    def __init__(self, *args, **kwargs):
        self.host = None
        extra_queueSystemConfigCount = int(kwargs.pop('queueSystemConfCont', 0))
        extra_queueConfigCount = int(kwargs.pop('queueConfCont', 0))
        super(HostForm, self).__init__(*args, auto_id=True, **kwargs)
        self.fields['queueSystemConfigCount'].initial = extra_queueSystemConfigCount
        self.fields['queueConfigCount'].initial = extra_queueConfigCount
        indexes = None
        if args is not None and len(args)>0:
            indexes = self.getQueueConfigIndexes(args[0], extra_queueConfigCount)
        self.createDynamicFields(extra_queueSystemConfigCount, extra_queueConfigCount, queueConfigIndexes=indexes)
        
                
    def createDynamicFields(self, extra_queueSystemConfigCount, extra_queueConfigCount, queueConfigIndexes = None):
        if extra_queueSystemConfigCount > 0:
            self.createQueueSystemConfigFields()            
            if queueConfigIndexes is None:
                for index in range(1, extra_queueConfigCount + 1):
                    self.createQueueConfigFields(index)
            else:
                for index in queueConfigIndexes:
                    self.createQueueConfigFields(index)   
                    
    def createQueueSystemConfigFields(self):
        self.fields['name'] = forms.CharField(label='Name', 
                                              required=True,
                                              error_messages={'required': 'Please, enter name for the queue system'},
                                              widget=forms.TextInput(attrs={'class' : 'generalInput', 'size' : 20}))
        self.fields['mandatory'] = forms.BooleanField(label='Mandatory',
                                                      required=True,
                                                      error_messages={'required': 'Please, select if queue is mandatory'})
        self.fields['submitTemplate'] = forms.CharField(label='Submit template', 
                                                        required=True,
                                                        error_messages={'required': 'Please, insert submit template'},
                                                        widget=forms.Textarea(attrs={'cols': 35, 'rows': 5}))
        self.fields['submitCommand'] = forms.CharField(label='Submit command', 
                                                       required=True,
                                                       error_messages={'required': 'Please, insert submit command'},
                                                       widget=forms.TextInput(attrs={'class' : 'generalInput', 'size' : 20}))
        self.fields['checkCommand'] = forms.CharField(label='Check command', 
                                                      required=True,
                                                      error_messages={'required': 'Please, insert check command'},
                                                      widget=forms.TextInput(attrs={'class' : 'generalInput', 'size' : 20}))
        self.fields['cancelCommand'] = forms.CharField(label='Cancel command', 
                                                       required=True,
                                                       error_messages={'required': 'Please, insert cancel command'},
                                                       widget=forms.TextInput(attrs={'class' : 'generalInput', 'size' : 20}))   
                    
    def createQueueConfigFields(self, index):
        self.fields['name_{index}'.format(index=index)] = forms.CharField(label='Name',
                                                                          required=True,
                                                                          error_messages={'required': 'Please, enter name for the queue configuration'},
                                                                          widget=forms.TextInput(attrs={'class' : 'generalInput', 'size' : 20}))
        self.fields['maxCores_{index}'.format(index=index)] = forms.IntegerField(label='Max. cores',
                                                                                 required=True,
                                                                                 error_messages={'required': 'Please, enter the maximum cores number'},
                                                                                 widget=forms.TextInput(attrs={'class' : 'generalInput', 'size' : 20}))
        self.fields['allowMPI_{index}'.format(index=index)] = forms.BooleanField(label='Allow MPI',
                                                                                 required=False)
        self.fields['allowThreads_{index}'.format(index=index)] = forms.BooleanField(label='Allow threads',
                                                                                     required=False)
        self.fields['maxHours_{index}'.format(index=index)] = forms.IntegerField(label='Max hours',
                                                                                 required=True,
                                                                                 error_messages={'required': 'Please, enter the maximum hours number'},
                                                                                 widget=forms.TextInput(attrs={'class' : 'generalInput', 'size' : 20}))                
        self.fields['queueConfigId_{index}'.format(index=index)] = forms.CharField(widget=forms.HiddenInput(), required = False)
    
    def getFormHost(self):
        if self.host is None:
            self.host = HostConfig()
            self.host.setObjId(None)
        self.host.setLabel(self.cleaned_data['label'])
        self.host.setHostName(self.cleaned_data['hostName'])
        self.host.setUserName(self.cleaned_data['userName'])
        self.host.setHostPath(self.cleaned_data['hostPath'])
        self.host.setPassword(self.cleaned_data['password'])
        self.host.setMpiCommand(self.cleaned_data['mpiCommand'])   
        if int(self.cleaned_data['queueSystemConfigCount']) > 0:
            if self.host.getQueueSystem() is None:
                queueSystemConfig = QueueSystemConfig()
            else:
                queueSystemConfig = self.host.getQueueSystem()
            queueSystemConfig.setName(self.cleaned_data['name'])
            queueSystemConfig.setMandatory(self.cleaned_data['mandatory'])
            queueSystemConfig.setSubmitTemplate(self.cleaned_data['submitTemplate'])
            queueSystemConfig.setSubmitCommand(self.cleaned_data['submitCommand'])
            queueSystemConfig.setCheckCommand(self.cleaned_data['checkCommand'])
            queueSystemConfig.setCancelCommand(self.cleaned_data['cancelCommand'])
            queueConfigCount = int(self.cleaned_data['queueConfigCount'])
            if queueConfigCount > 0:
                queuesList = List()
                index = 1
                cont = 0
                while cont < queueConfigCount:
#                 for index in range(1, int(self.cleaned_data['queueConfigCount'])+1):
                    if ('queueConfigId_{index}'.format(index=index)) in self.cleaned_data:
                        queueConfig = None
                        if self.cleaned_data['queueConfigId_{index}'.format(index=index)] == '':
                            queueConfig = QueueConfig()
                            queueConfig.setObjId(None)
                        else:                        
                            objId = int(self.cleaned_data['queueConfigId_{index}'.format(index=index)])
                            queueConfig = queueSystemConfig.getQueueConfig(objId)                            
                        queueConfig.setName(self.cleaned_data['name_{index}'.format(index=index)])
                        queueConfig.setMaxCores(self.cleaned_data['maxCores_{index}'.format(index=index)])
                        queueConfig.setAllowMPI(self.cleaned_data['allowMPI_{index}'.format(index=index)])
                        queueConfig.setAllowThreads(self.cleaned_data['allowThreads_{index}'.format(index=index)])
                        queueConfig.setMaxHours(self.cleaned_data['maxHours_{index}'.format(index=index)])
                        queuesList.append(queueConfig)
                        cont += 1
                    index += 1
                        
                if queueSystemConfig.getQueues() is not None:
                    queueSystemConfig.getQueues().clear()
                    for queue in queuesList:
                        queueSystemConfig.getQueues().append(queue)
                else:
                    queueSystemConfig.setQueues(queuesList)
            self.host.setQueueSystem(queueSystemConfig)                    
        return self.host
    
    def setFormHost(self, host):
        self.host = host
        self.fields['objId'].initial = host.getObjId()
        self.fields['label'].initial = host.getLabel()
        self.fields['hostName'].initial = host.getHostName()
        self.fields['userName'].initial = host.getUserName()
        self.fields['hostPath'].initial = host.getHostPath()
        self.fields['password'].initial = host.getPassword()
        self.fields['mpiCommand'].initial = host.getMpiCommand() 
        self.fields['queueSystemConfigCount'].initial = 0
        queueSystem = host.getQueueSystem()
        if queueSystem is not None:
            self.fields['queueSystemConfigCount'].initial = 1
            lenQueues = 0;
            if queueSystem.getQueues() is not None:
                lenQueues = len(queueSystem.getQueues())
            self.createDynamicFields(1, lenQueues)
            queueSystem = host.getQueueSystem()
            self.fields['name'].initial = queueSystem.getName()
            self.fields['mandatory'].initial = queueSystem.getMandatory()
            self.fields['submitTemplate'].initial = queueSystem.getSubmitTemplate()
            self.fields['submitCommand'].initial = queueSystem.getSubmitCommand()
            self.fields['checkCommand'].initial = queueSystem.getCheckCommand()
            self.fields['cancelCommand'].initial = queueSystem.getCancelCommand() 
            index = 1
            if queueSystem.getQueues() is None:
                self.fields['queueConfigCount'].initial = 0
            else:
                self.fields['queueConfigCount'].initial = len(queueSystem.getQueues())
            if self.fields['queueConfigCount'].initial > 0:
                for queueConfig in queueSystem.getQueues():
                    self.fields['queueConfigId_{index}'.format(index=index)].initial = queueConfig.getObjId()
                    self.fields['name_{index}'.format(index=index)].initial = queueConfig.getName()
                    self.fields['maxCores_{index}'.format(index=index)].initial = queueConfig.getMaxCores()
                    self.fields['allowMPI_{index}'.format(index=index)].initial = queueConfig.getAllowMPI()
                    self.fields['allowThreads_{index}'.format(index=index)].initial = queueConfig.getAllowThreads()
                    self.fields['maxHours_{index}'.format(index=index)].initial = queueConfig.getMaxHours()
                    index += 1
        
    def getBasicHostFields(self):
        list = []
        list.append(self['label'])
        list.append(self['hostName'])
        list.append(self['userName'])
        list.append(self['password'])
        list.append(self['hostPath'])
        list.append(self['mpiCommand'])
        return list
    
    def getQueueSystemConfigFields(self):
        list = []
        if 'queueSystemConfigCount' in self.fields and self.fields['queueSystemConfigCount'].initial != '':
            queueSystemConfigCount = int(self.fields['queueSystemConfigCount'].initial)
            if queueSystemConfigCount > 0:
                list.append(BoundField(self,self.fields['name'], 'name'))
                list.append(BoundField(self,self.fields['mandatory'], 'mandatory'))
                list.append(BoundField(self,self.fields['submitTemplate'], 'submitTemplate'))
                list.append(BoundField(self,self.fields['submitCommand'], 'submitCommand'))
                list.append(BoundField(self,self.fields['checkCommand'], 'checkCommand'))
                list.append(BoundField(self,self.fields['cancelCommand'], 'cancelCommand'))
        return list
    
    def getQueueConfigFields(self):
        resultList = []
        if 'queueConfigCount' in self.fields and self.fields['queueConfigCount'].initial != '':
            queueConfigCount = int(self.fields['queueConfigCount'].initial)
            if queueConfigCount > 0:
                for index in range(1, queueConfigCount+1):
                    list = []
                    list.append(BoundField(self,self.fields['name_{index}'.format(index=index)], 'name_{index}'.format(index=index)))
                    list.append(BoundField(self,self.fields['maxCores_{index}'.format(index=index)], 'maxCores_{index}'.format(index=index)))
                    list.append(BoundField(self,self.fields['allowMPI_{index}'.format(index=index)], 'allowMPI_{index}'.format(index=index)))
                    list.append(BoundField(self,self.fields['allowThreads_{index}'.format(index=index)], 'allowThreads_{index}'.format(index=index)))
                    list.append(BoundField(self,self.fields['maxHours_{index}'.format(index=index)], 'maxHours_{index}'.format(index=index)))
                    resultList.append(list)
        return resultList
    
    def getQueueConfigIndexes(self, post, number):
        indexes = []
        index = 1
        cont = 0        
        while cont < number:
            field = post.get('queueConfigId_{index}'.format(index=index))
            if field is not None:
                indexes.append(index)
                cont += 1
            index +=1
        return indexes
        
        
######################### Initialize Showj Form (Button toolbar) #####################        
class ShowjForm(forms.Form):
    # Init graphical components
    zoom = forms.CharField(required=True,
                            widget=forms.TextInput(attrs={'class' : 'menuInputNumber'}))
    
    goto = forms.IntegerField(required=True,
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
    
    
    # Init hidden fields
    path = forms.CharField(widget=forms.HiddenInput())
    allowRender = forms.BooleanField(widget=forms.HiddenInput())
    mode = forms.CharField(widget=forms.HiddenInput())
    colRowMode = forms.CharField(widget=forms.HiddenInput())
    
    mirrorY = forms.BooleanField(label='Invert Y axis')

#    imageWidth = forms.IntegerField(widget=forms.HiddenInput(), initial=0)
#    imageHeight = forms.IntegerField(widget=forms.HiddenInput(), initial=0)
    
     
    
    
    def __init__(self, dataset, tableLayoutConfiguration, *args, **kwargs):
        super(ShowjForm, self).__init__(*args, **kwargs)
        
        blockComboBoxValues = tuple(zip(dataset.listTables(), dataset.listTables()))
        self.fields['blockComboBox'] = forms.ChoiceField(label='Select Block',
                                                         required=False,
                                                         choices = blockComboBoxValues)

        
        labelsToRenderComboBoxValues = getLabelsToRenderComboBoxValues(tableLayoutConfiguration.columnsLayout)
        if len(labelsToRenderComboBoxValues) > 0:
            self.fields['labelsToRenderComboBox'] = forms.ChoiceField(label='Select Label',
                                                            required=False,
                                                            choices = labelsToRenderComboBoxValues)
            if self.data['mode'] != 'gallery':
                self.fields['labelsToRenderComboBox'].widget=forms.HiddenInput()
        else:
            self.fields['zoom'].widget.attrs['readonly'] = True
            
        
        if self.data['mode'] != 'gallery': 
            self.fields['cols'].widget=forms.HiddenInput()
            self.fields['rows'].widget=forms.HiddenInput()
        
        if self.data['colRowMode'] == 'Off':
            self.fields['cols'].widget.attrs['readonly'] = True
            self.fields['rows'].widget.attrs['readonly'] = True
                    
def getLabelsToRenderComboBoxValues(columnsLayout):
    labelsToRender = [columnLayout.label for columnLayout in columnsLayout.values() if (columnLayout.typeOfColumn == 'image')]
    return tuple(zip(labelsToRender,labelsToRender))

    

class VolVisualizationForm(forms.Form):   
    volPath = forms.CharField(label='Volume path', 
                            required=True,
                            error_messages={'required': 'Please, enter a volume path'})
    volumeTypes = ((0, "Byte"),(1, "Integer"),(2, "Float"))
    volType = forms.ChoiceField(label='Volume data type',
                             required=True,
                             choices = volumeTypes)
    threshold = forms.DecimalField(label='Threshold', required=True,
                                widget=forms.TextInput(attrs={'size' : 10})) 
