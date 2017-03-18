# **************************************************************************
# *
# * Authors:    Adrian Quintana (aquintana@cnb.csic.es)
# *                
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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

from django import forms
from django.forms.forms import BoundField
from pyworkflow.object import List
import xmipp
from pyworkflow.utils.properties import Message
import pyworkflow.em.showj as sj

######################### Initialize Showj Form (Button toolbar) #####################        
class ShowjForm(forms.Form):
    #Init message properties file
    messagesForm = Message()
    
    # Init graphical components
    zoom = forms.CharField(required=True,
                            widget=forms.TextInput(attrs={'class' : 'menuInputNumber'}))
    
    goto = forms.IntegerField(required=True,
                              max_value=100,
                              min_value=1,
                              localize=False,
                              widget=forms.TextInput(attrs={'class' : 'menuInputNumber'}))
    
    
    cols = forms.IntegerField(label=messagesForm.LABEL_COLS,
                              required=False,
                              max_value=100,
                              min_value=1,
                              localize=False,
                              widget=forms.TextInput(attrs={'class' : 'menuInputNumber'}))

    rows = forms.IntegerField(label=messagesForm.LABEL_ROWS,
                              required=False,
                              max_value=100,
                              min_value=1,
                              localize=False,
                              widget=forms.TextInput(attrs={'class' : 'menuInputNumber'}))
    
    resliceComboBox = forms.ChoiceField(label=messagesForm.LABEL_RESLICE,
                                required=False,
                                choices=((xmipp.VIEW_Z_NEG, messagesForm.RESLICE_Z_NEGATIVE),
                                           (xmipp.VIEW_Y_NEG, messagesForm.RESLICE_Y_NEGATIVE),
                                           (xmipp.VIEW_X_NEG, messagesForm.RESLICE_X_NEGATIVE),
                                           (xmipp.VIEW_Y_POS, messagesForm.RESLICE_Y_POSITIVE),
                                           (xmipp.VIEW_X_POS, messagesForm.RESLICE_X_POSITIVE))
                                )

    # Init hidden fields
    path = forms.CharField(widget=forms.HiddenInput())
    allowRender = forms.BooleanField(widget=forms.HiddenInput(), required=False)
    mode = forms.CharField(widget=forms.HiddenInput())
    colRowMode = forms.CharField(widget=forms.HiddenInput())
#    dims = forms.CharField(widget=forms.HiddenInput())
    
    imageMaxWidth = forms.CharField(widget=forms.HiddenInput())
    imageMinWidth = forms.CharField(widget=forms.HiddenInput())
    imageMaxHeight = forms.CharField(widget=forms.HiddenInput())
    imageMinHeight = forms.CharField(widget=forms.HiddenInput())
    
    mirrorY = forms.BooleanField(label=messagesForm.LABEL_MIRRORY, required=False)
#    applyTransformation = forms.BooleanField(label='Apply Transform', required=False)
    
#    CHOICES = (('transformMatrix', 'Apply Transform Matrix',), (IMG_ONLY_SHIFTS, 'Only Shifts',), ('wrap', 'Wrap',))
#    transformationChoice = forms.ChoiceField(widget=forms.RadioSelect, choices=(('transformMatrix', 'Apply Transform Matrix',), (IMG_ONLY_SHIFTS, 'Only Shifts',), ('wrap', 'Wrap',))
#                                             )
    applyTransformMatrix = forms.BooleanField(label=messagesForm.LABEL_APPLY_TRANSFORM, required=False)
    onlyShifts = forms.BooleanField(label=messagesForm.LABEL_ONLY_SHIFTS, required=False)
    wrap = forms.BooleanField(label=messagesForm.LABEL_WRAP, required=False)

    typeVolume = forms.CharField(widget=forms.HiddenInput())

    
    def __init__(self, dataset=None, tableLayoutConfiguration=None, labelsToRender = None, *args, **kwargs):
        #Init message properties file
        messagesForm = Message()
        
        super(ShowjForm, self).__init__(*args, **kwargs)
        
        if dataset is not None:
            self.fields['blockComboBox'] = forms.ChoiceField(label=messagesForm.LABEL_BLOCK_SELECTION,
                                                             required=False,
                                                             choices=tuple(zip(dataset.listTables(), dataset.listTables())))
            
            # This rendarable columns cannot have to access to the extra params
            namesToRender, labelsToRender = tableLayoutConfiguration.getRenderableColumns(extra=labelsToRender)
            
            if namesToRender:
                labelChoices = tuple(zip(namesToRender, labelsToRender))
                self.fields[sj.LABEL_SELECTED] = forms.ChoiceField(label=messagesForm.LABEL_LABEL_SELECTION,
                                                                required=False,
                                                                choices=labelChoices)
                if self.data[sj.MODE] != sj.MODE_GALLERY:
                    self.fields[sj.LABEL_SELECTED].widget = forms.HiddenInput()
            else:
                self.fields[sj.IMG_ZOOM].widget.attrs['readonly'] = True
                if self.data[sj.MODE] == sj.MODE_GALLERY:
                    self.fields[sj.GOTO].widget.attrs['readonly'] = True
                
            if dataset.getNumberSlices() > 1:
                labelColumn = dataset.getTable().getColumnValues(self.data[sj.LABEL_SELECTED])
                volumesToRenderComboBoxValues = tuple(zip(labelColumn, labelColumn))
                
                self.fields[sj.VOL_SELECTED] = forms.ChoiceField(label=messagesForm.LABEL_VOLUME_SELECTION,
                                                                required=False,
                                                                choices=volumesToRenderComboBoxValues)
                if self.data[sj.MODE] == sj.MODE_TABLE:
                    self.fields[sj.VOL_SELECTED].widget = forms.HiddenInput()
                if self.data[sj.MODE] != sj.MODE_GALLERY:                
                    self.fields[sj.VOL_VIEW].widget = forms.HiddenInput()
            else:
                self.fields[sj.VOL_VIEW].widget = forms.HiddenInput()
                
        else:
            self.fields[sj.IMG_ZOOM].widget = forms.HiddenInput()
            self.fields[sj.GOTO].widget = forms.HiddenInput()
            
            self.fields[sj.VOL_VIEW].widget = forms.HiddenInput()

            self.fields[sj.IMG_MIRRORY].widget = forms.HiddenInput()
            self.fields[sj.IMG_APPLY_TRANSFORM].widget = forms.HiddenInput()
            self.fields[sj.IMG_ONLY_SHIFTS].widget = forms.HiddenInput()
            self.fields[sj.IMG_WRAP].widget = forms.HiddenInput()
               
        
        if self.data[sj.MODE] != sj.MODE_GALLERY: 
            self.fields[sj.COLS].widget = forms.HiddenInput()
            self.fields[sj.ROWS].widget = forms.HiddenInput()
        
        if self.data[sj.MANUAL_ADJUST] == 'Off':
            self.fields[sj.COLS].widget.attrs['readonly'] = True
            self.fields[sj.ROWS].widget.attrs['readonly'] = True
            
        if self.data[sj.MODE] == 'column':    
            self.fields[sj.GOTO].widget.attrs['readonly'] = True
        
                                  
class DocumentForm(forms.Form):
    docfile = forms.FileField(
        label='Select a file'
    )
