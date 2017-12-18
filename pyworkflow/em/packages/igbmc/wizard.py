# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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

import os

import pyworkflow as pw
import pyworkflow.em.wizard as emwiz
import pyworkflow.utils as pwutils
from pyworkflow.em.viewer import CoordinatesObjectView

from protocol_gempicker import ProtGemPicker



#===============================================================================
# MASKS
#===============================================================================

class GemPickerMaskWizard(emwiz.ParticleMaskRadiusWizard):
    _targets = [(ProtGemPicker, ['maskRadius'])]
    
    def _getParameters(self, protocol):
        
        label, value = self._getInputProtocol(self._targets, protocol)
        
        protParams = {}
        protParams['input']= protocol.inputReferences
        protParams['label']= label
        protParams['value']= value
        return protParams
    
    def _getProvider(self, protocol):
        _objs = self._getParameters(protocol)['input'] 
        return emwiz.ParticleMaskRadiusWizard._getListProvider(self, _objs)
    
    def show(self, form):
        params = self._getParameters(form.protocol)
        _value = params['value']
        _label = params['label']
        emwiz.ParticleMaskRadiusWizard.show(self, form, _value, _label,
                                            emwiz.UNIT_PIXEL)
        
    
#===============================================================================
# PICKER
#===============================================================================

class GemPickerWizard(emwiz.EmWizard):
    _targets = [(ProtGemPicker, ['diameter', 'thresholdLow', 'thresholdHigh'])]

    def show(self, form):
        prot = form.protocol
        micSet = prot.getInputMicrographs()

        if not micSet:
            print 'must specify input micrographs'
            return

        project = prot.getProject()
        micfn = micSet.getFileName()

        # Prepare a temporary folder to convert some input files
        # and put some of the intermediate result files
        coordsDir = project.getTmpPath(micSet.getName())
        pwutils.cleanPath(coordsDir)
        pwutils.makePath(coordsDir)
        prot.convertInputs(coordsDir)

        pickerConfig = os.path.join(coordsDir, 'picker.conf')
        f = open(pickerConfig, "w")

        pickScript = pw.join('em', 'packages', 'igbmc', 'run_gempicker.py')

        pickCmd = prot._getPickArgs(threshold=False, workingDir=coordsDir)[0]
        convertCmd = pw.join('apps', 'pw_convert.py')

        args = {
                "pickScript": pickScript,
                "pickCmd": pickCmd,
                "convertCmd": convertCmd,
                'coordsDir': coordsDir,
                'micsSqlite': micSet.getFileName(),
                'thresholdLow': prot.thresholdLow,
                'thresholdHigh': prot.thresholdHigh,
                "useGPU": prot.useGPU
          }

        f.write("""
        parameters = thresholdLow,thresholdHigh
        thresholdLow.value =  %(thresholdLow)s
        thresholdLow.label = Threshold Low
        thresholdLow.help = Low value cut-off
        thresholdHigh.value =  %(thresholdHigh)s
        thresholdHigh.label = Threshold High
        thresholdHigh.help = High value cut-off
        autopickCommand = %(pickScript)s %%(micrograph) %(coordsDir)s %(useGPU)s %(pickCmd)s --thresh=%%(thresholdLow) --threshHigh=%%(thresholdHigh)
        convertCommand = %(convertCmd)s --coordinates --from gempicker --to xmipp --input  %(micsSqlite)s --output %(coordsDir)s
        """ % args)

        f.close()

        process = CoordinatesObjectView(project, micfn, coordsDir, prot,
                                        pickerProps=pickerConfig).show()
        process.wait()
        myprops = pwutils.readProperties(pickerConfig)

        if myprops['applyChanges'] == 'true':
            form.setVar('thresholdLow', myprops['thresholdLow.value'])
            form.setVar('thresholdHigh', myprops['thresholdHigh.value'])
