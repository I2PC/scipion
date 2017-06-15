# **************************************************************************
# *
# * Authors:     Grigory Sharov (sharov@igbmc.fr)
# *
# * L'Institut de genetique et de biologie moleculaire et cellulaire (IGBMC)
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
import pyworkflow.gui.dialog as dialog
from pyworkflow.em.viewer import CoordinatesObjectView
from pyworkflow.em.constants import *

from protocol_gautomatch import ProtGautomatch



#===============================================================================
# MASKS
#===============================================================================

class GautomatchParticleWizard(emwiz.ParticleMaskRadiusWizard):
    _targets = [(ProtGautomatch, ['particleSize'])]
    
    def _getParameters(self, protocol):
        
        label, value = self._getInputProtocol(self._targets, protocol)
        
        protParams = {}
        protParams['input'] = protocol.inputReferences
        protParams['label'] = label
        protParams['value'] = value
        return protParams
    
    def _getProvider(self, protocol):
        _objs = self._getParameters(protocol)['input'] 
        return emwiz.ParticleMaskRadiusWizard._getListProvider(self, _objs)
    
    def show(self, form):
        params = self._getParameters(form.protocol)
        _value = params['value']
        _label = params['label']
        emwiz.ParticleMaskRadiusWizard.show(self, form, _value, _label,
                                            emwiz.UNIT_ANGSTROM)

# ===============================================================================
# FILTERS
# ===============================================================================

class GautomatchBandpassWizard(emwiz.FilterParticlesWizard):
    _targets = [(ProtGautomatch, ['lowPass', 'highPass'])]

    def _getParameters(self, protocol):

        label, value = self._getInputProtocol(self._targets, protocol)

        protParams = {}
        protParams['input'] = protocol.inputMicrographs
        protParams['label'] = label
        protParams['value'] = value
        protParams['mode'] = FILTER_NO_DECAY
        return protParams

    def _getProvider(self, protocol):
        _objs = self._getParameters(protocol)['input']
        return emwiz.FilterParticlesWizard._getListProvider(self, _objs)

    def show(self, form):
        protocol = form.protocol
        provider = self._getProvider(protocol)
        params = self._getParameters(protocol)

        if provider is not None:

            args = {'mode': params['mode'],
                    'lowFreq': params['value'][1],
                    'highFreq': params['value'][0],
                    'unit': UNIT_ANGSTROM
                    }

            args['showDecay'] = False

            d = emwiz.BandPassFilterDialog(form.root, provider, **args)

            if d.resultYes():
                form.setVar('lowPass', d.getHighFreq())
                form.setVar('highPass', d.getLowFreq())

        else:
            dialog.showWarning("Input micrographs", "Select micrographs first", form.root)
        
    
#===============================================================================
# PICKER
#===============================================================================

class GautomatchPickerWizard(emwiz.EmWizard):
    _targets = [(ProtGautomatch, ['threshold'])]

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
        refStack = os.path.join(coordsDir, 'references.mrcs')
        prot.convertReferences(refStack)

        # Get current values of the properties
#         micfn = os.path.join(coordsDir, 'micrographs.xmd')
#         writeSetOfMicrographs(micSet, micfn)
        pickerConfig = os.path.join(coordsDir, 'picker.conf')
        f = open(pickerConfig, "w")

        pickScript = pw.join('em', 'packages', 'gautomatch',
                             'run_gautomatch.py')

        pickCmd = prot.getArgs(threshold=False, mindist=False)
        convertCmd = pw.join('apps', 'pw_convert.py')

        args = {
                "pickScript": pickScript,
                "pickCmd": pickCmd,
                "convertCmd": convertCmd,
                'coordsDir': coordsDir,
                'micsSqlite': micSet.getFileName(),
                'threshold': prot.threshold.get(),
                "mindist": prot.minDist.get(),
                "refStack": refStack
          }

        # If Gautomatch will guess advanced parameter we don't need to send
        # the min distance to the wizard.
        if prot.advanced:
            f.write("""
            parameters = threshold
            threshold.value =  %(threshold)s
            threshold.label = Threshold
            threshold.help = Particles with CCC above the threshold will be picked
            autopickCommand = %(pickScript)s %%(micrograph) %(refStack)s %(coordsDir)s %(pickCmd)s --cc_cutoff %%(threshold)
            convertCommand = %(convertCmd)s --coordinates --from gautomatch --to xmipp --input  %(micsSqlite)s --output %(coordsDir)s
            """ % args)

        else:
            f.write("""
            parameters = threshold,mindist
            threshold.value =  %(threshold)s
            threshold.label = Threshold
            threshold.help = Particles with CCC above the threshold will be picked
            mindist.value = %(mindist)s
            mindist.label = Min search distance
            mindist.help = Use value of 0.9~1.1X diameter
            autopickCommand = %(pickScript)s %%(micrograph) %(refStack)s %(coordsDir)s %(pickCmd)s --cc_cutoff %%(threshold) --min_dist %%(mindist)
            convertCommand = %(convertCmd)s --coordinates --from gautomatch --to xmipp --input  %(micsSqlite)s --output %(coordsDir)s
            """ % args)

        f.close()

        process = CoordinatesObjectView(project, micfn, coordsDir, prot,
                                        pickerProps=pickerConfig).show()
        process.wait()
        myprops = pwutils.readProperties(pickerConfig)

        if myprops['applyChanges'] == 'true':
            form.setVar('threshold', myprops['threshold.value'])
            if not prot.advanced:
                form.setVar('minDist', myprops['mindist.value'])
            else:
                pass  # TODO: We could even in future parse the 'guessed' params
