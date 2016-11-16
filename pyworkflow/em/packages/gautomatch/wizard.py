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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************

import os

import pyworkflow as pw
import pyworkflow.em.wizard as emwiz
import pyworkflow.utils as pwutils
from pyworkflow.em.viewer import CoordinatesObjectView

from protocol_gautomatch import ProtGautomatch



#===============================================================================
# MASKS
#===============================================================================

class GautomatchParticleWizard(emwiz.ParticleMaskRadiusWizard):
    _targets = [(ProtGautomatch, ['particleSize'])]
    
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
                                            emwiz.UNIT_ANGSTROM)
        
    
#===============================================================================
# PICKER
#===============================================================================

class GautomatchPickerWizard(emwiz.EmWizard):
    _targets = [(ProtGautomatch, ['diameter', 'threshold'])]

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
                "mindist": prot.maxDist.get(),
                "refStack": refStack
          }

        # If Gautomatch will guess advaced parameter we don't need to send
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
        #form.setVar('diameter', myprops['diameter.value'])
        form.setVar('threshold', myprops['threshold.value'])
        if not prot.advanced:
            form.setVar('maxDist', myprops['mindist.value'])
        else:
            pass # TODO: We could even in a future to parse the 'guessed' params