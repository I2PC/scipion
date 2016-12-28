# **************************************************************************
# *
# * Authors:     Jose Gutierrez (jose.gutierrez@cnb.csic.es)
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

from pyworkflow.em.packages.xmipp3.constants import *

from constants import *
from pyworkflow.em import *
from pyworkflow.em.wizard import *
from protocol_classify3d import ProtRelionClassify3D
from protocol_refine3d import ProtRelionRefine3D
from protocol_classify2d import ProtRelionClassify2D
from protocol_preprocess import ProtRelionPreprocessParticles
from protocol_autopick import ProtRelionAutopickFom, ProtRelionAutopick
from protocol_sort import ProtRelionSortParticles
from pyworkflow.utils.utils import readProperties

#===============================================================================
# MASKS
#===============================================================================

class RelionBackRadiusWizard(ParticleMaskRadiusWizard):
    _targets = [(ProtRelionPreprocessParticles, ['backRadius'])]
    _unit = UNIT_PIXEL
    
    def _getProtocolImages(self, protocol):
        return protocol.inputParticles
    
    def _getParameters(self, protocol):
        
        label, value = self._getInputProtocol(self._targets, protocol)
        
        protParams = {}
        protParams['input']= self._getProtocolImages(protocol)
        protParams['label']= label
        protParams['value']= value
        return protParams  
    
    def _getProvider(self, protocol):
        _objs = self._getParameters(protocol)['input']    
        return ParticleMaskRadiusWizard._getListProvider(self, _objs)
    
    def show(self, form):
        params = self._getParameters(form.protocol)
        _value = params['value']
        _label = params['label']
        ParticleMaskRadiusWizard.show(self, form, _value, _label, units=self._unit)


class RelionPartMaskDiameterWizard(RelionBackRadiusWizard):
    _targets = [(ProtRelionClassify2D, ['maskDiameterA']),
                (ProtRelionRefine3D, ['maskDiameterA']),
                (ProtRelionClassify3D, ['maskDiameterA']),
                (ProtRelionClassify2D, ['maskDiameterA'])]
    _unit = UNIT_ANGSTROM

    def _getParameters(self, protocol):
        protParams = RelionBackRadiusWizard._getParameters(self, protocol)
        # adjust to from diameter to radius
        protParams['value'] = protParams['value'] / 2

        return protParams

    def setVar(self, form, label, value):
        # adjust again from radius to diameter
        form.setVar(label, value * 2)

# We need this specific wizard for the sort protocol because
# this protocol have a particular way to grab the input images
class RelionSortMaskWizard(RelionPartMaskDiameterWizard):
    _targets = [(ProtRelionSortParticles, ['maskDiameterA'])]

    def _getProvider(self, protocol):
        if protocol.isInputClasses():
            images = [cls.clone()
                      for cls in protocol.inputSet.get().iterRepresentatives()]
        else:
            images = self._getParticles(protocol._allParticles(iterate=True))
        return ListTreeProvider(images)

    def _getProtocolImages(self, protocol):
        return None
        #return protocol._allParticles(iterate=True)


#===============================================================================
# FILTER
#===============================================================================

class RelionVolFilterWizard(FilterVolumesWizard):
    _targets = [(ProtRelionClassify3D, ['initialLowPassFilterA']),
                (ProtRelionRefine3D, ['initialLowPassFilterA'])]
    
    def _getParameters(self, protocol):
        
        label, value = self._getInputProtocol(self._targets, protocol)
        
        protParams = {}
        protParams['input']= protocol.referenceVolume
        protParams['label']= label
        protParams['value']= value
        protParams['mode'] = FILTER_LOW_PASS_NO_DECAY
        return protParams  
    
    def _getProvider(self, protocol):
        _objs = self._getParameters(protocol)['input']    
        return FilterVolumesWizard._getListProvider(self, _objs)
    
    # def show(self, form):
    #     params = self._getParameters(form.protocol)
    #     # Value should be LowFreq=0, HighFreq and Decay for Low pass filter
    #     _value = params['value']
    #     _label = params['label']
    #     FilterVolumesWizard.show(self, form, _value, _label,
    #                              mode=FILTER_LOW_PASS,
    #                              unit=UNIT_ANGSTROM,
    #                              showDecay=False)

    def show(self, form):
        params = self._getParameters(form.protocol)
        protocol = form.protocol
        provider = self._getProvider(protocol)

        if provider is not None:

            args = {'mode': params['mode'],
                    'highFreq': params['value'],
                    'unit': UNIT_ANGSTROM
                    }

            args['showLowFreq'] = False
            args['showDecay'] = False

            d = BandPassFilterDialog(form.root, provider, **args)

            if d.resultYes():
                form.setVar('initialLowPassFilterA', d.samplingRate/d.getHighFreq())

        else:
            dialog.showWarning("Input volumes", "Select volumes first", form.root)
            

#===============================================================================
# PICKING
#===============================================================================

class RelionPartDiameter(RelionPartMaskDiameterWizard):  
    _targets = [(ProtRelionAutopickFom, ['particleDiameter'])]
    
    def _getProtocolImages(self, protocol):
        return protocol.inputReferences 


class RelionAutopickParams(EmWizard):  
    _targets = [(ProtRelionAutopick, ['pickingThreshold',
                                      'interParticleDistance'])]
    
    def show(self, form):
        autopickProt = form.protocol
        autopickFomProt = autopickProt.getInputAutopick()
        # Get current values of the properties
#         _, values = self._getInputProtocol(self._targets, autopickProt)
#         threshold, distance = values
#         autopickFomProt.setStepsExecutor() # allow to use runJob
#         autopickFomProt.autopickStep(threshold, distance, '--read_fom_maps')
#         print "Writing Xmipp coordinate files."
#         micFn, coordsDir = autopickFomProt.writeXmippCoords()
        project = autopickProt.getProject()
        micSet = autopickFomProt.getInputMicrographs()
        micfn = micSet.getFileName()
        coordsDir = project.getTmpPath(micSet.getName())
        cleanPath(coordsDir)
        makePath(coordsDir)
        pickerProps = os.path.join(coordsDir, 'picker.conf')
        f = open(pickerProps, "w")
        
        
        
        args = {
          "picker" : "%s relion_autopick" % pw.getScipionScript(),
          "convert" : pw.join('apps', 'pw_convert.py'),
          'coordsDir':coordsDir,
          'micsSqlite': micSet.getFileName(),
          "diameter": autopickFomProt.particleDiameter,
          "threshold": autopickProt.pickingThreshold,
          "apix": micSet.getSamplingRate(),
          'ang': autopickFomProt.angularSampling,
          'lowpass':autopickFomProt.lowpassFilterRefs,
          'ref': 'input_references.star',
          'min_distance': autopickProt.interParticleDistance,
          'protDir': autopickFomProt.getWorkingDir()
          }

        autopickCommand = '%(picker)s --i extra/%%(micrographName).star --o autopick --particle_diameter %(diameter)s --angpix %(apix)s --ref %(ref)s --ang %(ang)s --lowpass %(lowpass)s --threshold %%(threshold) --min_distance %%(ipd) --read_fom_maps'%args
        if autopickFomProt.refsHaveInvertedContrast:
            autopickCommand += ' --invert'
        
        if autopickFomProt.refsCtfCorrected:
            autopickCommand += ' --ctf'
        args['autopickCommand'] = autopickCommand
        f.write("""
        parameters = ipd,threshold
        ipd.value = %(min_distance)s
        ipd.label = Minimum inter-particles distance
        ipd.help = some help
        threshold.value =  %(threshold)s
        threshold.label = Threshold
        threshold.help = some help
        runDir = %(protDir)s
        autopickCommand = %(autopickCommand)s 
        convertCommand = %(convert)s --coordinates --from relion --to xmipp --input  %(micsSqlite)s --output %(coordsDir)s --extra %(protDir)s/extra
        """ % args)
        f.close()
        process = CoordinatesObjectView(autopickProt.getProject(), micfn, coordsDir, autopickFomProt, pickerProps=pickerProps).show()
        process.wait()
        myprops = readProperties(pickerProps)
        form.setVar('pickingThreshold', myprops['threshold.value'])
        form.setVar('interParticleDistance', myprops['ipd.value'])
    