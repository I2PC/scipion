# **************************************************************************
# *
# * Authors:     Josue Gomez Blanco (jgomez@cnb.csic.es)
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

import os
import sys

import pyworkflow.em.metadata as md
from pyworkflow import VERSION_1_1
from pyworkflow.protocol.params import PointerParam, BooleanParam, LabelParam
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pyworkflow.em.protocol import ProtOperateParticles

from convert import writeSetOfParticles, convertBinaryVol, relionToLocation


class ProtRelionSubtract(ProtOperateParticles):
    """
    Subtract volume projections from the experimental particles.
    The particles must have projection alignment in order to
    properly generate volume projections.
    """
    _label = 'subtract projection'
    _lastUpdateVersion = VERSION_1_1
    
    def _initialize(self):
        self._createFilenameTemplates()
    
    def _createFilenameTemplates(self):
        """ Centralize how files are called. """
        myDict = {
                  'input_star': self._getPath('input_particles.star'),
                  'output': self._getExtraPath('output_particles'),
                  'output_star': self._getExtraPath('output_particles.star'),
                  'volume_masked': self._getTmpPath('volume_masked.mrc'),
                  }
        self._updateFilenamesDict(myDict)
    
    #--------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputParticles', PointerParam,
                      pointerClass='SetOfParticles',
                      pointerCondition='hasAlignmentProj',
                      label="Input particles", important=True,
                      help='Select the experimental particles.')
        form.addParam('inputVolume', PointerParam, pointerClass='Volume',
                      label="Input volume",
                      help='Provide the input volume that will be used to '
                           'calculate projections, which will be subtracted '
                           'from the experimental particles. Make sure this '
                           'map was calculated by RELION from the same '
                           'particles as above, and preferably with those '
                           'orientations, as it is crucial that the absolute '
                           'greyscale is the same as in the experimental '
                           'particles.')
        form.addParam('refMask', PointerParam, pointerClass='VolumeMask',
                      label='Reference mask (optional)', allowsNull=True,
                      help="Provide a soft mask where the protein density "
                           "you wish to subtract from the experimental "
                           "particles is white (1) and the rest of the "
                           "protein and the solvent is black (0). "
                           "That is: *the mask should INCLUDE the part of the "
                           "volume that you wish to SUBTRACT.*")
        form.addSection(label='CTF')
        form.addParam('doCTF', BooleanParam, default=True,
                      label='Do CTF-correction?',
                      help='If set to Yes, CTFs will be corrected inside the '
                           'MAP refinement. The resulting algorithm '
                           'intrinsically implements the optimal linear, or '
                           'Wiener filter. Note that input particles should '
                           'contains CTF parameters.')
        form.addParam('haveDataBeenPhaseFlipped', LabelParam,
                      label='Have data been phase-flipped?      '
                            '(Don\'t answer, see help)',
                      help='The phase-flip status is recorded and managed by '
                           'Scipion. \n In other words, when you import or '
                           'extract particles, \nScipion will record whether '
                           'or not phase flipping has been done.\n\n'
                           'Note that CTF-phase flipping is NOT a necessary '
                           'pre-processing step \nfor MAP-refinement in '
                           'RELION, as this can be done inside the internal\n'
                           'CTF-correction. However, if the phases have been '
                           'flipped, the program will handle it.')
        form.addParam('ignoreCTFUntilFirstPeak', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label='Ignore CTFs until first peak?',
                      help='If set to Yes, then CTF-amplitude correction will '
                           'only be performed from the first peak '
                           'of each CTF onward. This can be useful if the CTF '
                           'model is inadequate at the lowest resolution. '
                           'Still, in general using higher amplitude contrast '
                           'on the CTFs (e.g. 10-20%) often yields better '
                           'results. Therefore, this option is not generally '
                           'recommended.')
        
        form.addParallelSection(threads=0, mpi=0)
    
    #--------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):
        self._initialize()
        
        imgSet = self.inputParticles.get()
        partSetId = imgSet.getObjId()
        
        self._insertFunctionStep('convertInputStep', partSetId)
        if self.refMask.get() is not None:
            self._insertFunctionStep('applyMaskStep')
        self._insertFunctionStep('removeStep')
        self._insertFunctionStep('createOutputStep')
    
    #--------------------------- STEPS functions -------------------------------
    def convertInputStep(self, particlesId):
        """ Write the input images as a Xmipp metadata file. 
        particlesId: is only need to detect changes in
        input particles and cause restart from here.
        """
        imgSet = self.inputParticles.get()
        writeSetOfParticles(imgSet, self._getFileName('input_star'),
                            self._getExtraPath())
    
    def applyMaskStep(self):
        import pyworkflow.em.packages.xmipp3 as xmipp3
        from pyworkflow.em.packages.xmipp3.convert import getImageLocation
        
        params = ' -i %s --mult %s -o %s' % (getImageLocation(self.inputVolume.get()),
                                             getImageLocation(self.refMask.get()),
                                             self._getFileName('volume_masked'))
        self.runJob('xmipp_image_operate', params, env=xmipp3.getEnviron())
    
    def removeStep(self):
        volume = self.inputVolume.get()
        if self.refMask.get() is not None:
            volFn = self._getFileName('volume_masked')
        else:
            volFn = convertBinaryVol(volume, self._getExtraPath())
        
        params = ' --i %s --subtract_exp --angpix %0.3f' % (volFn,
                                                            volume.getSamplingRate())
        if self._getInputParticles().isPhaseFlipped():
            params += ' --ctf_phase_flipped'

        if self.doCTF:
            params += ' --ctf'
            if self.ignoreCTFUntilFirstPeak:
                params += ' --ctf_intact_first_peak'

        params += ' --ang %s  --o %s ' % (
            self._getFileName('input_star'),
            self._getFileName('output'))

        try:
            self.runJob('relion_project', params)
        except Exception, ex:
            fn = self._getFileName('output_star')
            if not os.path.exists(fn):
                sys.stderr.write('The file %s was not produced\n' % fn)
                raise ex
            else:
                sys.stderr.write('----Everything OK-----\n')

    def createOutputStep(self):
        imgSet = self._getInputParticles()
        outImgSet = self._createSetOfParticles()
        outImgsFn = self._getFileName('output_star')
         
        outImgSet.copyInfo(imgSet)
        outImgSet.setAlignmentProj()
        outImgSet.copyItems(imgSet,
                            updateItemCallback=self._updateItem,
                            itemDataIterator=md.iterRows(outImgsFn))
        self._defineOutputs(outputParticles=outImgSet)
        self._defineTransformRelation(imgSet, outImgSet)
    
    #--------------------------- INFO functions --------------------------------
    def _validate(self):
        """ Should be overwritten in subclasses to
        return summary message for NORMAL EXECUTION. 
        """
        errors = []
        self.validatePackageVersion('RELION_HOME', errors)
        self._validateDim(self._getInputParticles(), self.inputVolume.get(),
                          errors, 'Input particles', 'Input volume')

        return errors
    
    def _summary(self):
        """ Should be overwritten in subclasses to
        return summary message for NORMAL EXECUTION. 
        """
        summary = []
        if not hasattr(self, 'outputParticles'):
            summary.append("Output is not ready yet.")
        else:
            summary.append('Projections of the masked input volume %s were '
                           'subtracted from original particles %s.' %
                           (self.getObjectTag('inputVolume'),
                            self.getObjectTag('inputParticles')))

        if self._getInputParticles().isPhaseFlipped():
            flipMsg = "Your input images are ctf-phase flipped"
            summary.append(flipMsg)

        return summary
    
    #--------------------------- UTILS functions -------------------------------
    def _updateItem(self, item, row):
        newFn = row.getValue(md.RLN_IMAGE_NAME)
        newLoc = relionToLocation(newFn)
        item.setLocation(newLoc)

    def _getInputParticles(self):
        return self.inputParticles.get()
