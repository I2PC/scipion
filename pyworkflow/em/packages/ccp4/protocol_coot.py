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

import pyworkflow.utils as pwutils
#import pyworkflow.protocol.params as params
#import pyworkflow.em as em
from pyworkflow import VERSION_1_2
from pyworkflow.utils.properties import Message
from pyworkflow.em import PointerParam
from pyworkflow.em.protocol import EMProtocol
from convert import (getProgram)



class CCP4ProtCoot(EMProtocol):
    """Coot is aninteractive graphical application for
macromolecular model building, model completion
and validation.
"""
    _label = 'refinement'
    _version = VERSION_1_2

    # --------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        #set of volumes
        form.addParam('inputVolumes', PointerParam, pointerClass="SetOfVolumes",
                      label='Input Volume/s',
                      help="Set of volumes to process")

        #set of PDF files
        form.addParam('setOfPdbFiles', PointerParam, pointerClass="SetOfPdbFiles",
                      label='Input PDB/s',
                      help="Set of PDB files to process")

        # --------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):
        #test loop over inputVol
        f = open("/tmp/kk","w")
        for vol in self.inputVolumes:
            f.write("vol",vol.get().getLocation())
        #test loop over pdbFiles
        for pdbFile in self.setOfPdbFiles:
            f.write("pdbFile",pdbFile.get().getFileName())

        convertId = self._insertFunctionStep('convertInputStep')
        runCootId = self._insertFunctionStep('runCootStep', prerequisites=[convertId])
        self._insertFunctionStep('createOutputStep', prerequisites=[runCootId])
        f.close()
    # --------------------------- STEPS functions ---------------------------------------------------

    def convertInputStep(self):
        """ convert 3Dmaps to MRC '.mrcs' format
        """
        #micDir = self.getMicrographsDir()  # put output and mics in extra dir
        #pwutils.makePath(micDir)
        pass

    def runCootStep(self, micName, refStack, args):
        # We convert the input micrograph on demand if not in .mrc
        #runCoot(micName, refStack, self.getMicrographsDir(), args, env=self._getEnviron())
        pass

    def createOutputStep(self):
        #micSet = self.getInputMicrographs()
        #ih = em.ImageHandler()
        #coordSet = self._createSetOfCoordinates(micSet)
        #if self.boxSize and self.boxSize > 0:
        #    coordSet.setBoxSize(self.boxSize.get())
        #else:
        #    coordSet.setBoxSize(self.inputReferences.get().getXDim or 100)
        pass

    # --------------------------- INFO functions --------------------------------------------
    def _validate(self):
        errors = []
        # Check that the program exists
        program = getProgram(os.environ['COOT'])
        if program is None:
            errors.append("Missing variables COOT and/or CCP4_HOME")
        elif not os.path.exists(program):
            errors.append("Binary '%s' does not exists.\n" % program)

        # If there is any error at this point it is related to config variables
        if errors:
            errors.append("Check configuration file: ~/.config/scipion/scipion.conf")
            errors.append("and set COOT and CCP4_HOME variables properly.")
            if program is not None:
                errors.append("Current values:")
                errors.append("CCP4_HOME = %s" % os.environ['CCP4_HOME'])
                errors.append("COOT = %s" % os.environ['COOT'])

        return errors

    def _summary(self):
        #Think on how to update this summary with created PDB
        summary = []
        summary.append("Number of input micrographs: %d"
                       % self.getInputMicrographs().getSize())
        if (self.getOutputsSize() > 0):
            summary.append("Number of particles picked: %d" % self.getCoords().getSize())
            summary.append("Particle size: %d px" % self.getCoords().getBoxSize())
            summary.append("Threshold min: %0.2f" % self.threshold)
        else:
            summary.append(Message.TEXT_NO_OUTPUT_CO)
        return summary

    def _methods(self):
        methodsMsgs = []
        if self.getInputMicrographs() is None:
            return ['Input micrographs not available yet.']
        methodsMsgs.append("Input micrographs %s."
                           % (self.getObjectTag(self.getInputMicrographs())))

        if self.getOutputsSize() > 0:
            output = self.getCoords()
            methodsMsgs.append("%s: User picked %d particles with a particle "
                               "size of %d px and minimal threshold %0.2f."
                               % (self.getObjectTag(output), output.getSize(),
                                  output.getBoxSize(), self.threshold))
        else:
            methodsMsgs.append(Message.TEXT_NO_OUTPUT_CO)

        return methodsMsgs

    def _citations(self):
        return ['Emsley_2004']

    # --------------------------- UTILS functions --------------------------------------------------
    def getArgs(self, threshold=True, mindist=True):
        """ Return the Gautomatch parameters for picking one micrograph.
         The command line will depends on the protocol selected parameters.
         Params:
            micFn: micrograph filename
            refStack: filename with the references stack (.mrcs)
        """
        args = ' --apixM %0.2f' % self.inputMicrographs.get().getSamplingRate()
        args += ' --ang_step %d' % self.angStep
        args += ' --diameter %d' % (2 * self.particleSize.get())
        args += ' --lp %d' % self.lowPass
        args += ' --hp %d' % self.highPass
        args += ' --gid %d' % self.GPUId

        if self.inputReferences.get():
            args += ' --apixT %0.2f' % self.inputReferences.get().getSamplingRate()

        if not self.invertTemplatesContrast:
            args += ' --dont_invertT'

        if threshold:
            args += ' --cc_cutoff %0.2f' % self.threshold

        if not self.advanced:
            args += ' --speed %d' % self.speed
            args += ' --boxsize %d' % self.boxSize
            if mindist:
                args += ' --min_dist %d' % self.minDist
            args += ' --lsigma_cutoff %0.2f' % self.localSigmaCutoff
            args += ' --lsigma_D %d' % self.localSigmaDiam
            args += ' --lave_max %0.2f' % self.localAvgMax
            args += ' --lave_min %0.2f' % self.localAvgMin
            args += ' --lave_D %d' % self.localAvgDiam

        if self.preFilt:
            args += ' --do_pre_filter --pre_lp %d' % self.prelowPass
            args += ' --pre_hp %d' % self.prehighPass

        if self.exclusive:
            if self.inputBadCoords.get():
                args += ' --exclusive_picking --excluded_suffix _rubbish.star'
            if self.inputDefects.get():
                args += ' --global_excluded_box %s' % self._getExtraPath('micrographs/defects.star')

        if self.writeCC:
            args += ' --write_ccmax_mic'
        if self.writeFilt:
            args += ' --write_pref_mic'
        if self.writeBg:
            args += ' --write_bg_mic'
        if self.writeBgSub:
            args += ' --write_bgfree_mic'
        if self.writeSigma:
            args += ' --write_lsigma_mic'
        if self.writeMsk:
            args += ' --write_mic_mask'

        return args


    def getOutputName(self, fn, key):
        """ Give a key, append the mrc extension
        and prefix the protocol working dir.
        """
        template = pwutils.removeBaseExt(fn) + key + '.mrc'

        return pwutils.join(self.getMicrographsDir(), template)
