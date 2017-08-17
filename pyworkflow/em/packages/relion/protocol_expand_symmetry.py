# **************************************************************************
# *
# * Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk)
# *
# * MRC Laboratory of Molecular Biology, MRC-LMB
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

from pyworkflow.protocol.params import StringParam
from pyworkflow.em.protocol import ProtProcessParticles
from pyworkflow.em import ALIGN_PROJ
from convert import (writeSetOfParticles, getVersion, V1_3,
                     V1_4, readSetOfParticles, setRelionAttributes)
import pyworkflow.em.metadata as md

 
class ProtRelionExpandSymmetry(ProtProcessParticles):
    """
    Given an input set of particles with angular assignment,
    expand the set by applying a pseudo-symmetry.

    Be aware that input symmetry values follow Xmipp conventions as described in:
    http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Symmetry
    """
    _label = 'expand symmetry'

    @classmethod
    def isDisabled(cls):
        return getVersion() in [V1_3, V1_4]

    #--------------------------- DEFINE param functions -------------------------------------
    def _defineProcessParams(self, form):
        form.addParam('symmetryGroup', StringParam, default="c1",
                      label='Symmetry group',
                      help="See http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Symmetry"
                           " for a description of the symmetry groups format in Xmipp.\n")
        form.addParallelSection(threads=0, mpi=0)

    #--------------------------- INSERT steps functions -------------------------------------
    def _insertAllSteps(self):
        imgsFn = self._getPath('input_particles.star')
        self._insertFunctionStep('convertInputStep', imgsFn)
        self._insertFunctionStep('expandSymmetryStep', imgsFn)
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions --------------------------------------------

    def convertInputStep(self, outputFn):
        """ Create a metadata with the images and geometrical information. """
        writeSetOfParticles(self.inputParticles.get(), outputFn, self._getPath())

    #--------------------------- STEPS functions --------------------------------------------
    def expandSymmetryStep(self, imgsFn):
        outImagesMd = self._getExtraPath('expanded_particles.star')
        args = " --i %s --sym %s --o %s" % (imgsFn, self.symmetryGroup.get(),
                                            outImagesMd)
        self.runJob("relion_particle_symmetry_expand", args)

    def createOutputStep(self):
        imgSet = self.inputParticles.get()
        partSet = self._createSetOfParticles()
        partSet.copyInfo(imgSet)
        outImagesMd = self._getExtraPath('expanded_particles.star')

        mdOut = md.MetaData(outImagesMd)
        mdOut.removeLabel(md.RLN_IMAGE_ID)  # remove repeating rlnImageId in mdOut
        mdOut.write(outImagesMd, md.MD_OVERWRITE)

        readSetOfParticles(outImagesMd, partSet, alignType=ALIGN_PROJ,
                           postprocessImageRow=self._postprocessImageRow)

        self._defineOutputs(outputParticles=partSet)
        self._defineSourceRelation(imgSet, partSet)

    #--------------------------- INFO functions ---------------------------------------------
    def _summary(self):
        summary = []
        if not hasattr(self, 'outputParticles'):
            summary.append("Output particles not ready yet.")
        else:
            summary.append("Symmetry: %s" % self.symmetryGroup.get())
        return summary
    
    def _validate(self):
        errors = []
        self.validatePackageVersion('RELION_HOME', errors)

        return errors
        
    def _citations(self):
        return []
    
    def _methods(self):
        methods = []

        return methods

    #--------------------------- Utils functions --------------------------------------------
    def _postprocessImageRow(self, img, imgRow):
        setRelionAttributes(img, imgRow, md.RLN_MLMODEL_GROUP_NAME)
