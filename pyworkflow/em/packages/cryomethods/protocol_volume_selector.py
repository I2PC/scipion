# **************************************************************************
# *
# * Authors:     Josue Gomez Blanco (josue.gomez-blanco@mcgill.ca)
# *              Javier Vargas Balbuena (javier.vargasbalbuena@mcgill.ca)
# *
# * Department of Anatomy and Cell Biology, McGill University
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
import pyworkflow.em as em
from pyworkflow.em.packages.cryomethods.protocol_base import ProtocolBase
from convert import writeSetOfParticles


class ProtInitialVolumeSelector(ProtocolBase):
    """
    Protocol to obtain a better initial volume using as input a set of
    volumes and particles. The protocol uses a small subset (usually 1000/2000)
    particles from the input set of particles to estimate a better and reliable
    volume(s) to use as initial volume in an automatic way.
    """

    _label = 'volume selector'

    IS_VOLSELECTOR = True

    def __init__(self, **args):
        ProtocolBase.__init__(self, **args)

    # -------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        self._defineInputParams(form)
        self._defineCTFParams(form)
        self._defineOptimizationParams(form)
        self._defineSamplingParams(form)
        self._defineAdditionalParams(form)

    # -------------------------- STEPS functions -------------------------------
    def convertInputStep(self, resetDeps):
        """ Create the input file in STAR format as expected by Relion.
        If the input particles comes from Relion, just link the file.
        Params:
            particlesId, volumesId: use this parameters just to force redo of
            convert if either the input particles and/or input volumes are
            changed.
        """
        self._imgFnList = []
        imgSet = self._getInputParticles()
        imgStar = self._getFileName('input_star')

        subset = em.SetOfParticles(filename=":memory:")

        newIndex = 1
        for img in imgSet.iterItems(orderBy='RANDOM()', direction='ASC'):
            self._scaleImages(newIndex, img)
            newIndex += 1
            subset.append(img)
            subsetSize = self.subsetSize.get()
            minSize = min(subsetSize, imgSet.getSize())
            if subsetSize   > 0 and subset.getSize() == minSize:
                break
        writeSetOfParticles(subset, imgStar, self._getExtraPath(),
                            alignType=em.ALIGN_NONE,
                            postprocessImageRow=self._postprocessParticleRow)
        self._convertInput(subset)
        self._convertRef()

    def runClassifyStep(self, params):
        """ Execute the relion steps with the give params. """
        params += ' --j %d' % self.numberOfThreads.get()
        self.runJob(self._getProgram(), params)

    def createOutputStep(self):
        # create a SetOfVolumes and define its relations
        volumes = self._createSetOfVolumes()
        self._fillVolSetFromIter(volumes, self._lastIter())

        self._defineOutputs(outputVolumes=volumes)
        self._defineSourceRelation(self.inputVolumes, volumes)

    # --------------------------- INFO functions -------------------------------
    def _validateNormal(self):
        """ Should be overwritten in subclasses to
        return summary message for NORMAL EXECUTION.
        """
        return []

    def _summaryNormal(self):
        """ Should be overwritten in subclasses to
        return summary message for NORMAL EXECUTION.
        """
        return []

    def _methods(self):
        """ Should be overwritten in each protocol.
        """
        return []

    # -------------------------- UTILS functions ------------------------------
    def _setSamplingArgs(self, args):
        """ Set sampling related params. """
        args['--healpix_order'] = self.angularSamplingDeg.get()
        args['--offset_range'] = self.offsetSearchRangePix.get()
        args['--offset_step'] = (self.offsetSearchStepPix.get() *
                                 self._getSamplingFactor())

    def _getResetDeps(self):
        return "%s, %s, %s" % (self._getInputParticles().getObjId(),
                               self.inputVolumes.get().getObjId(),
                               self.targetResol.get())

    def _getClassId(self, volFile):
        result = None
        s = self._classRegex.search(volFile)
        if s:
            result = int(s.group(1)) # group 1 is 2 digits class number
        return self.volDict[result]
