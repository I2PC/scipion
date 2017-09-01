# ******************************************************************************
# *
# * Authors:     Amaya Jimenez Moreno (ajimenez@cnb.csic.es)
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
# ******************************************************************************

from pyworkflow.em import ALIGN_NONE, ALIGN_2D
from pyworkflow.em.protocol import ProtAlign2D
import pyworkflow.em.metadata as md
import pyworkflow.protocol.params as params

from convert import createItemMatrix, writeSetOfParticles


class XmippProtGpuCrrCL2D(ProtAlign2D):
    """ Aligns a set of particles using the GPU Correlation algorithm. """
    _label = 'align with GPU Correlation'


    # --------------------------- DEFINE param functions -----------------------
    def _defineAlignParams(self, form):
        form.addParam('useReferenceImages', params.BooleanParam, default=False,
                      label='Use a Set of Reference Images ?',
                      help='If you set to *Yes*, you should provide a '
                           'set of reference images.\n'
                           'If *No*, the default generation is done by '
                           'averaging subsets of the input images.')
        form.addParam('referenceImages', params.PointerParam,
                      condition='useReferenceImages',
                      pointerClass='SetOfParticles', allowsNull=True,
                      label="Reference images",
                      help='Set of images that will serve as class reference')
        form.addParam('maximumShift', params.IntParam, default=10,
                      label='Maximum shift (px):')
        form.addParam('keepBest', params.IntParam, default=1,
                      label='Number of best images:',
                      help='Number of the best images to keep for every class')
        form.addParam('numberOfIterations', params.IntParam, default=10,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Number of iterations:',
                      help='Maximum number of iterations')


    # --------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        """" Mainly prepare the command line for call cuda corrrelation program"""

        # Convert input images if necessary
        if self.useReferenceImages:
            self.refSet = self._getExtraPath('imagesRef.xmd')
            self._insertFunctionStep('convertSetStep', self.refSet, False)

        self.imgsExp = self._getExtraPath('imagesExp.xmd')
        self._insertFunctionStep('convertSetStep', self.imgsExp, True)

        if not self.useReferenceImages:
            #First step: divide the metadata input file to generate
            # a couple of references
            self._params = {'imgsExp': self.imgsExp}
            args = ('-i %(imgsExp)s -n 2')
            self._insertRunJobStep("xmipp_metadata_split", args % self._params)

            #Second step: calculate the means of the previous metadata
            expSet1 = self.imgsExp[0:-4] + '000001.xmd'
            avg1 = self.imgsExp[0:-4] + '_000001'
            expSet2 = self.imgsExp[0:-4] + '000002.xmd'
            avg2 = self.imgsExp[0:-4] + '_000002'
            self._params = {'imgsSet': expSet1,
                            'outputAvg': avg1}
            args = ('-i %(imgsSet)s --save_image_stats %(outputAvg)s -v 0')
            self._insertRunJobStep("xmipp_image_statistics", args % self._params)

            self._params = {'imgsSet': expSet2,
                            'outputAvg': avg2}
            args = ('-i %(imgsSet)s --save_image_stats %(outputAvg)s -v 0')
            self._insertRunJobStep("xmipp_image_statistics", args % self._params)


            #Third step: generate a single metadata with the two previous averages
            self.refSet = self._getExtraPath('refSet.xmd')
            self._params = {'avg1': avg1+'average.xmp',
                            'avg2': avg2+'average.xmp',
                            'outputMd': self.refSet}
            args = ('-i %(avg1)s --set union %(avg2)s -o %(outputMd)s')
            self._insertRunJobStep("xmipp_metadata_utilities", args % self._params)

        #Fouth step: calling program xmipp_cuda_correlation
        self._params = {'imgsRef': self.refSet,
                        'imgsExp': self.imgsExp,
                        'outputFile': 'test.xmd',
                        'extraDir': self._getExtraPath(),
                        'keepBest': self.keepBest.get(),
                        'maxshift': self.maximumShift.get(),
                        'iter': self.numberOfIterations.get(),
                        'outputClassesFile': 'classes.xmd',
                        }
        args = ('-i_ref %(imgsRef)s -i_exp %(imgsExp)s -o %(outputFile)s '
                '--odir %(extraDir)s --keep_best %(keepBest)d '
                '--maxShift %(maxshift)d --classify %(outputClassesFile)s')

        self._insertRunJobStep("xmipp_cuda_correlation", args % self._params)

        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions --------------------------
    def convertSetStep(self, imgs, exp):
        if exp:
            writeSetOfParticles(self.inputParticles.get(), imgs,
                                alignType=ALIGN_NONE)
        else:
            writeSetOfParticles(self.referenceImages.get(), imgs,
                                alignType=ALIGN_NONE)


    def createOutputStep(self):
        """ Store the setOfParticles object
        as result of the protocol.
        """
        return

        particles = self.inputParticles.get()
        # Define the output average
        avgFile = self._getExtraPath('level_00', 'class_classes.stk')
        avg = em.Particle()
        avg.setLocation(1, avgFile)
        avg.copyInfo(particles)
        self._defineOutputs(outputAverage=avg)
        self._defineSourceRelation(self.inputParticles, avg)

        # Generate the Set of Particles with alignment
        alignedSet = self._createSetOfParticles()
        alignedSet.copyInfo(particles)
        alignedSet.setRepresentative(avg)
        alignedSet.copyItems(particles,
                             updateItemCallback=self._createItemMatrix,
                             itemDataIterator=md.iterRows(self.imgsFn,
                                                          sortByLabel=md.MDL_ITEM_ID))
        alignedSet.setAlignment(em.ALIGN_2D)
        self._defineOutputs(outputParticles=alignedSet)
        self._defineSourceRelation(self.inputParticles, alignedSet)

    # --------------------------- INFO functions -------------------------------
    def _validate(self):
        return

        errors = []
        if self.numberOfMpi <= 1:
            errors.append('Mpi needs to be greater than 1.')
        if self.useReferenceImages:
            if self.referenceImage.hasValue():
                refImage = self.referenceImage.get()
                [x1, y1, z1] = refImage.getDim()
                [x2, y2, z2] = self.inputParticles.get().getDim()
                if x1 != x2 or y1 != y2 or z1 != z2:
                    errors.append('The input images and the reference image '
                                  'have different sizes')
            else:
                errors.append("Please, enter a reference image")
        return errors

    def _summary(self):
        return
        summary = []
        if not hasattr(self, 'outputParticles'):
            summary.append("Output alignment not ready yet.")
        else:
            summary.append("Input Particles: %s"
                           % self.inputParticles.get().getSize())
            if self.useReferenceImages:
                summary.append("Aligned with reference image: %s"
                               % self.referenceImage.get().getNameId())
            else:
                summary.append("Aligned with no reference image.")
        return summary

    def _citations(self):
        return ['Sorzano2010a']

    def _methods(self):
        return
        methods = []
        if not hasattr(self, 'outputParticles'):
            methods.append("Output alignment not ready yet.")
        else:
            if self.useReferenceImages:
                methods.append(
                    "We aligned images %s with respect to the reference image "
                    "%s using CL2D [Sorzano2010a]"
                    % (self.getObjectTag('inputParticles'),
                       self.getObjectTag('referenceImage')))
            else:
                methods.append(
                    "We aligned images %s with no reference using CL2D "
                    "[Sorzano2010a]" % self.getObjectTag('inputParticles'))
            methods.append(" and produced %s images."
                           % self.getObjectTag('outputParticles'))
        return methods

    def _createItemMatrix(self, item, row):
        createItemMatrix(item, row, align=ALIGN_2D)
