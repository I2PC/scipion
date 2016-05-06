# **************************************************************************
# *
# * Authors:     Grigory Sharov (sharov@igbmc.fr)
# *              J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************

from os.path import join, exists

from pyworkflow.protocol.params import IntParam, PointerParam, EnumParam, FloatParam
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pyworkflow.em.convert import ImageHandler
import pyworkflow.utils as pwutils

from ..constants import MODULATION
from protocol_base import ImagicProtocol


class ImagicProtMSA(ImagicProtocol):
    """Multivariate Statistical Analysis module of IMAGIC.

    It calculates eigenimages (eigenvectors) and eigenvalues of
    a set of input aligned images using an iterative eigenvector
    algorithm optimized for (extremely) large data sets.

    """
    _label = 'msa'
    MSA_DIR = 'MSA'

    def __init__(self, **kwargs):
        ImagicProtocol.__init__(self, **kwargs)

        self._params = {'eigen_img': join(self.MSA_DIR, 'eigen_img'),
                        'msa_pixvec_coord': join(self.MSA_DIR, 'msa_pixvec_coord'),
                        'msa_eigen_pixel': join(self.MSA_DIR, 'msa_eigen_pixel')
                        }

# --------------------------- DEFINE param functions --------------------------------------------

    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputParticles', PointerParam, label="Input particles", important=True,
                      pointerClass='SetOfParticles',
                      help='Select the input particles to perform MSA.')
        form.addParam('distanceType', EnumParam, default=MODULATION, choices=['EUCLIDIAN', 'CHISQUARE', 'MODULATION'],
                      label='MSA distance type',  expertLevel=LEVEL_ADVANCED,
                      help='Select general metric for square distance between images:\n\n'
                           'a. Euclidian metric (Principal Components Analysis: PCA)\n'
                           'b. Chi-square metric (Correspondence Analysis: CA)\n'
                           'c. Modulation metric (Modulation Analysis: MA, recommended)')
        form.addParam('numberOfFactors', IntParam, default=25,
                      label='Number of factors (eigenimages)',
                      help='A 64x64 image can be expressed as a vector of 4096 dimensions. '
                           'In this step, we will reduce this number of dimensions to the number of factors '
                           'specified here. These factors will represent the largest systematic variations '
                           'in the data.\n\nThe number of eigenimages that should be used depends on the '
                           'complexity of the input data.')
        form.addParam('numberOfIterations', IntParam, default=25,
                      label='Number of iterations',
                      help='The calculation of the eigenimages and the related eigenvalues '
                           'are determined in an iterative process.\n'
                           'Please give the number of iterations wanted.\n\n'
                           'NOTE: The iterations will stop automatically if the eigenimage'
                           '(eigenvector eigenvalue) calculations are converging')
        form.addParam('overcorrectionFactor', FloatParam, default=0.8, expertLevel=LEVEL_ADVANCED,
                      label='Overcorrection factor [0 - 0.9]',
                      help='The overcorrection factor is a very important parameter '
                           'in the MSA program. It determines the convergence speed '
                           'of the Eigenvector Eigenvalue algorithm. However, if '
                           'a too large overcorrection is chosen, the algorithm may '
                           'start oscillating. Oscillations of the algorithm may be '
                           'observed in the plot of the sum of the eigenvalues versus '
                           'iteration number which is part of the output of this '
                           'program. Divergence may thus only be detected a posteriori.\n\n'
                           'The accepted values for OVER_CORRECTION lie between 0 and 0.9.')
        form.addParam('maskType', EnumParam,
                      choices=['circular', 'object'], default=0,
                      display=EnumParam.DISPLAY_HLIST,
                      label='Mask type',
                      help='Select which type of mask do you want to apply. '
                           'Only the pixels beneath this mask will be analyzed. '
                           'In the simplest case, a circular mask can be used. '
                           'Alternatively, a custom mask can be used '
                           'which follows the contour of the particle (but not too tightly).')
        form.addParam('radius', IntParam, default=-1,
                      label='Mask radius (px)', condition='maskType==0',
                      help='If -1, the entire image (in pixels) will be considered.')
        form.addParam('maskImage', PointerParam, label="Mask image", condition='maskType==1',
                      pointerClass='Mask',
                      help="Select a mask file")

        form.addParallelSection(threads=0, mpi=1)

    # --------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('convertInputStep')

        if self.maskType > 0:
            self._insertFunctionStep('convertMaskStep', self.maskImage.get().getObjId())
        else:
            self._insertFunctionStep('createMaskStep')

        self._insertFunctionStep('msaStep')

    # --------------------------- STEPS functions --------------------------------------------

    def convertInputStep(self):
        # we need to put all images into a single stack to ease the call of imagic programs
        # TODO: skip writeStack, convert directly via e2proc2d.py
        inputParticles = self.inputParticles.get()
        tmpStack = self._getTmpPath('input_particles.stk')
        inputParticles.writeStack(tmpStack, applyTransform=True)
        ImageHandler().convert(tmpStack, self.getParticlesStack())

    def convertMaskStep(self, maskType):
        """ Convert the input mask to Imagic. """
        if maskType > 0:  # mask from file
            maskFn = self._getTmpPath('mask.img')
            ImageHandler().convert(self.maskImage.get(), maskFn)

    def createMaskStep(self):
        """ Create a circular mask in Imagic format. """
        inputParticles = self.inputParticles.get()
        radius = self.radius.get()

        if self.maskType.get() == 0:
            if radius < 0:  # usually -1
                radiusMask = inputParticles.getDim()[0] / 2  # use half of input dim
            else:
                radiusMask = radius
            outMask = self._getTmpPath('mask.img')
            ih = ImageHandler()
            ih.createCircularMask(radiusMask, inputParticles.getFirstItem(), outMask)

    def msaStep(self):
        """ Run MSA on input particles. """

        distances = ['EUCLIDIAN', 'CHISQUARE', 'MODULATION']
        distance_name = distances[self.distanceType.get()]

        self._params.update({'msa_dir': self.MSA_DIR,
                             'msa_distance': distance_name,
                             'num_factors': self.numberOfFactors.get(),
                             'num_iter': self.numberOfIterations.get(),
                             'overcorrectionFactor': self.overcorrectionFactor.get(),
                             'mpi_procs': self.numberOfMpi.get()
                             })

        msaDir = self._getPath(self.MSA_DIR)
        if exists(msaDir):
            pwutils.cleanPath(msaDir)
        pwutils.makePath(msaDir)

        self.runTemplate('msa/msa-run.b', self._params)

    # --------------------------- INFO functions --------------------------------------------

    def _validate(self):
        errors = []

        if self.maskImage.get():
            # check pixel size and image size
            pixel_inp = self.inputParticles.get().getSamplingRate()
            pixel_mask = self.maskImage.get().getSamplingRate()
            if pixel_inp != pixel_mask:
                errors.append('Pixel sizes of input images and mask should be the same!')

            if self.maskImage.get().getDim()[0] != self.inputParticles.get().getDim()[0]:
                errors.append('Image size of input images and mask is not the same!')

        # check overcorrection factor value
        value = round(self.overcorrectionFactor.get(), 1)
        if not 0.0 <= value <= 0.9:
            errors.append('Overcorrection factor value should be in [0 - 0.9] range!')

        # for compatibility with old IMAGIC versions
        if not (self.numberOfFactors.get() <= 64 and self.numberOfIterations.get() <= 64):
            errors.append('For compatibility with old IMAGIC versions, '
                          'number of eigenimages and iterations should be <= 64')

        # check radius size vs input particles
        radiusmax = self.inputParticles.get().getDim()[0] / 2
        if not self.radius.get() <= radiusmax:
            errors.append('Radius cannot be bigger than half-size of input images!')

        return errors

    def _citations(self):
        return ['Borland1990']

    def _summary(self):
        summary = []
        summary.append('This protocol generates only eigenimages (factors), that will be used later for'
                       ' MSA-based classification.')

        if self.distanceType == 0:
            summary.append('Distance type: *Euclidian*')
        if self.distanceType == 1:
            summary.append('Distance type: *ChiSquare*')
        if self.distanceType == 2:
            summary.append('Distance type: *Modulation*')

        summary.append('Number of factors: *%s*' % self.numberOfFactors)

        if self.maskType == 0:  # circular mask
            if self.radius == -1:
                summary.append('Mask: *Circular, of radius 1/2 image dimension*')
            else:
                summary.append('Mask: *Circular, of radius: %s*' % self.radius)
        else:  # custom mask
            summary.append('Mask: *Custom file*')

        return summary

    def _methods(self):

        msg = "\nInput particles %s were subjected to MSA using " % self.getObjectTag('inputParticles')

        if self.distanceType == 0:
            msg += "Euclidian metric, "
        if self.distanceType == 1:
            msg += "ChiSquare metric, "
        if self.distanceType == 2:
            msg += "Modulation metric, "

        msg += "computing %s factors, and using a " % self.numberOfFactors

        if self.maskType == 0:  # circular mask
            if self.radius == -1:
                msg += "circular mask of radius half the image dimension."
            else:
                msg += "circular mask of radius %s pixels." % self.radius
        else:  # custom mask
            msg += "custom mask %s." % self.getObjectTag('maskImage')

        return [msg]

    # --------------------------- UTILS functions --------------------------------------------

    def getParticlesStack(self):
        return self._getPath('input_particles.img')

    def _getOutputPath(self, fn):
        """ Return the output file from the run directory and the MSA dir. """
        return self._getPath(self.MSA_DIR, fn)

    def getOutputEigenImages(self):
        return self._getOutputPath('eigen_img.img')

    def getOutputLis(self):
        return self._getOutputPath('msa.lis')

    def getOutputPlt(self):
        return self._getOutputPath('msa.plt')
