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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
from pyworkflow import VERSION_1_1
from pyworkflow.em import ProtClassify2D, Float
from pyworkflow.protocol.params import PointerParam, IntParam, BooleanParam
from pyworkflow.protocol.constants import LEVEL_ADVANCED
import pyworkflow.utils as pwutils

from ..imagic import ImagicPltFile, ImagicLisFile
from protocol_base import ImagicProtocol


class ImagicProtMSAClassify(ProtClassify2D, ImagicProtocol):
    """ This classification protocol is a post-processor to the MSA.

        It is based on variance-oriented hierarchical ascendant classification program
        (an enhanced Ward-type algorithm).
    """
    _label = 'msa-classify'
    _version = VERSION_1_1
    CLASS_DIR = 'MSA-cls'

    def __init__(self, **kwargs):
        ImagicProtocol.__init__(self, **kwargs)

        self._params = {'cls_dir': self.CLASS_DIR,
                        'msa_cls_img': 'classes'}


# --------------------------- DEFINE param functions --------------------------------------------

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputMSA', PointerParam, label="Input particles", important=True,
                      pointerClass='ImagicProtMSA',
                      help='Input images after MSA')
        form.addParam('numberOfFactors', IntParam, default=15,
                      label='Number of eigenimages to use',
                      help='Select the first N eigenimages to use for classification.\n'
                           'Typically all but the first few are noisy.')
        form.addParam('numberOfClasses', IntParam, default=10,
                      label='Number of classes',
                      help='Desired number of classes.')
        form.addParam('percentIgnore', IntParam, default=0, expertLevel=LEVEL_ADVANCED,
                      label='Percent of images to ignore',
                      help='This option allows for a percentage of the original images to be ignored. '
                           'The last individual images to be merged into a class are set inactive '
                           'in the HAC algorithm. For noisy raw data a value of 15% could be tried, '
                           'for example (this is STATISTICS, remember?).')
        form.addParam('doDownweight', BooleanParam, default=False, expertLevel=LEVEL_ADVANCED,
                      label='Downweight small classes?',
                      help='A consequence of downweighting small classes is that classes with only '
                           'one member will contain only zeroes')
        form.addParam('percentIgnoreBad', IntParam, default=0, expertLevel=LEVEL_ADVANCED,
                      label='Percent of worst class members to ignore',
                      help='Here you get a final chance to polish your classes. Since in the CLS file '
                           'the sequence of images in a class is sorted by their contribution to the '
                           'internal variance of that class, then  we can enhance the class qualities by '
                           'ignoring the last images of each class. The fraction of the images to be ignored '
                           'is what you are supposed to specify here.')

    # --------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):

        self._insertFunctionStep('classifyStep')

        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions --------------------------------------------

    def classifyStep(self):
        """ Run MSA-CL and MSA-SUM from IMAGIC. """

        inputFile = self.inputMSA.get().getParticlesStack()
        inputFileBase = pwutils.removeExt(inputFile)
        inputFileImg = inputFileBase + '.img'
        inputFileHed = inputFileBase + '.hed'

        pwutils.createLink(inputFileImg, self._getTmpPath("particles.img"))
        pwutils.createLink(inputFileHed, self._getTmpPath("particles.hed"))
        inputFn = "tmp/particles"

        if self.doDownweight.get():
            downweight = 'YES'
        else:
            downweight = 'NO'

        self._params.update({'particles': inputFn,
                             'eigs_num': self.numberOfFactors.get(),
                             'cls_num': self.numberOfClasses.get(),
                             'perc_ign': self.percentIgnore.get(),
                             'downweight': downweight,
                             'perc_ign_bad': self.percentIgnoreBad.get()
                             })

        classDir = self._getPath(self.CLASS_DIR)
        pwutils.cleanPath(classDir)
        pwutils.makePath(classDir)

        self.runTemplate('msa/msa-cls.b', self._params)

    def createOutputStep(self):
        """ Create the SetOfClass from the cls file with the images-class
        assignment and the averages for each class.
        """
        particles = self.inputMSA.get().inputParticles.get()
        classes2D = self._createSetOfClasses2D(particles)
        # Load the class assignment file from results
        plt = ImagicPltFile(self._getPath(self.CLASS_DIR, 'class_assignment.plt'))
        self._loadClassInfo(self.numberOfClasses.get())

        # Here we are assuming that the order of the class assignment rows
        # is the same for the input particles and the generated img stack
        classes2D.classifyItems(updateItemCallback=self._updateParticle,
                                updateClassCallback=self._updateClass,
                                itemDataIterator=plt.iterRows())

        self._defineOutputs(outputClasses=classes2D)
        self._defineSourceRelation(particles, classes2D)

    # --------------------------- INFO functions --------------------------------------------

    def _validate(self):
        errors = []
        return errors

    def _citations(self):
        return ['vanHeel1984', 'vanHeel1989', 'Borland1990']

    def _summary(self):
        summary = []
        summary.append('Number of classes: *%s*' % self.numberOfClasses.get())
        summary.append('Number of eigenimages: *%s*' % self.numberOfFactors.get())
        return summary

    def _methods(self):
        msg = "\nInput particles after MSA run were divided into "
        msg += "%s classes by hierarchical ascendant classification (HAC) " % self.numberOfClasses.get()
        msg += "using first %s eigenimages." % self.numberOfFactors.get()
        return [msg]

    # --------------------------- UTILS functions --------------------------------------------

    def _getOutputPath(self, fn):
        """ Return the output file from the run directory and the CLASS dir. """
        return self._getPath(self.CLASS_DIR, fn)

    def getOutputLis(self):
        return self._getOutputPath('classes.lis')

    def _updateParticle(self, item, row):
        _, classNum = row
        item.setClassId(classNum)

    def _updateClass(self, item):
        classId = item.getObjId()
        avgFile = self._getPath(self.CLASS_DIR,
                                self._params['msa_cls_img'] + '_avg.img')
        rep = item.getRepresentative()
        rep.setSamplingRate(item.getSamplingRate())
        rep.setLocation(classId, avgFile)

        item._intraClassVariance = Float(self.varianceDict[classId])
        item._representationQuality = Float(self.quality1Dict[classId])
        item._overallQuality = Float(self.quality2Dict[classId])

    def _loadClassInfo(self, cls):
        fn = self.getOutputLis()
        self.varianceDict, self.quality1Dict, self.quality2Dict = ImagicLisFile(fn, cls).getParams()
