# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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

from os.path import exists, join

from pyworkflow.em import ProtClassify2D, Particle, Class2D, ImageHandler
from pyworkflow.protocol.params import PointerParam, IntParam, BooleanParam
from pyworkflow.protocol.constants import LEVEL_ADVANCED
import pyworkflow.utils as pwutils

from ..imagic import ImagicPltFile, MsaFile
from protocol_base import ImagicProtocol


class ImagicProtMSAClassify(ProtClassify2D, ImagicProtocol):
    """ This classification protocol is a post-processor to the MSA.

        It is based on variance-oriented hierarchical ascendant classification program
        (an enhanced Ward-type algorithm).
    """
    _label = 'msa-classify'

    def __init__(self, **kwargs):
        ImagicProtocol.__init__(self, **kwargs)
        self._classDir = 'MSA-cls'

        self._params = {'cls_dir': self._classDir,
                        'msa_cls_img': 'classes'}

    def getClassDir(self):
        return self._classDir

# --------------------------- DEFINE param functions --------------------------------------------

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputParticles', PointerParam, label="Input particles", important=True,
                      pointerClass='SetOfParticles',
                      pointerCondition='hasMSA',
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

        # self._insertFunctionStep('checkInputStep', self.inputParticles.get().strId())

        self._insertFunctionStep('classifyStep')

        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions --------------------------------------------

    # def checkInputStep(self, strId):
    #    pass

    def classifyStep(self):
        """ Run MSA-CL and MSA-SUM from IMAGIC. """

        inputFile = self.inputParticles.get().getFirstItem().getFileName()
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

        classDir = self._getPath(self._classDir)
        if exists(classDir):
            pwutils.cleanPath(classDir)
        pwutils.makePath(classDir)

        self.runTemplate('msa/msa-cls.b', self._params)

    def createOutputStep(self):
        """ Create the SetOfClass from the cls file with the images-class
        assignment and the averages for each class.
        """
        particles = self.inputParticles.get()
        sampling = particles.getSamplingRate()
        classes2D = self._createSetOfClasses2D(particles)

        for classId in range(1, self.numberOfClasses.get() + 1):
            class2D = Class2D()
            class2D.setObjId(classId)

            classDir = self._getPath(self._classDir)
            outputRootname = join(classDir, self._params.get('msa_cls_img'))

            avgImg = Particle()
            avgImg.setSamplingRate(sampling)

            avgFile = outputRootname + '_avg.img'
            avgImg.setLocation(classId, avgFile)

            class2D.setRepresentative(avgImg)
            class2D.copyInfo(particles)
            classes2D.append(class2D)

            pltClassFile = join(classDir, 'class_assignment.plt')
            plt = ImagicPltFile(pltClassFile)

            for imgId in plt.Values(classId):
                img = particles[imgId]
                if img is None:
                    print ">>> WARNING: Particle with id '%d' found in plt file, but not in input particles." % imgId
                else:
                    class2D.append(img)

            # Update information about number of particles in the class
            classes2D.update(class2D)

        classDir = self._getPath(self._classDir)
        lis = MsaFile()
        lis.filename.set(classDir + '/' + self._params.get('msa_cls_img') + '.lis')

        self._defineOutputs(outputClasses=classes2D, lisFile=lis)
        self._defineSourceRelation(self.inputParticles, classes2D)
        self._defineSourceRelation(self.inputParticles, lis)

    # --------------------------- INFO functions --------------------------------------------

    def _validate(self):
        errors = []

        particlesFn = self.inputParticles.get().getFirstItem().getFileName()
        eigName = 'MSA/eigen_img.img'
        eigPath = pwutils.findRootFrom(particlesFn, eigName)

        if eigPath is None:
            errors.append('Eigenimages file not found! Did you run MSA protocol on input particles?')
        else:
            eigNum = self.numberOfFactors.get()
            eigFn = join(eigPath, eigName)
            ih = ImageHandler()
            if ih.isImageFile(eigFn):
                _, _, _, n = ih.getDimensions(eigFn)
                if not eigNum <= n:
                    errors.append('Maximum number of eigenimages in MSA run was %d' % n)
        return errors

    def _citations(self):
        return ['vanHeel1984', 'vanHeel1989', 'Borland1990']

    def _summary(self):
        summary = []
        summary.append('Number of classes: *%s*' % self.numberOfClasses.get())
        summary.append('Number of eigenimages: *%s*' % self.numberOfFactors.get())
        return summary

    def _methods(self):
        msg = "\nInput particles %s " % self.getObjectTag('inputParticles')
        msg += "were divided into %s classes by MSA-based classification " % self.numberOfClasses.get()
        msg += "using first %s eigenimages." % self.numberOfFactors.get()
        return [msg]
