from pyworkflow.protocol.params import PointerParam, IntParam
from pyworkflow.em.protocol import ProtParticles
from pyworkflow.em.data import SetOfParticles
from convert import writeSetOfParticles
from pyworkflow.utils import removeExt
import pyworkflow.em as em
import pyworkflow.em.metadata as md
import pyworkflow.object as pwobj

class ProtAlignPartsLmbfgs(ProtParticles):
    """
    """
    _label = 'alignparts_lmbfgs'

    def __init__(self, **kwargs):
        ProtParticles.__init__(self, **kwargs)

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputSet', PointerParam,
                      pointerClass='SetOfParticles',
                      label='Input particles', important=True,
                      help='Select the particles to perform local motion correction on.')

        form.addParam('boxSize', IntParam, label='Box Size',
                      important=True,
                      help= 'Size of particle boxes extracted')

        line = form.addLine('Movie Dimension')
        line.addParam('nx', IntParam, label='x direction',
                      important=True)
        line.addParam('ny', IntParam, label='y direction',
                      important=True)
    # --------------------------- INSERT steps functions ------------------------

    def isInputAutoPicking(self):
        inputSet = self.inputSet.get()
        return (isinstance(inputSet, SetOfParticles) and
                not inputSet.hasAlignmentProj())

    def _insertAllSteps(self):
        self._createFilenameTemplates()
        self._insertFunctionStep('convertInputStep',
                                 self.inputSet.get().getObjId())
        self._insertLocalMotionCorStep()
        self._insertFunctionStep('createOutputStep')

    def _insertLocalMotionCorStep(self):
        """ Prepare command line arguments before calling alignparts_lmbfgs"""
        # Join all keys in a single line, value pairs of the args dict
        args = {}
        self._setArgs(args)
        params = ' '.join(['%s %s' % (k, str(v)) for k, v in args.iteritems()])

        self._insertFunctionStep('runLocalMotionCorStep', params)

    # --------------------------- STEPS functions -------------------------------

    def _createFilenameTemplates(self):
        """ Centralize how files are called"""
        myDict = {
            'input_particles':self._getExtraPath('input_particles.star'),
            'output_aligned':self._getExtraPath('output_aligned.star')
        }
        self._updateFilenamesDict(myDict)

    def convertInputStep(self, inputId):
        """ Input file in star format."""
        inputSet = self.inputSet.get()
        imgStar = self._getFileName('input_particles')

        self.info("Converting set from '%s' into '%s'" %
                  (inputSet.getFileName(), imgStar))

        #write input_particles star file
        writeSetOfParticles(inputSet, imgStar, self._getPath())

    def createOutputStep(self): #output in terms of coordinates / star file
        outputStar = self._getFileName('output_aligned')

    def runLocalMotionCorStep(self, params):
        """ Execute local motion correction with given params. """
        try:
            self.runJob("/nethome/shalini/Desktop/cluster/Downloads/alignparts_starfilehandler.py", params)
            print("blah")
        except:
            print("Not found")

    # --------------------------- INFO functions --------------------------------

    def _validate(self):
        return []

    def _citations(self):
        return []

    def _summary(self):
        return []

    def _methods(self):
        return []

    # --------------------------- UTILS functions --------------------------------------------

    def _setArgs(self, args):
        args.update({'--istar': self._getFileName('input_particles'),
                     '--ostar': self._getFileName('output_aligned'),
                     '--coord': 'coord',
                     '--boxsize': self.boxSize.get(),
                     '--nx': self.nx.get(),
                     '--ny': self.ny.get()
                     })

