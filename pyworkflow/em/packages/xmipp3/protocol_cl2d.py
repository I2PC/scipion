# ******************************************************************************
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
# ******************************************************************************

import pyworkflow.em as em
import pyworkflow.em.metadata as md
from os.path import join, dirname, exists
from glob import glob

import pyworkflow.protocol.params as param
import pyworkflow.protocol.constants as const
from pyworkflow.em.protocol import ProtClassify2D, SetOfClasses2D
from pyworkflow.utils.path import cleanPath, makePath
from pyworkflow.em.packages.xmipp3.convert import writeSetOfParticles, \
    createItemMatrix, writeSetOfClasses2D, xmippToLocation, rowToAlignment

# Comparison methods enum
CMP_CORRELATION = 0
CMP_CORRENTROPY = 1

# Clustering methods enum
CL_CLASSICAL = 0
CL_ROBUST = 1

# Classes keys
CLASSES = ''
CLASSES_CORE = '_core'
CLASSES_STABLE_CORE = '_stable_core'

# Suggested number of images per class
IMAGES_PER_CLASS = 200
        
class XmippProtCL2D(ProtClassify2D):
    """ Classifies a set of images using a clustering algorithm 
            to subdivide the original dataset
             into a given number of classes """
    
    _label = 'cl2d'
    
    def __init__(self, **args):
        ProtClassify2D.__init__(self, **args) 
        if self.numberOfMpi.get() < 2:
            self.numberOfMpi.set(2)       

    #--------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputParticles', param.PointerParam,
                      label="Input images",
                      important=True, pointerClass='SetOfParticles',
                      help='Select the input images to be classified.')        
        form.addParam('numberOfClasses', param.IntParam, default=64,
                      label='Number of classes:',
                      help='Number of classes (or references) to be generated.')
        form.addParam('randomInitialization', param.BooleanParam, default=True,
                      expertLevel=const.LEVEL_ADVANCED,
                      label='Random initialization of classes:',
                      help="Initialize randomly the first classes. If you "
                           "don't initialize randomly, you must supply a set "
                           "of initial classes")
        form.addParam('initialClasses', param.PointerParam,
                      label="Initial classes",
                      condition="not randomInitialization",
                      pointerClass='SetOfClasses2D, SetOfAverages',
                      help='Set of initial classes to start the classification')
        form.addParam('numberOfInitialClasses', param.IntParam, default=4,
                      expertLevel=const.LEVEL_ADVANCED,
                      label='Number of initial classes:',
                      condition="randomInitialization",
                      help='Initial number of classes used in the first level.')
        form.addParam('numberOfIterations', param.IntParam, default=10,
                      expertLevel=const.LEVEL_ADVANCED,
                      label='Number of iterations:',
                      help='Maximum number of iterations within each level.')
        form.addParam('comparisonMethod', param.EnumParam,
                      choices=['correlation', 'correntropy'],
                      label="Comparison method", default=CMP_CORRELATION,
                      expertLevel=const.LEVEL_ADVANCED,
                      display=param.EnumParam.DISPLAY_COMBO,
                      help='Use correlation or correntropy')
        form.addParam('clusteringMethod', param.EnumParam,
                      choices=['classical', 'robust'],
                      label="Clustering method", default=CL_CLASSICAL,
                      expertLevel=const.LEVEL_ADVANCED,
                      display=param.EnumParam.DISPLAY_COMBO,
                      help='Use the classical clustering criterion or the '
                           'robust')
        form.addParam('extraParams', param.StringParam,
                      expertLevel=const.LEVEL_ADVANCED,
                      label='Additional parameters',
                      help='Additional parameters for classify_CL2D: \n'
                           ' --verbose, --corrSplit, ...')

        form.addSection(label='Core analysis')
        form.addParam('doCore', param.BooleanParam, default=True,
                      label='Perform core analysis',
                      help='An image belongs to the core if it is close (see '
                           'Junk Zscore and PCA Zscore) to the class center')
        form.addParam('thZscore', param.FloatParam, default=3,
                      label='Junk Zscore', expertLevel=const.LEVEL_ADVANCED,
                      condition='doCore',
                      help='Which is the average Z-score to be considered as '
                           'junk. Typical values go from 1.5 to 3. For the '
                           'Gaussian distribution 99.5% of the data is '
                           'within a Z-score of 3. Lower Z-scores reject more '
                           'images. Higher Z-scores accept more images.')
        form.addParam('thPCAZscore', param.FloatParam, default=3,
                      condition='doCore', expertLevel=const.LEVEL_ADVANCED,
                      label='PCA Zscore',
                      help='Which is the PCA Z-score to be considered as junk. '
                           'Typical values go from 1.5 to 3. For the Gaussian '
                           'distribution 99.5% of the data is within a '
                           'Z-score of 3. Lower Z-scores reject more images. '
                           'Higher Z-scores accept more images.')
        form.addParam('doStableCore', param.BooleanParam, default=True,
                      condition='doCore', label='Perform stable core analysis',
                      help='Two images belong to the stable core if they have '
                           'been essentially together along the classification '
                           'process')
        form.addParam('tolerance', param.IntParam, default=1, label='Tolerance',
                      expertLevel=const.LEVEL_ADVANCED,
                      condition='doCore and doStableCore',
                      help='An image belongs to the stable core if it has been '
                           'with other images in the same class in all the '
                           'previous levels except possibly a few of them. '
                           'Tolerance defines how few is few. Tolerance=0 '
                           'means that an image must be in all previous levels '
                           'with the rest of images in the core.',)
        form.addParam("computeHierarchy", param.BooleanParam, default=False,
                      label="Compute class hierarchy",
                      expertLevel=const.LEVEL_ADVANCED)
        form.addParam("analyzeRejected", param.BooleanParam, default=False,
                      label="Analyze rejected particles",
                      expertLevel=const.LEVEL_ADVANCED,
                      help='To see the analysis you need to browse the '
                           'execution directory and go into the different '
                           'levels')
        form.addParallelSection(threads=0, mpi=4)

    #--------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):
        """ Mainly prepare the command line for call cl2d program"""

        # Convert input images if necessary
        self.imgsFn = self._getExtraPath('images.xmd')
        self.initialClassesFn = self._getExtraPath('initialClasses.xmd')
        self._insertFunctionStep('convertInputStep',
                                 self.inputParticles.get().getObjId(), 
                                 self.initialClasses.get().getObjId() if self.initialClasses.get() else None)

        # Prepare arguments to call program: xmipp_classify_CL2D
        self._params = {'imgsFn': self.imgsFn,
                        'extraDir': self._getExtraPath(),
                        'nref': self.numberOfClasses.get(),
                        'nref0': self.numberOfInitialClasses.get(),
                        'iter': self.numberOfIterations.get(),
                        'extraParams': self.extraParams.get(''),
                        'thZscore': self.thZscore.get(),
                        'thPCAZscore': self.thPCAZscore.get(),
                        'tolerance': self.tolerance.get(),
                        'initialClassesFn': self.initialClassesFn
        }
        args = '-i %(imgsFn)s --odir %(extraDir)s --oroot level --nref ' \
               '%(nref)d --iter %(iter)d %(extraParams)s'
        if self.comparisonMethod == CMP_CORRELATION:
            args += ' --distance correlation'
        if self.clusteringMethod == CL_CLASSICAL:
            args += ' --classicalMultiref'
        if self.randomInitialization:
            args += ' --nref0 %(nref0)d'
        else:
            args += ' --ref0 %(initialClassesFn)s'

        self._insertClassifySteps("xmipp_classify_CL2D", args, subset=CLASSES)

        #TODO: Added this If. Check with COSS error if makes sense.
        #Also, if conditions below are enough to validate that classes core
        # and stable core are not empty

        if not self.randomInitialization:
            self.numberOfInitialClasses.set(self.initialClasses.get().getSize())

        # Analyze cores and stable cores
        if self.numberOfClasses > self.numberOfInitialClasses and self.doCore:
            program = "xmipp_classify_CL2D_core_analysis"
            args = "--dir %(extraDir)s --root level "
            # core analysis
            self._insertClassifySteps(program, args + "--computeCore %(thZscore)f %(thPCAZscore)f", subset=CLASSES_CORE)
            if self.analyzeRejected:
                self._insertFunctionStep('analyzeOutOfCores', CLASSES_CORE)

            if self.numberOfClasses > (2 * self.numberOfInitialClasses.get()) and self.doStableCore: # Number of levels should be > 2
                # stable core analysis
                self._insertClassifySteps(program, args + "--computeStableCore %(tolerance)d", subset=CLASSES_STABLE_CORE)
                if self.analyzeRejected:
                    self._insertFunctionStep('analyzeOutOfCores', CLASSES_STABLE_CORE)

    def _insertClassifySteps(self, program, args, subset=CLASSES):
        """ Defines four steps for the subset:
        1. Run the main program.
        2. Evaluate classes
        3. Sort the classes.
        4. And create output
        """
        self._insertRunJobStep(program, args % self._params)
        self._insertFunctionStep('evaluateClassesStep', subset)
        self._insertFunctionStep('sortClassesStep', subset)
        self._insertFunctionStep('createOutputStep', subset)


    #--------------------------- STEPS functions -------------------------------
    def convertInputStep(self, particlesId, classesId):
        writeSetOfParticles(self.inputParticles.get(),self.imgsFn,alignType=em.ALIGN_NONE)
        if not self.randomInitialization:
            if isinstance(self.initialClasses.get(), SetOfClasses2D):
                writeSetOfClasses2D(self.initialClasses.get(),self.initialClassesFn, writeParticles=False)
            else:
                writeSetOfParticles(self.initialClasses.get(),self.initialClassesFn)

    def sortClassesStep(self, subset=''):
        """ Sort the classes and provided a quality criterion. """
        nproc = self.numberOfMpi.get()
        if nproc < 2:
            nproc = 2 # Force at leat two processor because only MPI version is available
        levelMdFiles = self._getLevelMdFiles(subset)
        for mdFn in levelMdFiles:
            fnRoot = join(dirname(mdFn), "classes%s_sorted" % subset)
            params = "-i classes@%s --oroot %s" % (mdFn, fnRoot)
            self.runJob("xmipp_image_sort", params, numberOfMpi=nproc)
            mdFnOut = fnRoot + ".xmd"
            mdOut = md.MetaData(mdFnOut)
            for objId in md:
                mdOut.setValue(md.MDL_ITEM_ID,long(md.getValue(md.MDL_REF,objId)),objId)
            mdOut.write("classes_sorted@" + mdFn, md.MD_APPEND)

    def evaluateClassesStep(self, subset=''):
        """ Calculate the FRC and output the hierarchy for
        each level of classes.
        """
        levelMdFiles = self._getLevelMdFiles(subset)
        hierarchyFnOut = self._getExtraPath("classes%s_hierarchy.txt" % subset)
        prevMdFn = None
        for mdFn in levelMdFiles:
            self.runJob("xmipp_classify_evaluate_classes", "-i " + mdFn, numberOfMpi=1)
            if self.computeHierarchy and prevMdFn is not None:
                args = "--i1 %s --i2 %s -o %s" % (prevMdFn, mdFn, hierarchyFnOut)
                if exists(hierarchyFnOut):
                    args += " --append"
                self.runJob("xmipp_classify_compare_classes",args, numberOfMpi=1)
            prevMdFn = mdFn

    def createOutputStep(self, subset=''):
        """ Store the SetOfClasses2D object
        resulting from the protocol execution.
        """
        levelMdFiles = self._getLevelMdFiles(subset)
        lastMdFn = levelMdFiles[-1]


        inputParticles = self.inputParticles.get()
        classes2DSet = self._createSetOfClasses2D(inputParticles, subset)

        classes2DSet.copyItems(inputParticles,
                               updateItemCallback=self._createItemMatrix,
                               itemDataIterator=md.iterRows(self.imgsFn,
                                                            sortByLabel=md.MDL_ITEM_ID))


        result = {'outputClasses' + subset: classes2DSet}
        self._defineOutputs(**result)
        self._defineSourceRelation(self.inputParticles, classes2DSet)
    
    def analyzeOutOfCores(self,subset):
        """ Analyze which images are out of cores """
        levelMdFiles = self._getLevelMdFiles(subset)
        for fn in levelMdFiles:
            mdAll=md.MetaData()
            blocks = md.getBlocksInMetaDataFile(fn)
            fnDir=dirname(fn)
            # Gather all images in block
            for block in blocks:
                if block.startswith('class0'):
                    mdClass=md.MetaData(block+"@"+fn)
                    mdAll.unionAll(mdClass)
            if mdAll.size()>0:
                # Compute difference to images
                fnSubset=join(fnDir,"images%s.xmd"%subset)
                mdAll.write(fnSubset)
                fnOutOfSubset=join(fnDir,"imagesOut.xmd")
                self.runJob("xmipp_metadata_utilities","-i %s --set subtraction %s -o %s"%(self.imgsFn,fnSubset,fnOutOfSubset),
                            numberOfMpi=1,numberOfThreads=1)
                
                # Remove disabled and intermediate files
                mdClass=md.MetaData(fnOutOfSubset)
                mdClass.removeDisabled()
                fnRejected="images_rejected@"+fn
                mdClass.write(fnRejected,md.MD_APPEND)
                cleanPath(fnOutOfSubset)
                cleanPath(fnSubset)
                
                # If enough images, make a small summary
                if mdClass.size()>100:
                    from math import ceil
                    fnRejectedDir=join(fnDir,"rejected%s"%subset)
                    makePath(fnRejectedDir)
                    Nclasses=int(ceil(mdClass.size()/300))
                    self.runJob("xmipp_classify_CL2D",
                                "-i %s --nref0 1 --nref %d --iter 5 --distance correlation --classicalMultiref --classifyAllImages --odir %s"\
                                %(fnRejected,Nclasses,fnRejectedDir))

    #--------------------------- INFO functions --------------------------------------------
    def _validate(self):
        validateMsgs = []
        if self.numberOfMpi <= 1:
            validateMsgs.append('Mpi needs to be greater than 1.')
        if self.numberOfInitialClasses > self.numberOfClasses:
            validateMsgs.append('The number of final classes cannot be smaller than the number of initial classes')
        if isinstance(self.initialClasses.get(), SetOfClasses2D):
            if not self.initialClasses.get().hasRepresentatives():
                validateMsgs.append("The input classes should have representatives.")
        return validateMsgs
    
    def _warnings(self):
        validateMsgs = []
        if self.inputParticles.get().getSamplingRate() < 3:
            validateMsgs.append("The sampling rate is smaller than 3 A/pix, consider downsampling the input images to speed-up the process. "\
                         "Probably you don't want such a precise 2D classification.")
        return validateMsgs       

    def _citations(self):
        citations=['Sorzano2010a']
        if self.doCore:
            citations.append('Sorzano2014')
        return citations

    def _summaryLevelFiles(self, summary, levelFiles, subset):
        if levelFiles:
            levels = [self._getLevelFromFile(fn) for fn in levelFiles]
            summary.append('Computed classes%s, levels: %s' % (subset, levels))

    def _summary(self):
        summary = []
        summary.append("Input Particles: *%d*\nClassified into *%d* classes\n" % (self.inputParticles.get().getSize(),
                                                                              self.numberOfClasses.get()))
        #summary.append('- Used a _clustering_ algorithm to subdivide the original dataset into the given number of classes')

        levelFiles = self._getLevelMdFiles()
        if levelFiles:
            self._summaryLevelFiles(summary, levelFiles, CLASSES)
            self._summaryLevelFiles(summary, self._getLevelMdFiles(CLASSES_CORE), CLASSES_CORE)
            self._summaryLevelFiles(summary, self._getLevelMdFiles(CLASSES_STABLE_CORE), CLASSES_STABLE_CORE)
        else:
            summary.append("Output classes not ready yet.")

        return summary

    def _methods(self):
        strline = ''
        if hasattr(self, 'outputClasses'):
            strline += 'We classified %d particles from %s ' % (self.inputParticles.get().getSize(), 
                                                                self.getObjectTag('inputParticles'))
            strline += 'into %d classes %s using CL2D [Sorzano2010a]. ' % (self.numberOfClasses, 
                                                                           self.getObjectTag('outputClasses'))
            strline += '%s method was used to compare images and %s clustering criterion. '%\
                           (self.getEnumText('comparisonMethod'), self.getEnumText('clusteringMethod'))
            if self.numberOfClasses > self.numberOfInitialClasses and self.doCore:
                strline+='We also calculated the class cores %s' % self.getObjectTag('outputClasses_core')
                if self.numberOfClasses > (2 * self.numberOfInitialClasses.get()) and self.doStableCore: # Number of levels should be > 2
                    strline += ' and the class stable cores %s' % self.getObjectTag('outputClasses_stable_core')
                strline+=' [Sorzano2014].'
        return [strline]
    
    #--------------------------- UTILS functions -------------------------------
    def _getLevelMdFiles(self, subset=''):
        """ Grab the metadata class files for each level. """
        levelMdFiles = glob(self._getExtraPath("level_??/level_classes%s.xmd" % subset))
        levelMdFiles.sort()
        return levelMdFiles
    
    def _getLevelFromFile(self, filename):
        """ Return the level number from the filename. """
        folder = dirname(filename) 
        # Folder should ends in level_?? where ?? is level number
        # so we split by '_' and take the last part as int
        return int(folder.split('_')[-1])
    
    def _createItemMatrix(self, item, row):
        createItemMatrix(item, row, align=em.ALIGN_2D)

    def _updateParticle(self, item, row):
        item.setClassId(row.getValue(md.MDL_REF))
        item.setTransform(rowToAlignment(row, em.ALIGN_2D))

    def _updateClass(self, item):
        classId = item.getObjId()

        if classId in self._classesInfo:
            index, fn, _ = self._classesInfo[classId]
            item.setAlignment2D()
            item.getRepresentative().setLocation(index, fn)

    def _loadClassesInfo(self, filename):
        """ Read some information about the produced 2D classes
        from the metadata file.
        """
        self._classesInfo = {}  # store classes info, indexed by class id

        mdClasses = md.MetaData(filename)

        for classNumber, row in enumerate(md.iterRows(mdClasses)):
            index, fn = xmippToLocation(row.getValue(md.MDL_IMAGE))
            # Store info indexed by id, we need to store the row.clone() since
            # the same reference is used for iteration
            self._classesInfo[classNumber + 1] = (index, fn, row.clone())

    def _fillClassesFromLevel(self, clsSet, level):
        """ Create the SetOfClasses2D from a given iteration. """
        self._loadClassesInfo(self._getIterMdClasses(level))
        dataXmd = self._getIterMdImages(level)
        clsSet.classifyItems(updateItemCallback=self._updateParticle,
                             updateClassCallback=self._updateClass,
                             itemDataIterator=md.iterRows(dataXmd,
                                                          sortByLabel=md.MDL_ITEM_ID))

