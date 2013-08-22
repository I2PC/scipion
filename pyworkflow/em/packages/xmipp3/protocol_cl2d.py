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
"""
This sub-package contains wrapper around CL2D Xmipp program
"""

from os.path import join, dirname, exists
from pyworkflow.em import *  
import xmipp
from data import *
from convert import createXmippInputImages, readSetOfClasses2D
#from xmipp3 import XmippProtocol
from glob import glob

# Comparison methods enum
CMP_CORRELATION = 0
CMP_CORRENTROPY = 1

# Clustering methods enum
CL_CLASSICAL = 0
CL_ROBUST = 1


class XmippDefCL2D(Form):
    """Create the definition of parameters for
    the XmippProtCL2D protocol"""
    def __init__(self):
        Form.__init__(self)
    
        self.addSection(label='Input')
        self.addParam('inputImages', PointerParam, label="Input images", important=True, 
                      pointerClass='SetOfParticles',
                      help='Select the input images from the project.'
                           'It should be a SetOfImages class')        
        self.addParam('numberOfReferences', IntParam, default=64,
                      label='Number of references:',
                      help='Number of references (or classes) to be generated.')
        self.addParam('numberOfInitialReferences', IntParam, default=4, expertLevel=LEVEL_ADVANCED,
                      label='Number of initial references:',
                      help='Initial number of references used in the first level.')
        self.addParam('numberOfIterations', IntParam, default=4, expertLevel=LEVEL_ADVANCED,
                      label='Number of iterations:',
                      help='Maximum number of iterations within each level.')         
        self.addParam('comparisonMethod', EnumParam, choices=['correlation', 'correntropy'],
                      label="Comparison method", default=CMP_CORRELATION,
                      display=EnumParam.DISPLAY_COMBO,
                      help='Use correlation or correntropy')
        self.addParam('clusteringMethod', EnumParam, choices=['classical', 'robust'],
                      label="Clustering method", default=CL_CLASSICAL,
                      display=EnumParam.DISPLAY_COMBO,
                      help='Use the classical clustering criterion or the robust')
        self.addParam('extraParams', StringParam, expertLevel=LEVEL_EXPERT,
              label='Additional parameters',
              help='Additional parameters for classify_CL2D:\n  --verbose, --corrSplit, ...')   
        
        self.addSection(label='Core analysis')        
        self.addParam('thZscore', FloatParam, default=3,
                      label='Junk Zscore',
                      help='Which is the average Z-score to be considered as junk. Typical values'
                           'go from 1.5 to 3. For the Gaussian distribution 99.5% of the data is'
                           'within a Z-score of 3. Lower Z-scores reject more images. Higher Z-scores'
                           'accept more images.')
        self.addParam('thPCAZscore', FloatParam, default=3,
                      label='PCA Zscore',
                      help='Which is the PCA Z-score to be considered as junk. Typical values'
                           'go from 1.5 to 3. For the Gaussian distribution 99.5% of the data is'
                           'within a Z-score of 3. Lower Z-scores reject more images. Higher Z-scores'
                           'accept more images.')        
        self.addParam('tolerance', IntParam, default=1,
                      label='Tolerance',
                      help='An image belongs to the stable core if it has been with other images in the same class'
                           'in all the previous levels except possibly a few of them. Tolerance defines how few is few.'
                           'Tolerance=0 means that an image must be in all previous levels with the rest of images in'
                           'the core.')          
        
        self.addParallelSection(threads=0, mpi=2)
        
        
class XmippProtCL2D(ProtAlign, ProtClassify):
    """ Protocol to preprocess a set of micrographs in the project. """
    _definition = XmippDefCL2D()
    _label = 'Xmipp CL2D'

    def _defineSteps(self):
        """ Mainly prepare the command line for call cl2d program"""
        
        # Convert input images if necessary
        imgsFn = createXmippInputImages(self, self.inputImages.get())
        
        # Prepare arguments to call program: xmipp_classify_CL2D
        self._params = {'imgsFn': imgsFn, 
                        'extraDir': self._getExtraPath(),
                        'nref': self.numberOfReferences.get(), 
                        'nref0': self.numberOfInitialReferences.get(),
                        'iter': self.numberOfIterations.get(),
                        'extraParams': self.extraParams.get(''),
                        'thZscore': self.thZscore.get(),
                        'thPCAZscore': self.thPCAZscore.get(),
                        'tolerance': self.tolerance.get()
                      }
        args = '-i %(imgsFn)s --odir %(extraDir)s --oroot level --nref %(nref)d --iter %(iter)d %(extraParams)s'
        if self.comparisonMethod == CMP_CORRELATION:
            args += ' --distance correlation'
        if self.clusteringMethod == CL_CLASSICAL:
            args += ' --classicalMultiref'
        if not self.extraParams.hasValue() or not '--ref0' in self.extraParams.get():
            args += ' --nref0 %(nref0)d'
    
        self._defineClassifySteps("xmipp_classify_CL2D", args)
        
        # Analyze cores and stable cores
        if self.numberOfReferences > self.numberOfInitialReferences:
            program = "xmipp_classify_CL2D_core_analysis"
            args = "--dir %(extraDir)s --root level "
            # core analysis
            self._defineClassifySteps(program, args + "--computeCore %(thZscore)f %(thPCAZscore)f", subset='_core')
            if self.numberOfReferences > (2 * self.numberOfInitialReferences.get()): # Number of levels should be > 2
                # stable core analysis
                self._defineClassifySteps(program, args + "--computeStableCore %(tolerance)d", subset='_stable_core')
        
    def _defineClassifySteps(self, program, args, subset=''):
        """ Defines four steps for the subset:
        1. Run the main program.
        2. Evaluate classes
        3. Sort the classes.
        4. And create output
        """
        self._insertRunJobStep(program, args % self._params)
        self._insertFunctionStep('evaluateClasses', subset)
        self._insertFunctionStep('sortClasses', subset)
        self._insertFunctionStep('createOutput', subset)        
        
    def _getLevelMdFiles(self, subset=''):
        """ Grab the metadata class files for each level. """
        levelMdFiles = glob(self._getExtraPath("level_??/level_classes%s.xmd" % subset))
        levelMdFiles.sort()
        return levelMdFiles        
    
    def sortClasses(self, subset=''):
        """ Sort the classes and provided a quality criterion. """
        nproc = self.numberOfMpi.get()
        if nproc < 2:
            nproc = 2 # Force at leat two processor because only MPI version is available
        levelMdFiles = self._getLevelMdFiles(subset)
        for mdFn in levelMdFiles:
            fnRoot = join(dirname(mdFn), "classes%s_sorted" % subset)
            params = "-i classes@%s --oroot %s" % (mdFn, fnRoot)
            self.runJob(None, "xmipp_image_sort", params, nproc)
            mdFnOut = fnRoot + ".xmd"
            md = xmipp.MetaData(mdFnOut)
            md.write("classes_sorted@" + mdFn, xmipp.MD_APPEND)
            #deleteFile(log,fnRoot+".xmd")
        
    def evaluateClasses(self, subset=''):
        """ Calculate the FRC and output the hierarchy for 
        each level of classes.
        """
        levelMdFiles = self._getLevelMdFiles(subset)
        hierarchyFnOut = self._getExtraPath("classes%s_hierarchy.txt" % subset)
        prevMdFn = None
        for mdFn in levelMdFiles:
            self.runJob(None, "xmipp_classify_evaluate_classes","-i " + mdFn)
            if prevMdFn is not None:
                args = "--i1 %s --i2 %s -o %s" % (prevMdFn, mdFn, hierarchyFnOut)
                if exists(hierarchyFnOut):
                    args += " --append"
                self.runJob(None, "xmipp_classify_compare_classes",args)
            prevMdFn = mdFn
            
    def createOutput(self, subset=''):
        """ Store the setOfClasses2D object 
        as result of the protocol. 
        """
        levelMdFiles = self._getLevelMdFiles(subset)
        lastMdFn = levelMdFiles[-1]
        result = {'outputClassification' + subset: readSetOfClasses2D(self._createSetOfClasses2D(subset), lastMdFn, 'classes_sorted')}
        self._defineOutputs(**result)

    def _summary(self):
        summary = []
        if not hasattr(self, 'outputClassification'):
            summary.append("Output classes not ready yet.")
        else:
            summary.append("Input Images: %s" % self.inputImages.get().getNameId())
            summary.append("Number of references: %d" % self.numberOfReferences.get())
            summary.append("Output classes: %s" % self.outputClassification.get())
        return summary