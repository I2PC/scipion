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
This module contains the protocol base class for Relion protocols
"""

import xmipp
from pyworkflow.em import *
from pyworkflow.protocol.params import BooleanParam, PointerParam, IntParam



class ProtRelionBase(EMProtocol):
    """ This class cointains the common functionalities for all Relion protocols.
    In subclasses there should be little changes about how to create the command line
    and the files produced.
    
    Most of the Relion protocols, have two modes: NORMAL or CONTINUE. That's why
    some of the function have a template pattern approach to define the behaivour
    depending on the case.
    """
    IS_CLASSIFY = True
    IS_2D = False
    
    def __init__(self, **args):        
        EMProtocol.__init__(self, **args)
        
    def _initialize(self):
        """ This function is mean to be called after the 
        working dir for the protocol have been set. (maybe after recovery from mapper)
        """
        self._createFilenameTemplates()
        self._createIterTemplates()
        
        self.FileKeys = ['data', 'optimiser', 'sampling'] 
        self.ClassLabel = xmipp.MDL_REF # by default 3d
        self.ClassFnTemplate = '%(rootDir)s/relion_it%(iter)03d_class%(ref)03d.mrc:mrc'
        self.outputClasses = 'classes_ref3D.xmd'
        self.outputVols = 'volumes.xmd'
        
        
    def _createFilenameTemplates(self):
        """ Centralize how files are called for iterations and references. """
        self.extraIter = self.extraPath('relion_it%(iter)03d_')
        myDict = {
                  'input_particles': self._getPath('input_particles.star'),
                  'data_sorted_xmipp': self.extraIter + 'data_sorted_xmipp.star',
                  'classes_xmipp': self.extraIter + 'classes_xmipp.xmd',
                  'angularDist_xmipp': self.extraIter + 'angularDist_xmipp.xmd',
                  'all_avgPmax_xmipp': self.tmpPath('iterations_avgPmax_xmipp.xmd'),
                  'all_changes_xmipp': self.tmpPath('iterations_changes_xmipp.xmd'),
                  'selected_volumes': self.tmpPath('selected_volumes_xmipp.xmd')
                  }
        # add to keys, data.star, optimiser.star and sampling.star
        for key in self.FileKeys:
            myDict[key] = self.extraIter + '%s.star' % key
            key_xmipp = key + '_xmipp'             
            myDict[key_xmipp] = self.extraIter + '%s.xmd' % key
        # add other keys that depends on prefixes
        for p in self._getPrefixes():            
            myDict['%smodel' % p] = self.extraIter + '%smodel.star' % p
            myDict['%svolume' % p] = self.extraIter + p + 'class%(ref3d)03d.mrc:mrc'

        self._fnDict = myDict
    
    def _createIterTemplates(self):
        """ Setup the regex on how to find iterations. """
        self._iterTemplate = self._getFileName('data', iter=0).replace('000','???')
        # Iterations will be identify by _itXXX_ where XXX is the iteration number
        # and is restricted to only 3 digits.
        self._iterRegex = re.compile('_it(\d{3,3})_')
        
        
    #--------------------------- DEFINE param functions --------------------------------------------   
    def _defineParams(self, form):
        form.addSection(label='Input')
        # Some hidden variables to be used for conditions
        form.addHidden('isClassify', BooleanParam, default=self.IS_CLASSIFY)
        form.addHidden('is2D', BooleanParam, default=self.IS_2D)
        
        form.addParam('doContinue', BooleanParam, default=False,
                      label='Continue from a previous run?',
                      help='If you set to *Yes*, you should select a previous\n'
                      'run of type *%s* class and most of the input parameters\n'
                      'will be taken from it.' % self.getClassName())
        form.addParam('inputParticles', PointerParam, pointerClass='SetOfParticles',
                      condition='not doContinue',
                      label="Input particles",  
                      help='Select the input images from the project.')   
        form.addParam('previousRun', PointerParam, pointerClass=self.getClassName(),
                      condition='doContinue',
                      label='Select previous run',
                      help='Select a previous run to continue from.')
        
        form.addParam('numberOfClasses', IntParam, default=3, 
                      condition='not doContinue and isClassify',
                      label='Number of classes:',
                      help='The number of classes (K) for a multi-reference refinement.\n'
                           'These classes will be made in an unsupervised manner from a single\n'
                           'reference by division of the data into random subsets during the\n'
                           'first iteration.')
        form.addParam('referenceVolume', PointerParam, pointerClass='Volume',
                      condition='not doContinue and not is2D',
                      label="Initial 3D map", 
                      help='Initial reference 3D map, it should have the same \n'
                           'dimensions and the same pixel size as your input particles.')
        form.addParam('isAbsoluteGreyScale', BooleanParam, default=False,
                      condition='not doContinue and not is2D',
                      label="Is initial 3D map on absolute greyscale?", 
                      help='The probabilities are based on squared differences, \n'
                           'so that the absolute grey scale is important.       \n'
                           'Probabilities are calculated based on a Gaussian noise model,\n'
                           'which contains a squared difference term between the reference and the experimental image.\n' 
                           'This has a consequence that the reference needs to be on the same absolute intensity \n'
                           'grey-scale as the experimental images. RELION and XMIPP reconstruct maps at their absolute\n'
                           'intensity grey-scale. Other packages may perform internal normalisations of the reference\n' 
                           'density, which will result in incorrect grey-scales. Therefore: if the map was reconstructed\n'
                           'in RELION or in XMIPP, set this option to Yes, otherwise set it to No. If set to No, RELION \n'
                           'will use a (grey-scale invariant) cross-correlation criterion in the first iteration, and \n'
                           'prior to the second iteration the map will be filtered again using the initial low-pass filter.\n'
                           'This procedure is relatively quick and typically does not negatively affect the outcome of the\n'
                           'subsequent MAP refinement. Therefore, if in doubt it is recommended to set this option to No.')        
        form.addParam('symmetryGroup', StringParam, default='c1',
                      condition='not doContinue and not is2D',
                      label="Symmetry group", 
                      help='See http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Symmetry \n'
                           'for a description of the symmetry groups format')                     
 
    #--------------------------- INSERT steps functions --------------------------------------------  
    def _insertAllSteps(self): 
        self._insertFunctionStep('convertInputStep')
        
        if self.doContinue:
            self._insertStepsContinue()
        else:
            self._insertStepsNormal()
            
        self._insertFunctionStep('createOutputStep')
 
   
    #--------------------------- STEPS functions --------------------------------------------       
    def convertInputStep(self):
        """ Create the input file in STAR format as expected by Relion.
        If the input particles comes from Relion, just link the file. 
        """
        from convert import createRelionInputParticles
        createRelionInputParticles(self.inputParticles.get(), 
                                   self._getFileName('input_particles'))
        
    def createOutputStep(self):
        pass # should be implemented in subclasses
        
    
    #--------------------------- INFO functions -------------------------------------------- 
    def _validate(self):
        errors = []
        if self.doContinue:
            errors += self._validateContinue()
        else:
            errors += self._validateNormal()
        return errors
            
    def _validateNormal(self):
        """ Should be overriden in subclasses to 
        return summary message for NORMAL EXECUTION. 
        """
        return []
    
    def _validateContinue(self):
        """ Should be overriden in subclasses to
        return summary messages for CONTINUE EXECUTION.
        """
        return []
    
    def _citations(self):
        cites = []
        return cites
    
    def _summary(self):
        if self.DoContinue:
            return self._summaryContinue()
        return self._summaryNormal()

    def _summaryNormal(self):
        """ Should be overriden in subclasses to 
        return summary message for NORMAL EXECUTION. 
        """
        return []
    
    def _summaryContinue(self):
        """ Should be overriden in subclasses to
        return summary messages for CONTINUE EXECUTION.
        """
        return []


    #--------------------------- UTILS functions --------------------------------------------
    def _getProgram(self):
        """ Get the program name depending on the MPI use or not. """
        program = 'relion_refine'
        if self.NumberOfMpi > 1:
            program += '_mpi'
        return program
    
    def _getFileName(self, key, **args):
        """ Retrieve a filename from the templates. """
        return self._fnDict[key] % args
    
    def _getIterNumber(self, index):
        """ Return the list of iteration files, give the iterTemplate. """
        result = None
        files = sorted(glob(self._iterTemplate))
        if files:
            f = files[index]
            s = self._iterRegex.search(f)
            if s:
                result = int(s.group(1)) # group 1 is 3 digits iteration number
        return result
        
    def _lastIter(self):
        return self._getIterNumber(-1)

    def firstIter(self):
        return self._getIterNumber(0) or 1