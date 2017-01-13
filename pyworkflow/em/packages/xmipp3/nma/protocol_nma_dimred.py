# **************************************************************************
# *
# * Authors:  J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es), Nov 2014
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


from os.path import basename

from pyworkflow.object import String
from pyworkflow.protocol.params import (PointerParam, StringParam, EnumParam, IntParam,
                                        LEVEL_ADVANCED)
from pyworkflow.em.protocol import ProtAnalysis3D
import xmipp 



DIMRED_PCA = 0
DIMRED_LTSA = 1
DIMRED_DM = 2
DIMRED_LLTSA = 3
DIMRED_LPP = 4
DIMRED_KPCA = 5
DIMRED_PPCA = 6
DIMRED_LE = 7
DIMRED_HLLE = 8
DIMRED_SPE = 9
DIMRED_NPE = 10

# Values to be passed to the program
DIMRED_VALUES = ['PCA', 'LTSA', 'DM', 'LLTSA', 'LPP', 'kPCA', 'pPCA', 'LE', 'HLLE', 'SPE', 'NPE']

# Methods that allows mapping
DIMRED_MAPPINGS = [DIMRED_PCA, DIMRED_LLTSA, DIMRED_LPP, DIMRED_PPCA, DIMRED_NPE]

       
class XmippProtDimredNMA(ProtAnalysis3D):
    """ This protocol will take the images with NMA deformations
    as points in a N-dimensional space (where N is the number
    of computed normal modes) and will project them in a reduced
    spaced (usually with less dimensions).
    """
    _label = 'nma dimred'
    
    def __init__(self, **kwargs):
        ProtAnalysis3D.__init__(self, **kwargs)
        self.mappingFile = String()
    
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputNMA', PointerParam, pointerClass='XmippProtAlignmentNMA',
                      label="Conformational distribution",                        
                      help='Select a previous run of the NMA alignment.')
        
        form.addParam('dimredMethod', EnumParam, default=DIMRED_PCA,
                      choices=['Principal Component Analysis (PCA)',
                               'Local Tangent Space Alignment',
                               'Diffusion map',
                               'Linear Local Tangent Space Alignment',
                               'Linearity Preserving Projection',
                               'Kernel PCA',
                               'Probabilistic PCA',
                               'Laplacian Eigenmap',
                               'Hessian Locally Linear Embedding',
                               'Stochastic Proximity Embedding',
                               'Neighborhood Preserving Embedding'],
                      label='Dim-Red method',
                      help=""" Dimensionality Reduction method.
    PCA
       Principal Component Analysis 
    LTSA <k=12>
       Local Tangent Space Alignment, k=number of nearest neighbours 
    DM <s=1> <t=1>
       Diffusion map, t=Markov random walk, s=kernel sigma 
    LLTSA <k=12>
       Linear Local Tangent Space Alignment, k=number of nearest neighbours 
    LPP <k=12> <s=1>
       Linearity Preserving Projection, k=number of nearest neighbours, s=kernel sigma 
    kPCA <s=1>
       Kernel PCA, s=kernel sigma 
    pPCA <n=200>
       Probabilistic PCA, n=number of iterations 
    LE <k=7> <s=1>
       Laplacian Eigenmap, k=number of nearest neighbours, s=kernel sigma 
    HLLE <k=12>
       Hessian Locally Linear Embedding, k=number of nearest neighbours 
    SPE <k=12> <global=1>
       Stochastic Proximity Embedding, k=number of nearest neighbours, global embedding or not 
    NPE <k=12>
       Neighborhood Preserving Embedding, k=number of nearest neighbours 
""")
        form.addParam('extraParams', StringParam, level=LEVEL_ADVANCED,
                      label="Extra params", 
                      help='This parameters will be passed to the program.')
                      
        form.addParam('reducedDim', IntParam, default=2,
                      label='Reduced dimension')
        form.addParallelSection(threads=0, mpi=0)    
    
    
    #--------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):
        # Take deforamtions text file and the number of images and modes
        inputSet = self.getInputParticles()
        rows = inputSet.getSize()
        reducedDim = self.reducedDim.get()
        method = self.dimredMethod.get()
        extraParams = self.extraParams.get('')
        
        deformationsFile = self.getDeformationFile()
        
        self._insertFunctionStep('convertInputStep', 
                                 deformationsFile, inputSet.getObjId())
        self._insertFunctionStep('performDimredStep', 
                                 deformationsFile, method, extraParams,
                                 rows, reducedDim) 
        self._insertFunctionStep('createOutputStep')
        
        
    #--------------------------- STEPS functions --------------------------------------------   
    
    def convertInputStep(self, deformationFile, inputId):
        """ Iterate through the images and write the 
        plain deformation.txt file that will serve as 
        input for dimensionality reduction.
        """
        inputSet = self.getInputParticles()
        f = open(deformationFile, 'w')
        
        for particle in inputSet:
            f.write(' '.join(particle._xmipp_nmaDisplacements))
            f.write('\n')
        f.close()
    
    def performDimredStep(self, deformationsFile, method, extraParams,
                          rows, reducedDim):
        outputMatrix = self.getOutputMatrixFile()
        methodName = DIMRED_VALUES[method]
        # Get number of columes in deformation files
        # it can be a subset of inputModes
        f = open(deformationsFile)
        columns = len(f.readline().split()) # count number of values in first line
        f.close()
        
        args = "-i %(deformationsFile)s -o %(outputMatrix)s -m %(methodName)s %(extraParams)s"
        args += "--din %(columns)d --samples %(rows)d --dout %(reducedDim)d"
        if method in DIMRED_MAPPINGS:
            mappingFile = self._getExtraPath('projector.txt')
            args += " --saveMapping %(mappingFile)s"
            self.mappingFile.set(mappingFile)
        self.runJob("xmipp_matrix_dimred", args % locals())
        
    def createOutputStep(self):
        pass

    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        summary = []
        return summary
    
    def _validate(self):
        errors = []
        return errors
    
    def _citations(self):
        return []
    
    def _methods(self):
        return []
    
    #--------------------------- UTILS functions --------------------------------------------

    def getInputParticles(self):
        """ Get the output particles of the input NMA protocol. """
        return self.inputNMA.get().outputParticles
    
    def getInputPdb(self):
        return self.inputNMA.get().getInputPdb()
    
    def getOutputMatrixFile(self):
        return self._getExtraPath('output_matrix.txt')
    
    def getDeformationFile(self):
        return self._getExtraPath('deformations.txt')
    
    def getProjectorFile(self):
        return self.mappingFile.get()
    
    def getMethodName(self):
        return DIMRED_VALUES[self.dimredMethod.get()]
