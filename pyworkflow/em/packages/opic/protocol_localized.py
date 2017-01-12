# *****************************************************************************
# *
# * Authors:     Josue Gomez Blanco (jgomez@cnb.csic.es)
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
# *****************************************************************************

import math

from os.path import join, exists
from pyworkflow.protocol.params import (PointerParam, BooleanParam, StringParam,
                                        EnumParam, NumericRangeParam,
                                        PathParam, FloatParam, LEVEL_ADVANCED)
from pyworkflow.em.protocol import ProtParticles, ProtParticlePicking
from pyworkflow.em.data import Coordinate, SetOfCoordinates

from pyworkflow.utils.path import replaceBaseExt

from convert import particleToRow, rowToSubcoordinate, setEnviron


CMM = 0
HAND = 1


class ProtLocalizedRecons(ProtParticlePicking, ProtParticles):
    """ This class cointains a re-implementation to a method for the
    localized three-dimensional reconstruction of such subunits. 
    After determining the particle orientations, local areas 
    corresponding to the subunits can be extracted and treated as 
    single particles.
    """
    _label = 'localized subparticles'
    
    def __init__(self, **args):
        ProtParticlePicking.__init__(self, **args)
        ProtParticles.__init__(self, **args)
        
    
    #--------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputParticles', PointerParam,
                      pointerClass='SetOfParticles',
                      pointerCondition='hasAlignmentProj',
                      important=True,
                      label="Input particles",  
                      help='Select the input images from the project.')

        group = form.addGroup('Symmetry')
        group.addParam('symmetryGroup', StringParam, default='c1',
                      label="Symmetry", 
                      help='If the molecule is asymmetric, set Symmetry group '
                           'to C1. Note their are multiple possibilities for '
                           'icosahedral symmetry: \n'
                           '* I1: No-Crowther 222 (standard in Heymann, '
                           'Chagoyen & Belnap, JSB, 151 (2005) 196-207)\n'
                           '* I2: Crowther 222 \n'
                           '* I3: 52-setting (as used in SPIDER?) \n'
                           '* I4: A different 52 setting \n')

        group.addParam('randomize', BooleanParam, default=False,
                      label='Randomize the order of the symmetry matrices?',
                      help='Useful for preventing preferred orientations.')
        group.addParam('relaxSym', BooleanParam, default=False,
                       expertLevel=LEVEL_ADVANCED,
                       label='Relax symmetry?',
                       help='Create one random subparticle for each particle ')

        group = form.addGroup('Vectors')
        group.addParam('defineVector',  EnumParam, default=CMM,
                      label='Is vector defined by?',
                      choices=['cmm file', 'string'],
                      display=EnumParam.DISPLAY_HLIST,
                      help='')
        group.addParam('vector', NumericRangeParam, default='0,0,1',
                      label='Location vectors', condition="defineVector==1",
                      help='Vector defining the location of the '
                           'subparticles. The vector is defined by 3 '
                           'values x,y,z separated by comma. \n'
                           'More than one vector can be specified separated by'
                           'semicolon. For example: \n'
                           '0,0,1            # Defines only one vector.\n'
                           '0,0,1; 1,0,0;    # Defines two vectors.'
                       )
        group.addParam('vectorFile', PathParam, default='',
                      condition="defineVector==0",
                      label='file obtained by Chimera: ',
                      help='CMM file defining the location(s) of the '
                           'sub-particle(s). Use instead of vector. ')
        group.addParam('length', StringParam, default=-1,
                      label='Alternative length of the vector (A)',
                      help='Use to adjust the sub-particle center. If it '
                           'is <= 0, the length of the given vector. '
                           'Different values must be separated by commas.')

        form.addSection('Sub-particles')
        form.addParam('alignSubparticles', BooleanParam, default=True,
                      label='Align the subparticles?',
                      help='Align sub-particles to the standard orientation. ')
        form.addParam('unique', FloatParam, default=-1,
                      label='Angle to keep unique sub-particles (deg)',
                      help='Keep only unique subparticles within angular '
                           'distance. It is useful to remove overlapping '
                           'sub-particles on symmetry axis.')
        form.addParam('mindist', FloatParam, default=-1,
                      label='Minimum distance between sub-particles (px)',
                      help='In pixels. Minimum distance between the '
                           'subparticles in the image. All overlapping ones '
                           'will be discarded.')
        form.addParam('side', FloatParam, default=-1,
                      label='Angle to keep sub-particles from side views (deg)',
                      help='Keep only particles within specified angular '
                           'distance from side views. All others will be '
                           'discarded. ')
        form.addParam('top', FloatParam, default=-1,
                      label='Angle to keep sub-particles from top views (deg)',
                      help='Keep only particles within specified angular '
                           'distance from top views. All others will be '
                           'discarded. ')

        form.addParallelSection(threads=0, mpi=0)
    
    #--------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('createOutputStep')
    
    #--------------------------- STEPS functions ------------------------------
    def createOutputStep(self):
        setEnviron() # Set the environment to access localrec modules
        import localrec
        import pyrelion

        inputSet = self._getInputParticles()
        outputSet = self._createSetOfCoordinates(inputSet)
        params = {"symmetryGroup" : self.symmetryGroup.get(),
                  "vector" : self.vector.get(),
                  "vectorFile" : self.vectorFile.get(),
                  "length" : self.length.get(),
                  "unique" : self.unique.get(),
                  "mindist" : self.mindist.get(),
                  "side" : self.side.get(),
                  "top" : self.top.get(),
                  "pxSize" : self.inputParticles.get().getSamplingRate(),
                  "dim" : self.inputParticles.get().getXDim()
                  }

        symMatrices = localrec.matrix_from_symmetry(self.symmetryGroup.get())

        if self.defineVector == CMM:
            cmmFn = params["vectorFile"]
            vector = " "
        else:
            cmmFn = " "
            vector = params["vector"]
            
        subpartVectorList = localrec.load_vectors(cmmFn, vector,
                                                  params["length"],
                                                  params["pxSize"])
        
        # Define some conditions to filter subparticles
        filters = localrec.load_filters(math.radians(params["side"]),
                                        math.radians(params["top"]),
                                        params["mindist"])
        
        coord = Coordinate()
        for part in inputSet:
            partItem = pyrelion.Item()
            particleToRow(part, partItem)
            
            subparticles, _ = localrec.create_subparticles(partItem,
                                                     symMatrices,
                                                     subpartVectorList,
                                                     params["dim"],
                                                     self.relaxSym,
                                                     self.randomize,
                                                     "subparticles",
                                                     params["unique"],
                                                     0,
                                                     self.alignSubparticles,
                                                     "",
                                                     True,
                                                     filters)
            for subpart in subparticles:
                rowToSubcoordinate(subpart, coord, part)
                coord.setObjId(None) # Force to insert as a new item
                outputSet.append(coord)
        
        self._defineOutputs(outputCoordinates=outputSet)
        self._defineSourceRelation(self.inputParticles, outputSet)
    
    #--------------------------- INFO functions -------------------------------
    def _validate(self):
        errors = []
        return errors
    
    def _citations(self):
        return ['Serban2015']

    def _summary(self):
        summary = []
        return summary
    
    def _methods(self):
        return []
    
    #--------------------------- UTILS functions ------------------------------
    def _getInputParticles(self):
        return self.getInputParticlesPointer().get()
    
    def _getSymMatrices(self):
        pass
        #matricesObjs = lr.matrix_from_symmetry(self.symmetryGroup.get())
        # We implement the binding to remove the dependency with relion. When
        # we test the new implementation (code below) and the results were
        # different. THe nunmber of particles removed are diverge in dependency of
        # how you obtain the symmetry matrices.
        # There aren't any obvious bug in the matrices.
        
        
#         matricesObjs = []
#         matrices = md.getSymmetryMatrices(self.symmetryGroup.get())
#         for matrix in matrices:
#             a = matrix[0]
#             a.extend(matrix[1])
#             a.extend(matrix[2])
#             print "Cadena Xmipp: ", a
#             matricesObjs.append(lr.Matrix3(a))
#        return matricesObjs

    def getInputParticlesPointer(self):
        return self.inputParticles

    def registerCoords(self, coordsDir):
        """ This methods is usually inherited from all Pickers
        and it is used from the Java picking GUI to register
        a new SetOfCoordinates when the user click on +Particles button.
        """
        suffix = self.__getOutputSuffix()
        outputName = self.OUTPUT_PREFIX + suffix

        from pyworkflow.em.packages.opic.convert import readSetOfCoordinates
        inputset = self._getInputParticles()
        # micrographs are the input set if protocol is not finished
        outputset = self._createSetOfCoordinates(inputset, suffix=suffix)
        readSetOfCoordinates(coordsDir, self._getInputParticles(), outputset)
        summary = self.getSummary(outputset)
        outputset.setObjComment(summary)
        outputs = {outputName: outputset}
        self._defineOutputs(**outputs)
        self._defineSourceRelation(self.getInputParticlesPointer(), outputset)
        self._store()

    def __getOutputSuffix(self):
        """ Get the name to be used for a new output.
        For example: outputCoordiantes7.
        It should take into account previous outputs
        and number with a higher value.
        """
        maxCounter = -1
        for attrName, _ in self.iterOutputAttributes(SetOfCoordinates):
            suffix = attrName.replace(self.OUTPUT_PREFIX, '')
            try:
                counter = int(suffix)
            except:
                counter = 1 # when there is not number assume 1
            maxCounter = max(counter, maxCounter)

        return str(maxCounter+1) if maxCounter > 0 else '' # empty if not outputs
