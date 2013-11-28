# **************************************************************************
# *
# * Authors:     Jose Gutierrez (jose.gutierrez@cnb.csic.es)
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
This sub-package contains wrapper around XmippProtPreprocessVolumes protocol
"""


from pyworkflow.em import *  
from pyworkflow.utils import *  
from xmipp3 import XmippProtocol
from convert import createXmippInputImages, readSetOfVolumes

#from xmipp3 import XmippProtocol
        
#class XmippProtPreprocessVolumes(ProtPreprocessVolumes, XmippProtocol):
class XmippProtPreprocessVolumes(XmippProtocol):        
    """ Protocol for Xmipp-based preprocess for volumes """
    _label = 'Xmipp Preprocess Volumes'
      
    # Aggregation constants
    AGG_AVERAGE=0
    AGG_SUM=1
    
    # Segmentation type
    SEG_VOXEL=0
    SEG_AMIN=1
    SEG_DALTON=2
    SEG_AUTO=3
    

    def _defineParams(self, form):
        form.addSection(label='Input')
        # Volumes to process
        form.addParam('inputVolumes', PointerParam, label="Volumes to process", important=True, 
                      pointerClass='SetOfVolumes',
                      help='Can be a density volume or a PDB file.')  
        form.addParam('inputVoxel', TextParam, label='Input voxel size (A/voxel)')
        
        form.addSection(label='Preprocess steps')
        # Change hand
        form.addParam('doChangeHand', BooleanParam, default=False,
                      label="Change hand", 
                      help='Change hand by applying a mirror along X.')
        form.addParam('doRandomize', BooleanParam, default=False,
                      label="Randomize phases", 
                      help='Randomize phases beyond a certain frequency.')
        # Maximum Resolution
        form.addParam('doMaxResRand', TextParam, default=40,
                      label="Maximum Resolution", condition='doRandomize == True',
                      help='Angstroms.')
        form.addParam('doLowPass', BooleanParam, default=False,
                      label="Low pass filter", 
                      help='Low pass filter the volume.')
        # Maximum Resolution
        form.addParam('doMaxResFilt', TextParam, default=40,
                      label="Maximum Resolution", condition='doLowPass == True',
                      help='Angstroms.')        
        form.addParam('doSymmetrize', BooleanParam, default=False,
                      label="Symmetrize", 
                      help='Symmetrize the input model.')
        # Symmetry group
        form.addParam('doSymmetryGroup', TextParam, default='c1',
                      label="Symmetry group", condition='doSymmetrize == True',
                      help='See [http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Symmetry]'
                      'for a description of the symmetry groups format, If no symmetry is present, give c1.')  
        form.addParam('doAggretion', EnumParam, choices=['Average', 'Sum'], 
                      display=EnumParam.DISPLAY_COMBO,
                      default=0, label='Aggregation mode', condition = 'DoSymmetrize==True',
                      help='Symmetrized volumes can be averaged or summed.')
        # Spherical Mask
        form.addParam('doMask', BooleanParam, default=False,
                      label="Spherical Mask", 
                      help='Remove all pixels beyond a certain radius.')
        form.addParam('maskRadius', TextParam, default=-1,
                      label="Maximum Resolution", condition='doMask == True',
                      help='In pixels. Set to -1 for half of the size of the volume.') 
        # Adjust gray values
        form.addParam('doAdjust', BooleanParam, default=False,
                      label="Adjust gray values", 
                      help='Adjust input gray values so that it is compatible with a set of projections.')
        form.addParam('setOfProjections', TextParam,
                      pointerClass='SetOfImages',
                      label="Set of projections", condition='doAdjust==True',
                      help='Metadata with a set of images to which the model should conform. The set of images should have the'
                      'final pixel size and the final size of the model.')
        # Normalize background
        form.addParam('doNormalize', BooleanParam, default=False,
                      label="Normalize background", 
                      help='Set background to have zero mean and standard deviation 1.')
        form.addParam('setOfProjections', TextParam,
                      pointerClass='SetOfImages', default=-1,
                      label="Mask Radius", condition='doNormalize==True',
                      help='In pixels. Set to -1 for half of the size of the volume.')
        # Threshold
        form.addParam('doThreshold', BooleanParam, default=False,
                      label="Threshold", 
                      help='Remove voxels below a certain value.')
        form.addParam('threshold', TextParam, default=0,
                      label="Threshold value", condition='doThreshold==True',
                      help='Grey value below which all voxels should be set to 0.')
        # Segment
        form.addParam('doSegment',
                      default=False, label="Segment", 
                      help='Separate the molecule from its background.')
        form.addParam('segmentationType', EnumParam, choices=['Voxel mass', 'Aminoacid mass','Dalton mass','Automatic'],
                      default=0, display=EnumParam.DISPLAY_COMBO,
                      label="Segmentation Type", condition='doSegment==True',
                      help='Type of segmentation.')
        
#        form.addParallelSection(threads=1, mpi=8)         
             
    def _defineSteps(self):
        pass
    

    