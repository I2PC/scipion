# **************************************************************************
# *
# * Authors:     Javier Vargas and Adrian Quintana (jvargas@cnb.csic.es aquintana@cnb.csic.es)
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
from pyworkflow.em.packages.xmipp3.convert import createXmippInputClasses2D,\
    createXmippInputVolumes
"""
This sub-package contains wrapper around Screen Classes Xmipp program
"""

from pyworkflow.em import *  
import xmipp
from xmipp3 import ProjMatcher
from convert import readSetOfClasses2D   

        
class XmippProtScreenClasses(ProtAlignClassify, ProjMatcher):
    """ Protocol to screen a set of classes in the project using a volume as reference """
    _label = 'screen classes'
    
    def __init__(self, **args):
        ProtAlignClassify.__init__(self, **args)

    def _defineParams(self, form):
        form.addSection(label='Input')
        
        form.addParam('inputClasses', PointerParam, label="Set of classes", important=True,
                      pointerClass='SetOfClasses2D', pointerCondition='hasAverages',
                      help='Provide a set of classes object')
        form.addParam('inputVolume', PointerParam, label="Volume to compare classes to", important=True,
#                      pointerClass='SetOfVolumes',
                      pointerClass='Volume',
                      help='Volume to be used for class comparison')
        form.addParam('symmetryGroup', StringParam, default="c1",
                      label='Symmetry group', 
                      help='See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Symmetry for a description of the symmetry groups format'
                        'If no symmetry is present, give c1')
        form.addParam('angularSampling', FloatParam, default=5, expertLevel=LEVEL_EXPERT,
                      label='Angular sampling rate',
                      help='In degrees.'
                      ' This sampling defines how fine the projection gallery from the volume is explored.')
        
        form.addParallelSection(mpi=8)
        
    def _insertAllSteps(self):
        """ Mainly prepare the command line for call cl2d program"""
        # Convert input images if necessary
        self.Xdim = self.inputClasses.get().getDimensions()[0]
        
        #TODO: This should be deleted when inputVolume was set to Volume type
        #self.volume = self.inputVolume.get().getFirstItem()
        self.volume = self.inputVolume.get()
        
        self.fn = createXmippInputClasses2D(self, self.inputClasses.get())
        self.visualizeInfoOutput = String("classes_aligned@%s" % self.fn)  
        
        self.fnAngles = self._getExtraPath('angles.xmd')
        self.images = "classes@%s" % self.fn
        self._insertFunctionStep("projMatchStep",\
                                 self.volume.getFileName(), self.angularSampling.get(),\
                                 self.symmetryGroup.get(), self.images,\
                                 self.fnAngles, self.Xdim)
        
        # Reorganize output and produce difference images 
        self._insertRunJobStep("xmipp_metadata_utilities", "-i classes@%s --set join %s --mode append" % (self.fn, self.fnAngles), numberOfMpi=1)
        self._insertFunctionStep("produceAlignedImagesStep", False, self.fn, self.images)
        self._insertRunJobStep("xmipp_metadata_utilities", "-i classes_aligned@%s --operate sort maxCC desc --mode append" % (self.fn), numberOfMpi=1)  
                
    def getVisualizeInfo(self):
        return self.visualizeInfoOutput
     
    def _summary(self):
        summary = []
        summary.append("Set of classes: [%s] " % self.inputClasses.get().getNameId())
        summary.append("Volume: [%s] " % self.inputVolume.getNameId())
        summary.append("Symmetry: %s " % self.symmetryGroup.get())
        return summary
    
    def _methods(self):
        
        methods = []
        methods.append("Set of classes: [%s] " % self.inputClasses.get().getNameId())
        methods.append("Volume: [%s] " % self.inputVolume.getNameId())
        methods.append("Symmetry: %s " % self.symmetryGroup.get())       
        methods.append("angularSampling: %s " % self.angularSampling.get())
        
        return methods
