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
from pyworkflow.em.packages.xmipp3.convert import createXmippInputClasses2D
"""
This sub-package contains wrapper around Screen Classes Xmipp program
"""

from pyworkflow.em import *  
import xmipp
from protlib_projmatch import projMatch,produceAlignedImages

        
        
class XmippProtScreenClasses(ProtAlignClassify):
    """ Protocol to screen a set of classes in the project. """
    _label = 'screen classes'
    
    def __init__(self, **args):
        ProtAlignClassify.__init__(self, **args)

    def _defineParams(self, form):
        form.addSection(label='Input')
        
        form.addParam('inputClasses', PointerParam, label="Set of classes", important=True,
                      pointerClass='SetOfClasses2D', pointerCondition='hasAverages',
                      help='Provide a set of classes object')
        form.addParam('inputVolume', PointerParam, label="Volume to compare classes to", important=True,
                      pointerClass='SetOfVolumes',
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
        self.inputClasses.get().printAll()
        print "taka",self.inputClasses.get().getDimensions()
        self.Xdim = self.inputClasses.get().getDimensions()[0]
        
        #TODO: This should be deleted when inputVolume was set to Volume type
        volume= self.inputVolume.get()[0]
        
        self.fn = createXmippInputClasses2D(self, self.inputClasses.get())
        
        fnAngles=self._getExtraPath('angles.xmd')
        self._insertFunctionStep("projMatch",
                                 Volume=volume, AngularSampling=self.angularSampling.get(), SymmetryGroup=self.symmetryGroup.get(),
                                 Images="classes@%s"%self.fn, ExtraDir=self._getExtraPath(),
                                 fnAngles=fnAngles, NumberOfMpi=self.numberOfMpi.get())
        
        # Reorganize output and produce difference images 
        self._insertRunJobStep("xmipp_metadata_utilities", "-i classes@%s --set join %s --mode append"%(fnOutputClass,fnAngles), numberOfMpi=1)
#        self._insertFunctionStep("runJob",programname="xmipp_metadata_utilities", params="-i classes@%s --set join %s --mode append"%(fnOutputClass,fnAngles),NumberOfMpi=1)  
        
        self._insertFunctionStep("produceAlignedImages",
                                 fnIn='classes@'+self.fn, fnOut='classes_aligned@'+self.fn, fnDiff=self._getExtraPath("diff.stk"),
                                 volumeIsCTFCorrected=False)
        
#        self.insertStep("produceAlignedImages",fnIn='classes@'+fnOutputClass, fnOut='classes_aligned@'+fnOutputClass, fnDiff=self.extraPath("diff.stk"),
#                        volumeIsCTFCorrected=False)
        self._insertRunJobStep("xmipp_metadata_utilities", "-i classes_aligned@%s --operate sort maxCC desc --mode append"%(fnOutputClass), numberOfMpi=1)  
   
                
        self._insertFunctionStep('createOutput')


                        
    def createOutput(self):
        print "output"


    def _summary(self):
        summary = []
        if not hasattr(self, 'outputClasses'):
            summary.append("Output classes not ready yet.")
        else:
            summary.append("Set of classes: [%s] " % self.Classes)
            summary.append("Volume: [%s] " % self.Volume)
            summary.append("Symmetry: %s " % self.SymmetryGroup)
        return summary
