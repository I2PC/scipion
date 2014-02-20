# **************************************************************************
# *
# * Authors:     Josue Gomez Blanco (jgomez@cnb.csic.es)
# *              Vahid Abrishami    (vabrishami@cnb.csic.es)
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
# *  e-mail address 'jgomez@cnb.csic.es'
# *
# **************************************************************************
"""
This sub-package contains wrapper around Initial Volume Simulated Annealing Xmipp program
"""


from pyworkflow.em import *  
from pyworkflow.utils import *  
import xmipp
from convert import createXmippInputClasses2D, createXmippInputVolumes
        
        
class XmippProtInitVolSimAnneal(ProtInitialVolume):
    """ Protocol for Xmipp-based Initial Volume Simulated Annealing. """
    _label = 'simulated annealing'


    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputClasses', PointerParam, label="Input classes", important=True, 
                      pointerClass='SetOfClasses2D', pointerCondition='hasAverages',
                      help='Select the input classes2D from the project.\n'
                           'It should be a SetOfClasses2D class with class representative')
        form.addParam('symmetryGroup', TextParam, default='c1',
                      label="Symmetry group",
                      help='See [[http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Symmetry][Symmetry]]'
                      'for a description of the symmetry groups format, If no symmetry is present, give c1.')  
        form.addParam('isRefVolume', BooleanParam, default=False,
                      label="There are a Reference Volume?", 
                       help='Set to _True or False_ if you want to use a Reference Volume')
        form.addParam('refVolume', PointerParam, label='Initial 3D reference volumes',
                      pointerClass='SetOfVolumes', condition="isRefVolume",
                      help='You may provide a very rough initial volume as a way to constraint the angular search. \n\n'
                      '_For instance_: when reconstructing a fiber, you may provide a cylinder so that side views \n'
                      'are assigned to the correct tilt angle, although the rotational angle may be completely wrong')
        form.addParam('angularSampling', IntParam, default=5,
                      label='Angular sampling',
                      help='Angular sampling in degrees for generating the projection gallery.')
        form.addParam('numberOfSimAnnealRef', IntParam, default=10,
                      label='Number of simulated annealing iterations',
                      help='During the simulated annealing iterations, all those particles \n'
                      'positively contributing to the improvement of the volume are considered. \n'
                      'In this way, the same image may participate several times from different \n'
                      'projection directions (but different weights) depending \n'
                      'on whether it improves the correlation with the volume or not')
        form.addParam('initTemperature', FloatParam, default=0.1,
                      label='Initial temperature',
                      help='The initial temperature determines whether wrong orientations are\n'
                      'considered or not. At the beginning of the simulated annealing\n'
                      'iterations it is important to take wrong directions in order to avoid\n'
                      'local minima.\n'
                      'You may increase this value and the number of simulated annealing iterations\n'
                      'as a way to have a slower cooling.')          
        form.addParam('DontApplyPositiveConstraint', BooleanParam, default=False,
                      label="Don't apply non-negative constraints?", 
                      help='In between simulated annealing iterations, the reconstructed \n'
                      'volume is constrained to have non-negative values. This helps \n'
                      'to reduce the search space, but it might be incorrect depending \n'
                      'on the normalization of the input images')
        form.addParam('numIterGreedy', IntParam, default=2,
                      label='Number of greedy iterations',
                      help='After simulating annealing, you may run some greedy \n'
                      'iterations by simply applying a projection matching approach.')
        form.addParam('percentRejection', IntParam, default=50,
                      label='Percentage of rejected particles',
                      help='At each iteration, the lowest correlated particles are\n'
                      'removed from the 3D reconstruction, although they may participate\n'
                      'in the next iteration')    
        
        form.addParallelSection(threads=2, mpi=1)    


    #--------------------------- INSERT steps functions --------------------------------------------
    def _prepareParams(self):
        """ Step to define the protocol params """
        self._params = {'symmetry': self.symmetryGroup.get(),
                        'angSampling': self.angularSampling.get(),
                        'numSimAnnealRef': self.numberOfSimAnnealRef.get(),
                        'temperature': self.initTemperature.get(),
                        'numIterGreedy':self.numIterGreedy.get(),
                        'reject':self.percentRejection.get()
                       }

        # Convert input images if necessary
        classFn = createXmippInputClasses2D(self, self.inputClasses.get())
        self._params['classFn'] = classFn
        
        outVolFn = self._getExtraPath("output_volume")
        self._params['volOutput'] = outVolFn
    
    def _insertAllSteps(self):
        self._prepareParams()
        self._insertSteps()
    
    def _insertSteps(self):
        self._insertInitialVolumeStep()
        self._insertFunctionStep('createOutputStep')
    
    def _insertInitialVolumeStep(self):
        """Prepare the command line and run job step"""
        args='-i %(classFn)s --oroot %(volOutput)s --sym %(symmetry)s --randomIter %(numSimAnnealRef)d --greedyIter %(numIterGreedy)d --rejection %(reject)f --keepIntermediateVolumes --T0 %(temperature)f --angularSampling %(angSampling)f'
        
        if self.DontApplyPositiveConstraint:
            args += " --dontApplyPositive"
        
        if self.isRefVolume.get():
            # convert the set of Volumes into Xmipp enviroment
            volSet = self.refVolume.get()
            refVolMd = createXmippInputVolumes(self, volSet)
            self._params['volume'] = refVolMd
            args += " --initial %(volume)s"
        self._insertRunJobStep("xmipp_volume_initial_simulated_annealing", args % self._params)
    
    #--------------------------- STEPS functions ---------------------------------------------------
    def createOutputStep(self):
        classes2DSet = self.inputClasses.get()
        
        # create a SetOfVolumes
        volumes = self._createSetOfVolumes()
        volumes.setSamplingRate(self.inputClasses.get().getAverages().getSamplingRate())
        
        # rename the projection file to reconstruct the last Volume.
        moveFile(self._getExtraPath('output_volume.xmd'), self._getExtraPath('volume_projections.xmd'))
        iters = self.numberOfSimAnnealRef.get() + self.numIterGreedy.get()
        
        # create a Volume object to store the final volume obtained
        finalVolFn = self._getExtraPath('output_volume.vol')
        finalVol = Volume()
        finalVol.setFileName(finalVolFn)
        
        # create a SetOfVolumes with volumes created in all iterations
        for k in range(1, iters + 1):
            volFn = self._getExtraPath('output_volume_iter%02d.vol' % k)
            vol = Volume()
            vol.setFileName(volFn)
            volumes.append(vol)
        
        self._defineOutputs(outputVolume=finalVol)
        self._defineOutputs(outputVolumes=volumes)
        self._defineSourceRelation(classes2DSet, finalVol)
        self._defineSourceRelation(classes2DSet, volumes)
    
    #--------------------------- INFO functions ----------------------------------------------------
    def _validate(self):
        validateMsgs = []
        # if there are Volume references, it cannot be empty.
        if self.isRefVolume.get() and not self.refVolume.hasValue():
            validateMsgs.append('Please provide an initial reference volume.')
        return validateMsgs
    
    def _citations(self):
        """ there aren't citations for this protocol """
        pass
    
    def _summary(self):
        summary = []
        summary.append("Input images:  %s" % self.inputClasses.get().getNameId())
        if self.isRefVolume.get():
            summary.append("Reference volumes(s): [%s]" % self.refVolume.get())
            
        if not hasattr(self, 'outputVolumes'):
            summary.append("Output volumes not ready yet.")
        else:
            summary.append("Output volumes: %s" % self.outputVolumes.getNameId())
        return summary
    
    def _methods(self):
        return self._summary()  # summary is quite explicit and serve as methods
