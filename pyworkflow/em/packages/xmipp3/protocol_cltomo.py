# **************************************************************************
# *
# * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
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
# *  e-mail address 'xmipp@cnb.csic.es'
# *
# **************************************************************************

from pyworkflow.em import *  
from constants import *
from convert import writeSetOfVolumes, readSetOfClassesVol, readSetOfVolumes
from xmipp import MetaData
from xmipp3 import getEnviron


class XmippProtCLTomo(ProtClassify3D):
    """ Averages a set of subtomograms taking into account the missing edge. """
    _label = 'cltomo'
    
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label='General parameters')
        form.addParam('inputVolumes', PointerParam, pointerClass="SetOfVolumes", label='Set of volumes',
                      help="Set of volumes to align")
        form.addParam('numberOfReferences',IntParam,label='Number of references', default=3,
                      help="How many references are computed at the end of the process")
        form.addParam('numberOfIterations',IntParam,label='Number of iterations', default=15,
                      expertLevel=LEVEL_ADVANCED,help="How many iterations at each of the Clustering levels")
        form.addParam('generateAligned',BooleanParam,default=True,label='Generate aligned subvolumes',
                      help="If set to true, it will be created a new set of volumes with all of them aligned")
        form.addParam('align',BooleanParam,default=True,label="Align",
                      help="Do not align if volumes are already aligned, only classify")
        
        form.addSection(label='Initial references')
        form.addParam('doGenerateInitial',BooleanParam,default=True,label='Generate initial volume',
                      help="Let CLTomo to automatically generate the initial references")
        form.addParam('numberOfReferences0',IntParam,label='Number of initial references', default=1, condition="doGenerateInitial",
                      help="How many initial volumes. If set to 1, all subvolumes are aligned to a single reference, "\
                           "and then they are classified")
        form.addParam('randomizeOrientation',BooleanParam,default=False,label='Randomize orientation', condition="doGenerateInitial",
                      help="Use this option if all the input volumes have the same missing wedge or if they have not been previously aligned.")
        form.addParam('referenceList', PointerParam, pointerClass="SetOfVolumes", label='Set of initial volumes',
                      condition="not doGenerateInitial", help="Set of initial volumes")
        
        form.addSection(label='Constraints')
        form.addParam('symmetry',StringParam,default='c1',label='Symmetry group',
                      help="See http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Symmetry for a description of the symmetry groups format."
                           "If no symmetry is present, give c1")
        form.addParam('inputMask', PointerParam, pointerClass="VolumeMask", label="Spatial mask", allowsNull=True)
        form.addParam('maximumResolution',FloatParam,default=0.25,label='Maximum resolution (pixels^-1)',
                      help="The maximum (Nyquist) resolution is 0.5. Use smaller values, e.g. 0.45, to prevent high-resolution artifacts.")
        form.addParam('sparsity',FloatParam,default=90,label='Sparsity in Fourier space',
                      help="A value of 90 drops 90% of the smallest Fourier coefficients")
        form.addParam('dwtSparsity',FloatParam,default=90,label='Sparsity in wavelet space',
                      help="A value of 95 drops 95% of the smallest wavelet coefficients")

        form.addSection(label='Search limits')
        form.addParam('maxRot',FloatParam,default=360,label='Maximum rotational angle',help="In degrees")
        form.addParam('maxTilt',FloatParam,default=360,label='Maximum tilt angle',help="In degrees")
        form.addParam('maxPsi',FloatParam,default=360,label='Maximum in-plane angle',help="In degrees")
        form.addParam('maxShiftX',FloatParam,default=10,label='Maximum shift X',help="In voxels")
        form.addParam('maxShiftY',FloatParam,default=10,label='Maximum shift Y',help="In voxels")
        form.addParam('maxShiftZ',FloatParam,default=10,label='Maximum shift Z',help="In voxels")
        form.addParallelSection(threads=0, mpi=4)

    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('runCLTomo')
        self._insertFunctionStep('createOutput')
    
    #--------------------------- STEPS functions --------------------------------------------

    def runCLTomo(self):
        fnVols=self._getPath('input_volumes.xmd')
        writeSetOfVolumes(self.inputVolumes.get(), fnVols)
        params= ' -i '            + fnVols + \
                ' --oroot '       + self._getExtraPath("results") + \
                ' --iter '        + str(self.numberOfIterations.get()) + \
                ' --nref '        + str(self.numberOfReferences.get()) + \
                ' --sym '         + self.symmetry.get() + \
                ' --maxFreq '     + str(self.maximumResolution.get()) + \
                ' --sparsity '    + str(self.sparsity.get()/100.0) + \
                ' --DWTsparsity ' + str(self.dwtSparsity.get()/100.0) + \
                ' --maxShiftX '   + str(self.maxShiftX.get()) + \
                ' --maxShiftY '   + str(self.maxShiftY.get()) + \
                ' --maxShiftZ '   + str(self.maxShiftZ.get()) + \
                ' --maxRot '      + str(self.maxRot.get()) + \
                ' --maxTilt '     + str(self.maxTilt.get()) + \
                ' --maxPsi '      + str(self.maxPsi.get())
        if self.doGenerateInitial.get():
            params+=' --nref0 '+str(self.numberOfReferences0.get())
            if self.randomizeOrientation.get():
                params+=' --randomizeStartingOrientation'
        else:
            fnInitialVols=self._getExtraPath('intial_volumes.xmd')
            writeSetOfVolumes(self.referenceList.get(), fnInitialVols)
            params+=' --ref0 '+fnInitialVols
        if self.inputMask.hasValue():
            params+=' --mask binary_file '+self.inputMask.get().getLocation()
        if self.generateAligned.get():
            params+=" --generateAlignedVolumes"
        if not self.align:
            params+=" --dontAlign"

        self.runJob('xmipp_mpi_classify_CLTomo','%d %s'%(self.numberOfMpi.get(),params),env=self.getCLTomoEnviron(),numberOfMpi=1)
    
    def createOutput(self):
        import glob
        levelFiles=glob.glob(self._getExtraPath("results_classes_level*.xmd"))
        if levelFiles:
            levelFiles.sort()
            lastLevelFile=levelFiles[-1]
            setOfClasses = self._createSetOfClassesVol()
            setOfClasses.setImages(self.inputVolumes.get())
            readSetOfClassesVol(setOfClasses,lastLevelFile)
            self._defineOutputs(outputClasses=setOfClasses)
            self._defineSourceRelation(self.inputVolumes, self.outputClasses)
        if self.generateAligned:
            setOfVolumes = self._createSetOfVolumes()
            fnAligned = self._getExtraPath('results_aligned.xmd')
            readSetOfVolumes(fnAligned,setOfVolumes)
            volumeList = self.inputVolumes.get()
            setOfVolumes.setSamplingRate(volumeList.getSamplingRate())
            self._defineOutputs(alignedVolumes=setOfVolumes)
            self._defineTransformRelation(self.inputVolumes, self.alignedVolumes)
    
    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        messages = []
        if self.doGenerateInitial.get():
            messages.append('Number of initial references: %d'%self.numberOfReferences0.get())
            if self.randomizeOrientation.get():
                messages.append('Input subvolume orientations were randomized')
        if self.dontAlign.get():
            messages.append('Input subvolumes were assumed to be already aligned')
        messages.append('Number of output references: %d'%self.numberOfReferences.get())
        return messages

    def _validate(self):
        errors=[]
        (Xdim1, Ydim1, Zdim1)=self.inputVolumes.get().getDimensions()
        if Xdim1!=Ydim1 or Ydim1!=Zdim1:
            errors.append("Input subvolumes are not cubic")
        N0=-1
        if not self.doGenerateInitial.get():
            if not self.referenceList.hasValue():
                errors.append("If references are not self generated, you have to provide a reference set of volumes")
            else:
                (Xdim2, Ydim2, Zdim2) = self.referenceList.get().getDimensions()
                if Xdim2!=Ydim2 or Ydim2!=Zdim2:
                    errors.append("Reference subvolumes are not cubic")
                if Xdim1!=Xdim2:
                    errors.append("Input and reference subvolumes are of different size")
        else:
            N0=self.numberOfReferences0.get()
        if N0>0 and N0>self.numberOfReferences.get():
            errors.append("The number of initial references have to be smaller or equal than the number of final references")            
        return errors

    def _citations(self):
        return ['Chen2013']

    #--------------------------- UTILS functions --------------------------------------------
    def getCLTomoEnviron(self):
        env = getEnviron()
        env.set('PYTHONPATH', os.path.join(os.environ['SCIPION_HOME'], 'software','lib','python2.7','site-packages','sh_alignment'),
                Environ.BEGIN)
        return env
