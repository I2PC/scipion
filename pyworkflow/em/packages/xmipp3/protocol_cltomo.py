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
"""
This sub-package contains protocols for performing subtomogram averaging.
"""

from pyworkflow.em import *  
from constants import *
from convert import createXmippInputVolumes, readSetOfClasses3D, readSetOfVolumes
from xmipp import MetaData


class XmippProtCLTomo(ProtClassify3D):
    """ Perform subtomogram averaging """
    _label = 'cltomo'
    
    def _defineParams(self, form):
        form.addSection(label='General parameters')
        form.addParam('volumelist', PointerParam, pointerClass="SetOfVolumes", label='Set of volumes',
                      help="Set of volumes to align")
        form.addParam('numberOfReferences',IntParam,label='Number of references', default=3,
                      help="How many classes are computed at the end of the process")
        form.addParam('numberOfIterations',IntParam,label='Number of iterations', default=15,
                      expertLevel=LEVEL_ADVANCED,help="How many iterations at each of the Clustering levels")
        form.addParam('generateAligned',BooleanParam,default=True,label='Generate aligned subvolumes',
                      help="If set to true, it will be created a new set of volumes with all of them aligned")
        form.addParam('dontAlign',BooleanParam,default=False,label="Do not align",
                      help="Volumes are already aligned, only classify")
        
        form.addSection(label='Initial classes')
        form.addParam('doGenerateInitial',BooleanParam,default=True,label='Generate initial volume',
                      help="Let CLTomo to automatically generate the initial classes")
        form.addParam('numberOfReferences0',IntParam,label='Number of initial', default=1, condition="doGenerateInitial",
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
        form.addParam('inputMask', PointerParam, pointerClass="VolumeMask", label="Spatial mask", allowNull=True)
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
        form.addParam('maxShiftX',FloatParam,default=360,label='Maximum shift X',help="In voxels")
        form.addParam('maxShiftY',FloatParam,default=360,label='Maximum shift Y',help="In voxels")
        form.addParam('maxShiftZ',FloatParam,default=360,label='Maximum shift Z',help="In voxels")

    def _insertAllSteps(self):
        self._insertFunctionStep('runCLTomo')
        self._insertFunctionStep('createOutput')
    
    def createOutput(self):
        import glob
        levelFiles=glob.glob(self._getExtraPath("results_classes_level*.xmd"))
        if levelFiles:
            levelFiles.sort()
            lastLevelFile=levelFiles[-1]
            setOfClasses = self._createSetOfClasses3D()
            setOfClasses.setVolumes(self.volumelist.get())
            readSetOfClasses3D(setOfClasses,lastLevelFile)
            self._defineOutputs(outputClasses=setOfClasses)
            self._defineSourceRelation(self.volumelist, self.outputClasses)
        if self.generateAligned.get():
            setOfVolumes = self._createSetOfVolumes()
            fnAligned = self._getExtraPath('results_aligned.xmd')
            self.runJob('xmipp_metadata_selfile_create', '-p %s -o %s -s'%(self._getExtraPath('results_aligned.stk'),fnAligned))
            md=MetaData(fnAligned)
            md.addItemId()
            md.write(fnAligned)
            readSetOfVolumes(fnAligned,setOfVolumes)
            volumeList=self.volumelist.get()
            setOfVolumes.setSamplingRate(volumeList.getSamplingRate())
            setOfVolumes.write()
            self._defineOutputs(alignedVolumes=setOfVolumes)
            self._defineTransformRelation(self.volumelist, self.alignedVolumes)
            
    def _summary(self):
        messages = []
        if self.doGenerateInitial.get():
            messages.append('Number of initial references: %d'%self.numberOfReferences0.get())
            if self.randomizeOrientation.get():
                messages.append('Input subvolume orientations were randomized')
        if self.dontAlign.get():
            messages.append('Input subvolumes were assumed to be already aligned')
        messages.append('Number of output classes: %d'%self.numberOfReferences.get())
        return messages

    def _validate(self):
        errors=[]
        (Xdim1, Ydim1, Zdim1, _)=self.volumelist.get().getDimensions()
        if Xdim1!=Ydim1 or Ydim1!=Zdim1:
            errors.append("Input subvolumes are not cubic")
        N0=-1
        if not self.doGenerateInitial.get():
            if not self.referenceList.hasValue():
                errors.append("If references are not self generated, you have to provide a reference set of volumes")
            else:
                (Xdim2, Ydim2, Zdim2, N0) = self.referenceList.get().getDimensions()
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

    def runCLTomo(self):
        params= ' -i '            + createXmippInputVolumes(self,self.volumelist.get()) + \
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
            params+=' --ref0 '+createXmippInputVolumes(self,self.referenceList.get(),self._getPath('references.xmd'))
        if self.inputMask.hasValue():
            params+=' --mask binary_file '+self.inputMask.get().getLocation()
        if self.generateAligned.get():
            params+=" --generateAlignedVolumes"
        if self.dontAlign.get():
            params+=" --dontAlign"

        self.runJob('xmipp_mpi_classify_CLTomo','%d %s'%(self.numberOfMpi.get(),params))
