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
from pyworkflow.em.packages.xmipp3.convert import locationToXmipp
"""
This sub-package contains protocols for performing subtomogram averaging.
"""

from pyworkflow.em import *  
from constants import *
from xmipp import MetaData, MDL_ANGLE_ROT, MDL_SHIFT_Z

class XmippProtHelicalParameters(ProtPreprocessVolumes):
    """ Estimate helical parameters and symmetrize.
    
         Helical symmetry is defined as V(r,rot,z)=V(r,rot+k*DeltaRot,z+k*Deltaz)."""
    _label = 'helical symmetry'
    
    def _defineParams(self, form):
        form.addSection(label='General parameters')
        form.addParam('inputVolume', PointerParam, pointerClass="Volume", label='Input volume')
        form.addParam('cylinderRadius',IntParam,label='Cylinder radius', default=-1,
                      help="The helix is supposed to occupy the this radius around the Z axis. Leave it as -1 for symmetrizing the whole volume")
        form.addParam('dihedral',BooleanParam,default=False,label='Apply dihedral symmetry')

        form.addSection(label='Search limits')
        form.addParam('rot0',FloatParam,default=0,label='Minimum rotational angle',help="In degrees")
        form.addParam('rotF',FloatParam,default=360,label='Maximum tilt angle',help="In degrees")
        form.addParam('rotStep',FloatParam,default=5,label='Angular step',help="In degrees")
        form.addParam('z0',FloatParam,default=0,label='Minimum shift Z',help="In voxels")
        form.addParam('zF',FloatParam,default=10,label='Maximum shift Z',help="In voxels")
        form.addParam('zStep',FloatParam,default=0.5,label='Shift step',help="In voxels")
        self.deltaZ=Float()
        self.deltaRot=Float()

    def _insertAllSteps(self):
        self._insertFunctionStep('coarseSearch')
        self._insertFunctionStep('fineSearch')
        self._insertFunctionStep('symmetrize')
        if self.dihedral.get():
            self._insertFunctionStep('applyDihedral')
        self._insertFunctionStep('createOutput')
    
    def createOutput(self):
        volume=Volume()
        volume.setFileName(self._getPath('volume_symmetrized.vol'))
        volume.setSamplingRate(self.inputVolume.get().getSamplingRate())
        self._defineOutputs(outputVolume=volume)
        self._defineTransformRelation(self.inputVolume, self.outputVolume)
        
        md=MetaData(self._getExtraPath('fineParams.xmd'))
        objId=md.firstObject()
        self.deltaRot.set(md.getValue(MDL_ANGLE_ROT,objId))
        self.deltaZ.set(md.getValue(MDL_SHIFT_Z,objId))
        
    def _summary(self):
        messages = []
        if self.deltaZ.hasValue():
            messages.append('DeltaZ=%f (voxels) %f (Angstroms)'%(self.deltaZ.get(),self.deltaZ.get()*self.inputVolume.get().getSamplingRate()))
            messages.append('DeltaRot=%f (degrees)'%self.deltaRot.get())      
        return messages

    def _citations(self):
        papers=[]
        return papers

    def _methods(self):
        messages = []      
        return messages

    def coarseSearch(self):
        fnVol=locationToXmipp(*self.inputVolume.get().getLocation())
        args="-i %s --sym helical -z %f %f %f --rotHelical %f %f %f --thr %d -o %s"%(fnVol,
                                                                                     float(self.z0.get()),float(self.zF.get()),
                                                                                     float(self.zStep.get()),float(self.rot0.get()),
                                                                                     float(self.rotF.get()),float(self.rotStep.get()),
                                                                                     self.numberOfThreads.get(),
                                                                                     self._getExtraPath('coarseParams.xmd'))
        if self.cylinderRadius.get()>0:
            [xdim,_,_,_]=self.inputVolume.get().getDim()
            args+=" --mask cylinder %d %d"%(int(-self.cylinderRadius.get()),int(-xdim))
        self.runJob('xmipp_volume_find_symmetry',args)

    def fineSearch(self):
        md=MetaData(self._getExtraPath('coarseParams.xmd'))
        objId=md.firstObject()
        rot0=md.getValue(MDL_ANGLE_ROT,objId)
        z0=md.getValue(MDL_SHIFT_Z,objId)
        fnVol=locationToXmipp(*self.inputVolume.get().getLocation())
        args="-i %s --sym helical --localHelical %f %f -o %s"%(fnVol,z0,rot0,self._getExtraPath('fineParams.xmd'))
        if self.cylinderRadius.get()>0:
            [xdim,_,_,_]=self.inputVolume.get().getDim()
            args+=" --mask cylinder %d %d"%(int(-self.cylinderRadius.get()),int(-xdim))
        self.runJob('xmipp_volume_find_symmetry',args)

    def symmetrize(self):
        md=MetaData(self._getExtraPath('fineParams.xmd'))
        objId=md.firstObject()
        rot0=md.getValue(MDL_ANGLE_ROT,objId)
        z0=md.getValue(MDL_SHIFT_Z,objId)
        fnOut=self._getPath('volume_symmetrized.vol')
        fnVol=locationToXmipp(*self.inputVolume.get().getLocation())
        args="-i %s --sym helical --helixParams %f %f -o %s"%(fnVol,z0,rot0,fnOut)
        self.runJob('xmipp_transform_symmetrize',args)
        if self.cylinderRadius.get()>0:
            [xdim,_,_,_]=self.inputVolume.get().getDim()
            args="-i %s --mask cylinder %d %d"%(fnOut,int(-self.cylinderRadius.get()),int(-xdim))
            self.runJob('xmipp_transform_mask',args)

    def applyDihedral(self):
        fnRotated=self._getTmpPath('volume_rotated.vol')
        fnOut=self._getPath('volume_symmetrized.vol')
        self.runJob("xmipp_transform_geometry","-i %s -o %s --rotate_volume axis 180 1 0 0"%(fnOut,fnRotated))
        if self.cylinderRadius.get()>0:
            [xdim,_,_,_]=self.inputVolume.get().getDim()
            maskArgs=" --mask cylinder %d %d"%(int(-self.cylinderRadius.get()),int(-xdim))
        else:
            maskArgs=""
        self.runJob("xmipp_volume_align","--i1 %s --i2 %s --rot 0 360 3 -z -2 2 0.5 --apply"%(fnOut,fnRotated)+maskArgs)
        self.runJob("xmipp_volume_align","--i1 %s --i2 %s --local --apply"%(fnOut,fnRotated)+maskArgs)
        self.runJob("xmipp_image_operate","-i %s --plus %s -o %s"%(fnOut,fnRotated,fnOut))
        self.runJob("xmipp_image_operate","-i %s --divide 2"%(fnOut))
        self.runJob("xmipp_transform_mask","-i %s "%(fnOut)+maskArgs)
