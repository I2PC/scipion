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

import pyworkflow
import pyworkflow.object as pwobj
from pyworkflow.em import *  
from xmipp import MetaData, MDL_ANGLE_ROT, MDL_ANGLE_TILT
from pyworkflow.em.packages.xmipp3.convert import getImageLocation
from pyworkflow.protocol.constants import LEVEL_ADVANCED


class XmippProtRotationalSymmetry(ProtPreprocessVolumes):
    """ Estimate the orientation of a rotational axis and symmetrize.
    """
    _label = 'rotational symmetry'
    
    GLOBAL_SEARCH = 0
    LOCAL_SEARCH = 1
    GLOBAL_LOCAL_SEARCH = 2

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label='General parameters')
        form.addParam('inputVolume', PointerParam, pointerClass="Volume", label='Input volume')
        form.addParam('symOrder',IntParam, default=3,label='Symmetry order', help="3 for a three-fold symmetry axis, 4 for a four-fold symmetry axis, ...")
        form.addParam('searchMode', EnumParam, choices=['Global','Local','Global+Local'], default=self.GLOBAL_LOCAL_SEARCH)
        
        form.addParam('rot',FloatParam,default=0,label='Initial rotational angle',condition='searchMode==1', help="In degrees")
        form.addParam('tilt',FloatParam,default=0,label='Initial tilt angle',condition='searchMode==1', help="In degrees")

        form.addParam('rot0',FloatParam,default=0,label='Minimum rotational angle',condition='searchMode!=1', help="In degrees")
        form.addParam('rotF',FloatParam,default=360,label='Maximum rotational angle',condition='searchMode!=1', help="In degrees")
        form.addParam('rotStep',FloatParam,default=5,label='Angular step',condition='searchMode!=1', help="In degrees")
        form.addParam('tilt0',FloatParam,default=0,label='Minimum tilt angle',condition='searchMode!=1', help="In degrees")
        form.addParam('tiltF',FloatParam,default=180,label='Maximum tilt angle',condition='searchMode!=1', help="In degrees")
        form.addParam('tiltStep',FloatParam,default=5,label='Tilt step',condition='searchMode!=1', help="In degrees")
        self.rotSym=Float()
        self.tiltSym=Float()

        form.addParallelSection(threads=4, mpi=0)

    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('copyInput')
        if self.searchMode.get()!=self.LOCAL_SEARCH:
            self._insertFunctionStep('coarseSearch')
        if self.searchMode.get()!=self.GLOBAL_SEARCH:
            self._insertFunctionStep('fineSearch')
        self._insertFunctionStep('symmetrize')
        self._insertFunctionStep('createOutput')
        self.fnVol = getImageLocation(self.inputVolume.get())
        self.fnVolSym=self._getPath('volume_symmetrized.vol')
        [self.height,_,_]=self.inputVolume.get().getDim()
    
    #--------------------------- STEPS functions --------------------------------------------
    def copyInput(self):
        ImageHandler().convert(self.inputVolume.get(), self.fnVolSym)
                        
    def coarseSearch(self):
        self.runJob("xmipp_volume_find_symmetry","-i %s -o %s --rot %f %f %f --tilt %f %f %f --sym rot %d --thr %d"%
                    (self.fnVolSym,self._getExtraPath('coarse.xmd'),self.rot0.get(), self.rotF.get(), self.rotStep.get(),
                     self.tilt0.get(), self.tiltF.get(), self.tiltStep.get(), self.symOrder.get(),self.numberOfThreads))

    def getAngles(self, fnAngles=""):
        if fnAngles=="":
            if self.searchMode.get()==self.GLOBAL_SEARCH:
                fnAngles = self._getExtraPath("coarse.xmd")
            else:
                fnAngles = self._getExtraPath("fine.xmd")
        md = MetaData(fnAngles)
        objId = md.firstObject()
        rot0 = md.getValue(MDL_ANGLE_ROT, objId)
        tilt0 = md.getValue(MDL_ANGLE_TILT, objId)
        return (rot0,tilt0)

    def fineSearch(self):
        if self.searchMode.get()==self.LOCAL_SEARCH:
            rot0 = self.rot.get()
            tilt0 = self.tilt.get()
        else:
            rot0, tilt0 = self.getAngles(self._getExtraPath('coarse.xmd'))
        self.runJob("xmipp_volume_find_symmetry","-i %s -o %s --localRot %f %f --sym rot %d"%
                    (self.fnVolSym,self._getExtraPath('fine.xmd'),rot0, tilt0, self.symOrder.get()))

    def symmetrize(self):
        rot0, tilt0 = self.getAngles()
        self.runJob("xmipp_transform_geometry","-i %s --rotate_volume euler %f %f 0 --dont_wrap"%(self.fnVolSym,rot0,tilt0))

    def createOutput(self):
        volume = Volume()
        volume.setFileName(self.fnVolSym)
        volume.copyInfo(self.inputVolume.get())
        self._defineOutputs(outputVolume=volume)
        self._defineTransformRelation(self.inputVolume, self.outputVolume)
        
        rot0, tilt0 = self.getAngles()
        self._defineOutputs(rotSym=pwobj.Float(rot0),tiltSym=pwobj.Float(tilt0))

    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        messages = []
        if self.rotSym.hasValue():
            messages.append('Rot. Angle of Symmetry axis=%f (degrees)'%self.rotSym.get())
            messages.append('Tilt.Angle of Symmetry axis=%f (degrees)'%self.tiltSym.get())
        return messages

    def _methods(self):
        messages = []      
        messages.append('We looked for the %d-fold rotational axis of the volume %s using Xmipp [delaRosaTrevin2013]. ' % (self.symOrder.get(),self.getObjectTag('inputVolume'))+
                        'We found it to be with an orientation given by a rotational angle of %f and a tilt angle of %f degrees.'
                        %(self.rotSym.get(),self.tiltSym.get()))
        return messages

