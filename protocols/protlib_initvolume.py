#!/usr/bin/env xmipp_python
'''
#/***************************************************************************
# * Authors:     C.O.S Sorzano (coss@cnb.csic.es)
# *
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
# ***************************************************************************
'''
# This library contains some common utilities 
# for all particles related protocols: Extract, Import
from protlib_base import *
from xmipp import MetaData, MetaDataInfo
from protlib_utils import runJob, runShowJ

class ProtInitVolumeBase(XmippProtocol):
    '''This class will serve as base for init volume related protocols'''
    def __init__(self, protocolName, scriptname, project):
        XmippProtocol.__init__(self, protocolName, scriptname, project)        
        self.Import = 'from protlib_initvolume import *; '
        self.Xdim=MetaDataInfo(self.Classes)[0]
        
    def defineSteps(self):
        self.insertStep('createDir',path=self.ExtraDir)
        self.insertStep("linkAcquisitionInfo",InputFile=self.Classes,dirDest=self.WorkingDir)
        fnOutputReducedClass = self.extraPath("reducedClasses.xmd")
        fnOutputReducedClassNoExt = os.path.splitext(fnOutputReducedClass)[0]
    
        # Low pass filter and resize        
        self.MaxFreq = float(self.MaxFreq)
        self.Ts = float(self.Ts)
        K = 0.25*(self.MaxFreq/self.Ts)
        Xdim2 = self.Xdim/K
        if (Xdim2 < 32):
            self.Xdim2 = 32
            K = self.Xdim/Xdim2
        else:
            self.Xdim2 = Xdim2
            
        freq = self.Ts/self.MaxFreq
        self.Ts = K*self.Ts

        self.insertRunJobStep("xmipp_transform_filter","-i %s -o %s.stk --save_metadata_stack %s.xmd --fourier low_pass %f"
                                                %(self.Classes,fnOutputReducedClassNoExt,fnOutputReducedClassNoExt,freq))
        self.insertRunJobStep("xmipp_image_resize","-i %s.stk --dim %d %d" %(fnOutputReducedClassNoExt,self.Xdim2,self.Xdim2))
        self.insertRunJobStep("rm", fnOutputReducedClassNoExt+".stk_tmp.xmd",NumberOfMpi=1)
        
    def summary(self):
        message=[]
        message.append("Input images: [%s]"%self.Classes)
        if self.InitialVolume!="":
            message.append("Initial volume: [%s]"%self.InitialVolume)
        message.append("Symmetry: "+self.SymmetryGroup)
        return message

    def validate(self):
        errors = []
        if float(self.MaxFreq)<2*float(self.Ts):
            errors.append("Maximum frequency cannot be smaller than twice the sampling rate")
        return errors
