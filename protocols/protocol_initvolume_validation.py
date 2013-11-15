#!/usr/bin/env xmipp_python
'''
#/***************************************************************************
# * Authors:     Javier Vargas (jvargas@cnb.csic.es)
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

class ProtInitVolValidate(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.initvolume_validation.name, scriptname, project)
        self.Import = 'from protocol_initvolume_validation import *'
        
    def defineSteps(self):
        self.insertStep('createDir',path=self.ExtraDir)
        self.insertStep("linkAcquisitionInfo",InputFile=self.fnClasses,dirDest=self.WorkingDir)
        self.insertStep("createNoisyImages",WorkingDir=self.WorkingDir, fnClasses=self.fnClasses, 
                        fnProjections=self.fnProjections,SymmetryGroup=self.SymmetryGroup)      
  
    def summary(self):
        message = []
        message.append("Set of classes: [%s] " % self.fnClasses)
        message.append("Volume: [%s] " % self.fnInitialVolume)
        message.append("Symmetry: %s " % self.SymmetryGroup)
        return message
    
    def validate(self):
        errors = []
        return errors    

    def visualize(self):
        #fnAligned = 'classes_aligned@' + self.workingDirPath('classes.xmd')
        #runShowJ(fnAligned, extraParams="--mode metadata --render first")
        os.system('xmipp_chimera_client -i '+self.fnInitialVolume+' --mode projector 256 &')

def createNoisyImages(log,WorkingDir,fnClasses,fnProjections,SymmetryGroup):

    fnOut=os.path.join(WorkingDir,'extra/noisyImages')
    runJob(log,'xmipp_image_operate','-i %s -o %s.stk --reset --save_metadata_stack'%(fnProjections,fnOut))
    runJob(log,'xmipp_transform_add_noise','-i %s.stk --type gaussian 1'%fnOut)
    runJob(log,'xmipp_reconstruct_fourier','-i %s.xmd -o %s.vol --sym %s --weight --max_resolution 0.25'%(fnOut,fnOut,SymmetryGroup))

#def createNoisyVolume(log,WorkingDir,fnClasses,fnProjections,fnInitialVolume):



