# **************************************************************************
# *
# * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
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

import os

from pyworkflow.object import Float, String
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pyworkflow.protocol.params import PointerParam, StringParam, FloatParam
from pyworkflow.em.protocol.protocol import EMProtocol
from pyworkflow.em.protocol.protocol_3d import ProtAnalysis3D
from pyworkflow.utils import cleanPath
from xmipp3 import getMatlabEnviron
# import xmipp

        
class XmippProtVolumeStrain(ProtAnalysis3D):
    """Compare two states of a volume to analyze the local strains and rotations"""
    _label = 'calculate strain'
    
    def __init__(self, **args):
        ProtAnalysis3D.__init__(self, **args)
    
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        
        form.addParam('inputVolume0', PointerParam, label="Initial state", important=True,
                      pointerClass='Volume',
                      help='Initial state of the structure, it will be deformed to fit into the final state')
        form.addParam('inputVolumeF', PointerParam, label="Final state", important=True,
                      pointerClass='Volume',
                      help='Initial state of the structure, it will be deformed to fit into the final state')
        form.addParam('inputMask', PointerParam, label="Mask for the final state", important=True,
                      pointerClass='VolumeMask',
                      help='Binary mask that defines where the strains and rotations will be calculated')
        form.addParam('symmetryGroup', StringParam, default="c1",
                      label='Symmetry group', 
                      help='See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Symmetry for a description of the symmetry groups format'
                        'If no symmetry is present, give c1')
    
    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        fnVol0 = self.inputVolume0.get().getFileName()
        fnVolF = self.inputVolumeF.get().getFileName()
        fnMask = self.inputMask.get().getFileName()
        self._insertFunctionStep("calculateStrain",fnVol0,fnVolF,fnMask)
        self._insertFunctionStep("prepareOutput")
        self._insertFunctionStep("createChimeraScript")
    
    #--------------------------- STEPS functions ---------------------------------------------------
    def calculateStrain(self, fnVol0, fnVolF, fnMask):
        fnRoot=self._getExtraPath('result')
        mirtDir = os.path.join(os.environ['XMIPP_HOME'], 'external', 'mirt')
        # -wait -nodesktop
        args='''-r "diary('%s'); xmipp_calculate_strain('%s','%s','%s','%s'); exit"'''%(fnRoot+"_matlab.log",fnVolF,fnVol0,fnMask,fnRoot)
        self.runJob("matlab", args, env=getMatlabEnviron(mirtDir))
    
    def prepareOutput(self):
        volDim = self.inputVolume0.get().getDim()[0]
        fnRoot=self._getExtraPath('result')
        self.runJob("xmipp_image_convert", "-i %s_initial.raw#%d,%d,%d,0,float -o %s_initial.vol"%
                    (fnRoot,volDim,volDim,volDim,fnRoot))
        self.runJob("xmipp_transform_mirror","-i %s_initial.vol --flipX"%fnRoot)
        self.runJob("xmipp_image_convert", "-i %s_final.raw#%d,%d,%d,0,float -o %s_final.vol"%
                    (fnRoot,volDim,volDim,volDim,fnRoot))
        self.runJob("xmipp_transform_mirror","-i %s_final.vol --flipX"%fnRoot)
        self.runJob("xmipp_image_convert", "-i %s_initialDeformedToFinal.raw#%d,%d,%d,0,float -o %s_initialDeformedToFinal.vol"%
                    (fnRoot,volDim,volDim,volDim,fnRoot))
        self.runJob("xmipp_transform_mirror","-i %s_initialDeformedToFinal.vol --flipX"%fnRoot)
        self.runJob("xmipp_image_convert", "-i %s_strain.raw#%d,%d,%d,0,float -o %s_strain.vol"%
                    (fnRoot,volDim,volDim,volDim,fnRoot))
        self.runJob("xmipp_transform_mirror","-i %s_strain.vol --flipX"%fnRoot)
        self.runJob("xmipp_image_convert", "-i %s_localrot.raw#%d,%d,%d,0,float -o %s_localrot.vol"%
                    (fnRoot,volDim,volDim,volDim,fnRoot))
        self.runJob("xmipp_transform_mirror","-i %s_localrot.vol --flipX"%fnRoot)
        self.runJob("rm","-f "+self._getExtraPath('result_*.raw'))
        if self.symmetryGroup!="c1":
            self.runJob("xmipp_transform_symmetrize","-i %s --sym %s"%(fnRoot+"_strain.vol",self.symmetryGroup.get()))
            self.runJob("xmipp_transform_symmetrize","-i %s --sym %s"%(fnRoot+"_localrot.vol",self.symmetryGroup.get()))
    
    def createChimeraScript(self):
        fnRoot = "extra/result"
        scriptFile = self._getPath('result') + '_strain_chimera.cmd'
        fhCmd = open(scriptFile, 'w')
        fhCmd.write("open %s\n" % (fnRoot+"_final.vol"))
        fhCmd.write("open %s\n" % (fnRoot+"_strain.vol"))
        fhCmd.write("vol #1 hide\n")
        fhCmd.write("scolor #0 volume #1 cmap rainbow reverseColors True\n")
        fhCmd.close()

        scriptFile = self._getPath('result') + '_localrot_chimera.cmd'
        fhCmd = open(scriptFile, 'w')
        fhCmd.write("open %s\n" % (fnRoot+"_final.vol"))
        fhCmd.write("open %s\n" % (fnRoot+"_localrot.vol"))
        fhCmd.write("vol #1 hide\n")
        fhCmd.write("scolor #0 volume #1 cmap rainbow reverseColors True\n")
        fhCmd.close()

        scriptFile = self._getPath('result') + '_morph_chimera.cmd'
        fhCmd = open(scriptFile, 'w')
        fhCmd.write("open %s\n" % (fnRoot+"_initial.vol"))
        fhCmd.write("open %s\n" % (fnRoot+"_final.vol"))
        fhCmd.write("vol #0 hide\n")
        fhCmd.write("vol #1 hide\n")
        fhCmd.write("vop morph #0,1 frames 50\n")
        fhCmd.close()

    #--------------------------- INFO functions --------------------------------------------
    def _validate(self):
        errors = []
        xdim0 = self.inputVolume0.get().getDim()[0]
        xdimF = self.inputVolumeF.get().getDim()[0]
        if xdim0 != xdimF:
            errors.append("Make sure that the two volumes have the same size")
        return errors    
